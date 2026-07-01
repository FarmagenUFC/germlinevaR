# Tests for the 0.99.2 fix of read.gvr.dual() multi-ALT indel collisions.
#
# The pre-0.99.2 SnpEff-side coord derivation was a verbatim copy of the VEP
# helper but without bcftools-norm-style REF/ALT trimming; on multi-ALT indels
# whose alleles compressed onto the same padded MAF key (e.g. chr1:6095864 in
# freebayes+snpEff+VEP dual output) the Phase 2 join threw:
#
#     "trying to assign N items to M items"
#
# 0.99.2 fixes this by:
#
#   (a) normalize_alleles = TRUE (default): bcftools-norm-style prefix/suffix
#       trimming inside .gvr_coords() before deriving MAF-like coords, so
#       distinct multi-ALT records land on distinct MAF keys.
#
#   (b) A dedupe layer inside .gvr_dual_attach_snpeff() that collapses any
#       residual same-key duplicates (e.g. legitimate two-event collisions in
#       repeat regions) by keeping the highest-impact SnpEff block, with
#       empty-blocks treated as strictly worst.
#
#   (c) A belt-and-suspenders `mult = "first"` on the data.table join that
#       still succeeds even if step (b) is somehow bypassed.
#
# Fixture: tests/testthat/fixtures/example.dual.collision.vcf.gz
#   Record 1 (chr1:6095864 pattern): separates under the fix.
#   Record 2 (chr10:68173832 pattern): STILL collides post-trim (repeat
#     region) so exercises the dedupe layer.

stage_collision_vcf <- function() {
    src <- test_path("fixtures", "example.dual.collision.vcf.gz")
    tmp <- tempfile("gvr_collision_")
    dir.create(tmp)
    file.copy(src, file.path(tmp, basename(src)))
    tmp
}


# -----------------------------------------------------------------------------
# Test 1 : .gvr_coords() with normalize_alleles = TRUE trims to distinct keys
# -----------------------------------------------------------------------------
# Verifies the fix mechanism at the unit level: given the chr1:6095864 REF and
# the two distinct ALTs, bcftools-norm-style trimming produces two DIFFERENT
# MAF-like (Start_Position, Reference_Allele, Tumor_Seq_Allele2) tuples.

test_that(".gvr_coords() with normalize_alleles = TRUE separates chr1:6095864 multi-ALT indel", {
    ref  <- "TCCCCCCCCCTGCCC"
    alt1 <- "TCCCCCCCCCCTGCCC"   # real 1-bp C insertion
    alt2 <- "TCCCCACCCCCTGCCC"   # compound event, same left-padding

    c1 <- germlinevaR:::.gvr_coords(6095864L, ref, alt1, normalize_alleles = TRUE)
    c2 <- germlinevaR:::.gvr_coords(6095864L, ref, alt2, normalize_alleles = TRUE)

    # Distinct MAF keys
    expect_false(identical(
        list(c1$start, c1$ref_allele, c1$tum_allele2),
        list(c2$start, c2$ref_allele, c2$tum_allele2)
    ))

    # Both are insertions
    expect_identical(c1$var_type, "INS")
    expect_identical(c2$var_type, "INS")

    # Concrete expected values: verified against the plan's canonical example
    expect_identical(c1$start,       6095874L)
    expect_identical(c1$ref_allele,  "-")
    expect_identical(c1$tum_allele2, "T")

    expect_identical(c2$start,       6095869L)
    expect_identical(c2$ref_allele,  "-")
    expect_identical(c2$tum_allele2, "C")
})


# -----------------------------------------------------------------------------
# Test 2 : .gvr_coords() with normalize_alleles = FALSE reproduces legacy coords
# -----------------------------------------------------------------------------
# Verifies the escape hatch: the pre-0.99.2 coord-derivation is preserved
# byte-identically when the user requests it (for reproducibility with an
# older analysis run). This is the "off switch" documented in NEWS.

test_that(".gvr_coords() with normalize_alleles = FALSE reproduces pre-0.99.2 collision coords", {
    ref  <- "TCCCCCCCCCTGCCC"
    alt1 <- "TCCCCCCCCCCTGCCC"
    alt2 <- "TCCCCACCCCCTGCCC"

    c1 <- germlinevaR:::.gvr_coords(6095864L, ref, alt1, normalize_alleles = FALSE)
    c2 <- germlinevaR:::.gvr_coords(6095864L, ref, alt2, normalize_alleles = FALSE)

    # Under legacy behaviour BOTH ALTs collapse to the same MAF key
    expect_identical(c1$start,       c2$start)
    expect_identical(c1$ref_allele,  c2$ref_allele)
    expect_identical(c1$tum_allele2, c2$tum_allele2)

    # Concrete legacy values (verified against the pre-fix diagnostics)
    expect_identical(c1$start,       6095878L)
    expect_identical(c1$ref_allele,  "-")
    expect_identical(c1$tum_allele2, "C")
})


# -----------------------------------------------------------------------------
# Test 3 : read.gvr.dual() on the collision fixture completes without error
# -----------------------------------------------------------------------------
# End-to-end: reads a hermetic 2-record dual-annotated VCF that reproduces
# BOTH the coord-separation case (record 1) and the residual-collision case
# (record 2). Pre-fix this call raised the Phase 2 join error.

test_that("read.gvr.dual() on collision fixture returns rows without erroring", {
    fix_dir <- stage_collision_vcf()

    res <- read.gvr.dual(
        folder         = fix_dir,
        add_abraom     = FALSE,
        canonical_only = FALSE,   # test fixture uses hand-written CSQ so CANONICAL may not match "YES"
        min_DP         = NULL,
        min_GQ         = NULL,
        verbose        = FALSE
    )

    expect_true(data.table::is.data.table(res))
    expect_identical(attr(res, "annotator"), "dual")

    # Record 1 contributes 2 rows (distinct keys); record 2 contributes 2 rows
    # (VEP splits per ALT and both attach to the deduped SnpEff row).
    expect_equal(nrow(res), 4L)

    # Record 1: two distinct MAF-normalized keys with distinct genes
    r1 <- res[Chromosome == "chr1"]
    expect_equal(nrow(r1), 2L)
    expect_true(!identical(r1$Start_Position[1], r1$Start_Position[2]) ||
                !identical(r1$Tumor_Seq_Allele2[1], r1$Tumor_Seq_Allele2[2]))
    expect_setequal(r1$Hugo_Symbol, c("RPS6KA1", "PLEKHG5"))

    # SnpEff annotations flowed through the coord-fix path (Phase 2 join)
    expect_true(all(!is.na(r1$snpeff_impact)))
    expect_setequal(r1$snpeff_gene, c("RPS6KA1", "PLEKHG5"))
})


# -----------------------------------------------------------------------------
# Test 4 : dedupe layer catches residual collisions in repeat regions
# -----------------------------------------------------------------------------
# Record 2 of the fixture is a chr10:68173832 pattern where BOTH multi-ALT
# indels normalize to the SAME MAF-like key even after bcftools-norm-style
# trimming (both are 1-bp insertions in a T-homopolymer that trim to
# (68173848, "-", "G")). This exercises .gvr_dual_dedupe_snpeff() and the
# `mult = "first"` join, and populates the audit-trail attribute.

test_that("read.gvr.dual() dedupe path resolves residual collision + populates audit-trail attribute", {
    fix_dir <- stage_collision_vcf()

    res <- read.gvr.dual(
        folder         = fix_dir,
        add_abraom     = FALSE,
        canonical_only = FALSE,
        min_DP         = NULL,
        min_GQ         = NULL,
        verbose        = FALSE
    )

    # Residual-collision record: both ALTs collapse to the same MAF key
    r2 <- res[Chromosome == "chr10"]
    expect_equal(nrow(r2), 2L)
    expect_true(all(r2$Start_Position == 68173848L))
    expect_true(all(r2$Reference_Allele == "-"))
    expect_true(all(r2$Tumor_Seq_Allele2 == "G"))

    # The dedupe path should keep the annotated CTNNA3 block over the empty
    # candidate, so BOTH r2 rows resolve to CTNNA3 on the SnpEff side.
    expect_true(all(r2$snpeff_gene == "CTNNA3"))

    # Audit trail: attr(<result>, "snpeff_collisions_discarded") should carry
    # the 1 discarded row for the chr10 key.
    disc <- attr(res, "snpeff_collisions_discarded")
    expect_true(!is.null(disc))
    expect_true(data.table::is.data.table(disc))
    expect_gte(nrow(disc), 1L)
    expect_true(any(disc$Chromosome == "chr10" &
                    disc$Start_Position == 68173848L))
})
