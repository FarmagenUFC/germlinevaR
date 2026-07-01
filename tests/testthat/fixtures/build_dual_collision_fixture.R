# Builds tests/testthat/fixtures/example.dual.collision.vcf.gz
#
# This script is checked into the repo alongside the fixture it produces
# so future maintainers can regenerate the fixture deterministically.
# Not called by the test suite at runtime.
#
# The fixture is a hermetic dual-annotated (VEP CSQ + SnpEff ANN) VCF that
# reproduces the two collision patterns exercised by test-dual-collision.R:
#
#   Record 1: chr1 multi-ALT indel that collides under legacy coord derivation
#   but SEPARATES to distinct MAF keys under bcftools-norm-style trimming.
#   This is the chr1:6095864 pattern from S4.freebayes_snpEff_VEP.ann.vcf.gz.
#
#   Record 2: chr10 multi-ALT indel in a homopolymer/repeat region where BOTH
#   ALTs still collide even after coord trimming (both events normalize to the
#   same 1-bp G insertion at the same start position). Exercises the dedupe
#   path. This is the chr10:68173832 pattern from S4.
#
# Both records carry CSQ and ANN blocks per ALT so the readers can populate
# the canonical VEP CSQ columns (80-field) and the SnpEff impact/gene columns.

# Resolve the fixture directory whether sourced (sys.frame(1L)$ofile is set)
# or run via `Rscript` (commandArgs() carries the file argument). Falls back
# to the current working directory so the script never errors on unusual
# invocations.
.this_dir <- local({
    of <- tryCatch(sys.frame(1L)$ofile, error = function(e) NULL)
    if (!is.null(of) && nzchar(of)) return(dirname(of))
    ca <- commandArgs(trailingOnly = FALSE)
    fa <- sub("^--file=", "", grep("^--file=", ca, value = TRUE))
    if (length(fa)) return(dirname(normalizePath(fa)))
    getwd()
})
out_path <- file.path(.this_dir, "example.dual.collision.vcf.gz")

# --- Header ------------------------------------------------------------------

hdr <- c(
    "##fileformat=VCFv4.2",
    "##contig=<ID=chr1,length=248956422>",
    "##contig=<ID=chr10,length=133797422>",
    '##INFO=<ID=DP,Number=1,Type=Integer,Description="Total depth">',
    '##INFO=<ID=AF,Number=A,Type=Float,Description="AF">',
    '##INFO=<ID=MQ,Number=1,Type=Float,Description="MQ">',
    '##INFO=<ID=QD,Number=1,Type=Float,Description="QD">',
    # Full 80-field CSQ (matches germlinevaR canonical spec)
    paste0(
        '##INFO=<ID=CSQ,Number=.,Type=String,Description="Consequence annotations from ',
        'Ensembl VEP. Format: Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|',
        'BIOTYPE|EXON|INTRON|HGVSc|HGVSp|cDNA_position|CDS_position|Protein_position|',
        'Amino_acids|Codons|Existing_variation|DISTANCE|STRAND|FLAGS|VARIANT_CLASS|',
        'SYMBOL_SOURCE|HGNC_ID|CANONICAL|MANE|MANE_SELECT|MANE_PLUS_CLINICAL|TSL|APPRIS|',
        'CCDS|ENSP|SWISSPROT|TREMBL|UNIPARC|UNIPROT_ISOFORM|GENE_PHENO|SIFT|PolyPhen|',
        'DOMAINS|miRNA|AF|AFR_AF|AMR_AF|EAS_AF|EUR_AF|SAS_AF|gnomADe_AF|gnomADe_AFR_AF|',
        'gnomADe_AMR_AF|gnomADe_ASJ_AF|gnomADe_EAS_AF|gnomADe_FIN_AF|gnomADe_MID_AF|',
        'gnomADe_NFE_AF|gnomADe_REMAINING_AF|gnomADe_SAS_AF|gnomADg_AF|gnomADg_AFR_AF|',
        'gnomADg_AMI_AF|gnomADg_AMR_AF|gnomADg_ASJ_AF|gnomADg_EAS_AF|gnomADg_FIN_AF|',
        'gnomADg_MID_AF|gnomADg_NFE_AF|gnomADg_REMAINING_AF|gnomADg_SAS_AF|MAX_AF|',
        'MAX_AF_POPS|FREQS|CLIN_SIG|SOMATIC|PHENO|PUBMED|MOTIF_NAME|MOTIF_POS|',
        'HIGH_INF_POS|MOTIF_SCORE_CHANGE|TRANSCRIPTION_FACTORS">'
    ),
    paste0(
        '##INFO=<ID=ANN,Number=.,Type=String,Description="Functional annotations: ',
        '\'Allele | Annotation | Annotation_Impact | Gene_Name | Gene_ID | Feature_Type | ',
        'Feature_ID | Transcript_BioType | Rank | HGVS.c | HGVS.p | cDNA.pos / cDNA.length | ',
        'CDS.pos / CDS.length | AA.pos / AA.length | Distance | ERRORS / WARNINGS / INFO\'">'
    ),
    '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
    '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Depth">',
    '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths">',
    '##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype quality">',
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tCOLLISION01"
)

# --- Record 1 : chr1:6095864 pattern (SEPARATES under normalize_alleles=TRUE) --
# REF=TCCCCCCCCCTGCCC (15 nt)
# ALT1=TCCCCCCCCCCTGCCC (16 nt) -> real 1-bp C insertion; post-trim = (6095874, "-", "T", INS)
# ALT2=TCCCCACCCCCTGCCC (16 nt) -> compound event; post-trim = (6095869, "-", "C", INS)
# Under legacy coord derivation BOTH collapse to (6095878, "-", "C", INS).
#
# CSQ blocks: distinct genes (RPS6KA1 for ALT1, PLEKHG5 for ALT2) so we can
# assert the join preserves separate annotations for each ALT.

rec1_csq <- paste(
    # ALT1 CSQ block (RPS6KA1 splice_region_variant on canonical transcript)
    "TCCCCCCCCCCTGCCC|splice_region_variant&intron_variant|LOW|RPS6KA1|ENSG00000117676|Transcript|ENST00000260060|protein_coding||1/22|c.100+8_100+9insC|||||||||1|||SNV|HGNC|10430|YES|MANE_Select|ENST00000260060.5||||CCDS287|ENSP00000260060|Q15418|||||||probably_damaging(0.9)|||0.001|||||||||||||||||||||||||0.001|EAS_AF||||||||||||",
    # ALT2 CSQ block (PLEKHG5 intron_variant)
    "TCCCCACCCCCTGCCC|intron_variant|MODIFIER|PLEKHG5|ENSG00000171680|Transcript|ENST00000377728|protein_coding||3/29|c.500+50_500+51insCTAG|||||||||1|||indel|HGNC|29105|YES|MANE_Select|ENST00000377728.5||||CCDS44088|ENSP00000366957|O94827||||||||||||||||||||||||||||||||||||||",
    sep = ","
)

rec1_ann <- paste(
    # ALT1 ANN block (RPS6KA1 splice_region_variant LOW)
    "TCCCCCCCCCCTGCCC|splice_region_variant&intron_variant|LOW|RPS6KA1|ENSG00000117676|transcript|ENST00000260060|protein_coding|1/22|c.100+8_100+9insC||||||",
    # ALT2 ANN block (PLEKHG5 intron_variant MODIFIER)
    "TCCCCACCCCCTGCCC|intron_variant|MODIFIER|PLEKHG5|ENSG00000171680|transcript|ENST00000377728|protein_coding|3/29|c.500+50_500+51insCTAG||||||",
    sep = ","
)

rec1_info <- sprintf(
    "DP=40;AF=0.4,0.4;MQ=60;QD=15;CSQ=%s;ANN=%s",
    rec1_csq, rec1_ann
)

rec1 <- paste(
    "chr1", "6095864", ".",
    "TCCCCCCCCCTGCCC",
    "TCCCCCCCCCCTGCCC,TCCCCACCCCCTGCCC",
    "500", "PASS", rec1_info,
    "GT:DP:AD:GQ", "0/1/2:40:10,15,15:60",
    sep = "\t"
)

# --- Record 2 : chr10:68173832 pattern (STILL COLLIDES post-trim) -------------
# REF=ATTTTTTTTTTTTTTGTTAGAG  (22 nt, T-homopolymer)
# ALT1=ATTTTTTTTTTTTTTTTGTAGAG (23 nt) -- 1-bp T insertion, trims to (68173848, "-", "G")
# ALT2=ATTTTTTTTTTTTTTTGGTAGAG (23 nt) -- 1-bp G insertion, ALSO trims to (68173848, "-", "G")
#
# Both events are genomically distinct but bcftools-norm-style trim in a repeat
# context collapses them to the same minimal (pos, ref, alt). The DEDUPE layer
# resolves the collision: keep the higher-impact ANN block; discard the empty/
# lower-impact one; log to attr(<result>, "snpeff_collisions_discarded").
#
# We give ALT1 an empty ANN candidate (0 blocks) and ALT2 a MODIFIER ANN block
# so the deduper unambiguously keeps ALT2.

rec2_csq <- paste(
    # ALT1 CSQ block: intergenic (MODIFIER)
    "ATTTTTTTTTTTTTTTTGTAGAG|intergenic_variant|MODIFIER||||||||||||||||1||||||||||||||||||||||||||||||||||||||||||||||||||||||",
    # ALT2 CSQ block: intron variant
    "ATTTTTTTTTTTTTTTGGTAGAG|intron_variant|MODIFIER|CTNNA3|ENSG00000183230|Transcript|ENST00000373544|protein_coding||5/17|c.900-100_900-99insG|||||||||1|||indel|HGNC|2510|YES|MANE_Select|ENST00000373544.5||||CCDS31198|ENSP00000362645|Q9UI47||||||||||||||||||||||||||||||||||||||",
    sep = ","
)

rec2_ann <- paste(
    # ALT1 ANN: empty candidate (no block matches ALT1)
    "ATTTTTTTTTTTTTTTGGTAGAG|intron_variant|MODIFIER|CTNNA3|ENSG00000183230|transcript|ENST00000373544|protein_coding|5/17|c.900-100_900-99insG||||||",
    sep = ","
)

rec2_info <- sprintf(
    "DP=35;AF=0.3,0.3;MQ=55;QD=12;CSQ=%s;ANN=%s",
    rec2_csq, rec2_ann
)

rec2 <- paste(
    "chr10", "68173832", ".",
    "ATTTTTTTTTTTTTTGTTAGAG",
    "ATTTTTTTTTTTTTTTTGTAGAG,ATTTTTTTTTTTTTTTGGTAGAG",
    "450", "PASS", rec2_info,
    "GT:DP:AD:GQ", "0/1/2:35:10,12,13:55",
    sep = "\t"
)

# --- Write -------------------------------------------------------------------

lines <- c(hdr, rec1, rec2)
con <- gzfile(out_path, "wb")
writeLines(lines, con)
close(con)
message("Wrote ", out_path, " (", file.info(out_path)$size, " bytes)")
