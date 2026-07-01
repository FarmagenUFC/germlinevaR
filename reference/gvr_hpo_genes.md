# Resolve HPO phenotype terms to associated genes

Retrieves genes associated with one or more Human Phenotype Ontology
(HPO) terms using the HPO phenotype-to-gene association table. This is
mainly used by
[`read.gvr()`](https://farmagenufc.github.io/germlinevaR/reference/read.gvr.md)
and
[`gvr_filter()`](https://farmagenufc.github.io/germlinevaR/reference/gvr_filter.md)
through their `hpo` argument, but can also be called directly to inspect
the gene set before filtering.

The HPO table is downloaded once and cached under
`tools::R_user_dir("germlinevaR", "cache")`. Cached files older than
`max_age_days` days trigger an automatic re-download; setting
`refresh_cache = TRUE` forces refresh regardless of age. Use
`hpo_path = <local file>` for offline / air-gapped operation.

Input HPO identifiers are lenient-parsed: `"HP:0003002"`,
`"hp:0003002"`, `"3002"`, and `"0003002"` are all accepted and
normalised to canonical `"HP:0003002"` before lookup. When normalisation
happens a [`message()`](https://rdrr.io/r/base/message.html) is emitted
so typos remain visible.

No ontology-descendant expansion is performed in this release: only
genes associated with the exact HPO terms supplied are returned. To
include descendants, resolve them yourself first (e.g. via
`ontologyIndex::get_descendants()`).

## Usage

``` r
gvr_hpo_genes(
  hpo,
  hpo_path = NULL,
  hpo_url = "https://purl.obolibrary.org/obo/hp/hpoa/phenotype_to_genes.txt",
  cache_dir = NULL,
  refresh_cache = FALSE,
  max_age_days = 30,
  verbose = TRUE
)
```

## Arguments

- hpo:

  Character vector of HPO identifiers, e.g. `"HP:0003002"` or
  `c("HP:0003002", "HP:0001939")`. Lenient input accepted (see
  Description).

- hpo_path:

  Optional path to a local `phenotype_to_genes.txt` file. If supplied,
  no download is attempted.

- hpo_url:

  URL used when `hpo_path = NULL`. Default is the canonical HPO/OBO
  phenotype-to-gene association file.

- cache_dir:

  Directory used to cache the downloaded HPO file. `NULL` uses
  `tools::R_user_dir("germlinevaR", "cache")`.

- refresh_cache:

  Logical. If `TRUE`, force re-download of the HPO file ignoring
  `max_age_days`. Default `FALSE`.

- max_age_days:

  Numeric. Automatic refresh threshold in days for the cached HPO file.
  Default 30. Set to `Inf` to disable age-based refresh.

- verbose:

  Logical. If `TRUE`, print progress messages (download, cache age,
  per-code resolution counts, normalisation events). Default TRUE.

## Value

A sorted, upper-cased, deduplicated character vector of HGNC gene
symbols associated with the supplied HPO term(s). If nothing resolves,
returns `character(0)` and emits a warning listing the unresolved terms.

## See also

Other germlinevaR:
[`gvr_filter()`](https://farmagenufc.github.io/germlinevaR/reference/gvr_filter.md),
[`gvr_list_panels()`](https://farmagenufc.github.io/germlinevaR/reference/gvr_list_panels.md),
[`gvr_panel_genes()`](https://farmagenufc.github.io/germlinevaR/reference/gvr_panel_genes.md),
[`gvr_plot()`](https://farmagenufc.github.io/germlinevaR/reference/gvr_plot.md),
[`gvr_summary()`](https://farmagenufc.github.io/germlinevaR/reference/gvr_summary.md),
[`read.gvr()`](https://farmagenufc.github.io/germlinevaR/reference/read.gvr.md),
[`read.gvr.snpeff()`](https://farmagenufc.github.io/germlinevaR/reference/read.gvr.snpeff.md)

## Examples

``` r
# Runnable example: use the tiny HPO fixture shipped in inst/extdata so
# the example needs no network access.
hpo_fx <- system.file("extdata", "hpo_phenotype_to_genes_mini.tsv",
                      package = "germlinevaR")
genes <- gvr_hpo_genes("HP:0003002", hpo_path = hpo_fx)
#> gvr_hpo_genes: resolved 1 HPO term(s) -> 10 gene rows (per-code: HP:0003002=10)
head(genes)
#> [1] "ATM"   "BARD1" "BRCA1" "BRCA2" "BRIP1" "CHEK2"

# Lenient input forms all normalise to the canonical HP:0003002:
identical(genes,
          gvr_hpo_genes("3002", hpo_path = hpo_fx))
#> gvr_hpo_genes: normalised 1 HPO input(s) to canonical form: '3002' -> 'HP:0003002'
#> gvr_hpo_genes: resolved 1 HPO term(s) -> 10 gene rows (per-code: HP:0003002=10)
#> [1] TRUE

# Multiple terms are unioned:
gvr_hpo_genes(c("HP:0003002", "HP:0025022"), hpo_path = hpo_fx)
#> gvr_hpo_genes: resolved 2 HPO term(s) -> 13 gene rows (per-code: HP:0003002=10, HP:0025022=3)
#>  [1] "ATM"    "BARD1"  "BRCA1"  "BRCA2"  "BRIP1"  "CHEK2"  "COL5A1" "COL5A2"
#>  [9] "NF1"    "PALB2"  "PMS2"   "TNXB"   "TP53"  

# \donttest{
# Network use (downloads and caches the full HPO table on first call):
# genes <- gvr_hpo_genes("HP:0003002")

# Force refresh of the cached HPO file:
# genes <- gvr_hpo_genes("HP:0003002", refresh_cache = TRUE)

# Use with read.gvr() to restrict to HPO-implicated genes:
# gvr <- read.gvr("path/to/vcfs", hpo = "HP:0003002")

# Or on an already-loaded table:
# gvr_flt <- gvr_filter(gvr, hpo = "HP:0003002")
# }
```
