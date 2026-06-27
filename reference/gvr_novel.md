# Filter a read.gvr table down to candidate novel variants

Returns only the rows of a
[`read.gvr()`](https://farmagenufc.github.io/germlinevaR/reference/read.gvr.md)
/
[`gvr_filter()`](https://farmagenufc.github.io/germlinevaR/reference/gvr_filter.md)
table that have NO database evidence in any of the standard catalogues:
no rsID in `dbSNP_RS`, AND no allele frequency in `gnomADe_AF`, `AF`, or
`ABraOM_AF`. This is the canonical "show me only the novel variants"
subsetter and complements
[`gvr_filter()`](https://farmagenufc.github.io/germlinevaR/reference/gvr_filter.md)
(which uses presence of AF as a continuous filter).

## Usage

``` r
gvr_novel(gvr, verbose = TRUE)
```

## Arguments

- gvr:

  data.table / data.frame from
  [`read.gvr()`](https://farmagenufc.github.io/germlinevaR/reference/read.gvr.md)
  (or compatible, e.g. the output of
  [`gvr_filter()`](https://farmagenufc.github.io/germlinevaR/reference/gvr_filter.md)).
  Filtered in a copy; the input object is not modified.

- verbose:

  logical(1). If `TRUE` (default), print a step-by-step audit line for
  each of the 4 cascaded filters showing rows-in -\> rows-out and
  percentage retained, ending with a one-line summary.

## Value

A `data.table` with the same columns as `gvr` but only the rows that
pass all four "no database evidence" checks. If none of the 4 columns is
present in `gvr`, the output equals the input and a warning is issued.

## Details

A column value counts as "missing/empty" iff `is.na(x) | x == ""`. This
is the same convention used throughout the package (see
[`gvr_filter()`](https://farmagenufc.github.io/germlinevaR/reference/gvr_filter.md)
and
[`gvr_summary()`](https://farmagenufc.github.io/germlinevaR/reference/gvr_summary.md)):
[`read.gvr()`](https://farmagenufc.github.io/germlinevaR/reference/read.gvr.md)
writes `""` for absent values, but NA can creep in if upstream code
re-coerced columns.

A row is kept iff ALL FOUR of the following columns are missing/empty:

- `dbSNP_RS` - rsID from VEP `Existing_variation`

- `gnomADe_AF` - gnomAD exome allele frequency

- `AF` - gnomAD genome / VCF allele frequency

- `ABraOM_AF` - ABraOM SABE-609 Brazilian-cohort allele frequency

If any of the four columns is absent from the input (e.g. ABraOM
disabled at read time), the function treats that column as "all
missing" - i.e. it does not exclude any row on the basis of the missing
column. The audit log will say `(column not present, skipped)` for
transparency.

The cascade order is fixed (dbSNP_RS -\> gnomADe_AF -\> AF -\>
ABraOM_AF) so the verbose audit is reproducible across runs.

The input is never modified: filtering operates on an internal
`data.table` copy and the output has the same column set and order as
the input.

## See also

[`read.gvr()`](https://farmagenufc.github.io/germlinevaR/reference/read.gvr.md),
[`gvr_filter()`](https://farmagenufc.github.io/germlinevaR/reference/gvr_filter.md),
[`gvr_summary()`](https://farmagenufc.github.io/germlinevaR/reference/gvr_summary.md)

## Examples

``` r
## Load the shipped example table and find candidate novel variants
gvr <- readRDS(system.file("extdata", "example_gvr.rds",
    package = "germlinevaR"))
nov <- gvr_novel(gvr, verbose = FALSE)
dim(nov)
#> [1]   3 115

## Sanity-check: every kept row really is novel
stopifnot(all(is.na(nov$dbSNP_RS)   | nov$dbSNP_RS   == ""))
stopifnot(all(is.na(nov$gnomADe_AF) | nov$gnomADe_AF == ""))

## Combine with gvr_filter() to restrict to filtered novel variants
filt <- gvr_filter(gvr, ABraOM_AF = NULL, verbose = FALSE)
nov_filt <- gvr_novel(filt, verbose = FALSE)
dim(nov_filt)
#> [1]   3 115
```
