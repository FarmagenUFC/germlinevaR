# Clear the auto-fetched protein-domain cache used by `gvr_lollipop`

Deletes cached InterPro domain `.rds` files written by
`gvr_lollipop(..., domains = "auto")`. Useful when InterPro releases a
new version (roughly quarterly) and you want fresh annotations on the
next call. Also clears the matching in-memory session cache so the next
call actually re-fetches instead of returning the in-memory copy.

## Usage

``` r
gvr_domain_cache_clear(
  gene = NULL,
  organism = NULL,
  cache_dir = NULL,
  verbose = TRUE
)
```

## Arguments

- gene:

  Character(1) or `NULL`. Gene symbol to delete from the cache. `NULL`
  (default) deletes every cached domain file in the cache directory
  (i.e., every file matching `^domains_interpro_.*\\.rds$`).

- organism:

  Integer or character or `NULL`. NCBI taxonomy id used as part of the
  cache key. When `gene` is given but `organism` is `NULL`, every
  taxonomy variant of that gene is deleted (i.e., the pattern is
  `^domains_interpro_<GENE>_.*\\.rds$`). Default `NULL`.

- cache_dir:

  Character(1) or `NULL`. Override the cache directory. `NULL` (default)
  triggers the precedence chain above.

- verbose:

  Logical(1). Print one line per deleted file. Default `TRUE`.

## Value

Invisibly, a character vector of the file paths that were deleted
(length 0 if none).

## Details

The cache directory is resolved with the same precedence chain used by
[`gvr_lollipop()`](https://farmagenufc.github.io/germlinevaR/reference/gvr_lollipop.md),
but here the chain only **finds** an existing directory; it never
creates one. The precedence (first hit wins) is:

1.  `cache_dir` argument (explicit override).

2.  Environment variable `GVR_CACHE_DIR`.

3.  R option `getOption("germlinevaR.cache_dir")`.

4.  `tools::R_user_dir("germlinevaR", "cache")`.

5.  `file.path(tempdir(), "germlinevaR_cache")`.

If the resolved directory does not exist, the helper returns invisibly
with a `verbose` message and no files deleted (this is not an error; a
missing cache directory simply means there was nothing to clear).

## See also

[`gvr_lollipop()`](https://farmagenufc.github.io/germlinevaR/reference/gvr_lollipop.md)

## Examples

``` r
## Clear the cache for a specific gene (safe to run; no-op if cache is empty)
gvr_domain_cache_clear(gene = "TP53")
#> gvr_domain_cache_clear: no cache directory found (nothing to clear).

# \donttest{
  ## Clear everything (all genes, all organisms)
  gvr_domain_cache_clear()
#> gvr_domain_cache_clear: no cache directory found (nothing to clear).

  ## Clear only TP53 across all organisms
  gvr_domain_cache_clear(gene = "TP53")
#> gvr_domain_cache_clear: no cache directory found (nothing to clear).

  ## Clear only human TP53
  gvr_domain_cache_clear(gene = "TP53", organism = 9606)
#> gvr_domain_cache_clear: no cache directory found (nothing to clear).
# }
```
