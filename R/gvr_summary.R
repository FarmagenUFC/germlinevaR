#' Multi-section summary of a germline MAF (read.gvr / gvr_filter output)
#'
#' @description
#' Produces a multi-section overview of a MAF-style table - the output of
#' [read.gvr()], or of [gvr_filter()] - covering variant burden, affected genes,
#' functional classes, clinical significance and predicted impact. Every section is
#' returned as a tidy `data.table` with one column per sample plus a `Total` column.
#' Optionally writes a multi-sheet Excel workbook and/or a multi-page PDF report
#' (both into a `gvr_summary/` subfolder of `out_dir`). For a cohort oncoplot, see
#' [gvr_plot()].
#'
#' @details
#' Sections returned (as a named list of `data.table`s):
#' \itemize{
#'   \item `overview` - cohort-level counts: total variants, variants per sample,
#'     distinct genes affected, and variants with no gene symbol.
#'   \item `top_genes` - the `top_n_genes` genes with the most variants
#'     (per-sample + `Total`); unknown/blank genes excluded.
#'   \item `variant_classification` - counts per `Variant_Classification` (maftools
#'     functional class), per-sample + `Total`, sorted by `Total` descending.
#'   \item `variant_type` - counts per `Variant_Type` (SNP/DEL/INS/ONP/DNP/TNP).
#'   \item `clin_sig` - counts per `CLIN_SIG` token (ClinVar categories). `CLIN_SIG`
#'     is split on `&` and `/`, so a variant annotated `"pathogenic&benign"`
#'     increments BOTH categories; category counts therefore sum to \eqn{\ge} the
#'     number of variants. A `missing/unclassified` row counts NA/"" `CLIN_SIG`.
#'   \item `top_genes_per_sample` - a named list of `data.table`s, one per sample,
#'     each containing the `top_n_genes` genes with the most variants in that
#'     sample (columns: `Hugo_Symbol`, `<sample_name>`). Unknown/blank genes
#'     excluded.
#'   \item `top_variants` - the `top_n_variants` most frequent variants by
#'     `dbSNP_RS` (rsID), with per-sample counts + `Total`. Rows with
#'     blank/missing `dbSNP_RS` are excluded. If the `dbSNP_RS` column is absent,
#'     this section is skipped with a warning.
#'   \item `impact` - counts per VEP `IMPACT` (HIGH/MODERATE/LOW/MODIFIER), in
#'     severity order rather than count order.
#' }
#'
#' The section tables are the core return value. By default the function also writes
#' an Excel workbook (`save_excel = TRUE`) and a PDF report (`save_pdf = TRUE`) into
#' `out_dir/gvr_summary/`; set either to `FALSE` to skip it. The cohort oncoplot lives
#' in [gvr_plot()].
#'
#' Data conventions:
#' \itemize{
#'   \item "Missing" means `NA` OR empty string `""`.
#'   \item Unknown/blank gene symbols are `Hugo_Symbol` in `c(".", "", "Unknown")`;
#'     these are excluded from the distinct-gene tally and from `top_genes`, but their
#'     variants are still counted in the totals (and reported as "variants with no
#'     gene symbol").
#'   \item Works on ANY MAF-shaped table; it makes no assumption about prior
#'     filtering. It is commonly run on [gvr_filter()] output to summarise the
#'     retained hits.
#' }
#'
#' @param maf A `data.table`/`data.frame` MAF from [read.gvr()] or [gvr_filter()].
#'   Required columns: `Hugo_Symbol`, `Variant_Classification`, `Variant_Type`,
#'   `IMPACT`, `CLIN_SIG`.
#' @param sample_col Name of the per-sample column. Default `"Tumor_Sample_Barcode"`.
#'   If absent, all rows are pooled into a single sample `"All"` (with a warning).
#' @param top_n_genes Integer; number of genes to report in `top_genes` (by total
#'   variant count). Default `20`.
#' @param top_n_variants Integer; number of variants to report in `top_variants`
#'   (by `dbSNP_RS` frequency). Default `20`. Ignored if the `dbSNP_RS` column is
#'   absent from the input.
#' @param save_excel Logical; if `TRUE` (default), write a multi-sheet `.xlsx`.
#'   The workbook is written into the `gvr_summary/` subfolder of `out_dir`
#'   (see `out_dir`). Pass `FALSE` for a compute-only run.
#' @param save_pdf Logical; if `TRUE` (default), write a multi-page PDF dashboard
#'   report into the `gvr_summary/` subfolder of `out_dir`. Page 1 is a hero page (a
#'   row of KPI cards above two grouped/faceted bar charts - top genes and variant
#'   classification); the following pages hold the section tables (packed
#'   two-per-row where they fit, else full-width) plus the functional-impact chart.
#'   The layout adapts to cohort size (faceting and column pagination for many
#'   samples; see the examples). Requires \pkg{gridExtra}, \pkg{ggplot2} and
#'   \pkg{scales}; if unavailable, the PDF is skipped with a warning and the sections
#'   are still returned. Pass `FALSE` for a compute-only run.
#' @param save_html Logical; if `TRUE` (default), write an interactive HTML dashboard
#'   (`<file_prefix>_report.html`) into the `gvr_summary/` subfolder of `out_dir`. It
#'   mirrors the PDF dashboard - a row of KPI cards, bar charts (top genes,
#'   variant classification, functional impact, top variants) as interactive
#'   \pkg{plotly} charts (grouped for \eqn{\le 6} samples, faceted small-multiples
#'   for \eqn{> 6}), and all section tables as sortable, searchable \pkg{DT} tables.
#'   The Clinical significance table is interactive: clicking a CLIN_SIG token
#'   (e.g. "pathogenic") expands a detail panel showing the individual variants
#'   with that annotation. By default a
#'   single self-contained file is produced (assets inlined via pandoc); if pandoc is
#'   unavailable the report is written as `<file_prefix>_report.html` plus a sibling
#'   `<file_prefix>_report_files/` asset folder (a `verbose` message notes this).
#'   Requires \pkg{plotly}, \pkg{DT}, \pkg{htmlwidgets} and \pkg{htmltools}; if any are
#'   unavailable the HTML is skipped with a warning and the sections are still
#'   returned. Pass `FALSE` for a compute-only run.
#' @param out_dir Parent output directory. All written outputs (Excel, PDF and/or
#'   HTML) are placed in a `gvr_summary/` subfolder created inside `out_dir`. The
#'   subfolder is created only when `save_excel`, `save_pdf` or `save_html` is `TRUE`.
#'   Default `"."` (current working directory), i.e. outputs go to `./gvr_summary/`.
#' @param file_prefix Base filename for written outputs. Default `"gvr_summary"`, giving
#'   `<file_prefix>.xlsx`, `<file_prefix>_report.pdf` and `<file_prefix>_report.html`
#'   (no timestamp). Filenames are fixed, so re-running into the same `out_dir`
#'   overwrites the previous files (a message is printed when `verbose = TRUE`).
#' @param verbose Logical; if `TRUE` (default) print a compact console digest and the
#'   path(s) of any file(s) written.
#'
#' @return Invisibly, a named list of `data.table`s: `overview`, `top_genes`,
#'   `top_genes_per_sample`, `variant_classification`, `variant_type`, `clin_sig`,
#'   `top_variants`, `impact`. The `top_genes_per_sample` element is itself a named
#'   list (one data.table per sample). The `top_variants` section is absent if
#'   `dbSNP_RS` is not in the input. The return value is identical regardless of
#'   whether the Excel/PDF/HTML files are written.
#'
#' @section Dependencies:
#' Core summary uses \pkg{data.table}. The optional Excel export uses \pkg{openxlsx};
#' the optional PDF dashboard uses \pkg{gridExtra} + \pkg{ggplot2} + \pkg{scales},
#' rendered via the \code{grDevices::cairo_pdf} device (full Unicode, so en-dashes,
#' multiplication signs and similar punctuation render correctly). The optional
#' interactive HTML dashboard uses \pkg{plotly} (charts), \pkg{DT} (tables) and
#' \pkg{htmlwidgets} + \pkg{htmltools} (assembly); a single self-contained file is
#' produced when \pkg{pandoc} is available (used only for the optional asset-inlining
#' step), otherwise a `.html` + `_files/` asset folder is written. Each optional
#' output degrades gracefully: if its package(s) are unavailable, that output is
#' skipped with a warning and the section tables are still returned.
#'
#' @seealso [read.gvr()] to build the MAF, [gvr_filter()] to filter it before
#'   summarising, [gvr_plot()] for a cohort oncoplot.
#' @family germlinevaR
#' @author germlinevaR authors
#'
#' @examples
#' \dontrun{
#' maf <- read.gvr("/path/to/vcf_folder")
#'
#' ## ---- Default run -------------------------------------------------------
#' ## Returns the six section tables AND writes three files into ./gvr_summary/:
#' ##   * gvr_summary.xlsx         - one sheet per section
#' ##   * gvr_summary_report.pdf   - the print dashboard report (see below)
#' ##   * gvr_summary_report.html  - the interactive dashboard (see below)
#' s <- gvr_summary(maf)
#' s$variant_classification          # inspect a section
#' s$impact                          # HIGH -> MODIFIER, in severity order
#'
#' ## ---- The PDF dashboard -------------------------------------------------
#' ## gvr_summary_report.pdf is a single A4-portrait dashboard:
#' ##   Page 1 (hero):  a row of four KPI cards - total variants, number of
#' ##                   samples, distinct known genes, and HIGH-impact variants -
#' ##                   above two grouped bar charts: "Top genes (top 10)" and
#' ##                   "Variant classification (top 10)".
#' ##   Page 2+ (reference): the six section tables, packed two-per-row when both
#' ##                   fit the page width and stacked full-width otherwise, plus
#' ##                   the Functional-impact (VEP IMPACT) bar chart.
#' ## Distinct-gene / HIGH-impact KPIs and chart values mirror the returned tables.
#' s <- gvr_summary(maf)                       # 2-sample cohort -> ~3 PDF pages
#'
#' ## ---- The interactive HTML dashboard ------------------------------------
#' ## gvr_summary_report.html is the same dashboard, interactive: the four KPI
#' ## cards, the three bar charts as plotly (hover for exact counts, click the
#' ## legend to toggle samples, drag to zoom), and all six section tables as DT
#' ## widgets (per-column sort, a global search box, and pagination). By default
#' ## it is a single self-contained file you can email or open offline; if pandoc
#' ## is unavailable it is written alongside a gvr_summary_report_files/ folder.
#' s <- gvr_summary(maf)                       # writes gvr_summary_report.html
#' gvr_summary(maf, save_html = FALSE)          # opt out of the HTML report
#'
#' ## ---- Multi-sample / cohort behaviour -----------------------------------
#' ## The same call scales to any number of samples; both reports adapt:
#' ##   * <= 6 samples : charts are GROUPED bars (one bar per sample).
#' ##   * >  6 samples : charts switch to small-multiple FACETS (one panel per
#' ##                    sample) so labels stay legible - in the PDF and, via
#' ##                    plotly::subplot, in the HTML.
#' ##   * Wide tables  : in the PDF, when per-sample columns no longer fit the
#' ##                    page width, reference tables are COLUMN-PAGINATED - the
#' ##                    category and Total columns repeat on each part, titled
#' ##                    e.g. "Top genes (samples 1-7 of 8)". The HTML DT tables
#' ##                    instead keep every sample column and scroll horizontally.
#' ## No argument controls this; layout is chosen automatically from sample count
#' ## and measured table widths. A 20-sample cohort typically spans ~8 PDF pages.
#' s <- gvr_summary(maf, sample_col = "Tumor_Sample_Barcode")
#'
#' ## ---- Other modes -------------------------------------------------------
#' ## Compute only (no files written); still prints a console digest:
#' s <- gvr_summary(maf, save_excel = FALSE, save_pdf = FALSE, save_html = FALSE)
#'
#' ## Summarise filtered hits, writing into results/summary/gvr_summary/:
#' gvr_summary(gvr_filter(maf), out_dir = "results/summary")
#'
#' ## Report more genes and silence the console digest:
#' gvr_summary(maf, top_n_genes = 30, verbose = FALSE)
#' }
#'
#' @importFrom data.table as.data.table data.table setnames setcolorder setorder
#'   uniqueN copy melt :=
#' @importFrom utils head
#' @importFrom openxlsx createWorkbook
#' @importFrom grDevices pdf dev.off
#' @export


gvr_summary <- function(maf,
                        sample_col     = "Tumor_Sample_Barcode",
                        top_n_genes    = 20,
                        top_n_variants = 20,
                        save_excel     = TRUE,
                        save_pdf       = TRUE,
                        save_html      = TRUE,
                        out_dir        = ".",
                        file_prefix    = "gvr_summary",
                        verbose        = TRUE) {

  if (!requireNamespace("data.table", quietly = TRUE)) {
    stop("gvr_summary requires the 'data.table' package.")
  }
  dt <- data.table::as.data.table(maf)
  n_total <- nrow(dt)

  .is_missing <- function(v) is.na(v) | v == ""
  UNKNOWN_GENE <- c(".", "", "Unknown")
  # --- Resolve sample column ---------------------------------------------------
  if (!sample_col %in% names(dt)) {
    warning(sprintf("gvr_summary: sample column '%s' not found; pooling all rows into 'All'.",
                    sample_col))
    dt[, .__sample__ := "All"]
  } else {
    dt[, .__sample__ := as.character(get(sample_col))]
    dt[.is_missing(.__sample__), .__sample__ := "NA_sample"]
  }
  samples <- sort(unique(dt$.__sample__))

  # --- Column-existence guard for the analytic columns -------------------------
  req <- c("Hugo_Symbol", "Variant_Classification", "Variant_Type", "IMPACT", "CLIN_SIG")
  miss_req <- req[!req %in% names(dt)]
  if (length(miss_req) > 0) {
    stop(sprintf("gvr_summary: required column(s) not found: %s",
                 paste(miss_req, collapse = ", ")))
  }

  # --- Helper: counts of `valuevec` x sample as a wide data.table --------------
  #     rows = categories, cols = samples + Total. `order_levels` optionally fixes
  #     row order; otherwise sorted by Total desc.
  .count_by_sample <- function(valuevec, samplevec, category_name,
                               order_levels = NULL, drop_values = NULL) {
    keep <- rep(TRUE, length(valuevec))
    if (!is.null(drop_values)) keep <- !(valuevec %in% drop_values)
    v <- valuevec[keep]; s <- samplevec[keep]
    if (length(v) == 0L) {
      out <- data.table::data.table(X = character(0))
      data.table::setnames(out, "X", category_name)
      for (sm in samples) out[, (sm) := integer(0)]
      out[, Total := integer(0)]
      return(out[])
    }
    tab <- table(factor(v), factor(s, levels = samples))
    m <- data.table::as.data.table(unclass(tab), keep.rownames = TRUE)
    data.table::setnames(m, "rn", category_name)
    # ensure all sample columns exist (table drops empty factor combos only if level absent)
    for (sm in samples) if (!sm %in% names(m)) m[, (sm) := 0L]
    data.table::setcolorder(m, c(category_name, samples))
    m[, Total := rowSums(as.matrix(.SD)), .SDcols = samples]
    if (!is.null(order_levels)) {
      m <- m[match(order_levels, m[[category_name]])]
      m <- m[!is.na(m[[category_name]])]
    } else {
      data.table::setorder(m, -Total)
    }
    m[]
  }
  # ============================================================================
  # SECTION 1: overview (variants & genes)
  # ============================================================================
  is_known_gene <- !(dt$Hugo_Symbol %in% UNKNOWN_GENE)
  ov_rows <- list()
  # total variants
  per_sample_var <- vapply(samples, function(sm) sum(dt$.__sample__ == sm), integer(1))
  ov_rows[["Total variants"]] <- c(per_sample_var, Total = n_total)
  # distinct genes (known only)
  per_sample_genes <- vapply(samples, function(sm)
    data.table::uniqueN(dt$Hugo_Symbol[dt$.__sample__ == sm & is_known_gene]), integer(1))
  ov_rows[["Distinct genes (known)"]] <- c(per_sample_genes,
                                           Total = data.table::uniqueN(dt$Hugo_Symbol[is_known_gene]))
  # variants with no gene symbol
  per_sample_nogene <- vapply(samples, function(sm)
    sum(dt$.__sample__ == sm & !is_known_gene), integer(1))
  ov_rows[["Variants with no gene symbol"]] <- c(per_sample_nogene, Total = sum(!is_known_gene))

  overview <- data.table::data.table(Metric = names(ov_rows))
  for (sm in c(samples, "Total")) {
    overview[, (sm) := vapply(ov_rows, function(r) r[[sm]], numeric(1))]
  }

  # ============================================================================
  # SECTION 1b: top genes (known only), by Total variant count
  # ============================================================================
  gene_tab <- .count_by_sample(dt$Hugo_Symbol, dt$.__sample__, "Hugo_Symbol",
                               drop_values = UNKNOWN_GENE)
  top_genes <- utils::head(gene_tab, top_n_genes)

  # ============================================================================
  # SECTION 1c: per-sample top genes (known only), one table per sample
  # ============================================================================
  top_genes_per_sample <- list()
  for (sm in samples) {
    dt_sm <- dt[.__sample__ == sm & !(Hugo_Symbol %in% UNKNOWN_GENE)]
    if (nrow(dt_sm) > 0L) {
      sm_gene_tab <- dt_sm[, .(N = .N), by = Hugo_Symbol]
      data.table::setorder(sm_gene_tab, -N)
      sm_gene_tab <- utils::head(sm_gene_tab, top_n_genes)
      data.table::setnames(sm_gene_tab, "N", sm)
      top_genes_per_sample[[sm]] <- sm_gene_tab
    }
  }

  # ============================================================================
  # SECTION 2: functional classes
  # ============================================================================
  variant_classification <- .count_by_sample(dt$Variant_Classification, dt$.__sample__,
                                              "Variant_Classification")
  variant_type <- .count_by_sample(dt$Variant_Type, dt$.__sample__, "Variant_Type")

  # ============================================================================
  # SECTION 3: clinical categories (token-split on & and /)
  # ============================================================================
  cs   <- dt$CLIN_SIG
  miss <- .is_missing(cs)
  # expand non-missing CLIN_SIG into (token, sample) pairs
  idx_nm <- which(!miss)
  if (length(idx_nm) > 0) {
    # dedup tokens per row so multi-token CLIN_SIG values
    # (e.g. "pathogenic&pathogenic/likely_pathogenic") count each token once per variant
    toks_list <- lapply(strsplit(cs[idx_nm], "[&/]"), function(x) unique(trimws(x)))
    ntok <- lengths(toks_list)
    tok_vec <- unlist(toks_list, use.names = FALSE)
    samp_vec <- rep(dt$.__sample__[idx_nm], ntok)
    # drop any empty tokens produced by stray delimiters
    keep_tok <- tok_vec != ""
    clin_counts <- .count_by_sample(tok_vec[keep_tok], samp_vec[keep_tok], "CLIN_SIG")
  } else {
    clin_counts <- .count_by_sample(character(0), character(0), "CLIN_SIG")
  }
  # append missing/unclassified row
  per_sample_miss <- vapply(samples, function(sm) sum(miss & dt$.__sample__ == sm), integer(1))
  miss_row <- data.table::as.data.table(c(list(CLIN_SIG = "missing/unclassified"),
                                          as.list(per_sample_miss),
                                          list(Total = sum(miss))))
  clin_sig <- rbind(clin_counts, miss_row, fill = TRUE)

  # ============================================================================
  # SECTION 3b: top variants by dbSNP_RS (rsID frequency)
  # ============================================================================
  top_variants <- NULL
  has_dbsnp <- "dbSNP_RS" %in% names(dt)
  if (has_dbsnp) {
    rs <- dt$dbSNP_RS
    rs_miss <- .is_missing(rs) | rs == "novel" | tolower(rs) == "unknown"
    dt_rs <- dt[!rs_miss]
    if (nrow(dt_rs) > 0L) {
      # Group by dbSNP_RS; pick representative metadata (first non-blank per field)
      dt_rs[, .__rs__ := dbSNP_RS]
      # Aggregate: count per rsID x sample, plus representative Hugo/Chr/Pos/VC
      rs_tab <- dt_rs[, .(Hugo_Symbol = {
        v <- Hugo_Symbol[!(Hugo_Symbol %in% UNKNOWN_GENE) & !.is_missing(Hugo_Symbol)]
        if (length(v)) v[1L] else ""
      }, Chromosome = {
        v <- Chromosome[!.is_missing(Chromosome)]; if (length(v)) v[1L] else ""
      }, Start_Position = {
        v <- Start_Position[!.is_missing(Start_Position)]; if (length(v)) v[1L] else ""
      }, Variant_Classification = {
        v <- Variant_Classification[!.is_missing(Variant_Classification)]; if (length(v)) v[1L] else ""
      }, Total = .N), by = .__rs__]
      # Per-sample counts
      for (sm in samples) {
        rs_tab[, (sm) := dt_rs[.__sample__ == sm, .N, by = .__rs__][match(rs_tab$.__rs__, .__rs__), N]]
        rs_tab[is.na(get(sm)), (sm) := 0L]
      }
      data.table::setcolorder(rs_tab, c(".__rs__", "Hugo_Symbol", "Chromosome",
                                         "Start_Position", "Variant_Classification",
                                         samples, "Total"))
      data.table::setorder(rs_tab, -Total)
      data.table::setnames(rs_tab, ".__rs__", "dbSNP_RS")
      top_variants <- utils::head(rs_tab, top_n_variants)
    }
  } else {
    if (isTRUE(verbose))
      message("  Note: 'dbSNP_RS' column not found; top_variants section skipped.")
  }

  # ============================================================================
  # SECTION 4: impact severity (fixed severity order)
  # ============================================================================
  impact_order <- c("HIGH", "MODERATE", "LOW", "MODIFIER")
  present_impact <- c(intersect(impact_order, unique(dt$IMPACT)),
                      setdiff(unique(dt$IMPACT[!.is_missing(dt$IMPACT)]), impact_order))
  impact <- .count_by_sample(dt$IMPACT, dt$.__sample__, "IMPACT",
                             order_levels = present_impact)

  sections <- list(overview = overview,
                   top_genes = top_genes,
                   top_genes_per_sample = top_genes_per_sample,
                   variant_classification = variant_classification,
                   variant_type = variant_type,
                   clin_sig = clin_sig,
                   impact = impact)
  if (!is.null(top_variants)) sections$top_variants <- top_variants

  # ============================================================================
  # Verbose console digest
  # ============================================================================
  if (isTRUE(verbose)) {
    message(sprintf("gvr_summary: %d variants across %d sample(s): %s",
                    n_total, length(samples), paste(samples, collapse = ", ")))
    message(sprintf("  Distinct genes (known): %d   |  variants with no gene symbol: %d",
                    data.table::uniqueN(dt$Hugo_Symbol[is_known_gene]), sum(!is_known_gene)))
    topfc <- variant_classification[1:min(3, nrow(variant_classification))]
    message(sprintf("  Top functional classes: %s",
                    paste(sprintf("%s=%d", topfc$Variant_Classification, topfc$Total), collapse = ", ")))
    message(sprintf("  IMPACT: %s",
                    paste(sprintf("%s=%d", impact$IMPACT, impact$Total), collapse = ", ")))
    nrare <- clin_sig[CLIN_SIG != "missing/unclassified"]
    if (nrow(nrare) > 0) {
      message(sprintf("  CLIN_SIG (top tokens): %s",
                      paste(sprintf("%s=%d", utils::head(nrare$CLIN_SIG, 4),
                                    utils::head(nrare$Total, 4)), collapse = ", ")))
    }
    if (!is.null(top_variants) && nrow(top_variants) > 0) {
      message(sprintf("  Top variant (by rsID): %s (%d occurrences)",
                      top_variants$dbSNP_RS[1L], top_variants$Total[1L]))
    }
  }

  # ============================================================================
  # Output folder: ALL written summary outputs go into <out_dir>/gvr_summary/.
  # The subfolder is created only when something is actually written.
  # ============================================================================
  out_subdir <- file.path(out_dir, "gvr_summary")
  if (isTRUE(save_excel) || isTRUE(save_pdf) || isTRUE(save_html)) {
    if (!dir.exists(out_subdir))
      dir.create(out_subdir, recursive = TRUE, showWarnings = FALSE)
  }

  # --------------------------------------------------------------------------
  # Shared report metadata (built once; consumed by BOTH the PDF and HTML
  # renderers so the two reports show identical cohort figures).
  # --------------------------------------------------------------------------
  meta <- list(
    out_dir   = normalizePath(out_dir, mustWork = FALSE),
    n_samples = length(samples), samples = samples,
    n_total   = n_total,
    n_genes   = data.table::uniqueN(dt$Hugo_Symbol[is_known_gene]),
    n_nogene  = sum(!is_known_gene),
    generated = format(Sys.time(), "%Y-%m-%d %H:%M:%S"))

  # --------------------------------------------------------------------------
  # Shared helpers + constants used by BOTH the PDF and HTML renderers.
  # Hoisted to the outer scope so the two report paths cannot drift apart.
  # --------------------------------------------------------------------------
  fmt <- function(x) format(x, big.mark = ",", trim = TRUE)

  # Colour-blind-safe shared palette (Okabe-Ito ordering for the first 9
  # entries: PDF dashboard uses exactly these). The HTML branch extends this
  # with Paul-Tol "muted" + extras for cohorts up to 24 distinct legend entries
  # (defined locally inside the HTML renderer).
  PHYLO_PAL_BASE <- c("#0279EE", "#E69F00", "#009E73", "#CC79A7", "#56B4E9",
                      "#75A025", "#D55E00", "#332288", "#AA4499")

  # Common-prefix stripper for facet/group sample labels: identifies a leading
  # substring shared by all sample names (length >= 3, e.g. "ACRO") and strips
  # it from chart labels so legends stay compact. Returns the identity for
  # samples with no shared >=3-char prefix or for cohorts of size 1.
  .lab_fun <- local({
    sn <- as.character(samples); pre <- ""
    if (length(sn) > 1L) {
      mn <- min(nchar(sn)); i <- 0L
      while (i < mn && length(unique(substr(sn, 1L, i + 1L))) == 1L) i <- i + 1L
      if (i >= 3L) pre <- substr(sn[1], 1L, i)
    }
    function(x) if (nzchar(pre)) sub(paste0("^", pre), "", x) else x
  })

  # Build the report table spec list (Overview / Top genes / per-sample top
  # genes / VC / VT / Clinical / Impact / Top variants if present). Consumed
  # identically by the PDF and HTML renderers, with two cosmetic divergences:
  #   * PDF appends "(table)" to the impact title because the PDF puts a
  #     separate IMPACT chart on the same reference page (so the table needs
  #     a label that disambiguates from the chart).
  #   * HTML tags VC / clin_sig / impact with a `special` key so the renderer
  #     knows to render those entries as drill-down (summary + detail) pairs.
  .build_tbl_specs <- function(sections, html = FALSE) {
    specs <- list(
      list(s = "overview",  t = "Overview"),
      list(s = "top_genes", t = "Top genes"))
    for (sm in names(sections$top_genes_per_sample)) {
      specs <- c(specs, list(
        list(s = paste0("top_genes_per_sample.", sm),
             t = sprintf("Top genes \u2013 %s", sm),
             per_sample_key = sm)))
    }
    impact_t <- if (isTRUE(html)) "Functional impact" else "Functional impact (table)"
    rest <- list(
      list(s = "variant_classification", t = "Variant classification"),
      list(s = "variant_type",           t = "Variant type"),
      list(s = "clin_sig",               t = "Clinical significance"),
      list(s = "impact",                 t = impact_t))
    if (isTRUE(html)) {
      rest[[1]]$special <- "vc"
      rest[[3]]$special <- "clin_sig"
      rest[[4]]$special <- "impact"
      # vN+7: Top genes becomes a 4th clickable section (drill-down = variants
      # in the clicked gene, capped at 1,000 by composite IMPACT > nsamp >
      # gnomADe_AF ranking). The specs$top_genes index is fixed: items 1=overview,
      # 2=top_genes, then per-sample tables.
      specs[[2]]$special <- "top_genes"
    }
    specs <- c(specs, rest)
    if (!is.null(sections$top_variants))
      specs <- c(specs, list(list(s = "top_variants", t = "Top variants (by rsID)")))
    specs
  }

  # ============================================================================
  # Optional Excel export  ->  <out_dir>/gvr_summary/<file_prefix>.xlsx
  # Fixed filename (no timestamp): re-running overwrites the previous workbook.
  # ============================================================================
  if (isTRUE(save_excel)) {
    if (!requireNamespace("openxlsx", quietly = TRUE)) {
      warning("gvr_summary: 'openxlsx' not installed; skipping Excel export.")
    } else {
      xlsx_name  <- sprintf("%s.xlsx", file_prefix)
      final_xlsx <- file.path(out_subdir, xlsx_name)
      if (file.exists(final_xlsx) && isTRUE(verbose))
        message(sprintf("  Overwriting existing Excel: %s", final_xlsx))

      sheet_map <- c(overview = "Overview", top_genes = "Top_genes",
                     variant_classification = "Variant_classification",
                     variant_type = "Variant_type", clin_sig = "Clinical",
                     top_variants = "Top_variants", impact = "Impact")
      wb <- openxlsx::createWorkbook()
      hs <- openxlsx::createStyle(textDecoration = "bold", halign = "center")
      for (nm in names(sections)) {
        if (nm == "top_genes_per_sample") {
          # Per-sample top genes: one sheet per sample
          for (sm in names(sections$top_genes_per_sample)) {
            sh_nm <- sprintf("TopGenes_%s", sm)
            # Excel sheet names max 31 chars; truncate sample name if needed
            if (nchar(sh_nm) > 31L) sh_nm <- paste0(substr(sh_nm, 1L, 28L), "...") 
            openxlsx::addWorksheet(wb, sh_nm)
            openxlsx::writeData(wb, sh_nm, sections$top_genes_per_sample[[sm]], headerStyle = hs)
            openxlsx::freezePane(wb, sh_nm, firstRow = TRUE)
            openxlsx::setColWidths(wb, sh_nm, cols = seq_len(ncol(sections$top_genes_per_sample[[sm]])),
                                   widths = "auto")
          }
          next
        }
        sh <- if (nm %in% names(sheet_map)) sheet_map[[nm]] else nm
        openxlsx::addWorksheet(wb, sh)
        openxlsx::writeData(wb, sh, sections[[nm]], headerStyle = hs)
        openxlsx::freezePane(wb, sh, firstRow = TRUE)
        openxlsx::setColWidths(wb, sh, cols = seq_len(ncol(sections[[nm]])), widths = "auto")
      }
      # ==============================================================
      # vN+7: per-clickable XLSX sheets -- one sheet per token of each
      # of the 4 clickable categories (Top genes / CLIN_SIG / IMPACT /
      # Variant_Classification). Each sheet holds the FULL uncapped
      # matching variants for that token, 7-column fixed projection,
      # sorted by the same composite ranking as the HTML drill-downs.
      # Excel sheet names are capped at 31 chars; tokens are truncated
      # with collision-suffix on conflict.
      # ==============================================================
      .gvr_impact_rank_xl <- function(v) {
        ord <- c("HIGH" = 1L, "MODERATE" = 2L, "LOW" = 3L, "MODIFIER" = 4L)
        r <- ord[as.character(v)]; r[is.na(r)] <- 5L; as.integer(r)
      }
      .gvr_rank_xl <- function(d) {
        ir <- .gvr_impact_rank_xl(d$IMPACT)
        nsamp <- if (".__nsamp__" %in% names(d)) d$.__nsamp__ else rep(0L, nrow(d))
        gaf <- if ("gnomADe_AF" %in% names(d)) suppressWarnings(as.numeric(d$gnomADe_AF))
                else rep(NA_real_, nrow(d))
        chr <- if ("Chromosome" %in% names(d)) as.character(d$Chromosome) else rep("", nrow(d))
        pos <- if ("Start_Position" %in% names(d)) suppressWarnings(as.integer(d$Start_Position))
                else rep(0L, nrow(d))
        order(ir, -nsamp, gaf, chr, pos, na.last = TRUE)
      }
      # Precompute .__nsamp__ on dt if not already present.
      if (!(".__nsamp__" %in% names(dt))) {
        key_cols_xl <- intersect(c("Chromosome", "Start_Position",
                                   "Reference_Allele", "Tumor_Seq_Allele2"),
                                 names(dt))
        if (length(key_cols_xl) < 2L) {
          dt[, .__nsamp__ := 1L]
        } else {
          dt[, .__nsamp__ := data.table::uniqueN(.__sample__), by = key_cols_xl]
        }
      }
      .XLSX_ROW_CAP <- 1048575L      # Excel data-row hard limit
      .gvr_safe_sheet <- function(pfx, tok, used) {
        # Replace Excel-illegal chars, build "<PFX>_<token>", truncate to 31, dedup.
        clean <- gsub("[\\\\/?*\\[\\]:]", "_", tok, perl = TRUE)
        nm <- paste0(pfx, clean)
        if (nchar(nm) > 31L) nm <- substr(nm, 1L, 31L)
        base <- nm; k <- 1L
        while (nm %in% used) {
          suf <- paste0("_", k)
          nm <- paste0(substr(base, 1L, 31L - nchar(suf)), suf)
          k <- k + 1L
        }
        nm
      }
      .gvr_xl_cols <- c("Hugo_Symbol", "dbSNP_RS", "CLIN_SIG", "IMPACT",
                        "gnomADe_AF", "Variant_Classification",
                        "Tumor_Sample_Barcode")
      .gvr_xl_cols <- intersect(.gvr_xl_cols, names(dt))
      used_names <- names(wb)
      # Category specs: (prefix, tokens, filter_fn returning logical mask on dt)
      # token order: same composite ranking applied to a one-row-per-token aggregate
      # for IMPACT/VC/Hugo_Symbol; CLIN_SIG keeps the summary-table order (already by Total desc).
      cat_specs <- list()
      if (!is.null(sections$top_genes) && nrow(sections$top_genes) > 0L) {
        cat_specs$top_genes <- list(
          pfx = "GEN_",
          tokens = as.character(sections$top_genes$Hugo_Symbol),
          filter_fn = function(tok) dt$Hugo_Symbol == tok)
      }
      if (!is.null(sections$clin_sig) && nrow(sections$clin_sig) > 0L) {
        cs_tokens <- as.character(sections$clin_sig$CLIN_SIG)
        cs_tokens <- cs_tokens[cs_tokens != "missing/unclassified"]
        if (length(cs_tokens) > 0L) {
          cat_specs$clin_sig <- list(
            pfx = "CS_",
            tokens = cs_tokens,
            filter_fn = function(tok) {
              cs <- dt$CLIN_SIG; m <- !.is_missing(cs)
              # Multi-token match: split on &/ and check membership.
              out <- logical(length(cs))
              if (any(m)) {
                toks_list <- strsplit(cs[m], "[&/]")
                hit <- vapply(toks_list, function(x) tok %in% trimws(x), logical(1L))
                out[m] <- hit
              }
              out
            })
        }
      }
      if (!is.null(sections$impact) && nrow(sections$impact) > 0L) {
        cat_specs$impact <- list(
          pfx = "IMP_",
          tokens = as.character(sections$impact$IMPACT),
          filter_fn = function(tok) !is.na(dt$IMPACT) & dt$IMPACT == tok)
      }
      if (!is.null(sections$variant_classification) && nrow(sections$variant_classification) > 0L) {
        cat_specs$vc <- list(
          pfx = "VC_",
          tokens = as.character(sections$variant_classification$Variant_Classification),
          filter_fn = function(tok) !is.na(dt$Variant_Classification) & dt$Variant_Classification == tok)
      }
      for (cat_nm in names(cat_specs)) {
        sp <- cat_specs[[cat_nm]]
        for (tok in sp$tokens) {
          mask <- sp$filter_fn(tok)
          if (!any(mask)) next                    # skip empty pools
          sub <- dt[mask]
          ord <- .gvr_rank_xl(sub)
          sub <- sub[ord]
          if (nrow(sub) > .XLSX_ROW_CAP) sub <- sub[seq_len(.XLSX_ROW_CAP)]
          out_df <- as.data.frame(sub[, .gvr_xl_cols, with = FALSE], stringsAsFactors = FALSE)
          sh_nm <- .gvr_safe_sheet(sp$pfx, tok, used_names)
          used_names <- c(used_names, sh_nm)
          openxlsx::addWorksheet(wb, sh_nm)
          openxlsx::writeData(wb, sh_nm, out_df, headerStyle = hs)
          openxlsx::freezePane(wb, sh_nm, firstRow = TRUE)
          openxlsx::setColWidths(wb, sh_nm, cols = seq_len(ncol(out_df)), widths = "auto")
        }
      }
      # Write to a local temp file first, then shell-cp to out_subdir (FUSE-safe:
      # openxlsx uses zip random-access writes that can fail / 0-byte on S3-backed mounts).
      tmp_xlsx <- file.path(tempdir(), xlsx_name)
      wrote_ok <- tryCatch({ openxlsx::saveWorkbook(wb, tmp_xlsx, overwrite = TRUE); TRUE },
                           error = function(e) { warning(sprintf("gvr_summary: Excel write failed: %s", conditionMessage(e))); FALSE })
      if (wrote_ok) {
        cp <- system2("cp", c(shQuote(tmp_xlsx), shQuote(final_xlsx)))
        if (!file.exists(final_xlsx) || file.info(final_xlsx)$size == 0) {
          warning(sprintf("gvr_summary: copy to '%s' may have failed; Excel left at '%s'.",
                          final_xlsx, tmp_xlsx))
          final_xlsx <- tmp_xlsx
        }
        if (isTRUE(verbose)) message(sprintf("  Excel written: %s", final_xlsx))
      }
    }
  }

  # ============================================================================
  # Optional PDF report -> <out_dir>/gvr_summary/<file_prefix>_report.pdf (fixed name).
  # Layout: title/metadata page; each CHARTED section (top genes, variant
  # classification, predicted impact) shown as its table together with its bar
  # chart (side-by-side when the table is short, table full-width + chart below
  # when the table is tall); then the remaining chart-less tables (overview,
  # variant type, clinical significance) grouped/stacked on their own page(s).
  # Requires gridExtra + ggplot2; if unavailable the PDF is skipped with a warning
  # (sections are still returned).
  # ============================================================================
  if (isTRUE(save_pdf)) {
    have_pkgs <- requireNamespace("gridExtra", quietly = TRUE) &&
                 requireNamespace("ggplot2",   quietly = TRUE) &&
                 requireNamespace("grid",      quietly = TRUE) &&
                 requireNamespace("scales",    quietly = TRUE)
    if (!have_pkgs) {
      warning("gvr_summary: 'gridExtra'/'ggplot2'/'scales' not installed; skipping PDF report.")
    } else {
      # Internal renderer for the multi-page PDF (nested: single-use, not exported).
      # ASCII-only text (base pdf() substitutes non-ASCII glyphs); writes to a POSIX
      # temp path then shell-cp to final_pdf (R file.copy() can 0-byte on S3 FUSE).
      .gvr_summary_pdf <- function(sections, samples, meta, final_pdf, file_prefix = "gvr_summary") {
        PHYLO_BLUE <- "#0279EE"; PHYLO_GREEN <- "#75A025"; CREAM <- "#FAF9F3"; STONE <- "#ECE9E2"
        INK <- "#000000"; ORANGE <- "#E69F00"; VERMIL <- "#D55E00"

        W <- 8.27; H <- 11.69                       # A4 portrait inches
        USABLE_W_IN <- 7.3                          # printable body width (matches shipped)
        CEX_CORE_FLOOR <- 0.50                      # min table body font scale before col-paginate
        CEX_HEAD_FLOOR <- 0.55
        FACET_THRESHOLD <- 6L                       # > this many samples -> facet charts

        # PDF palette: shared 9-colour base + one neutral grey for a 10th sample
        # before cycling. (HTML extends further; see the HTML renderer.)
        pal <- c(PHYLO_PAL_BASE, "#888888")
        fill_vals <- stats::setNames(pal[((seq_along(samples) - 1L) %% length(pal)) + 1L], samples)

        # ---- measured grob size in inches (drive fit by measurement, not guesses) --------
        .grob_w_in <- function(g) grid::convertWidth(sum(g$widths),  "in", valueOnly = TRUE)
        .grob_h_in <- function(g) grid::convertHeight(sum(g$heights), "in", valueOnly = TRUE)

        # ---- dashboard table theme at a given core font scale ----------------------------
        .mk_theme <- function(cex_core = 0.72) {
          cex_head <- max(CEX_HEAD_FLOOR, cex_core + 0.06)
          gridExtra::ttheme_minimal(
            core    = list(fg_params = list(cex = cex_core, hjust = 1, x = 0.95),
                           bg_params = list(fill = c(CREAM, STONE), col = NA)),
            colhead = list(fg_params = list(col = "white", fontface = "bold", cex = cex_head),
                           bg_params = list(fill = PHYLO_BLUE, col = NA)))
        }

        # ---- robust table-grob factory (PORTED VERBATIM from shipped renderer): -----------
        #   row pagination -> font shrink -> column pagination. Returns flat list of grobs,
        #   each with attr "colspan" (e.g. "samples 1-6 of 20" or NULL) and "n_rows".
        .mk_table_grobs <- function(dt, rows_per_page = 34L) {
          df <- as.data.frame(dt, stringsAsFactors = FALSE)
          for (j in seq_len(ncol(df)))
            if (is.numeric(df[[j]])) df[[j]] <- format(df[[j]], big.mark = ",", trim = TRUE)
          # Dashboard tables are DENSE by default (start at 0.60, not gridExtra's larger
          # default) so the common few-sample report packs tightly; the ladder still steps
          # DOWN toward the floor when a wide many-sample table overflows the page width
          # (after which .mk_table_grobs column-paginates).
          cex_try <- c(0.60, 0.55, CEX_CORE_FLOOR)
          cex_use <- CEX_CORE_FLOOR
          for (cx in cex_try) {
            g0 <- gridExtra::tableGrob(df[1, , drop = FALSE], rows = NULL, theme = .mk_theme(cx))
            if (.grob_w_in(g0) <= USABLE_W_IN) { cex_use <- cx; break }
            cex_use <- cx
          }
          th <- .mk_theme(cex_use)
          ncol_df <- ncol(df)
          probe   <- gridExtra::tableGrob(df[1, , drop = FALSE], rows = NULL, theme = th)
          col_groups <- list(seq_len(ncol_df))
          if (ncol_df >= 3L && .grob_w_in(probe) > USABLE_W_IN) {
            cat_i <- 1L; total_i <- ncol_df
            mid   <- setdiff(seq_len(ncol_df), c(cat_i, total_i))
            wcol  <- vapply(seq_len(ncol_df), function(j)
              .grob_w_in(gridExtra::tableGrob(df[1, j, drop = FALSE], rows = NULL, theme = th)),
              numeric(1))
            fixed_w <- wcol[cat_i] + wcol[total_i]
            budget  <- USABLE_W_IN - fixed_w
            col_groups <- list(); cur <- integer(0); used <- 0
            for (j in mid) {
              if (length(cur) > 0L && used + wcol[j] > budget) {
                col_groups[[length(col_groups) + 1L]] <- c(cat_i, cur, total_i)
                cur <- integer(0); used <- 0
              }
              cur <- c(cur, j); used <- used + wcol[j]
            }
            if (length(cur)) col_groups[[length(col_groups) + 1L]] <- c(cat_i, cur, total_i)
          }
          n_mid_total <- if (length(col_groups) > 1L) ncol_df - 2L else NA_integer_
          npg_rows <- max(1L, ceiling(nrow(df) / rows_per_page))
          out <- list()
          for (cg_idx in seq_along(col_groups)) {
            cg <- col_groups[[cg_idx]]
            colspan_lbl <- if (length(col_groups) > 1L) {
              mids <- setdiff(cg, c(1L, ncol_df))
              sprintf("samples %d-%d of %d", mids[1] - 1L, mids[length(mids)] - 1L, n_mid_total)
            } else NULL
            for (pg in seq_len(npg_rows)) {
              rs <- ((pg - 1L) * rows_per_page + 1L):min(pg * rows_per_page, nrow(df))
              g  <- gridExtra::tableGrob(df[rs, cg, drop = FALSE], rows = NULL, theme = th)
              attr(g, "colspan") <- colspan_lbl
              attr(g, "n_rows")  <- length(rs)
              out[[length(out) + 1L]] <- g
            }
          }
          out
        }

        # ---- bar chart grob: grouped (<=6 samples) or faceted (>6) — PORTED -------------
        .bar_grob <- function(dt, cat_col, top = NULL, title = NULL) {
          d <- as.data.frame(dt, stringsAsFactors = FALSE)
          d <- d[d[[cat_col]] != "" & !is.na(d[[cat_col]]), , drop = FALSE]
          if (!is.null(top) && nrow(d) > top) d <- d[order(-d$Total)[seq_len(top)], , drop = FALSE]
          long <- do.call(rbind, lapply(samples, function(s)
            data.frame(Category = d[[cat_col]], Sample = s, n = d[[s]], stringsAsFactors = FALSE)))
          long$Category <- factor(long$Category, levels = d[[cat_col]][order(d$Total)])
          many <- length(samples) > FACET_THRESHOLD
          gg <- ggplot2::ggplot(long, ggplot2::aes(x = Category, y = n, fill = Sample)) +
            ggplot2::coord_flip() +
            ggplot2::scale_fill_manual(values = fill_vals) +
            ggplot2::labs(title = title, x = NULL, y = NULL, fill = NULL) +
            ggplot2::theme_minimal(base_size = 10) +
            ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", colour = PHYLO_BLUE, size = 13),
                           panel.grid.major.y = ggplot2::element_blank())
          if (many) {
            # Compact facet labels via the shared common-prefix stripper (.lab_fun).
            gg <- gg + ggplot2::geom_col(show.legend = FALSE) +
              ggplot2::facet_wrap(~ Sample, ncol = ceiling(sqrt(length(samples))),
                                  labeller = ggplot2::as_labeller(.lab_fun)) +
              ggplot2::scale_y_continuous(breaks = scales::breaks_extended(3),
                                          labels = scales::label_number(scale_cut = scales::cut_short_scale())) +
              ggplot2::theme(legend.position = "none",
                             strip.text = ggplot2::element_text(size = 7, face = "bold"),
                             axis.text  = ggplot2::element_text(size = 5.5),
                             panel.spacing = grid::unit(4, "pt"))
          } else {
            gg <- gg + ggplot2::geom_col(position = ggplot2::position_dodge(width = 0.75), width = 0.68) +
              ggplot2::scale_y_continuous(breaks = scales::breaks_extended(4),
                                          labels = scales::label_number(scale_cut = scales::cut_short_scale())) +
              ggplot2::theme(legend.position = "top", legend.justification = "left",
                             legend.key.size = grid::unit(10, "pt"),
                             legend.text = ggplot2::element_text(size = 8),
                             axis.text = ggplot2::element_text(size = 9))
          }
          ggplot2::ggplotGrob(gg)
        }

        # ---- left-justify a table grob within a full-width row ---------------------------
        .left_just <- function(g) {
          w_in <- .grob_w_in(g)
          gridExtra::arrangeGrob(g, grid::nullGrob(), ncol = 2,
                                 widths = grid::unit.c(grid::unit(w_in, "in"), grid::unit(1, "null")))
        }

        # ---- KPI card. Big-number font auto-shrinks for long values (e.g. 7-digit cohort
        #      totals) so the number never kisses the card edge. ----------------------------
        .mk_kpi <- function(value, label, fill = PHYLO_BLUE, fg = "white") {
          nchar_v <- nchar(value)
          num_fs  <- if (nchar_v <= 6L) 30 else if (nchar_v <= 8L) 24 else 20
          grid::grobTree(
            grid::roundrectGrob(gp = grid::gpar(fill = fill, col = NA), r = grid::unit(6, "pt")),
            grid::textGrob(value, y = 0.60, gp = grid::gpar(col = fg, fontface = "bold", fontsize = num_fs)),
            grid::textGrob(label, y = 0.22, gp = grid::gpar(col = fg, fontsize = 11)))
        }

        # ---- cairo_pdf blank-first-page guard --------------------------------------------
        .make_new_page <- function() {
          first <- TRUE
          function() { if (first) first <<- FALSE else grid::grid.newpage() }
        }

        # ============================ HERO PAGE ===========================================
        .render_hero <- function(np) {
          np()
          title_grob <- grid::grobTree(
            grid::textGrob("germlinevaR \u2013 Cohort Summary", x = 0.02, y = 0.70, hjust = 0,
                           gp = grid::gpar(fontface = "bold", fontsize = 24, col = PHYLO_BLUE)),
            grid::textGrob(sprintf("%d sample%s  \u00b7  %s total variants  \u00b7  %s distinct genes",
                                   meta$n_samples, if (meta$n_samples == 1L) "" else "s",
                                   fmt(meta$n_total), fmt(meta$n_genes)),
                           x = 0.02, y = 0.26, hjust = 0, gp = grid::gpar(fontsize = 11, col = INK)))
          hi <- sections$impact$Total[sections$impact$IMPACT == "HIGH"]; if (!length(hi)) hi <- 0L
          kpis <- list(
            .mk_kpi(fmt(meta$n_total),    "Total variants",         PHYLO_BLUE),
            .mk_kpi(fmt(meta$n_samples),  "Samples",                PHYLO_GREEN),
            .mk_kpi(fmt(meta$n_genes),    "Distinct genes (known)", ORANGE, fg = INK),
            .mk_kpi(fmt(hi),              "HIGH-impact variants",   VERMIL))
          cards <- lapply(kpis, function(g)
            gridExtra::arrangeGrob(g, vp = grid::viewport(width = 0.94, height = 0.86)))
          cards_row <- gridExtra::arrangeGrob(grobs = cards, nrow = 1)
          chart_top <- .bar_grob(sections$top_genes, "Hugo_Symbol",
                                 title = "Top genes by variant count (top 10)", top = 10L)
          chart_cls <- .bar_grob(sections$variant_classification, "Variant_Classification",
                                 title = "Variant classification (top 10)", top = 10L)
          hero <- gridExtra::arrangeGrob(
            grobs = list(title_grob, cards_row, chart_top, chart_cls), ncol = 1,
            heights = grid::unit.c(grid::unit(0.95, "in"), grid::unit(1.5, "in"),
                                   grid::unit(1, "null"), grid::unit(1, "null")))
          grid::pushViewport(grid::viewport(width = grid::unit(W - 1, "in"),
                                            height = grid::unit(H - 1, "in")))
          grid::grid.draw(hero); grid::popViewport()
        }

        # ===================== AUTO SIDE-BY-SIDE REFERENCE SECTION ========================
        # Each section -> .mk_table_grobs (already font-shrunk / column-paginated). A
        # section's grob is SIDE-BY-SIDE eligible only when it produced a SINGLE grob whose
        # measured width fits half the content width (so two fit with a gap). Otherwise it
        # is FULL-WIDTH and stacks. This degrades safely for many-sample (wide / paginated)
        # tables, matching the user's "side-by-side if it fits, else stack" rule.
        .render_reference <- function(np, target_dev,
                                      title_in = 0.30, tgap_in = 0.10, bgap_in = 0.18,
                                      col_gap_in = 0.30, draw_frac = 0.97) {
          content_w <- W - 1
          budget    <- (H - 1) * draw_frac
          tbl_specs <- .build_tbl_specs(sections, html = FALSE)
          # MEASUREMENT must happen on a throwaway device that is opened AND closed HERE,
          # so the cairo report device is current during all subsequent drawing. (Leaving
          # pdf(NULL) open across the draw loop sends every page to the throwaway device
          # and the report comes out empty/0-byte.)
          items <- local({
            grDevices::pdf(NULL, width = W, height = H)
            on.exit(grDevices::dev.off())
            grid::pushViewport(grid::viewport(width = grid::unit(content_w, "in"),
                                              height = grid::unit(budget, "in")))
            it <- list()
            for (sp in tbl_specs) {
              sec_dt <- if (!is.null(sp$per_sample_key))
                          sections$top_genes_per_sample[[sp$per_sample_key]]
                        else sections[[sp$s]]
              gl <- .mk_table_grobs(sec_dt)
              multi <- length(gl) > 1L
              for (gi in seq_along(gl)) {
                g  <- gl[[gi]]
                cs <- attr(g, "colspan")
                lbl <- if (multi && !is.null(cs)) sprintf("%s (%s)", sp$t, cs) else sp$t
                wi <- .grob_w_in(g); hi <- .grob_h_in(g)
                # Side-by-side eligible iff this section produced a SINGLE grob (not
                # column-paginated). Whether two eligible tables actually fit a row is
                # decided later by the real pairwise width sum (a$w + gap + b$w <= content_w),
                # which is more permissive than a rigid half-width gate (a wide table can
                # still pair with a narrow one). Column-paginated (multi) tables stay
                # full-width and stack -- the many-sample safe fallback.
                eligible <- !multi
                it[[length(it) + 1L]] <- list(type = "table", lbl = lbl, g = g,
                                              w = wi, h = hi, pair_ok = eligible)
              }
            }
            grid::popViewport()
            it
          })
          # Closing the throwaway measurement device above can leave a DIFFERENT device
          # current (R falls back to the most-recent surviving device, not necessarily our
          # report device). Force the report device current before any drawing.
          if (target_dev %in% grDevices::dev.list()) grDevices::dev.set(target_dev)
          imp_chart <- list(type = "chart", lbl = "Functional impact (VEP IMPACT)",
                            g = .bar_grob(sections$impact, "IMPACT", title = NULL),
                            w = content_w, h = 2.4, pair_ok = FALSE)

          # ---- height-balanced pairing among pair_ok items; others stay single --------
          ord  <- order(vapply(items, function(u) u$h, numeric(1)), decreasing = TRUE)
          done <- rep(FALSE, length(items)); rows <- list()
          for (ii in seq_along(ord)) {
            i <- ord[ii]; if (done[i]) next
            a <- items[[i]]
            if (!isTRUE(a$pair_ok)) { done[i] <- TRUE
              rows[[length(rows)+1]] <- list(kind = "single", a = a, h = a$h); next }
            partner <- NA_integer_
            for (jj in seq_along(ord)) {
              j <- ord[jj]
              if (j == i || done[j] || !isTRUE(items[[j]]$pair_ok)) next
              if (a$w + col_gap_in + items[[j]]$w <= content_w) { partner <- j; break }
            }
            if (!is.na(partner)) {
              b <- items[[partner]]; done[i] <- TRUE; done[partner] <- TRUE
              rows[[length(rows)+1]] <- list(kind = "pair", a = a, b = b, h = max(a$h, b$h))
            } else {
              done[i] <- TRUE
              rows[[length(rows)+1]] <- list(kind = "single", a = a, h = a$h)
            }
          }
          rows[[length(rows)+1]] <- list(kind = "single", a = imp_chart, h = imp_chart$h)

          per_row_overhead <- title_in + tgap_in + bgap_in
          pages <- list(); cur <- list(); used <- 0
          for (r in rows) {
            cost <- r$h + per_row_overhead
            if (length(cur) > 0 && used + cost > budget) { pages[[length(pages)+1]] <- cur; cur <- list(); used <- 0 }
            cur[[length(cur)+1]] <- r; used <- used + cost
          }
          if (length(cur)) pages[[length(pages)+1]] <- cur

          block_for <- function(item, full_width, row_h) {
            body <- if (item$type == "table") { if (full_width) .left_just(item$g) else item$g } else item$g
            # Top-anchor the title within its band (vjust = 1 at y = 1) so its
            # baseline sits high and there is full clearance to the table body that
            # follows the tgap spacer -- prevents the body riding up under the title.
            ttl <- grid::textGrob(item$lbl, x = 0.02, y = 1, hjust = 0, vjust = 1,
                                  gp = grid::gpar(fontsize = 10.5, fontface = "bold", col = PHYLO_BLUE))
            pad <- max(row_h - item$h, 0)
            gridExtra::arrangeGrob(ttl, grid::nullGrob(), body, grid::nullGrob(), ncol = 1,
                                   heights = grid::unit.c(grid::unit(title_in, "in"), grid::unit(tgap_in, "in"),
                                                          grid::unit(item$h, "in"), grid::unit(pad, "in")))
          }

          for (pg in pages) {
            np()
            grobs <- list(grid::textGrob("Reference tables", x = 0.02, hjust = 0, vjust = 1, y = 1,
                                         gp = grid::gpar(fontface = "bold", fontsize = 10, col = PHYLO_GREEN)))
            rel <- c(0.18)
            for (r in pg) {
              if (r$kind == "pair") {
                ga <- block_for(r$a, FALSE, r$h); gb <- block_for(r$b, FALSE, r$h)
                row_grob <- gridExtra::arrangeGrob(ga, grid::nullGrob(), gb, ncol = 3,
                              widths = grid::unit.c(grid::unit(1, "null"), grid::unit(col_gap_in, "in"),
                                                    grid::unit(1, "null")))
              } else {
                row_grob <- block_for(r$a, TRUE, r$h)
              }
              grobs[[length(grobs)+1]] <- row_grob
              rel <- c(rel, r$h + title_in + tgap_in)
              grobs[[length(grobs)+1]] <- grid::nullGrob(); rel <- c(rel, bgap_in)
            }
            if (sum(rel) < budget) { grobs[[length(grobs)+1]] <- grid::nullGrob(); rel <- c(rel, budget - sum(rel)) }
            grid::pushViewport(grid::viewport(y = 0.5, height = grid::unit(budget, "in"),
                                              width = grid::unit(content_w, "in")))
            grid::grid.draw(gridExtra::arrangeGrob(grobs = grobs, ncol = 1, heights = grid::unit(rel, "in")))
            grid::popViewport()
          }
        }

        # ============================ DRIVE THE DEVICE ====================================
        tmp_pdf <- file.path(tempdir(), basename(final_pdf))
        grDevices::cairo_pdf(tmp_pdf, width = W, height = H, onefile = TRUE)
        cairo_dev <- grDevices::dev.cur()
        on.exit(if (cairo_dev %in% grDevices::dev.list()) grDevices::dev.off(cairo_dev), add = TRUE)
        np <- .make_new_page()
        .render_hero(np)
        .render_reference(np, target_dev = cairo_dev)
        if (cairo_dev %in% grDevices::dev.list()) grDevices::dev.off(cairo_dev)
        on.exit()
        if (!file.exists(tmp_pdf) || file.info(tmp_pdf)$size == 0)
          stop("PDF device produced no/zero-byte file.")
        system2("cp", c("-f", shQuote(tmp_pdf), shQuote(final_pdf)))
        invisible(final_pdf)
      }

      pdf_name  <- sprintf("%s_report.pdf", file_prefix)
      final_pdf <- file.path(out_subdir, pdf_name)
      if (file.exists(final_pdf) && isTRUE(verbose))
        message(sprintf("  Overwriting existing PDF report: %s", final_pdf))
      ok_pdf <- tryCatch({
        .gvr_summary_pdf(sections, samples, meta, final_pdf, file_prefix)
        file.exists(final_pdf) && file.info(final_pdf)$size > 0
      }, error = function(e) {
        warning(sprintf("gvr_summary: PDF report failed: %s", conditionMessage(e))); FALSE })
      if (isTRUE(ok_pdf) && isTRUE(verbose)) message(sprintf("  PDF report written: %s", final_pdf))
    }
  }

  # ============================================================================
  # Optional interactive HTML report
  #   -> <out_dir>/gvr_summary/<file_prefix>_report.html (fixed name).
  # An interactive dashboard mirroring the PDF: KPI cards, the same three bar
  # charts (Top genes / Variant classification / Functional impact) as plotly
  # (grouped <=6 samples, faceted small-multiples >6), and all six section
  # tables as searchable/sortable DT widgets. Self-contained single file when
  # pandoc is available; otherwise a <prefix>_report.html + <prefix>_report_files/
  # asset folder (fallback) with a note. Requires plotly + DT + htmlwidgets +
  # htmltools; if unavailable the HTML is skipped with a warning (sections,
  # Excel and PDF are still produced).
  # ============================================================================
  if (isTRUE(save_html)) {
    have_html <- requireNamespace("plotly",      quietly = TRUE) &&
                 requireNamespace("DT",          quietly = TRUE) &&
                 requireNamespace("htmlwidgets", quietly = TRUE) &&
                 requireNamespace("htmltools",   quietly = TRUE)
    if (!have_html) {
      warning("gvr_summary: 'plotly'/'DT'/'htmlwidgets' not installed; skipping HTML report.")
    } else {
      # Pandoc availability probe (kept local so a test can mask it to force the
      # sidecar fallback). Needs both 'rmarkdown' (the pandoc wrapper used for
      # asset-inlining) and a working pandoc; if either is missing we return FALSE
      # and the renderer writes the sidecar (.html + _files/) form instead.
      .gvr_pandoc_ok <- function() {
        requireNamespace("rmarkdown", quietly = TRUE) &&
          isTRUE(tryCatch(rmarkdown::pandoc_available(), error = function(e) FALSE))
      }

      # --- Adaptive size-band policy (HTML drill-down tables) ----------------
      # The HTML report self-contains every detail-table row as inline JSON. At
      # the pilot (~22k variants) the file is ~8.6 MB; linear projection to
      # 300k+ variants overflows pandoc's --self-contained pass on Windows
      # (VirtualAlloc / pagefile). To keep the report self-containable at scale
      # we trim the inline payload as a function of cohort size:
      #
      #   small   (< .HTML_SMALL_MAX)   : full drill-downs (current behaviour)
      #   medium  (< .HTML_MEDIUM_MAX)  : detail tables pre-filtered to the
      #                                   biologically informative subset
      #                                   (no MODIFIER impact / no Intron
      #                                   classification / no missing CLIN_SIG).
      #                                   VC additionally capped at 50k rows.
      #   large   (>= .HTML_MEDIUM_MAX) : drill-downs disabled entirely. The
      #                                   CLIN_SIG / VC / IMPACT sections
      #                                   render as plain summary tables only;
      #                                   per-variant data lives in the Excel.
      #
      # Tunable in one place. Boundaries are conservative against the pilot's
      # MODIFIER/Intron fractions and pandoc's observed working set.
      .HTML_SMALL_MAX      <- 50000L
      .HTML_PAYLOAD_BUDGET <- 100000L   # asymptotic row budget for inline JSON

      .gvr_size_band <- function(n_total) {
        if (n_total < .HTML_SMALL_MAX) "small" else "medium"
      }


      # Internal renderer for the interactive HTML (nested: single-use, not exported).
      # Reuses the PDF palette / KPI / chart logic so both reports stay consistent.
      # Returns the path actually written, with attr "sidecar" (TRUE when pandoc was
      # unavailable and assets live in a sibling <prefix>_report_files/ folder) and
      # attr "files_dir" (that folder, or NA when self-contained).
      .gvr_summary_html <- function(sections, samples, meta, final_html, file_prefix = "gvr_summary") {

        PHYLO_BLUE <- "#0279EE"; PHYLO_GREEN <- "#75A025"; CREAM <- "#FAF9F3"; STONE <- "#ECE9E2"
        INK <- "#000000"; ORANGE <- "#E69F00"; VERMIL <- "#D55E00"

        # 24-colour qualitative palette: shared 9-colour PHYLO_PAL_BASE (Okabe-Ito
        # ordering, matches the PDF dashboard exactly) extended with Paul-Tol
        # "muted" + extras for cohorts up to 24 (cycles only beyond that --
        # an acceptable edge case).
        pal <- c(PHYLO_PAL_BASE,
                 "#117733", "#882255", "#88CCEE", "#999933", "#661100",
                 "#6699CC", "#DDCC77", "#44AA99", "#AA4466", "#4477AA",
                 "#228833", "#EE6677", "#BBBBBB", "#000000", "#EE3377")
        fill_vals <- stats::setNames(pal[((seq_along(samples) - 1L) %% length(pal)) + 1L], samples)

        # one interactive plotly bar chart: GROUPED horizontal bars (one bar per
        # sample) at every cohort size. Unlike the static PDF (which facets >6
        # samples because print cannot toggle), the HTML keeps a single grouped
        # plot and leans on plotly's interactivity - click a legend entry to
        # isolate/hide a sample, hover for exact counts - which stays legible and
        # useful even for large cohorts. Sample legend labels reuse the PDF's
        # common-prefix stripping so they stay compact.
        .plt_bar <- function(dt, cat_col, top = NULL, title = NULL) {
          d <- as.data.frame(dt, stringsAsFactors = FALSE)
          d <- d[d[[cat_col]] != "" & !is.na(d[[cat_col]]), , drop = FALSE]
          if (!is.null(top) && nrow(d) > top) d <- d[order(-d$Total)[seq_len(top)], , drop = FALSE]
          lev <- d[[cat_col]][order(d$Total)]          # ascending so biggest sits on top
          slab <- .lab_fun(samples)                    # compact (prefix-stripped) labels

          long <- do.call(rbind, lapply(seq_along(samples), function(k)
            data.frame(Category = factor(d[[cat_col]], levels = lev),
                       Sample   = factor(slab[k], levels = slab),
                       n        = d[[samples[k]]], stringsAsFactors = FALSE)))

          # Plot height grows with the number of categories (so grouped bars stay
          # readable) AND with the number of samples (so the per-sample legend shows
          # every entry instead of collapsing into a short scrollable box). Measured
          # in-browser, plotly's vertical legend needs ~20 px/entry of *content*, and
          # it only avoids engaging its internal scrollbox when the plotting area is
          # comfortably taller than that content. We therefore budget ~22 px/entry plus
          # generous top/bottom padding so all samples stay visible without scrolling.
          h_cat <- 90 + 26 * length(lev)
          h_leg <- 120 + 22 * length(samples)
          h_px  <- max(360, h_cat, h_leg)
          p <- plotly::plot_ly(long, x = ~n, y = ~Category, color = ~Sample,
                               colors = unname(fill_vals), type = "bar", orientation = "h",
                               height = h_px,
                               hovertemplate = "%{y}: %{x:,}<extra>%{fullData.name}</extra>")
          plotly::layout(p, barmode = "group",
                         title  = list(text = title, x = 0, font = list(color = PHYLO_BLUE, size = 16)),
                         xaxis  = list(title = "", tickformat = ","),
                         yaxis  = list(title = ""),
                         legend = list(title = list(text = "Sample"), font = list(size = 9)),
                         margin = list(l = 8, r = 8, t = 50, b = 8))
        }

        # one DT table for a section (sort/search/paginate; comma-formatted numerics)
        .dt_tbl <- function(dt, caption = NULL) {
          df <- as.data.frame(dt, stringsAsFactors = FALSE)
          num_cols <- names(df)[vapply(df, is.numeric, logical(1))]
          dtbl <- DT::datatable(
            df, rownames = FALSE, caption = caption,
            class = "stripe hover compact", filter = "none",
            options = list(pageLength = 10, lengthMenu = c(5, 10, 25, 50),
                           dom = "ftip", scrollX = TRUE,
                           columnDefs = list(list(className = "dt-right",
                                                  targets = which(names(df) %in% num_cols) - 1L))))
          if (length(num_cols))
            dtbl <- DT::formatCurrency(dtbl, num_cols, currency = "", interval = 3, mark = ",", digits = 0)
          dtbl
        }

        # KPI card (styled div mirroring the PDF cards)
        .kpi_card <- function(value, label, fill = PHYLO_BLUE, fg = "white") {
          htmltools::tags$div(
            style = sprintf(paste0("flex:1; min-width:150px; background:%s; color:%s; border-radius:6px;",
                                   "padding:14px 16px; box-shadow:0 1px 3px rgba(0,0,0,.12);"), fill, fg),
            htmltools::tags$div(value, style = "font-size:30px; font-weight:700; line-height:1.1;"),
            htmltools::tags$div(label, style = "font-size:12px; opacity:.92; margin-top:4px;"))
        }

        hi <- sections$impact$Total[sections$impact$IMPACT == "HIGH"]; if (!length(hi)) hi <- 0L
        cards <- htmltools::tags$div(
          style = "display:flex; gap:12px; flex-wrap:wrap; margin:14px 0 22px 0;",
          .kpi_card(fmt(meta$n_total),   "Total variants",         PHYLO_BLUE),
          .kpi_card(fmt(meta$n_samples), "Samples",                PHYLO_GREEN),
          .kpi_card(fmt(meta$n_genes),   "Distinct genes (known)", ORANGE, fg = INK),
          .kpi_card(fmt(hi),             "HIGH-impact variants",   VERMIL))

        header <- htmltools::tags$div(
          htmltools::tags$h1("germlinevaR \u2013 Cohort Summary",
                             style = sprintf("color:%s; margin:0; font-size:26px;", PHYLO_BLUE)),
          htmltools::tags$div(
            sprintf("%d sample%s \u00b7 %s total variants \u00b7 %s distinct genes \u00b7 generated %s",
                    meta$n_samples, if (meta$n_samples == 1L) "" else "s",
                    fmt(meta$n_total), fmt(meta$n_genes), meta$generated),
            style = sprintf("color:%s; font-size:13px; margin-top:4px;", INK)))

        sec_h <- function(txt) htmltools::tags$h2(
          txt, style = sprintf("color:%s; border-bottom:2px solid %s; padding-bottom:4px; margin-top:30px;",
                               PHYLO_GREEN, STONE))

        # Per-sample top genes: simple horizontal bar chart (single sample, no grouping)
        .plt_bar_single <- function(dt, cat_col, val_col, top = 10L, title = NULL) {
          d <- as.data.frame(dt, stringsAsFactors = FALSE)
          d <- d[d[[cat_col]] != "" & !is.na(d[[cat_col]]), , drop = FALSE]
          if (!is.null(top) && nrow(d) > top) d <- d[order(-d[[val_col]])[seq_len(top)], , drop = FALSE]
          lev <- d[[cat_col]][order(d[[val_col]])]   # ascending so biggest on top
          d[[cat_col]] <- factor(d[[cat_col]], levels = lev)
          h_px <- max(280, 90 + 26 * length(lev))
          p <- plotly::plot_ly(d, x = stats::as.formula(paste0("~", val_col)),
                               y = stats::as.formula(paste0("~", cat_col)),
                               type = "bar", orientation = "h", height = h_px,
                               hovertemplate = "%{y}: %{x:,}<extra></extra>",
                               marker = list(color = PHYLO_BLUE))
          plotly::layout(p, barmode = "group",
                         title  = list(text = title, x = 0, font = list(color = PHYLO_BLUE, size = 14)),
                         xaxis  = list(title = "", tickformat = ","),
                         yaxis  = list(title = ""),
                         margin = list(l = 8, r = 8, t = 44, b = 8),
                         showlegend = FALSE)
        }

        charts <- htmltools::tagList(
          sec_h("Charts"),
          .plt_bar(sections$top_genes, "Hugo_Symbol", top = 10L, title = "Top genes by variant count (top 10)"),
          .plt_bar(sections$variant_classification, "Variant_Classification", top = 10L,
                   title = "Variant classification (top 10)"),
          .plt_bar(sections$impact, "IMPACT", top = NULL, title = "Functional impact (VEP IMPACT)"),
          if (!is.null(sections$top_variants))
            .plt_bar(sections$top_variants, "dbSNP_RS", top = 10L,
                     title = "Top variants by rsID frequency (top 10)"),
          # Per-sample top genes charts (top 10 genes per sample)
          if (length(sections$top_genes_per_sample))
            lapply(names(sections$top_genes_per_sample), function(sm) {
              sm_dt <- sections$top_genes_per_sample[[sm]]
              .plt_bar_single(sm_dt, "Hugo_Symbol", sm, top = 10L,
                              title = sprintf("Top genes \u2013 %s (top 10)", sm))
            }))

        # --- Drill-down detail tables (CLIN_SIG, IMPACT, Variant_Classification) ---
        # Clicking a category token in the summary table filters the detail table
        # on the relevant column and toggles its visibility.
        # IMPORTANT: the detail container must NOT use display:none at build time,
        # because DT cannot initialise inside a hidden element (causes TN/2 warning).
        # Instead, the detail DT's own initComplete hides the container after it
        # has fully rendered.
        #
        # CRITICAL: the click handler uses dtApi.column(idx).search(val) for
        # COLUMN-SPECIFIC filtering (not dtApi.search(val) which does a global
        # search across all columns and returns wrong results).

        # --- Helper: build a drill-down pair (summary DT + detail DT) ---
        # Returns a named list with: summary_dtbl, detail_dtbl, container_id
        .mk_drilldown <- function(summary_section, token_map, container_id,
                                  filter_col, caption_summary, header_col_label,
                                  cap_note = NULL) {
          # ------------------------------------------------------------------
          # vN+7.1 redesign: per-token detail tables via JSON lookup.
          # The detail DT is now an EMPTY shell initialised with column names
          # only; clicking a summary cell looks up `token_map[<val>]` in a
          # hidden <script type="application/json"> block, clears the table,
          # and adds the per-token ranked rows. This guarantees that the row
          # order shown after a click is the SAME order written to the XLSX
          # `<PFX>_<token>` sheet (both built from the same per-token slice),
          # closing the CS7 acceptance gate that the prior dedup-by-locus
          # design failed for 11/18 CS_ tokens.
          # ------------------------------------------------------------------
          # Compose the "empty shell" for the detail DT: a 0-row data frame
          # carrying just the column names. The actual rows arrive client-side.
          if (length(token_map) > 0L) {
            shell_cols <- names(token_map[[1L]])
          } else {
            shell_cols <- character(0)
          }
          shell_df <- as.data.frame(
            stats::setNames(rep(list(character(0)), length(shell_cols)),
                            shell_cols),
            stringsAsFactors = FALSE)

          detail_dtbl <- suppressWarnings(DT::datatable(
            shell_df, rownames = FALSE,
            class = "stripe hover compact", filter = "none",
            options = list(pageLength = 5, lengthMenu = c(5, 10, 25),
                           dom = "ftip", scrollX = TRUE, searchDelay = 300,
                           initComplete = htmlwidgets::JS(
                             "function() {",
                             sprintf("  document.getElementById('%s').style.display = 'none';", container_id),
                             # Wire the Close button. Listener attached once;
                             # button persists across redraws.
                             sprintf("  var btn = document.getElementById('%s_close');", container_id),
                             "  if (btn) {",
                             "    btn.addEventListener('click', function() {",
                             sprintf("      var c = document.getElementById('%s');", container_id),
                             "      if (!c) return;",
                             "      c.style.display = 'none';",
                             "    });",
                             "  }",
                             "}"))))
          # Defensive eager init (see prior comment block in pre-N+7.1 code
          # for the htmlwidgets lazyRender gate rationale).
          detail_dtbl$x$lazyRender <- FALSE

          # Encode the per-token map as a hidden JSON <script> block. Each
          # token maps to an array of row arrays (dataframe="values" -> rows
          # as nested arrays, matching what dtApi.rows.add() consumes).
          # NA encoded as JSON null; DT renders null as an empty cell.
          json_id <- paste0(container_id, "_data")
          json_str <- if (length(token_map) > 0L) {
            as.character(jsonlite::toJSON(token_map, dataframe = "values",
                                          na = "null", auto_unbox = FALSE))
          } else {
            "{}"
          }
          # JSON-encoded strings cannot contain a literal "</script" sequence
          # without breaking the closing tag parser; escape the forward slash.
          json_str <- gsub("</", "<\\\\/", json_str, fixed = TRUE)
          data_script <- htmltools::tags$script(
            type = "application/json",
            id   = json_id,
            htmltools::HTML(json_str))

          # Summary DT: clickable first-column cells.
          summary_dtbl <- .dt_tbl(summary_section, caption = caption_summary)
          summary_dtbl$x$options$initComplete <- htmlwidgets::JS(
            "function() {",
            "  var tbl = this.api().table().body();",
            "  var cells = tbl.querySelectorAll('td:first-child');",
            sprintf("  var dataEl = document.getElementById('%s');", json_id),
            "  var tokenMap = {};",
            "  if (dataEl) {",
            "    try { tokenMap = JSON.parse(dataEl.textContent || '{}'); }",
            "    catch (e) { tokenMap = {}; }",
            "  }",
            "  cells.forEach(function(c) {",
            "    var val0 = c.textContent.trim();",
            "    if (val0 === 'missing/unclassified') return;",
            # Hide non-clickable behaviour for tokens not present in tokenMap
            # (defensive; zero-count rows are pre-filtered at section build).
            "    if (!Object.prototype.hasOwnProperty.call(tokenMap, val0)) return;",
            "    c.style.cursor = 'pointer';",
            "    c.style.textDecoration = 'underline';",
            "    c.style.color = '#0279EE';",
            "    c.addEventListener('click', function() {",
            "      var val = this.textContent.trim();",
            "      if (!Object.prototype.hasOwnProperty.call(tokenMap, val)) return;",
            sprintf("      var container = document.getElementById('%s');", container_id),
            "      if (!container) return;",
            "      var bodyTable = container.querySelector('.dataTables_scrollBody table');",
            "      if (!bodyTable) {",
            "        var allTables = container.querySelectorAll('table');",
            "        bodyTable = allTables[allTables.length - 1];",
            "      }",
            "      var dtApi = bodyTable ? new $.fn.dataTable.Api(bodyTable) : null;",
            "      if (dtApi && dtApi.context.length === 0) dtApi = null;",
            "      if (!dtApi) return;",
            "      var rows = tokenMap[val] || [];",
            "      dtApi.clear();",
            "      if (rows.length) dtApi.rows.add(rows);",
            # Reset to page 0 (top) on every new token click -- predictable
            # behaviour matching user expectation that a new token starts
            # at the top of its ranking.
            "      dtApi.page(0).draw();",
            "      container.style.display = 'block';",
            "      dtApi.columns.adjust();",
            "      var n = rows.length;",
            "      var nFmt = n.toLocaleString();",
            sprintf("      var v = document.getElementById('%s_value');", container_id),
            sprintf("      var ct = document.getElementById('%s_count');", container_id),
            "      if (v)  v.textContent  = val;",
            "      if (ct) ct.textContent = nFmt + ' variant' + (n === 1 ? '' : 's');",
            "    });",
            "  });",
            "}")

          # ------------------------------------------------------------------
          # Panel wrapper: cream background + 4px yellow left edge accent +
          # Close button. The panel is hidden at build time via the detail
          # DT's initComplete; clicking a summary token shows it. The
          # per-token JSON <script> block lives INSIDE the panel so it's
          # carried along with the panel in the final HTML.
          # ------------------------------------------------------------------
          panel_tag <- htmltools::tags$div(
            id    = container_id,
            style = paste0("display:none;",
                           " margin-top:12px;",
                           " background:#FAF9F3;",
                           " border-left:4px solid #E9ED4C;",
                           " border-radius:6px;",
                           " padding:14px 16px;"),
            # Header row: title (left, dynamic) + Close button (right).
            htmltools::tags$div(
              style = paste0("display:flex; align-items:center;",
                             " justify-content:space-between;",
                             " margin-bottom:10px;"),
              htmltools::tags$div(
                style = sprintf("font-size:15px; font-weight:bold; color:%s;", "#0279EE"),
                sprintf("%s = ", header_col_label),
                htmltools::tags$span(id = paste0(container_id, "_value"),
                                     htmltools::HTML("&hellip;")),
                htmltools::tags$span("  |  ",
                                     style = "color:#888; font-weight:normal;"),
                htmltools::tags$span(id = paste0(container_id, "_count"),
                                     style = "font-weight:bold;",
                                     htmltools::HTML("&hellip;"))),
              htmltools::tags$button(
                id    = paste0(container_id, "_close"),
                type  = "button",
                style = paste0("background:transparent; border:1px solid #ccc;",
                               " border-radius:4px; padding:4px 10px;",
                               " cursor:pointer; font-size:12px; color:#555;"),
                htmltools::HTML("&times; Close"))),
            if (!is.null(cap_note))
              htmltools::tags$div(
                style = paste0("font-size:11.5px; color:#666;",
                               " font-style:italic; margin-bottom:8px;"),
                cap_note),
            data_script,
            detail_dtbl)

          list(summary_dtbl = summary_dtbl,
               panel_tag    = panel_tag,
               container_id = container_id)
        }

        # --- Drill-down build (CLIN_SIG / IMPACT / Variant_Classification) ---
        # Each entry below drives one drill-down pair: a clickable summary DT (the
        # category section) plus a hidden detail DT showing the underlying variant
        # rows (curated column set). The detail rows are filtered to non-missing
        # values of the click-target column; rows where the filter column is
        # NA/"" can never be matched by the regex search anyway, so dropping them
        # keeps the detail table compact. Adding a new drill-down = adding one
        # entry here.
        dd_specs <- list(
          clin_sig = list(
            section_key = "clin_sig", filter_col = "CLIN_SIG",
            container = "gvr_clin_detail",
            cols = c("Hugo_Symbol", "dbSNP_RS", "CLIN_SIG", "IMPACT",
                     "gnomADe_AF", "Variant_Classification",
                     "Tumor_Sample_Barcode"),
            cap_summary  = "Clinical significance (click a token to see variants)",
            header_label = "CLIN_SIG"),
          impact = list(
            section_key = "impact", filter_col = "IMPACT",
            container = "gvr_impact_detail",
            cols = c("Hugo_Symbol", "dbSNP_RS", "CLIN_SIG", "IMPACT",
                     "gnomADe_AF", "Variant_Classification",
                     "Tumor_Sample_Barcode"),
            cap_summary  = "Functional impact (click a level to see variants)",
            header_label = "IMPACT"),
          vc = list(
            section_key = "variant_classification",
            filter_col = "Variant_Classification",
            container = "gvr_vc_detail",
            cols = c("Hugo_Symbol", "dbSNP_RS", "CLIN_SIG", "IMPACT",
                     "gnomADe_AF", "Variant_Classification",
                     "Tumor_Sample_Barcode"),
            cap_summary  = "Variant classification (click a class to see variants)",
            header_label = "Variant_Classification"),
          # vN+7: Top genes becomes a 4th clickable. Click handler filters detail
          # by exact Hugo_Symbol match.
          top_genes = list(
            section_key = "top_genes", filter_col = "Hugo_Symbol",
            container = "gvr_topgenes_detail",
            cols = c("Hugo_Symbol", "dbSNP_RS", "CLIN_SIG", "IMPACT",
                     "gnomADe_AF", "Variant_Classification",
                     "Tumor_Sample_Barcode"),
            cap_summary  = "Top genes (click a symbol to see variants)",
            header_label = "Hugo_Symbol"))

        # ----------------------------------------------------------------
        # Adaptive build: classify cohort size, then build the per-section
        # detail data accordingly.
        #   small  -> full detail (current behaviour, after non-missing filter)
        #   medium -> drop low-information rows per section
        #             (MODIFIER for IMPACT; Intron for VC; missing for CLIN_SIG
        #             which is already filtered out by the non-missing step);
        #             additionally cap VC detail at .HTML_SMALL_MAX (50k)
        #             rows ordered by IMPACT severity, since after dropping
        #             Intron the remaining classes can still exceed the cap.
        #   large  -> no detail tables at all; dd_map stays NULL and the
        #             tables branch falls through to plain summary tables.
        # ----------------------------------------------------------------
        # vN+7: composite ranking helper -- IMPACT > nsamp > gnomADe_AF >
        # genomic position. Returns integer order over rows of `d`. Sample
        # count is precomputed per locus (Chromosome+Start_Position+
        # Reference_Allele+Tumor_Seq_Allele2) so all drill-downs share it.
        .impact_rank <- function(v) {
          ord <- c("HIGH" = 1L, "MODERATE" = 2L, "LOW" = 3L, "MODIFIER" = 4L)
          r <- ord[as.character(v)]
          r[is.na(r)] <- 5L   # unknown impact sorts last
          as.integer(r)
        }
        .gvr_rank_variants <- function(d) {
          ir <- .impact_rank(d$IMPACT)
          nsamp <- if (".__nsamp__" %in% names(d)) d$.__nsamp__ else rep(0L, nrow(d))
          gaf <- if ("gnomADe_AF" %in% names(d)) suppressWarnings(as.numeric(d$gnomADe_AF))
                  else rep(NA_real_, nrow(d))
          chr <- if ("Chromosome" %in% names(d)) as.character(d$Chromosome) else rep("", nrow(d))
          pos <- if ("Start_Position" %in% names(d)) suppressWarnings(as.integer(d$Start_Position))
                  else rep(0L, nrow(d))
          order(ir, -nsamp, gaf, chr, pos, na.last = TRUE)
        }

        # vN+7: precompute distinct-sample count per locus on `dt` (once).
        # Falls back to (Chromosome,Start_Position) if allele cols missing.
        if (!(".__nsamp__" %in% names(dt))) {
          key_cols <- intersect(c("Chromosome", "Start_Position",
                                  "Reference_Allele", "Tumor_Seq_Allele2"),
                                names(dt))
          if (length(key_cols) < 2L) {
            dt[, .__nsamp__ := 1L]   # cannot key variants; tiebreak inert
          } else {
            dt[, .__nsamp__ := data.table::uniqueN(.__sample__),
               by = key_cols]
          }
        }

        # vN+7 (refined): per-TOKEN cap. Each clickable category emits one
        # detail JSON formed by concatenating per-token slices, each slice
        # capped at .HTML_TOPK_CAP rows. The click handler still uses the
        # existing token-bounded regex search to surface only the clicked
        # token's rows, so no JS changes are needed.
        #
        # Payload budget: total embedded rows across all 4 categories is
        # bounded by sum_over_tokens(min(pool_size, cap)). If projection
        # exceeds .HTML_PAYLOAD_BUDGET the drill-downs are dropped and the
        # tables fall back to plain summary view -- same behaviour as the
        # legacy 'large' band but driven by actual payload size, not
        # cohort size. With the standard token vocabularies (CLIN_SIG ~18,
        # IMPACT 4, VC ~17, top genes 20) the asymptote is ~60k rows,
        # well under the budget, so the fallback only triggers for unusual
        # parameter choices (e.g. top_n_genes >> 20).
        .HTML_TOPK_CAP <- 1000L

        # Per-category token enumeration -- exact same logic as the XLSX
        # writer above so click and sheet semantics stay aligned.
        .dd_token_specs <- list(
          clin_sig = list(
            tokens = {
              cs_t <- as.character(sections$clin_sig$CLIN_SIG)
              cs_t[cs_t != "missing/unclassified"]
            },
            filter_fn = function(tok) {
              cs <- dt$CLIN_SIG; m <- !.is_missing(cs)
              out <- logical(length(cs))
              if (any(m)) {
                toks_list <- strsplit(cs[m], "[&/]")
                hit <- vapply(toks_list, function(x) tok %in% trimws(x), logical(1L))
                out[m] <- hit
              }
              out
            }),
          impact = list(
            tokens = as.character(sections$impact$IMPACT),
            filter_fn = function(tok) !is.na(dt$IMPACT) & dt$IMPACT == tok),
          vc = list(
            tokens = as.character(sections$variant_classification$Variant_Classification),
            filter_fn = function(tok) !is.na(dt$Variant_Classification) &
                                       dt$Variant_Classification == tok),
          top_genes = list(
            tokens = as.character(sections$top_genes$Hugo_Symbol),
            filter_fn = function(tok) !is.na(dt$Hugo_Symbol) & dt$Hugo_Symbol == tok))

        # Project the embed payload first to decide on fallback.
        .dd_project <- function() {
          tot <- 0L
          for (cat_nm in names(dd_specs)) {
            spec <- .dd_token_specs[[cat_nm]]
            if (is.null(spec) || length(spec$tokens) == 0L) next
            for (tok in spec$tokens) {
              n_tok <- sum(spec$filter_fn(tok))
              tot <- tot + min(n_tok, .HTML_TOPK_CAP)
            }
          }
          tot
        }
        payload_proj <- .dd_project()
        budget_exceeded <- payload_proj > .HTML_PAYLOAD_BUDGET

        # Build per-category detail_df: concatenation of per-token top-K slices.
        .dd_build_cat <- function(cat_nm) {
          # vN+7.1: return a per-TOKEN map (named list, token -> data.frame
          # of top-K ranked rows). The XLSX writer above produces the
          # equivalent slice with `cap = Inf`; the per-token order is
          # identical because both call sites apply .gvr_rank_variants() to
          # the same per-token mask. This guarantees CS7: XLSX `<PFX>_<tok>`
          # row 1 == HTML token_map[tok] row 1.
          sp <- dd_specs[[cat_nm]]
          spec <- .dd_token_specs[[cat_nm]]
          cols <- intersect(sp$cols, names(dt))
          if (is.null(spec) || length(spec$tokens) == 0L) {
            return(list(token_map = list(), n_pre_total = 0L,
                        any_capped = FALSE, n_tokens = 0L,
                        n_nonempty = 0L))
          }
          token_map <- list()
          n_pre_total <- 0L
          any_capped  <- FALSE
          n_nonempty  <- 0L
          for (i in seq_along(spec$tokens)) {
            tok <- spec$tokens[[i]]
            mask <- spec$filter_fn(tok)
            n_tok <- sum(mask)
            n_pre_total <- n_pre_total + n_tok
            if (n_tok == 0L) next
            sub <- dt[mask]
            ord <- .gvr_rank_variants(sub)
            if (n_tok > .HTML_TOPK_CAP) {
              sub <- sub[ord[seq_len(.HTML_TOPK_CAP)]]
              any_capped <- TRUE
            } else {
              sub <- sub[ord]
            }
            # Project to the rendered 7 columns. No dedup, no concat: each
            # token's slice is independent and carries its own canonical
            # ranking.
            slice_df <- as.data.frame(sub[, ..cols], stringsAsFactors = FALSE)
            token_map[[as.character(tok)]] <- slice_df
            n_nonempty <- n_nonempty + 1L
          }
          list(token_map = token_map, n_pre_total = n_pre_total,
               any_capped = any_capped, n_tokens = length(spec$tokens),
               n_nonempty = n_nonempty)
        }

        dd_map <- if (budget_exceeded) NULL else lapply(names(dd_specs), function(cat_nm) {
          sp <- dd_specs[[cat_nm]]
          built <- .dd_build_cat(cat_nm)
          # If no token has any matching variants (e.g. degenerate cohort),
          # skip the drill-down entirely -- the section falls back to its
          # plain summary table via the dd_map[[special]] is.null branch.
          if (built$n_nonempty == 0L) return(NULL)
          cap_note <- if (built$any_capped) sprintf(
            "Each row group is capped at %s variants, ranked by IMPACT severity, sample count, then rarity (gnomADe_AF). The Excel workbook has the full uncapped data per group.",
            format(.HTML_TOPK_CAP, big.mark = ",")) else NULL
          .mk_drilldown(sections[[sp$section_key]], built$token_map,
                        sp$container, sp$filter_col,
                        sp$cap_summary, sp$header_label,
                        cap_note = cap_note)
        })
        if (!is.null(dd_map)) names(dd_map) <- names(dd_specs)

        tbl_specs <- .build_tbl_specs(sections, html = TRUE)

          .tables_head <- function() {
            head <- list(sec_h("Tables"))
            if (is.null(dd_map)) {
              banner <- htmltools::tags$div(
                style = paste0("margin:10px 0 14px 0; padding:10px 14px;",
                               " background:#FAF9F3; border-left:4px solid #FF9400;",
                               " border-radius:6px; font-size:12.5px; color:#333;"),
                htmltools::tags$strong(sprintf(
                  "Cohort has %s variants. ",
                  format(meta$n_total, big.mark = ","))),
                "Drill-down detail tables are disabled because the projected ",
                "inline payload would exceed the HTML self-contain budget. ",
                "The full per-variant data is available in the Excel workbook ",
                "(one sheet per token).")
              head <- c(head, list(banner))
            }
            head
          }

        tables <- htmltools::tagList(c(
          .tables_head(),
          lapply(tbl_specs, function(sp) {
            sec_dt <- if (!is.null(sp$per_sample_key))
                        sections$top_genes_per_sample[[sp$per_sample_key]]
                      else sections[[sp$s]]
            if (!is.null(sp$special) && !is.null(dd_map) && !is.null(dd_map[[sp$special]])) {
              # Drill-down section (small + medium bands): clickable summary
              # table on top, cream panel with dynamic header + detail table
              # below. The panel is hidden until a summary cell is clicked.
              dd <- dd_map[[sp$special]]
              htmltools::tags$div(
                style = "margin:14px 0 26px 0;",
                htmltools::tags$h3(sp$t,
                  style = sprintf("color:%s; font-size:15px; margin:0 0 6px 0;", PHYLO_BLUE)),
                dd$summary_dtbl,
                dd$panel_tag)
            } else {
              # Plain section (non-drill-down sections, OR drill-down sections
              # falling through here in the large band: dd_map is NULL).
              htmltools::tags$div(style = "margin:14px 0 26px 0;",
                                  htmltools::tags$h3(sp$t,
                                    style = sprintf("color:%s; font-size:15px; margin:0 0 6px 0;", PHYLO_BLUE)),
                                  .dt_tbl(sec_dt))
            }
          })))

        page <- htmltools::tags$div(
          style = paste0("max-width:1100px; margin:24px auto; padding:0 18px;",
                         " font-family:-apple-system,Segoe UI,Roboto,Helvetica,Arial,sans-serif;"),
          header, cards, charts, tables,
          htmltools::tags$div(style = "margin-top:40px; color:#888; font-size:11px;",
                              "Generated by germlinevaR \u00b7 gvr_summary()"))

        # ---- write ---------------------------------------------------------------
        # A composite page (HTML tags + several plotly/DT widgets) cannot go through
        # htmlwidgets::saveWidget() (it expects ONE widget; dependency resolution
        # fails on the plain tags). htmltools::save_html() handles arbitrary tag
        # trees and writes the widget JS/CSS into a sidecar lib folder. For a single
        # portable file we then inline that folder with pandoc
        # (rmarkdown::pandoc_self_contained_html). If pandoc is unavailable we keep
        # the sidecar form (<prefix>_report_files/) and flag it. All work happens
        # under tempdir() (the S3-backed mount has no random-access / file locking);
        # the caller shell-copies the finished artifact(s) out.
        widget   <- htmltools::browsable(page)
        base_nm  <- sub("\\.html$", "", basename(final_html))
        files_nm <- paste0(base_nm, "_files")
        work_dir <- file.path(tempdir(), paste0("gvrhtml_", as.integer(Sys.time())))
        if (dir.exists(work_dir)) unlink(work_dir, recursive = TRUE)
        dir.create(work_dir, showWarnings = FALSE, recursive = TRUE)

        staged_html <- file.path(work_dir, basename(final_html))
        # Suppress DT "data too big for client-side DataTables" warning via the
        # DT.warn.size option (checked in DT::datatable's preRenderHook).
        ow_dt <- options(DT.warn.size = FALSE)
        htmltools::save_html(widget, staged_html, libdir = files_nm)
        options(ow_dt)

        # Strip the leading "<!DOCTYPE html>" before pandoc: pandoc treats the input
        # as an HTML fragment and would otherwise escape the doctype into a visible
        # "<p>&lt;!DOCTYPE html&gt;</p>" text node atop the self-contained page.
        # Removing the line lets pandoc's --standalone wrapper emit one clean
        # document. (Sidecar mode keeps the original file, so its doctype is fine.)
        # Also inject a <title> tag (pandoc requires one; without it, a verbose
        # warning is emitted: "This document format requires a nonempty <title>").
        sc_html  <- file.path(work_dir, paste0(base_nm, ".selfcontained.html"))
        sidecar  <- FALSE
        ok_sc <- tryCatch({
          if (isTRUE(.gvr_pandoc_ok())) {
            pre_html <- file.path(work_dir, paste0(base_nm, ".prepandoc.html"))
            raw <- paste(readLines(staged_html, warn = FALSE), collapse = "\n")
            raw <- sub("(?is)^\\s*<!DOCTYPE[^>]*>\\s*", "", raw, perl = TRUE)
            writeLines(raw, pre_html)
            # Bypass rmarkdown::pandoc_self_contained_html(): that helper
            # reads the input as `markdown_strict`, in which pandoc's parser
            # ignores HTML <head> entirely.  As a result any <title> we
            # injected into <head> is never seen, and pandoc emits the
            # "[WARNING] This document format requires a nonempty <title>"
            # message to stderr — which suppressWarnings / suppressMessages /
            # capture.output(type="message") cannot reliably catch because
            # the subprocess writes directly to fd 2.
            #
            # Call pandoc_convert() directly with:
            #   --metadata title=...  : satisfies pandoc's title requirement
            #                            via the metadata layer (read in any
            #                            format, including markdown_strict)
            #   --quiet               : pandoc's own flag for suppressing
            #                            WARNING / INFO messages at the source
            # The minimal template emits only $body$ so the body HTML is
            # passed through unmodified into the self-contained wrapper.
            tpl <- tempfile(fileext = ".html")
            writeLines(c(
              "<!DOCTYPE html>",
              "<html lang=\"en\">",
              "<head>",
              "<meta charset=\"utf-8\" />",
              "<title>$title$</title>",
              "</head>",
              "<body>",
              "$body$",
              "</body>",
              "</html>"
            ), tpl)
            on.exit(unlink(tpl), add = TRUE)
            from_fmt <- if (rmarkdown::pandoc_available("1.17")) "markdown_strict" else "markdown"
            rmarkdown::pandoc_convert(
              input   = pre_html,
              from    = from_fmt,
              output  = sc_html,
              options = c(if (rmarkdown::pandoc_available("2.19"))
                            c("--embed-resources", "--standalone")
                          else "--self-contained",
                          "--template", tpl,
                          "--metadata", "title=germlinevaR Cohort Summary",
                          "--quiet"))
            file.exists(sc_html) && file.info(sc_html)$size > 0
          } else FALSE
        }, error = function(e) FALSE)

        if (isTRUE(ok_sc)) {
          out_path <- sc_html
        } else {
          sidecar  <- TRUE
          out_path <- staged_html
        }
        attr(out_path, "sidecar")   <- sidecar
        attr(out_path, "files_dir") <- if (sidecar) file.path(work_dir, files_nm) else NA_character_
        out_path
      }

      html_name  <- sprintf("%s_report.html", file_prefix)
      final_html <- file.path(out_subdir, html_name)
      files_dest <- file.path(out_subdir, sprintf("%s_report_files", file_prefix))
      if (file.exists(final_html) && isTRUE(verbose))
        message(sprintf("  Overwriting existing HTML report: %s", final_html))

      ok_html <- tryCatch({
        rendered <- .gvr_summary_html(sections, samples, meta, final_html, file_prefix)
        # FUSE-safe: copy the finished .html (and, in sidecar mode, its _files/ folder)
        # from tempdir() to out_subdir via shell cp (R file.copy() can 0-byte on S3).
        system2("cp", c("-f", shQuote(rendered), shQuote(final_html)))
        if (isTRUE(attr(rendered, "sidecar"))) {
          # refresh any stale assets from a previous run, then copy the folder
          if (dir.exists(files_dest)) unlink(files_dest, recursive = TRUE)
          system2("cp", c("-r", shQuote(attr(rendered, "files_dir")), shQuote(files_dest)))
          if (isTRUE(verbose))
            message(sprintf("  HTML report is NOT self-contained; assets in: %s", files_dest))
        }
        file.exists(final_html) && file.info(final_html)$size > 0
      }, error = function(e) {
        warning(sprintf("gvr_summary: HTML report failed: %s", conditionMessage(e))); FALSE })
      if (isTRUE(ok_html) && isTRUE(verbose)) message(sprintf("  HTML report written: %s", final_html))
    }
  }

  invisible(sections)
}

# NOTE: globalVariables() declarations for this package are consolidated in
# R/globals.R (one package-scoped block covering all functions).
