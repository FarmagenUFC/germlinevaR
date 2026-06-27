#' Cohort top-genes variant matrix from a germline gvr table (read.gvr / gvr_filter output)
#'
#' @description
#' Draws a top-genes variant matrix from an MAF-like table - the output of
#' [read.gvr()], or of [gvr_filter()] - and writes it to a PNG file. Rows are the
#' top-`top_n` genes (ranked by number of distinct samples mutated, then by variant
#' count); columns are samples. Each gene x sample cell shows the single MOST-SEVERE
#' `Variant_Classification` observed for that gene in that sample. The plot is drawn
#' with \pkg{ComplexHeatmap}.
#'
#' @details
#' This is the standalone top-genes variant matrix previously produced inside [gvr_summary()]. It
#' needs only the `Hugo_Symbol` and `Variant_Classification` columns (plus the
#' per-sample column).
#'
#' Cell collapse: when a gene has several variant classes in one sample, the cell is
#' coloured by the most severe class, using this ranking (high to low):
#' Translation_Start_Site, Nonsense_Mutation, Nonstop_Mutation, Splice_Site,
#' Frame_Shift_Del, Frame_Shift_Ins, In_Frame_Del, In_Frame_Ins, Missense_Mutation,
#' Splice_Region, Protein_altering_variant, Silent, 5'UTR, 3'UTR, 5'Flank, 3'Flank,
#' RNA, Intron, IGR, Targeted_Region. Any class outside this list ranks last and is
#' coloured grey. Colours follow a colourblind-safe (Okabe-Ito) palette.
#'
#' Annotations: a right-side bar shows each gene's total variant burden; a top bar
#' shows each sample's total variant burden (axis labelled in thousands). Empty cells
#' (gene not mutated in that sample) are light grey.
#'
#' Data conventions:
#' \itemize{
#'   \item "Missing" means `NA` OR empty string `""`.
#'   \item Unknown/blank gene symbols are `Hugo_Symbol` in `c(".", "", "Unknown")`;
#'     these are excluded from the top-genes variant matrix.
#'   \item Works on ANY MAF-like table; it makes no assumption about prior
#'     filtering. It is commonly run on [gvr_filter()] output.
#' }
#'
#' @param gvr An MAF-like `data.table`/`data.frame` from [read.gvr()] or [gvr_filter()].
#'   Required columns: `Hugo_Symbol`, `Variant_Classification`.
#' @param top_n Integer; number of genes (rows) shown, ranked by number of distinct
#'   samples mutated then by variant count. Default `20`.
#' @param sample_col Name of the per-sample column. Default `"Tumor_Sample_Barcode"`.
#'   If absent, all rows are pooled into a single sample `"All"` (with a warning).
#' @param out_dir Output directory for the PNG. Created if it does not exist.
#'   Default `"."` (current working directory).
#' @param file_prefix Base filename for the written PNG. Default `"gvr_plot"`;
#'   the file is written as `<file_prefix>.png` (fixed name, no timestamp), e.g.
#'   `gvr_plot.png`. An existing file at that path is overwritten (a message is
#'   emitted when `verbose = TRUE`).
#' @param sample_name_rot Numeric; rotation angle (degrees) for the sample-name
#'   labels at the top of the heatmap. Default `45`. Common alternatives are
#'   `0` (horizontal) and `90` (vertical). Must be a single finite numeric.
#' @param gene_name_size Numeric; font size in points for the gene-name row
#'   labels (left side of the heatmap). Default `10`. Must be a single
#'   finite positive numeric.
#' @param sample_name_size Numeric; font size in points for the sample-name
#'   column labels (top of the heatmap). Default `10`. Must be a single
#'   finite positive numeric.
#' @param axis_tick_size Numeric; font size in points for the y-axis tick
#'   labels of the top "Variant impact" bar and the right gene-burden bar.
#'   Default `7`. Must be a single finite positive numeric.
#' @param legend_label_size Numeric; font size in points for the label text
#'   of both side legends ("Impact" and "Most severe class"). Default `10`.
#'   Must be a single finite positive numeric.
#' @param legend_title_size Numeric; font size in points for the title text
#'   of both side legends ("Impact" and "Most severe class"). The bold
#'   face is preserved. Default `10`. Must be a single finite positive
#'   numeric.
#' @param legend_label_wrap_chars Numeric; if a (prettified) legend label is
#'   longer than this many characters, it is wrapped onto two lines at the
#'   space closest to the middle. Default `Inf` disables wrapping (so all
#'   labels stay on one line, matching the legacy behaviour). Use a finite
#'   integer (e.g. `14`) when long labels would clip the right edge at
#'   large `legend_label_size`. Labels without internal spaces (e.g.
#'   `5'UTR`, `RNA`) are never wrapped regardless of length.
#' @param impact_title_side One of `"left"` (default) or `"right"`; controls
#'   where the `"Variant impact"` annotation title is drawn relative to the
#'   top stacked-bar panel. `"left"` renders the title vertically (acting as a
#'   y-axis title for the impact panel); `"right"` renders it horizontally
#'   on the right of the panel (the previous default).
#' @param verbose Logical; if `TRUE` (default) print the path of the file written.
#'
#' @return Invisibly, the path of the written PNG (character), or `NA_character_` if
#'   no known-gene variants are present in the table.
#'
#' @section Dependencies:
#' Uses \pkg{ComplexHeatmap} (a \pkg{Bioconductor} package, listed in `Imports`).
#'
#' @seealso [gvr_summary()] for the tabular summary, [read.gvr()] to build the table,
#'   [gvr_filter()] to filter it before plotting.
#' @family germlinevaR
#' @author germlinevaR authors
#'
#' @examples
#' ## Load the shipped example table; write plot to a temp directory
#' gvr <- readRDS(system.file("extdata", "example_gvr.rds",
#'     package = "germlinevaR"))
#' p <- gvr_plot(gvr, out_dir = tempdir(), verbose = FALSE)
#' class(p)
#'
#' ## Smaller top-genes variant matrix of filtered hits to a temp folder
#' filt <- gvr_filter(gvr, ABraOM_AF = NULL, verbose = FALSE)
#' if (nrow(filt) > 0L) {
#'     gvr_plot(filt, top_n = 15, out_dir = tempdir(), verbose = FALSE)
#' }
#' @importFrom ComplexHeatmap rowAnnotation anno_barplot HeatmapAnnotation Heatmap Legend packLegend draw
#' @importFrom data.table as.data.table data.table setorder uniqueN :=
#' @importFrom grDevices png dev.off
#' @importFrom grid gpar grid.rect unit
#' @importFrom utils head
#' @export
gvr_plot <- function(gvr,
                     top_n             = 20,
                     sample_col        = "Tumor_Sample_Barcode",
                     out_dir           = ".",
                     file_prefix       = "gvr_plot",
                     sample_name_rot   = 45,
                     gene_name_size    = 10,
                     sample_name_size  = 10,
                     axis_tick_size    = 7,
                     legend_label_size = 10,
                     legend_title_size = 10,
                     legend_label_wrap_chars = Inf,
                     impact_title_side = c("left", "right"),
                     verbose           = TRUE) {

    if (!requireNamespace("data.table", quietly = TRUE)) {
        stop("gvr_plot requires the 'data.table' package.")
    }
    dt <- data.table::as.data.table(gvr)

    # --- Validate new layout args -----------------------------------------------
    impact_title_side <- match.arg(impact_title_side)
    if (!is.numeric(sample_name_rot) || length(sample_name_rot) != 1L ||
        !is.finite(sample_name_rot))
        stop("gvr_plot: 'sample_name_rot' must be a single finite numeric (degrees).",
            call. = FALSE)
    .check_pos_num <- function(x, nm) {
        if (!is.numeric(x) || length(x) != 1L || !is.finite(x) || x <= 0)
            stop(sprintf("gvr_plot: '%s' must be a single finite positive numeric (points).", nm),
                call. = FALSE)
    }
    .check_pos_num(gene_name_size,   "gene_name_size")
    .check_pos_num(sample_name_size, "sample_name_size")
    .check_pos_num(axis_tick_size,   "axis_tick_size")
    .check_pos_num(legend_label_size, "legend_label_size")
    .check_pos_num(legend_title_size, "legend_title_size")
    # legend_label_wrap_chars accepts Inf (the default no-wrap sentinel) as well
    # as a finite positive numeric. NA, zero, negative, character, or non-scalar
    # values are rejected with a clear message.
    if (!is.numeric(legend_label_wrap_chars) ||
        length(legend_label_wrap_chars) != 1L ||
        is.na(legend_label_wrap_chars) ||
        legend_label_wrap_chars <= 0)
        stop("gvr_plot: 'legend_label_wrap_chars' must be a single positive numeric (Inf is allowed and disables wrapping).",
            call. = FALSE)

    # --- Soft guard for IMPACT column (used by top annotation) ----------------
    has_impact <- "IMPACT" %in% names(dt)
    if (!has_impact)
        warning("gvr_plot: 'IMPACT' column not found; falling back to total-burden bar.")

    .is_missing <- function(v) is.na(v) | v == ""
    UNKNOWN_GENE <- c(".", "", "Unknown")

    # --- Variant_Classification severity order (high -> low). ----
    #     Used to collapse multi-class gene x sample cells to a single most-severe
    #     class. Any class not listed sorts LAST and renders grey.
    GVR_SEVERITY <- c("Translation_Start_Site", "Nonsense_Mutation", "Nonstop_Mutation",
        "Splice_Site", "Frame_Shift_Del", "Frame_Shift_Ins",
        "In_Frame_Del", "In_Frame_Ins", "Missense_Mutation",
        "Splice_Region", "Protein_altering_variant", "Silent",
        "5'UTR", "3'UTR", "5'Flank", "3'Flank", "RNA", "Intron",
        "IGR", "Targeted_Region")
    .sev_rank <- function(x) {
        r <- match(x, GVR_SEVERITY)
        r[is.na(r)] <- length(GVR_SEVERITY) + 1L
        r
    }
    # Display-only label prettifier: 'Missense_Mutation' -> 'Missense Mutation',
    # 'HIGH' -> 'High', etc. Preserves region tokens (5'UTR, 3'Flank) and
    # short uppercase acronyms (RNA, IGR). Used ONLY for legend labels; data
    # matching and palette keys remain in the original case throughout.
    .pretty_label <- function(x) {
        is_acronym <- function(tok)
            grepl("^[A-Z]{2,5}$", tok) &&
                !tok %in% c("HIGH", "MODERATE", "LOW", "MODIFIER")
        cap_first <- function(tok)
            if (nchar(tok) == 0L) tok
            else paste0(toupper(substr(tok, 1, 1)), tolower(substr(tok, 2, nchar(tok))))
        handle_token <- function(tok) {
            m <- regmatches(tok, regexec("^([0-9]+')(.+)$", tok))[[1]]
            if (length(m) == 3L) {
                rest <- m[3]
                if (is_acronym(rest)) return(paste0(m[2], rest))
                return(paste0(m[2], cap_first(rest)))
            }
            if (is_acronym(tok)) return(tok)
            cap_first(tok)
        }
        vapply(x, function(s) {
            parts <- strsplit(s, "_", fixed = TRUE)[[1]]
            paste(vapply(parts, handle_token, character(1), USE.NAMES = FALSE),
                collapse = " ")
        }, character(1), USE.NAMES = FALSE)
    }
    # Label-wrap helper: if a (display) label is longer than max_chars, split it
    # onto two lines at the space whose position is closest to nchar(s)/2. Labels
    # without spaces (acronyms like 5'UTR, RNA, IGR) are returned unchanged. If
    # max_chars is Inf (default) every label passes through untouched.
    .wrap_balanced <- function(x, max_chars) {
        if (!is.finite(max_chars)) return(x)
        vapply(x, function(s) {
            if (nchar(s) <= max_chars) return(s)
            sp <- gregexpr(" ", s, fixed = TRUE)[[1L]]
            if (length(sp) == 1L && sp[1L] == -1L) return(s)
            mid <- nchar(s) / 2
            pick <- sp[which.min(abs(sp - mid))]
            paste0(substr(s, 1L, pick - 1L), "\n", substr(s, pick + 1L, nchar(s)))
        }, character(1L), USE.NAMES = FALSE)
    }
    # Colorblind-safe palette (Okabe-Ito + extensions), keyed by class; "Other" = grey.
    GVR_CLASS_COLORS <- c(
        "Translation_Start_Site" = "#000000", "Nonsense_Mutation" = "#D55E00",
        "Nonstop_Mutation" = "#882255", "Splice_Site" = "#CC79A7",
        "Frame_Shift_Del" = "#E69F00", "Frame_Shift_Ins" = "#F0E442",
        "In_Frame_Del" = "#56B4E9", "In_Frame_Ins" = "#0072B2",
        "Missense_Mutation" = "#009E73", "Splice_Region" = "#44AA99",
        "Protein_altering_variant" = "#117733", "Silent" = "#999933",
        "5'UTR" = "#AA4499", "3'UTR" = "#DDCC77", "5'Flank" = "#88CCEE",
        "3'Flank" = "#332288", "RNA" = "#BBBBBB", "Intron" = "#DDDDDD",
        "IGR" = "#777777", "Targeted_Region" = "#666666", "Other" = "#CCCCCC")

    # --- Resolve sample column ---------------------------------------------------
    if (!sample_col %in% names(dt)) {
        warning(sprintf("gvr_plot: sample column '%s' not found; pooling all rows into 'All'.",
            sample_col))
        dt[, .__sample__ := "All"]
    } else {
        dt[, .__sample__ := as.character(get(sample_col))]
        dt[.is_missing(.__sample__), .__sample__ := "NA_sample"]
    }
    samples <- sort(unique(dt$.__sample__))

    # --- Column-existence guard --------------------------------------------------
    req <- c("Hugo_Symbol", "Variant_Classification")
    miss_req <- req[!req %in% names(dt)]
    if (length(miss_req) > 0) {
        stop(sprintf("gvr_plot: required column(s) not found: %s",
            paste(miss_req, collapse = ", ")))
    }

    # --- FUSE-safe PNG writer: render via `draw_fun(tmp)` to a tempdir path, then ---
    #     shell-cp to the final out_dir path (S3-backed mounts can 0-byte direct
    #     random-access writes). Returns the final path on success, NA on failure.
    .fuse_save_png <- function(final_path, draw_fun) {
        tmp <- file.path(tempdir(), basename(final_path))
        ok <- tryCatch(
            {
                draw_fun(tmp)
                file.exists(tmp) && file.info(tmp)$size > 0
            },
            error = function(e) {
                warning(sprintf("gvr_plot: plot render failed (%s): %s",
                    basename(final_path), conditionMessage(e)))
                FALSE
            })
        if (!ok) return(NA_character_)
        system2("cp", c("-f", shQuote(tmp), shQuote(final_path)))
        if (!file.exists(final_path) || file.info(final_path)$size == 0) {
            warning(sprintf("gvr_plot: copy to '%s' failed; left at '%s'.", final_path, tmp))
            return(tmp)
        }
        final_path
    }

    # --- Top-genes variant matrix builder (ComplexHeatmap): top-N genes x samples, each cell the ---
    #     single MOST-SEVERE Variant_Classification. Returns final path or NA.
    # ComplexHeatmap is a hard dependency (declared in DESCRIPTION Imports), so no
    # runtime presence check is needed; users always have it installed via Bioconductor.
    m <- dt[!(Hugo_Symbol %in% UNKNOWN_GENE)]
    if (nrow(m) == 0L) {
        warning("gvr_plot: no known-gene variants; skipping plot.")
        return(invisible(NA_character_))
    }
    gstat <- m[, .(n_var = .N, n_samp = data.table::uniqueN(.__sample__)), by = Hugo_Symbol]
    data.table::setorder(gstat, -n_samp, -n_var)
    top_g <- utils::head(gstat$Hugo_Symbol, top_n)
    sub <- m[Hugo_Symbol %in% top_g]
    sub[, .sr := .sev_rank(Variant_Classification)]
    cell <- sub[, .SD[which.min(.sr)], by = .(Hugo_Symbol, .__sample__)]
    mat <- matrix("", nrow = length(top_g), ncol = length(samples),
        dimnames = list(top_g, samples))
    for (i in seq_len(nrow(cell)))
        mat[cell$Hugo_Symbol[i], cell$.__sample__[i]] <- cell$Variant_Classification[i]
    # gene order: most-mutated at top (already ordered in top_g)
    mat <- mat[rev(top_g), , drop = FALSE]
    classes_present <- setdiff(unique(as.vector(mat)), "")
    col_map <- GVR_CLASS_COLORS[classes_present]
    col_map[is.na(col_map)] <- GVR_CLASS_COLORS[["Other"]]
    names(col_map) <- classes_present
    # per-gene total variant burden (right annotation) and per-sample burden (top)
    gene_burden <- m[Hugo_Symbol %in% top_g, .N, by = Hugo_Symbol]
    gb <- gene_burden$N[match(rownames(mat), gene_burden$Hugo_Symbol)]
    gb[is.na(gb)] <- 0
    samp_burden <- vapply(samples, function(s) sum(m$.__sample__ == s), integer(1))

    # IMPACT palette and per-sample counts matrix (rows = levels, cols = samples).
    # anno_barplot stacks rows bottom-to-top, so to render HIGH on top the matrix
    # row order must be MODIFIER (bottom) -> LOW -> MODERATE -> HIGH (top).
    IMPACT_LEVELS  <- c("HIGH", "MODERATE", "LOW", "MODIFIER")
    IMPACT_STACK   <- rev(IMPACT_LEVELS)   # bottom-to-top stacking order
    IMPACT_COLORS  <- c(HIGH = "#D55E00", MODERATE = "#E69F00",
        LOW  = "#009E73", MODIFIER  = "#BBBBBB")
    if (has_impact) {
        imp_mat <- do.call(rbind, lapply(IMPACT_STACK, function(lv) {
            vapply(samples, function(s)
                sum(m$.__sample__ == s & !is.na(m$IMPACT) & m$IMPACT == lv),
            integer(1))
        }))
        rownames(imp_mat) <- IMPACT_STACK
    }

    # Right annotation: per-gene total variant burden bar
    ra <- ComplexHeatmap::rowAnnotation(
        `Variants` = ComplexHeatmap::anno_barplot(
            gb, border = FALSE,
            gp = grid::gpar(fill = "#0279EE", col = NA),
            axis_param = list(gp = grid::gpar(fontsize = axis_tick_size))),
        width = grid::unit(1.8, "cm"),
        annotation_name_gp = grid::gpar(fontsize = 9))

    # Map impact_title_side -> annotation_name_side/_rot. Left = vertical title
    # acting as a y-axis label; right = horizontal (the previous default).
    .impact_name_side <- impact_title_side
    .impact_name_rot  <- if (impact_title_side == "left") 90 else 0

    # Top annotation: stacked IMPACT bar (or fallback total-burden bar).
    # Compute pretty y-axis ticks once and reuse for at, labels, and ylim so the
    # axis line extends all the way to the topmost tick (not just to the data max).
    ta <- if (has_impact) {
        .imp_at <- pretty(c(0, colSums(imp_mat)), n = 3)
        .imp_labels <- if (max(.imp_at) >= 1000)
            paste0(round(.imp_at / 1000), "k")
        else
            as.character(as.integer(.imp_at))
        ComplexHeatmap::HeatmapAnnotation(
            `Variant\nimpact` = ComplexHeatmap::anno_barplot(
                t(imp_mat), border = FALSE, beside = FALSE,
                gp = grid::gpar(fill = IMPACT_COLORS[IMPACT_STACK], col = NA),
                ylim = c(0, max(.imp_at)),
                axis_param = list(
                    at     = .imp_at,
                    labels = .imp_labels,
                    gp     = grid::gpar(fontsize = axis_tick_size))),
            height = grid::unit(1.6, "cm"),
            annotation_name_gp   = grid::gpar(fontsize = 9),
            annotation_name_side = .impact_name_side,
            annotation_name_rot  = .impact_name_rot)
    } else {
        # fallback: original total-burden bar
        .bur_at <- pretty(c(0, samp_burden), n = 3)
        .bur_labels <- if (max(.bur_at) >= 1000)
            paste0(round(.bur_at / 1000), "k")
        else
            as.character(as.integer(.bur_at))
        ComplexHeatmap::HeatmapAnnotation(
            `Burden` = ComplexHeatmap::anno_barplot(
                samp_burden, border = FALSE,
                gp = grid::gpar(fill = "#0279EE", col = NA),
                ylim = c(0, max(.bur_at)),
                axis_param = list(
                    at     = .bur_at,
                    labels = .bur_labels,
                    gp     = grid::gpar(fontsize = axis_tick_size))),
            height = grid::unit(1.6, "cm"),
            annotation_name_gp = grid::gpar(fontsize = 9))
    }
    cell_fun <- function(j, i, x, y, width, height, fill) {
        v <- mat[i, j]
        grid::grid.rect(x, y, width = width * 0.95, height = height * 0.95,
            gp = grid::gpar(fill = if (v == "") "#F2F2F2" else col_map[[v]], col = "white", lwd = 1))
    }
    ht <- ComplexHeatmap::Heatmap(
        mat, name = "Most severe\nclass", col = col_map, rect_gp = grid::gpar(type = "none"),
        cell_fun = cell_fun, na_col = "#F2F2F2",
        cluster_rows = FALSE, cluster_columns = FALSE,
        show_heatmap_legend = FALSE, row_names_side = "left", column_names_side = "top",
        column_names_rot = sample_name_rot, column_names_centered = TRUE,
        row_names_gp = grid::gpar(fontsize = gene_name_size),
        column_names_gp = grid::gpar(fontsize = sample_name_size),
        right_annotation = ra, top_annotation = ta,
        column_title = sprintf("Top %d genes \u00d7 %d sample(s) \u2014 cells show most-severe class",
            length(top_g), length(samples)),
        column_title_gp = grid::gpar(fontsize = 11))

    if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    # Fixed output name (no timestamp), mirroring gvr_summary; re-runs overwrite in
    # place and announce the overwrite when verbose.
    final_path <- file.path(out_dir, sprintf("%s.png", file_prefix))
    if (file.exists(final_path) && isTRUE(verbose))
        message(sprintf("gvr_plot: overwriting existing %s", final_path))
    # Swatch size: when label wrap is OFF (Inf) or no label actually wraps,
    # use the ComplexHeatmap default (4mm) so output stays bit-identical to
    # the unwrapped baseline. When at least one label wraps to two lines, pad
    # ALL swatches up to the wrapped-row height so the colored squares remain
    # uniformly sized across the legend (otherwise ComplexHeatmap only inflates
    # the wrapped row's swatch, making them visually inconsistent). The pad
    # height was measured empirically on ComplexHeatmap 2.22 as a tight linear
    # function of fontsize: h_mm(fs) ~= 0.6343 * fs + 1.41 for a 2-line label.
    .all_labels   <- c(if (has_impact) .pretty_label(IMPACT_LEVELS),
        .pretty_label(names(col_map)))
    .wrap_labels  <- .wrap_balanced(.all_labels, legend_label_wrap_chars)
    .any_wrapped  <- any(grepl("\n", .wrap_labels, fixed = TRUE))
    .swatch_side_mm <- if (.any_wrapped) {
        ceiling((0.6343 * legend_label_size + 1.41) * 10) / 10
    } else 4

    # IMPACT legend (anno_barplot does not auto-generate one).
    impact_lgd <- if (has_impact)
        ComplexHeatmap::Legend(
            labels      = .wrap_balanced(.pretty_label(IMPACT_LEVELS), legend_label_wrap_chars),
            title       = "Impact",
            legend_gp   = grid::gpar(fill = IMPACT_COLORS[IMPACT_LEVELS]),
            labels_gp   = grid::gpar(fontsize = legend_label_size),
            title_gp    = grid::gpar(fontsize = legend_title_size, fontface = "bold"),
            grid_height = grid::unit(.swatch_side_mm, "mm"),
            grid_width  = grid::unit(.swatch_side_mm, "mm"))
    else NULL

    # "Most severe class" legend, built manually from the heatmap palette. The
    # heatmap's own auto-legend is suppressed (show_heatmap_legend = FALSE) so
    # we can place IMPACT on top of "Most severe class" via heatmap_legend_list
    # ordering in draw() below.
    class_lgd <- ComplexHeatmap::Legend(
        labels      = .wrap_balanced(.pretty_label(names(col_map)), legend_label_wrap_chars),
        title       = "Most severe\nclass",
        legend_gp   = grid::gpar(fill = unname(col_map)),
        labels_gp   = grid::gpar(fontsize = legend_label_size),
        title_gp    = grid::gpar(fontsize = legend_title_size, fontface = "bold"),
        grid_height = grid::unit(.swatch_side_mm, "mm"),
        grid_width  = grid::unit(.swatch_side_mm, "mm"))

    path <- .fuse_save_png(final_path, function(tmp) {
        grDevices::png(tmp, width = max(1100, 360 + 150 * length(samples)),
            height = max(720, 110 + 34 * length(top_g)), res = 150)
        # Stack IMPACT (if present) on top of "Most severe class". The whole legend
        # column is positioned manually so it begins at the TOP of the column
        # annotation (i.e., next to the impact barplot), not at the heatmap body
        # top. ComplexHeatmap's align_heatmap_legend = "heatmap_top" anchors to
        # the body top, which leaves IMPACT visually next to the heatmap rather
        # than next to the impact barplot it describes; the manual overlay fixes
        # that.
        .pl <- if (!is.null(impact_lgd))
            ComplexHeatmap::packLegend(impact_lgd, class_lgd, direction = "vertical",
                row_gap = grid::unit(4, "mm"))
        else class_lgd
        ComplexHeatmap::draw(ht, show_heatmap_legend = FALSE,
            padding = grid::unit(c(2, 6, 2, 45), "mm"))
        # Compute device-absolute coordinates of the top column annotation and the
        # right annotation, then push a viewport for the legend column starting at
        # the top of the column annotation.
        .anno_name <- if (has_impact) "annotation_Variant\nimpact_1"
        else "annotation_Burden_1"
        grid::seekViewport(.anno_name)
        .anno_top <- grid::deviceLoc(grid::unit(0, "npc"), grid::unit(1, "npc"),
            valueOnly = TRUE)
        grid::upViewport(0)
        grid::seekViewport("annotation_Variants_1")
        .ra_right <- grid::deviceLoc(grid::unit(1, "npc"), grid::unit(0, "npc"),
            valueOnly = TRUE)
        grid::upViewport(0)
        .dev_w <- grDevices::dev.size("in")[1]
        .lgd_x <- .ra_right$x + 0.25
        .lgd_w <- .dev_w - .lgd_x - 0.1
        grid::pushViewport(grid::viewport(
            x = grid::unit(.lgd_x, "in"),
            y = grid::unit(.anno_top$y, "in"),
            width = grid::unit(.lgd_w, "in"),
            height = grid::unit(.anno_top$y - 0.1, "in"),
            just = c("left", "top")))
        ComplexHeatmap::draw(.pl, x = grid::unit(0, "npc"), y = grid::unit(1, "npc"),
            just = c("left", "top"))
        grid::upViewport()
        grDevices::dev.off()
    })
    if (!is.na(path) && isTRUE(verbose)) message(sprintf("gvr_plot: written %s", path))
    invisible(path)
}

# NOTE: globalVariables() declarations for this package are consolidated in
# R/globals.R (one package-scoped block covering all functions).
