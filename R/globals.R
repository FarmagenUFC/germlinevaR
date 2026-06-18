# germlinevaR: package-level global-variable declarations
#
# data.table's non-standard evaluation and ggplot2's aes() reference column
# names as bare symbols, which `R CMD check` flags as "no visible binding for
# global variable". Declaring them here (package-scoped, registered at build
# time) silences those NOTEs for every function in the package. This is the
# idiomatic CRAN placement: one dedicated R/globals.R rather than per-file blocks.

#' @importFrom utils globalVariables
#' @keywords internal
utils::globalVariables(c(
  # --- data.table special symbols (used inside dt[...] across functions) ----
  ".", ".N", ".SD", ".__sample__", ".sr", ".__nsamp__",
  # --- data.table column-vector prefix (used as `dt[, ..cols]` in gvr_summary)
  "..cols", "..final_cols",

  # --- read.gvr(): temp/internal columns created via := ---------------------
  ".rs", ".ref", ".alt",
  # --- read.gvr(): MAF columns referenced bare inside data.table calls ------
  "Genotype", "Tumor_Seq_Allele1", "Tumor_Seq_Allele2", "ABraOM_AF",
  "dbSNP_RS", "Reference_Allele",
  # --- read.gvr(): ABraOM lookup-table columns (built then keyed/joined) ----
  "avsnp147", "Ref", "Alt", "Frequencies", "rs", "ref", "alt", "af", "x.af",

  # --- gvr_summary(): section columns + nested-renderer aes() symbols -------
  "Hugo_Symbol", "Total", "Variant_Classification", "CLIN_SIG",
  "Category", "Sample", "n",
  # --- gvr_summary(): top-rsID aggregation + drill-down detail columns ------
  ".__rs__", ".__n__", "Chromosome", "Start_Position", "IMPACT",

  # --- gvr_plot(): matrix/summary columns ------------------------------
  "n_var", "n_samp", "N",

  # --- gvr_lollipop(): per-dot/stem aes() symbols + internal data.table cols ----
  "HGVSp_Short", "Tumor_Sample_Barcode",
  ".__aa_pos__", ".__hgvsp__", ".__vc__", ".__y__", ".n",
  "aa_pos", "y", "top", "vc", "hgvsp", "label",
  # --- gvr_lollipop() Phase G: geom_rect / geom_text aes() symbols ----------
  # bar_df (xmin/xmax/ymin/ymax), domain rects (also xmid + name), stems (xend)
  "xmin", "xmax", "ymin", "ymax", "xmid", "name", "xend",
  # --- gvr_lollipop() Phase H+: domain-label NSE columns added by aes() ----
  ".lbl", ".legend", "fill", ".lbl_resolved", ".font_mult",
  # --- gvr_lollipop() Phase I: hotspot_df aes() / column symbols -----------
  # n_in_window: hotspot row count (optional label/tooltip in future versions)
  "n_in_window",
  # --- gvr_summary() Phase N+9 Stage D: per-token rank-key precompute -------
  # .__ir__/.__gaf__/.__chr__/.__pos__ are columns built on .xl_proj via := ;
  # gnomADe_AF is referenced bare inside the := RHS for as.numeric(gnomADe_AF).
  ".__ir__", ".__gaf__", ".__chr__", ".__pos__", "gnomADe_AF",
  # --- gvr_genepos.plot(): exon/intron/UTR + parser column names ---------
  # Layout (.gvr_genepos_layout): exon_idx, region_kind, g_start, g_end,
  # length_bp, x_start, x_end. Parser (.gvr_parse_hgvsc): kind, pos, valid.
  # MAF subset internal cols added via := during plot pipeline.
  "exon_idx", "region_kind", "g_start", "g_end", "length_bp",
  "x_start", "x_end", "kind", "pos", "valid",
  ".__enst__", ".__kind__", ".__valid__", ".__x__", ".__top__",
  # --- gvr_genepos.plot(): VEP-annotated MAF cols used in transcript pick
  # MANE_SELECT/CANONICAL/Transcript_ID are MAF columns from read.gvr(),
  # referenced bare in dt[...] filters during transcript auto-resolution.
  "MANE_SELECT", "CANONICAL", "Transcript_ID",
  # ggplot aes() reuses 'x' as a column name in dot_df/label_df/stem_df.
  "x",
  # GVR_CLASS_COLORS is a package-level named-vector constant defined in
  # gvr_lollipop.R; declaring it here keeps R CMD check quiet when other
  # files reference it via the same-file scoping rule.
  "GVR_CLASS_COLORS",

  # --- read.gvr.dual(): VEP CSQ extension + SnpEff-derived columns ---------
  # FREQS is the 81st canonical CSQ field (added by VEP --everything in v113+):
  # names the gnomAD population where MAX_AF was observed (e.g. "gnomADg_ASJ").
  # The four LOF_*/NMD_* columns come from SnpEff's INFO LOF= and NMD= keys
  # (loss-of-function and nonsense-mediated-decay predictions). The four
  # snpeff_* columns hold SnpEff's parallel pick (consequence/impact/gene/hgvsc)
  # for side-by-side comparison with VEP's pick in the same row.
  "FREQS",
  "LOF_Gene", "LOF_Pct_Transcripts", "NMD_Gene", "NMD_Pct_Transcripts",
  "snpeff_consequence", "snpeff_impact", "snpeff_gene", "snpeff_hgvsc"
))
