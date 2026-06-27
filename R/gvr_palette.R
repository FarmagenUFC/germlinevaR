# =============================================================================
# gvr_palette.R - 52-palette color system for germlinevaR
# =============================================================================
# Ported 1:1 from the sibling package's plot_multiple_compounds.R
# (generate_colors_from_palette / list_available_palettes).
# Hex codes preserved verbatim so figures stay consistent across the two
# packages.  Exports: gvr_color_palette(), gvr_list_palettes().
# Optional dependencies (RColorBrewer, viridisLite) are loaded only when
# available; otherwise we fall back to hardcoded base colors via
# colorRampPalette().
# =============================================================================


#' Generate a Vector of Colors From a Named Palette
#'
#' `gvr_color_palette()` returns a character vector of `n` hex codes drawn
#' from one of 52 named palettes covering ColorBrewer, viridis, scientific
#' journals, publishers, corporate brands, colorblind-friendly schemes,
#' gradients, and classic R palettes.  This is the helper that
#' [gvr_lollipop()] uses internally to resolve the `variant_palette` and
#' `domain_palette` arguments, but it is also exported so users can build
#' their own color vectors for any plot.
#'
#' @param palette Character string. Name of the palette to use.  See
#'   [gvr_list_palettes()] for the full catalog of 52 names.  Defaults to
#'   `"okabe_ito"`, the colorblind-friendly palette recommended by
#'   Okabe & Ito (2008).
#' @param n Integer.  Number of colors to return.  When `n` exceeds the
#'   base length of a discrete palette, intermediate colors are interpolated
#'   via [grDevices::colorRampPalette()].
#'
#' @return A character vector of length `n` containing hex color codes.
#'
#' @details
#' The 52 palettes are grouped into the following families.  Names are
#' identical to those exposed by the sibling `plot_multiple_compounds()`
#' function so the two packages share the same color vocabulary.
#'
#' \emph{Base:} `"hue"`, `"ggplot2"`, `"default"`  (the ggplot2 default).
#'
#' \emph{ColorBrewer qualitative (good for distinct categories):}
#'   `"set1"`, `"set2"`, `"set3"`, `"dark2"`, `"paired"`, `"accent"`,
#'   `"pastel1"`, `"pastel2"`.
#'
#' \emph{ColorBrewer sequential (ordered data):}
#'   `"blues"`, `"reds"`, `"greens"`, `"purples"`, `"oranges"`, `"greys"`.
#'
#' \emph{ColorBrewer diverging:} `"spectral"`, `"rdylbu"`, `"rdylgn"`,
#'   `"piyg"`, `"prgn"`.
#'
#' \emph{Viridis (perceptually uniform):} `"viridis"`, `"magma"`,
#'   `"inferno"`, `"plasma"`.
#'
#' \emph{Journal palettes:} `"nature"`, `"science"`, `"cell"`, `"plos"`,
#'   `"elife"`.
#'
#' \emph{Publisher palettes:} `"bmc"`, `"frontiers"`, `"wiley"`, `"elsevier"`,
#'   `"oxford"`, `"springer"`, `"acs"`, `"rsc"`.
#'
#' \emph{Corporate brand palettes:} `"ibm"`, `"google"`, `"microsoft"`,
#'   `"twitter"`.
#'
#' \emph{Colorblind-friendly:} `"okabe_ito"`, `"colorblind"`, `"cud"`,
#'   `"tol"`.
#'
#' \emph{Gradient palettes:} `"blue_red"`, `"green_red"`, `"purple_orange"`,
#'   `"cool_warm"`, `"blue_yellow"`.
#'
#' \emph{Classic R:} `"rainbow"`, `"heat"`, `"terrain"`, `"topo"`, `"cm"`.
#'
#' \strong{Optional dependencies:}
#' ColorBrewer palettes use \pkg{RColorBrewer} when available and fall back
#' to interpolated base colors otherwise.  Viridis palettes use
#' \pkg{viridisLite} when available and fall back to manual hex sequences
#' otherwise.  Both packages are in `Suggests:`.
#'
#' @examples
#' # Get 8 Nature-journal colors
#' gvr_color_palette("nature", 8)
#'
#' # Default Okabe-Ito palette for 6 categories
#' gvr_color_palette(n = 6)
#'
#' # List all available palette names
#' gvr_list_palettes()
#'
#' if (requireNamespace("ggplot2", quietly = TRUE)) {
#'     # Use directly in a ggplot2 plot
#'     ggplot2::ggplot(iris,
#'         ggplot2::aes(Sepal.Length, Sepal.Width, color = Species)) +
#'         ggplot2::geom_point() +
#'         ggplot2::scale_color_manual(
#'             values = gvr_color_palette("viridis", 3))
#' }
#' @seealso [gvr_list_palettes()], [gvr_lollipop()].
#' @export
gvr_color_palette <- function(palette = "okabe_ito", n = 8L) {

    if (!is.character(palette) || length(palette) != 1L)
        stop("'palette' must be a single character string")
    n <- as.integer(n)
    if (is.na(n) || n < 1L)
        stop("'n' must be a positive integer")

    # ---------- Scientific-journal palettes ----------
    nature_palettes <- list(
        nature = function(n) {
            base <- c("#0072B2", "#D55E00", "#009E73", "#CC79A7",
                "#F0E442", "#56B4E9")
            if (n <= length(base)) base[seq_len(n)]
            else grDevices::colorRampPalette(base)(n)
        },
        science = function(n) {
            base <- c("#3070B0", "#B03070", "#30B070", "#B07030",
                "#7030B0", "#70B030")
            if (n <= length(base)) base[seq_len(n)]
            else grDevices::colorRampPalette(base)(n)
        },
        cell = function(n) {
            base <- c("#DC143C", "#4682B4", "#2E8B57", "#FF8C00",
                "#9370DB", "#20B2AA")
            if (n <= length(base)) base[seq_len(n)]
            else grDevices::colorRampPalette(base)(n)
        },
        plos = function(n) {
            base <- c("#3498DB", "#E74C3C", "#2ECC71", "#F39C12",
                "#9B59B6", "#1ABC9C")
            if (n <= length(base)) base[seq_len(n)]
            else grDevices::colorRampPalette(base)(n)
        },
        elife = function(n) {
            base <- c("#F04E4E", "#4EA5F0", "#4EF0A5", "#F0A54E",
                "#A54EF0", "#F04EA5")
            if (n <= length(base)) base[seq_len(n)]
            else grDevices::colorRampPalette(base)(n)
        }
    )

    # ---------- Corporate-brand palettes ----------
    corporate_palettes <- list(
        ibm = function(n) {
            base <- c("#648FFF", "#785EF0", "#DC267F", "#FE6100",
                "#FFB000", "#000000")
            if (n <= length(base)) base[seq_len(n)]
            else grDevices::colorRampPalette(base)(n)
        },
        google = function(n) {
            base <- c("#4285F4", "#EA4335", "#FBBC05", "#34A853",
                "#FF6D00", "#46BDC6")
            if (n <= length(base)) base[seq_len(n)]
            else grDevices::colorRampPalette(base)(n)
        },
        microsoft = function(n) {
            base <- c("#F65314", "#7CBB00", "#00A1F1", "#FFBB00",
                "#A0A0A0", "#505050")
            if (n <= length(base)) base[seq_len(n)]
            else grDevices::colorRampPalette(base)(n)
        },
        twitter = function(n) {
            base <- c("#1DA1F2", "#14171A", "#657786", "#AAB8C2",
                "#E1E8ED", "#F5F8FA")
            if (n <= length(base)) base[seq_len(n)]
            else grDevices::colorRampPalette(base)(n)
        }
    )

    # ---------- Colorblind-friendly palettes ----------
    colorblind_palettes <- list(
        okabe_ito = function(n) {
            base <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442",
                "#0072B2", "#D55E00", "#CC79A7", "#999999")
            if (n <= length(base)) base[seq_len(n)]
            else grDevices::colorRampPalette(base)(n)
        },
        colorblind = function(n) {
            base <- c("#0072B2", "#D55E00", "#009E73", "#CC79A7",
                "#F0E442", "#56B4E9")
            if (n <= length(base)) base[seq_len(n)]
            else grDevices::colorRampPalette(base)(n)
        },
        cud = function(n) {
            base <- c("#0072B2", "#E69F00", "#009E73", "#F0E442",
                "#56B4E9", "#D55E00", "#CC79A7")
            if (n <= length(base)) base[seq_len(n)]
            else grDevices::colorRampPalette(base)(n)
        },
        tol = function(n) {
            if (n <= 1L) return("#4477AA")
            if (n == 2L) return(c("#4477AA", "#EE6677"))
            if (n == 3L) return(c("#4477AA", "#EE6677", "#228833"))
            if (n == 4L) return(c("#4477AA", "#EE6677", "#228833", "#CCBB44"))
            if (n == 5L) return(c("#4477AA", "#EE6677", "#228833", "#CCBB44",
                "#66CCEE"))
            if (n == 6L) return(c("#4477AA", "#EE6677", "#228833", "#CCBB44",
                "#66CCEE", "#AA3377"))
            base <- c("#4477AA", "#EE6677", "#228833", "#CCBB44",
                "#66CCEE", "#AA3377")
            grDevices::colorRampPalette(base)(n)
        }
    )

    # ---------- Publisher palettes ----------
    journal_palettes <- list(
        bmc = function(n) {
            base <- c("#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C",
                "#FB9A99", "#E31A1C")
            if (n <= length(base)) base[seq_len(n)]
            else grDevices::colorRampPalette(base)(n)
        },
        frontiers = function(n) {
            base <- c("#4DBBD5", "#E64B35", "#00A087", "#3C5488",
                "#F39B7F", "#8491B4")
            if (n <= length(base)) base[seq_len(n)]
            else grDevices::colorRampPalette(base)(n)
        },
        wiley = function(n) {
            base <- c("#5A9BD5", "#ED7D31", "#A5A5A5", "#FFC000",
                "#4472C4", "#70AD47")
            if (n <= length(base)) base[seq_len(n)]
            else grDevices::colorRampPalette(base)(n)
        },
        elsevier = function(n) {
            base <- c("#F39800", "#DC143C", "#004080", "#009944",
                "#8B4513", "#4B0082")
            if (n <= length(base)) base[seq_len(n)]
            else grDevices::colorRampPalette(base)(n)
        },
        oxford = function(n) {
            base <- c("#002147", "#8B0000", "#006A4E", "#FF6B35",
                "#5E2A84", "#008080")
            if (n <= length(base)) base[seq_len(n)]
            else grDevices::colorRampPalette(base)(n)
        },
        springer = function(n) {
            base <- c("#B22222", "#0066CC", "#228B22", "#FF8C00",
                "#9400D3", "#20B2AA")
            if (n <= length(base)) base[seq_len(n)]
            else grDevices::colorRampPalette(base)(n)
        },
        acs = function(n) {
            base <- c("#0066CC", "#CC0000", "#009966", "#FF9900",
                "#660099", "#FF6600")
            if (n <= length(base)) base[seq_len(n)]
            else grDevices::colorRampPalette(base)(n)
        },
        rsc = function(n) {
            base <- c("#B31B1B", "#005F8C", "#2E8B57", "#FF7F0E",
                "#9467BD", "#17BECF")
            if (n <= length(base)) base[seq_len(n)]
            else grDevices::colorRampPalette(base)(n)
        }
    )

    # ---------- Gradient palettes ----------
    gradient_palettes <- list(
        blue_red = function(n) grDevices::colorRampPalette(
            c("#313695", "#4575B4", "#74ADD1", "#ABD9E9", "#E0F3F8",
                "#FFFFBF", "#FEE090", "#FDAE61", "#F46D43", "#D73027",
                "#A50026"))(n),
        green_red = function(n) grDevices::colorRampPalette(
            c("#1A9850", "#91CF60", "#D9EF8B", "#FEE08B", "#FC8D59",
                "#D73027"))(n),
        purple_orange = function(n) grDevices::colorRampPalette(
            c("#542788", "#998EC3", "#D8DAEB", "#FEE0B6", "#F1A340",
                "#B35806"))(n),
        cool_warm = function(n) grDevices::colorRampPalette(
            c("#3A5F8F", "#7BA0C0", "#C7D8E8", "#F1E3CF", "#EAAA7D",
                "#CC673B"))(n),
        blue_yellow = function(n) grDevices::colorRampPalette(
            c("#053061", "#2166AC", "#4393C3", "#92C5DE", "#D1E5F0",
                "#F7F7F7", "#FDDBC7", "#F4A582", "#D6604D", "#B2182B",
                "#67001F"))(n)
    )

    # ---------- ColorBrewer + Viridis + base palettes ----------
    palettes <- c(
        # Base
        list(
            hue     = function(n) scales::hue_pal()(n),
            ggplot2 = function(n) scales::hue_pal()(n),
            default = function(n) scales::hue_pal()(n)
        ),

        # ColorBrewer qualitative
        list(
            set1 = function(n) {
                if (requireNamespace("RColorBrewer", quietly = TRUE))
                    grDevices::colorRampPalette(RColorBrewer::brewer.pal(9, "Set1"))(n)
                else scales::hue_pal()(n)
            },
            set2 = function(n) {
                if (requireNamespace("RColorBrewer", quietly = TRUE))
                    grDevices::colorRampPalette(RColorBrewer::brewer.pal(8, "Set2"))(n)
                else scales::hue_pal()(n)
            },
            set3 = function(n) {
                if (requireNamespace("RColorBrewer", quietly = TRUE))
                    grDevices::colorRampPalette(RColorBrewer::brewer.pal(12, "Set3"))(n)
                else scales::hue_pal()(n)
            },
            dark2 = function(n) {
                if (requireNamespace("RColorBrewer", quietly = TRUE))
                    grDevices::colorRampPalette(RColorBrewer::brewer.pal(8, "Dark2"))(n)
                else scales::hue_pal()(n)
            },
            paired = function(n) {
                if (requireNamespace("RColorBrewer", quietly = TRUE))
                    grDevices::colorRampPalette(RColorBrewer::brewer.pal(12, "Paired"))(n)
                else scales::hue_pal()(n)
            },
            accent = function(n) {
                if (requireNamespace("RColorBrewer", quietly = TRUE))
                    grDevices::colorRampPalette(RColorBrewer::brewer.pal(8, "Accent"))(n)
                else scales::hue_pal()(n)
            },
            pastel1 = function(n) {
                if (requireNamespace("RColorBrewer", quietly = TRUE))
                    grDevices::colorRampPalette(RColorBrewer::brewer.pal(9, "Pastel1"))(n)
                else scales::hue_pal()(n)
            },
            pastel2 = function(n) {
                if (requireNamespace("RColorBrewer", quietly = TRUE))
                    grDevices::colorRampPalette(RColorBrewer::brewer.pal(8, "Pastel2"))(n)
                else scales::hue_pal()(n)
            }
        ),

        # ColorBrewer sequential
        list(
            blues = function(n) {
                if (requireNamespace("RColorBrewer", quietly = TRUE))
                    grDevices::colorRampPalette(RColorBrewer::brewer.pal(9, "Blues"))(n)
                else grDevices::colorRampPalette(c("#F7FBFF", "#08306B"))(n)
            },
            reds = function(n) {
                if (requireNamespace("RColorBrewer", quietly = TRUE))
                    grDevices::colorRampPalette(RColorBrewer::brewer.pal(9, "Reds"))(n)
                else grDevices::colorRampPalette(c("#FFF5F0", "#67000D"))(n)
            },
            greens = function(n) {
                if (requireNamespace("RColorBrewer", quietly = TRUE))
                    grDevices::colorRampPalette(RColorBrewer::brewer.pal(9, "Greens"))(n)
                else grDevices::colorRampPalette(c("#F7FCF5", "#00441B"))(n)
            },
            purples = function(n) {
                if (requireNamespace("RColorBrewer", quietly = TRUE))
                    grDevices::colorRampPalette(RColorBrewer::brewer.pal(9, "Purples"))(n)
                else grDevices::colorRampPalette(c("#FCFBFD", "#3F007D"))(n)
            },
            oranges = function(n) {
                if (requireNamespace("RColorBrewer", quietly = TRUE))
                    grDevices::colorRampPalette(RColorBrewer::brewer.pal(9, "Oranges"))(n)
                else grDevices::colorRampPalette(c("#FFF5EB", "#7F2704"))(n)
            },
            greys = function(n) {
                if (requireNamespace("RColorBrewer", quietly = TRUE))
                    grDevices::colorRampPalette(RColorBrewer::brewer.pal(9, "Greys"))(n)
                else grDevices::colorRampPalette(c("#FFFFFF", "#000000"))(n)
            }
        ),

        # ColorBrewer diverging
        list(
            spectral = function(n) {
                if (requireNamespace("RColorBrewer", quietly = TRUE))
                    grDevices::colorRampPalette(RColorBrewer::brewer.pal(11, "Spectral"))(n)
                else grDevices::colorRampPalette(
                    c("#D53E4F", "#F46D43", "#FDAE61", "#FEE08B",
                        "#E6F598", "#ABDDA4", "#66C2A5", "#3288BD"))(n)
            },
            rdylbu = function(n) {
                if (requireNamespace("RColorBrewer", quietly = TRUE))
                    grDevices::colorRampPalette(RColorBrewer::brewer.pal(11, "RdYlBu"))(n)
                else grDevices::colorRampPalette(
                    c("#D73027", "#FC8D59", "#FEE090",
                        "#E0F3F8", "#91BFDB", "#4575B4"))(n)
            },
            rdylgn = function(n) {
                if (requireNamespace("RColorBrewer", quietly = TRUE))
                    grDevices::colorRampPalette(RColorBrewer::brewer.pal(11, "RdYlGn"))(n)
                else grDevices::colorRampPalette(
                    c("#D73027", "#FC8D59", "#FEE08B",
                        "#D9EF8B", "#91CF60", "#1A9850"))(n)
            },
            piyg = function(n) {
                if (requireNamespace("RColorBrewer", quietly = TRUE))
                    grDevices::colorRampPalette(RColorBrewer::brewer.pal(11, "PiYG"))(n)
                else grDevices::colorRampPalette(
                    c("#C51B7D", "#E9A3C9", "#FDE0EF",
                        "#E6F5D0", "#A6DBA0", "#008837"))(n)
            },
            prgn = function(n) {
                if (requireNamespace("RColorBrewer", quietly = TRUE))
                    grDevices::colorRampPalette(RColorBrewer::brewer.pal(11, "PRGn"))(n)
                else grDevices::colorRampPalette(
                    c("#762A83", "#9970AB", "#C2A5CF",
                        "#E7D4E8", "#D9F0D3", "#ACD39E", "#5AAE61", "#1B7837"))(n)
            }
        ),

        # Viridis
        list(
            viridis = function(n) {
                if (requireNamespace("viridisLite", quietly = TRUE))
                    viridisLite::viridis(n)
                else grDevices::colorRampPalette(
                    c("#440154", "#3B528B", "#21908C", "#5DC863", "#FDE725"))(n)
            },
            magma = function(n) {
                if (requireNamespace("viridisLite", quietly = TRUE))
                    viridisLite::magma(n)
                else grDevices::colorRampPalette(
                    c("#000004", "#2D1263", "#721F81", "#B63679",
                        "#F8765C", "#FCFDBF"))(n)
            },
            inferno = function(n) {
                if (requireNamespace("viridisLite", quietly = TRUE))
                    viridisLite::inferno(n)
                else grDevices::colorRampPalette(
                    c("#000004", "#1F0C48", "#550F6D", "#A52C60",
                        "#E7683A", "#FCFDBF"))(n)
            },
            plasma = function(n) {
                if (requireNamespace("viridisLite", quietly = TRUE))
                    viridisLite::plasma(n)
                else grDevices::colorRampPalette(
                    c("#0D0887", "#46039F", "#7201A8", "#9C179E", "#BD3786",
                        "#D8576B", "#ED7953", "#FA9E3B", "#FDC926", "#F0F921"))(n)
            }
        ),

        # Classic R
        list(
            rainbow = function(n) grDevices::rainbow(n),
            heat    = function(n) grDevices::heat.colors(n),
            terrain = function(n) grDevices::terrain.colors(n),
            topo    = function(n) grDevices::topo.colors(n),
            cm      = function(n) grDevices::cm.colors(n)
        ),

        # Add journal, corporate, colorblind, publisher, gradient
        nature_palettes,
        corporate_palettes,
        colorblind_palettes,
        journal_palettes,
        gradient_palettes
    )

    if (palette %in% names(palettes)) {
        return(palettes[[palette]](n))
    }

    # Final fallback: try arbitrary RColorBrewer palette name
    if (requireNamespace("RColorBrewer", quietly = TRUE)) {
        if (palette %in% rownames(RColorBrewer::brewer.pal.info)) {
            max_n <- RColorBrewer::brewer.pal.info[palette, "maxcolors"]
            return(grDevices::colorRampPalette(
                RColorBrewer::brewer.pal(max_n, palette))(n))
        }
    }

    warning("Palette '", palette,
        "' not recognized. Using default hue palette.")
    scales::hue_pal()(n)
}


#' List All Available Palettes
#'
#' Returns the sorted character vector of every palette name accepted by
#' [gvr_color_palette()].  Use this for documentation, autocomplete, or
#' programmatic iteration over palettes.
#'
#' @return Character vector of palette names.
#'
#' @details
#' The 56 palette names returned by this function are organized into the
#' following 11 families.  The same taxonomy is used by
#' [gvr_color_palette()] and the sibling `plot_multiple_compounds()` plotting
#' function so the two packages share the same color vocabulary.
#'
#' \emph{Base:} `"hue"`, `"ggplot2"`.
#'
#' \emph{ColorBrewer qualitative (good for distinct categories):}
#'   `"set1"`, `"set2"`, `"set3"`, `"dark2"`, `"paired"`, `"accent"`,
#'   `"pastel1"`, `"pastel2"`.
#'
#' \emph{ColorBrewer sequential (ordered data):}
#'   `"blues"`, `"reds"`, `"greens"`, `"purples"`, `"oranges"`, `"greys"`.
#'
#' \emph{ColorBrewer diverging:} `"spectral"`, `"rdylbu"`, `"rdylgn"`,
#'   `"piyg"`, `"prgn"`.
#'
#' \emph{Viridis (perceptually uniform):} `"viridis"`, `"magma"`,
#'   `"inferno"`, `"plasma"`.
#'
#' \emph{Journal palettes:} `"nature"`, `"science"`, `"cell"`, `"plos"`,
#'   `"elife"`.
#'
#' \emph{Publisher palettes:} `"bmc"`, `"frontiers"`, `"wiley"`, `"elsevier"`,
#'   `"oxford"`, `"springer"`, `"acs"`, `"rsc"`.
#'
#' \emph{Corporate brand palettes:} `"ibm"`, `"google"`, `"microsoft"`,
#'   `"twitter"`.
#'
#' \emph{Colorblind-friendly:} `"okabe_ito"`, `"colorblind"`, `"cud"`,
#'   `"tol"`.
#'
#' \emph{Gradient palettes:} `"blue_red"`, `"green_red"`, `"purple_orange"`,
#'   `"cool_warm"`, `"blue_yellow"`.
#'
#' \emph{Classic R:} `"rainbow"`, `"heat"`, `"terrain"`, `"topo"`, `"cm"`.
#'
#' See [gvr_color_palette()] for example hex codes from each family.
#'
#' @examples
#' gvr_list_palettes()
#' # Generate a 5-color sample of every palette
#' lapply(gvr_list_palettes(), gvr_color_palette, n = 5)
#' @seealso [gvr_color_palette()], [gvr_lollipop()].
#' @export
gvr_list_palettes <- function() {
    c(
        # Base
        "hue", "ggplot2",
        # ColorBrewer qualitative
        "set1", "set2", "set3", "dark2", "paired", "accent",
        "pastel1", "pastel2",
        # ColorBrewer sequential
        "blues", "reds", "greens", "purples", "oranges", "greys",
        # ColorBrewer diverging
        "spectral", "rdylbu", "rdylgn", "piyg", "prgn",
        # Viridis
        "viridis", "magma", "inferno", "plasma",
        # Journal
        "nature", "science", "cell", "plos", "elife",
        # Publisher
        "bmc", "frontiers", "wiley", "elsevier", "oxford",
        "springer", "acs", "rsc",
        # Corporate
        "ibm", "google", "microsoft", "twitter",
        # Colorblind-friendly
        "okabe_ito", "colorblind", "cud", "tol",
        # Gradient
        "blue_red", "green_red", "purple_orange", "cool_warm", "blue_yellow",
        # Classic R
        "rainbow", "heat", "terrain", "topo", "cm"
    )
}
