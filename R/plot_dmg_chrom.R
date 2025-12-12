#' @title Plot differentially methylated regions (DMRs) hyper/hypo distributions by chromosome
#' @description Plot \bold{\emph{differentially methylated regions (DMRs)}} hyper/hypo distributions by chromosome.
#'
#' @return A \bold{\emph{ggplot object}} of chromosome-wise DMR distributions.
#' @param dmr_file DEG table from \bold{\emph{MethylKit}} analysis.
#' @param violin_scale Violin scale mode. (\bold{\emph{"count"}}, "area", "width").
#' @param violin_border Violin border width. (\bold{\emph{0.5}}).
#' @param point_shape Points shape (0-25). (\bold{\emph{8}}).
#' @param point_size Point size. (\bold{\emph{2}}).
#' @param jitter_width Horizontal jitter width. (\bold{\emph{0.2}}).
#' @param hyper_color Color for hyper-methylated points. (\bold{\emph{"#ff000088"}}).
#' @param hypo_color Color for hypo-methylated points. (\bold{\emph{"#00880088"}}).
#'
#' @importFrom magrittr %>%
#' @export
#'
#' @examples
#' # DMR results from Methylkit
#' dmr_file <- system.file(
#'     "extdata",
#'     "example.dmr",
#'     package = "GAnnoViz")
#'
#' # Plot
#' plot_dmg_chrom(
#'   dmr_file = dmr_file,
#'   violin_scale = "count",
#'   violin_border = 0.5,
#'   point_shape = 8,
#'   point_size = 2,
#'   jitter_width = 0.2,
#'   hyper_color = "#ff880088",
#'   hypo_color = "#0088ff88"
#' )
#'
plot_dmg_chrom <- function(dmr_file,
						   violin_scale = "count",
						   violin_border = 0.5,
						   point_shape = 8,
						   point_size = 2,
						   jitter_width = 0.2,
						   hyper_color = "#ff880088",
						   hypo_color = "#0088ff88") {
	# DMR results
  dmr <- utils::read.table(
    dmr_file,
    header = TRUE,
    sep = "\t",
    fill = TRUE,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

	# Chrom
	chrom_raw <- as.character(dmr$chr)
	chr_prefix <- any(stringr::str_detect(chrom_raw, "^(?i)chr"))
	keep_mask <- if (chr_prefix)
		stringr::str_detect(chrom_raw, "^(?i)chr")
	else
		rep(TRUE, length(chrom_raw))
	dmr <- dmr[keep_mask, , drop = FALSE]
	chrom_chr <- as.character(dmr$chr)
	chrom_num <- stringr::str_extract(chrom_chr, "\\d+")
	suppressWarnings({
		chrom_num <- as.numeric(chrom_num)
	})
	order_index <- order(is.na(chrom_num), chrom_num, chrom_chr)
	chrs <- unique(chrom_chr[order_index])

	# Data
	data <- dmr %>% dplyr::mutate(chr = factor(chr, levels = chrs))
	hyper_hypo <- data %>%
		dplyr::group_by(chr) %>%
		dplyr::summarise(
			n_positive = sum(meth.diff > 0, na.rm = TRUE),
			n_negative = sum(meth.diff < 0, na.rm = TRUE),
			max_y = suppressWarnings(max(meth.diff, na.rm = TRUE)),
			.groups = 'drop'
		) %>%
		dplyr::mutate(
			chr = factor(chr, levels = chrs),
			y_pos = max_y + (
				max(data$meth.diff, na.rm = TRUE) - min(data$meth.diff, na.rm = TRUE)
			) * 0.05,
			label = sprintf("+%d / -%d", n_positive, n_negative)
		)

	# Plot
	p <- ggplot2::ggplot(data, ggplot2::aes(x = chr, y = meth.diff)) +
		ggplot2::geom_violin(
			stat = "ydensity",
			position = "dodge",
			trim = TRUE,
			bounds = c(-Inf, Inf),
			scale = violin_scale,
			na.rm = FALSE,
			orientation = NA,
			show.legend = NA,
			inherit.aes = TRUE,
			width = 1,
			fill = "white",
			color = "black",
			linewidth = violin_border
		) +
		ggplot2::geom_point(
			ggplot2::aes(color = meth.diff > 0),
			position = ggplot2::position_jitter(width = jitter_width, height = 0),
			size = point_size,
			shape = point_shape
		) +
		ggplot2::scale_color_manual(
			values = c(`TRUE` = hyper_color, `FALSE` = hypo_color),
			name = "Diff",
			labels = c("Hypo", "Hyper")
		) +
		ggplot2::geom_text(
			data = hyper_hypo,
			ggplot2::aes(x = chr, y = 0, label = label),
			size = 3,
			hjust = 0.5,
			vjust = -1,
			angle = -90,
			color = "black"
		) +
		ggplot2::scale_x_discrete(drop = FALSE) +
		ggplot2::labs(x = "Chromosome", y = "Meth Diff (Hyper/Hypo)") +
		my_theme()
	p
}
