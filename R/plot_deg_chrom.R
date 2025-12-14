#' @title Plot differentially expressed genes (DEGs) hyper/hypo distributions by chromosome
#' @description Plot \bold{\emph{differentially expressed genes (DEGs)}} hyper/hypo distributions by chromosomes.
#'
#' @return A \bold{\emph{ggplot object}} of chromosome-wise DEG distributions.
#' @param deg_file DEG table from \bold{\emph{DESeq2}} analysis.
#' @param gff_file Genomic structural annotation \bold{\file{GFF3/GTF}} file path.
#' @param format Format of GFF3/GTF file. (\bold{\emph{"auto"}}, "gff3", "gtf").
#' @param id_col Gene IDs column name. (\bold{\emph{"GeneID"}}).
#' @param fc_col Log2(fold change) column name. (\bold{\emph{"log2FoldChange"}}).
#' @param violin_scale Violin scale mode. (\bold{\emph{"count"}}, "area", "width").
#' @param violin_border Violin border width. (\bold{\emph{0.5}}).
#' @param point_shape Points shape (0-25). (\bold{\emph{16}}).
#' @param point_size Point size. (\bold{\emph{2}}).
#' @param jitter_width Horizontal jitter width. (\bold{\emph{0.2}}).
#' @param hyper_color Color for up-regulated points. (\bold{\emph{"#ff000088"}}).
#' @param hypo_color Color for down-regulated points. (\bold{\emph{"#00880088"}}).
#'
#' @importFrom magrittr %>%
#' @export
#'
#' @examples
#' # DEG results from DESeq2
#' deg_file <- system.file(
#'     "extdata",
#'     "example.deg",
#'     package = "GAnnoViz")
#'
#' # Genomic structure annotation
#' gff_file <- system.file(
#'     "extdata",
#'     "example.gff3.gz",
#'     package = "GAnnoViz")
#'
#' # Plot
#' plot_deg_chrom(
#'   deg_file = deg_file,
#'   gff_file = gff_file,
#'   format = "auto",
#'   id_col = "GeneID",
#'   fc_col = "log2FoldChange",
#'   violin_scale = "count",
#'   violin_border = 0.5,
#'   point_shape = 16,
#'   point_size = 2,
#'   jitter_width = 0.2,
#'   hyper_color = "#ff000088",
#'   hypo_color = "#00880088"
#' )
#'
plot_deg_chrom <- function(deg_file,
						   gff_file,
						   format = "auto",
						   id_col = "GeneID",
						   fc_col = "log2FoldChange",
						   violin_scale = "count",
						   violin_border = 0.5,
						   point_shape = 16,
						   point_size = 2,
						   jitter_width = 0.2,
						   hyper_color = "#ff000088",
						   hypo_color = "#00880088") {
	# DEG results
	df <- anno_deg_chrom(
		deg_file = deg_file,
		gff_file = gff_file,
		format = format,
		id_col = id_col,
		fc_col = fc_col,
		use_strand = FALSE,
		drop_unmapped = TRUE
	)
	df <- df %>%
		dplyr::rename(chr = chrom,
					  log2fc = score,
					  geneid = gene_id)
	chrom_raw <- as.character(df$chr)
	chr_prefix <- any(stringr::str_detect(chrom_raw, "^(?i)chr"))
	keep_mask <- if (chr_prefix)
		stringr::str_detect(chrom_raw, "^(?i)chr")
	else
		rep(TRUE, length(chrom_raw))
	df <- df[keep_mask, , drop = FALSE]
	chrom_chr <- as.character(df$chr)
	chrom_num <- stringr::str_extract(chrom_chr, "\\d+")
	suppressWarnings({
		chrom_num <- as.numeric(chrom_num)
	})
	order_index <- order(is.na(chrom_num), chrom_num, chrom_chr)
	chrs <- unique(chrom_chr[order_index])

	# Data
	data <- df %>%
		dplyr::mutate(chr = factor(chr, levels = chrs))
	hyper_hypo <- data %>%
		dplyr::group_by(chr) %>%
		dplyr::summarise(
			n_positive = sum(log2fc > 0, na.rm = TRUE),
			n_negative = sum(log2fc < 0, na.rm = TRUE),
			max_y = suppressWarnings(max(log2fc, na.rm = TRUE)),
			.groups = 'drop'
		) %>%
		dplyr::mutate(
			chr = factor(chr, levels = chrs),
			y_pos = max_y + (
				max(data$log2fc, na.rm = TRUE) - min(data$log2fc, na.rm = TRUE)
			) * 0.05,
			label = sprintf("+%d / -%d", n_positive, n_negative)
		)

	# Plot
	p <- ggplot2::ggplot(data, ggplot2::aes(x = chr, y = log2fc)) +
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
			ggplot2::aes(color = log2fc > 0),
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
		ggplot2::labs(x = "Chromosome", y = "Log2FoldChange (Hyper/Hypo)") +
		my_theme()
	p
}
