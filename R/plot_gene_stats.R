#' @title Plot gene stats for chromosomes
#' @description \bold{\emph{Plot}} gene stats for chromosomes.
#' @author benben-miao
#'
#' @return A plot object of gene stats for chromosomes.
#' @param gff_file Path to \bold{\file{GFF3/GTF}} file as input.
#' @param format Format of GFF3/GTF file. (\bold{\emph{"auto"}}, "gff3", "gtf").
#' @param bar_width Bar width percent. (\bold{\emph{0.7}}).
#' @param bar_color Bar color with name or hex code. (\bold{\emph{"#0055ff55"}}).
#' @param lable_size Lable text size. (\bold{\emph{3}}).
#'
#' @importFrom magrittr %>%
#' @export
#'
#' @examples
#' # Example GFF3 file in GAnnoViz
#' gff_file <- system.file(
#'     "extdata",
#'     "example.gff",
#'     package = "GAnnoViz")
#'
#' # Plot gene stats
#' plot_gene_stats(
#'     gff_file = gff_file,
#'     format = "auto",
#'     bar_width = 0.7,
#'     bar_color = "#0055ff55",
#'     lable_size = 3)
#'
plot_gene_stats <- function(gff_file,
							format = "auto",
							bar_width = 0.7,
							bar_color = "#0055ff55",
							lable_size = 3) {
	# Extract genes
	genes <- extract_genes(gff_file = gff_file,
						   format = format,
						   gene_info = "all")

	# Chrom <- Seqnames
	chrom_raw <- as.character(GenomicRanges::seqnames(genes))

	chrom_chr <- chrom_raw[stringr::str_detect(chrom_raw, "^(?i)chr")]
	if (length(chrom_chr) == 0) {
		warning("No chromosomes starting with 'chr|Chr|CHR' were detected!")
	}

	chrom_num <- stringr::str_extract(chrom_chr, "\\d+")
	suppressWarnings({
		chrom_num <- as.numeric(chrom_num)
	})
	order_index <- order(is.na(chrom_num), chrom_num, chrom_chr)
	chrom_levels <- unique(chrom_chr[order_index])
	chrom_factor <- factor(chrom_chr, levels = chrom_levels)

	df <- data.frame(chrom = chrom_factor) %>%
		dplyr::count(chrom)

	# Plot
	ggplot2::ggplot(df, ggplot2::aes(x = chrom, y = n)) +
		ggplot2::geom_col(na.rm = FALSE,
						  width = bar_width,
						  fill = bar_color) +
		ggplot2::geom_text(
			ggplot2::aes(label = n),
			size = lable_size,
			color = "#000000",
			vjust = -0.5
		) +
		ggplot2::labs(x = "Chromosome", y = "Genes density") +
		my_theme()
}
