#' @title Plot gene structure (Promoter, 3'UTR, Exon, Intron, 5'UTR)
#' @description Plot gene structure (\bold{\emph{Promoter, 3'UTR, Exon, Intron, 5'UTR}}).
#' @author benben-miao
#'
#' @return Plot of gene structure.
#' @param gff_file Genomic structural annotation \bold{\file{GFF3/GTF}} file path.
#' @param format Format of GFF3/GTF file. (\bold{\emph{"auto"}}, "gff3", "gtf").
#' @param gene_id Gene id same as GFF3/GTF. (\bold{\emph{necessary}}).
#' @param x_breaks X axis breaks number. (\bold{\emph{10}}).
#' @param upstream Promoter upstream (bp). (\bold{\emph{2000}}).
#' @param downstream Promoter downstream (bp). (\bold{\emph{200}}).
#' @param feature_alpha Elements alpha. (\bold{\emph{0.8}}).
#' @param intron_width Intron line width. (\bold{\emph{1}}).
#' @param arrow_count Intron arrow number bold. (\bold{\emph{1}}).
#' @param arrow_length Intron arrows length（pt). (\bold{\emph{1}}).
#' @param arrow_unit Intron arrow length unit. (\bold{\emph{"pt"}}, "mm").
#' @param promoter_color Promoter color. (\bold{\emph{"#ff8800"}}).
#' @param utr5_color 5'UTR color. (\bold{\emph{"#008833"}}).
#' @param utr3_color 3'UTR color. (\bold{\emph{"#ff0033"}}).
#' @param exon_color Exon color. (\bold{\emph{"#0033ff"}}).
#' @param intron_color Intron color. (\bold{\emph{"#333333"}}).
#'
#' @importFrom magrittr %>%
#' @export
#'
#' @examples
#' # Example GFF3 file in SlideAnno
#' gff_file <- system.file(
#'   "extdata",
#'   "example.gff",
#'   package = "SlideAnno")
#'
#' # Plot gene structure
#' plot_gene_structure(
#'   gff_file = gff_file,
#'   format = "auto",
#'   gene_id = "HdF029609",
#'   upstream = 2000,
#'   downstream = 200,
#'   feature_alpha = 0.8,
#'   intron_width = 1,
#'   x_breaks = 10,
#'   arrow_length = 5,
#'   arrow_count = 1,
#'   arrow_unit = "pt",
#'   promoter_color = "#ff8800",
#'   utr5_color = "#008833",
#'   utr3_color = "#ff0033",
#'   exon_color = "#0033ff",
#'   intron_color = "#333333"
#' )
#'
plot_gene_structure <- function(gff_file,
								format = "auto",
								gene_id,
								x_breaks = 10,
								upstream = 2000,
								downstream = 200,
								feature_alpha = 0.8,
								intron_width = 1,
								arrow_count = 1,
								arrow_length = 5,
								arrow_unit = "pt",
								promoter_color = "#ff8800",
								utr5_color = "#008833",
								utr3_color = "#ff0033",
								exon_color = "#0033ff",
								intron_color = "#333333") {
	# TXDB
	txdb <- suppressWarnings(txdbmaker::makeTxDbFromGFF(file = gff_file, format = format))

	# Genes
	genes <- suppressWarnings(GenomicFeatures::genes(txdb))
	if (missing(gene_id)) {
		stop("Please input a correct gene_id!")
	}
	if (!gene_id %in% genes$gene_id) {
		stop("Input gene_id no included in GFF3/GTF!", gene_id)
	}
	gene_gr <- genes[genes$gene_id == gene_id]

	# Gene range
	tx_by_gene <- suppressWarnings(GenomicFeatures::transcriptsBy(txdb, by = "gene"))
	if (!gene_id %in% names(tx_by_gene)) {
		stop("Input gene_id not include transcript!", gene_id)
	}
	tx_gr <- tx_by_gene[[gene_id]]
	tx_names <- tx_gr$tx_name
	tx_strand <- as.character(tx_gr@strand)

	# Promoters, Exons, Introns
	exons_by_tx <- suppressWarnings(GenomicFeatures::exonsBy(txdb, by = "tx", use.names = TRUE))
	introns_by_tx <- suppressWarnings(GenomicFeatures::intronsByTranscript(txdb, use.names = TRUE))
	utr5_by_tx <- suppressWarnings(GenomicFeatures::fiveUTRsByTranscript(txdb, use.names = TRUE))
	utr3_by_tx <- suppressWarnings(GenomicFeatures::threeUTRsByTranscript(txdb, use.names = TRUE))
	promoters_tx <- suppressWarnings(
		GenomicFeatures::promoters(
			tx_gr,
			upstream = upstream,
			downstream = downstream,
			use.names = TRUE
		)
	)
	keep_tx <- tx_names

	# Elements data.frame
	to_df <- function(gr, feature, tx_label = NA_character_) {
		if (length(gr) == 0)
			return(data.frame())
		data.frame(
			chrom = as.character(GenomicRanges::seqnames(gr)),
			start = BiocGenerics::start(gr),
			end = BiocGenerics::end(gr),
			strand = as.character(GenomicRanges::strand(gr)),
			tx_name = if (!is.na(tx_label))
				tx_label
			else if (!is.null(gr$tx_name))
				gr$tx_name
			else
				NA_character_,
			feature = feature,
			stringsAsFactors = FALSE
		)
	}

	build_list_df <- function(grl, feature) {
		if (is.null(grl) || length(grl) == 0)
			return(data.frame())
		lst <- grl[names(grl) %in% keep_tx]
		if (length(lst) == 0)
			return(data.frame())
		dfs <- mapply(
			function(gr, nm)
				to_df(gr, feature, tx_label = nm),
			gr = as.list(lst),
			nm = names(lst),
			SIMPLIFY = FALSE
		)
		do.call(rbind, dfs)
	}

	df_exon <- build_list_df(exons_by_tx, "exon")
	df_intron <- build_list_df(introns_by_tx, "intron")
	df_utr5 <- build_list_df(utr5_by_tx, "utr5")
	df_utr3 <- build_list_df(utr3_by_tx, "utr3")

	promoters_list <- split(promoters_tx, promoters_tx$tx_name)
	df_promoter <- build_list_df(promoters_list, "promoter")

	df_all <- dplyr::bind_rows(df_promoter, df_utr5, df_utr3, df_exon, df_intron)
	df_all <- df_all[stats::complete.cases(df_all), , drop = FALSE]

	tx_levels <- unique(df_all$tx_name)
	df_all$tx_factor <- factor(df_all$tx_name, levels = tx_levels)
	df_all$y <- as.numeric(df_all$tx_factor)

	# Start, End
	xmin <- BiocGenerics::start(gene_gr)
	xmax <- BiocGenerics::end(gene_gr)
	xpad <- max(upstream, downstream)

	# Stats
	chrom_id <- as.character(GenomicRanges::seqnames(gene_gr))
	exons_gene <- suppressWarnings(GenomicFeatures::exonsBy(txdb, by = "gene"))
	exon_count <- if (gene_id %in% names(exons_gene))
		length(exons_gene[[gene_id]])
	else
		sum(df_all$feature == "exon")

	feature_order <- c("promoter", "utr5", "utr3", "exon", "intron")
	df_all$feature <- factor(df_all$feature, levels = feature_order)

	# Colors
	col_map <- c(
		promoter = promoter_color,
		utr5 = utr5_color,
		utr3 = utr3_color,
		exon = exon_color,
		intron = intron_color
	)

	df_intron_plot <- df_all[df_all$feature == "intron", , drop = FALSE]
	if (nrow(df_intron_plot) > 0) {
		df_intron_plot$x <- ifelse(df_intron_plot$strand == "-",
								   df_intron_plot$end,
								   df_intron_plot$start)
		df_intron_plot$xend <- ifelse(df_intron_plot$strand == "-",
									  df_intron_plot$start,
									  df_intron_plot$end)
	}

	# Arrows
	df_intron_arrows <- data.frame()
	if (nrow(df_intron_plot) > 0) {
		for (i in seq_len(nrow(df_intron_plot))) {
			s <- df_intron_plot$start[i]
			e <- df_intron_plot$end[i]
			y <- df_intron_plot$y[i]
			st <- df_intron_plot$strand[i]
			len <- e - s
			if (len <= 0)
				next
			cnt <- if (is.na(arrow_count)) {
				max(1, min(20, floor(len / 1000)))
			} else {
				max(1, as.integer(arrow_count))
			}
			seg_len <- max(1, len / (cnt * 6))
			if (st == "-") {
				pos <- seq(
					from = e,
					to = s + seg_len,
					length.out = cnt
				)
				x <- pos
				xend <- pos - seg_len
			} else {
				pos <- seq(
					from = s,
					to = e - seg_len,
					length.out = cnt
				)
				x <- pos
				xend <- pos + seg_len
			}
			df_i <- data.frame(
				x = x,
				xend = xend,
				y = y,
				yend = y,
				feature = "intron"
			)
			df_intron_arrows <- rbind(df_intron_arrows, df_i)
		}
	}

	# Plot
	p <- ggplot2::ggplot(df_all) +
		ggplot2::geom_rect(
			data = df_all[df_all$feature %in% c("promoter", "utr5", "utr3", "exon"), , drop = FALSE],
			ggplot2::aes(
				xmin = start,
				xmax = end,
				ymin = y - 0.3,
				ymax = y + 0.3,
				fill = feature
			),
			color = NA,
			alpha = feature_alpha
		) +
		ggplot2::geom_segment(
			data = df_intron_plot,
			ggplot2::aes(
				x = x,
				xend = xend,
				y = y,
				yend = y,
				color = feature
			),
			linewidth = intron_width,
			lineend = "round"
		) +
		ggplot2::geom_segment(
			data = df_intron_arrows,
			ggplot2::aes(
				x = x,
				xend = xend,
				y = y,
				yend = y,
				color = feature
			),
			linewidth = intron_width,
			lineend = "round",
			arrow = grid::arrow(
				length = grid::unit(arrow_length, arrow_unit),
				type = "closed"
			)
		) +
		ggplot2::scale_fill_manual(values = col_map, drop = FALSE) +
		ggplot2::scale_color_manual(values = col_map, drop = FALSE) +
		ggplot2::scale_y_continuous(
			breaks = seq_along(tx_levels),
			labels = tx_levels,
			expand = ggplot2::expansion(mult = c(0.05, 0.05))
		) +
		ggplot2::labs(
			x = "Genomic position",
			y = "Transcript",
			title = paste0(
				"Chrom: ",
				chrom_id,
				"; Gene: ",
				gene_id,
				"; Exons: ",
				exon_count,
				"; Introns: ",
				sum(df_all$feature == "intron")
			)
		) +
		ggplot2::scale_x_continuous(breaks = scales::pretty_breaks(n = x_breaks)) +
		ggplot2::coord_cartesian(xlim = c(xmin - xpad, xmax + xpad)) +
		my_theme()
	p
}
