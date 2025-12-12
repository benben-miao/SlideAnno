#' @title Annotate FST/DMR slide windows with genomic features
#' @description Annotate \bold{\emph{FST/DMR}} slide windows with \bold{\emph{genomic features}}.
#' @author benben-miao
#'
#' @return Data.frame with input intervals and \bold{\emph{annotation results}}.
#' @param gff_file Genomic structural annotation \bold{\file{GFF3/GTF}} file path.
#' @param format Format of GFF3/GTF file. (\bold{\emph{"auto"}}, "gff3", "gtf").
#' @param genomic_ranges Genomic ranges file (e.g., FST, DMR). (\bold{\emph{chromosome (prefix: chr), start, end}}).
#' @param chrom_col Chromosome column name. (\bold{\emph{FST: "CHROM", DMR: "chr"}}).
#' @param start_col Start column name. (\bold{\emph{FST: "BIN_START", DMR: "start"}}).
#' @param end_col End column name. (\bold{\emph{FST: "BIN_END", DMR: "end"}}).
#' @param upstream Promoter upstream (bp). (\bold{\emph{2000}}).
#' @param downstream Promoter downstream (bp). (\bold{\emph{200}}).
#' @param ignore_strand Whether to ignore strand when computing overlaps. (\bold{\emph{TRUE}}).
#' @param features Which feature types to annotate. c(\bold{\emph{"promoter", "UTR5", "gene", "exon", "intron", "CDS", "UTR3", "intergenic"}}).
#'
#' @export
#'
#' @examples
#' # Example GFF3 file in GAnnoViz
#' gff_file <- system.file(
#'   "extdata",
#'   "example.gff",
#'   package = "GAnnoViz")
#'
#' # Annotate FST
#' fst_table <- system.file(
#'     "extdata",
#'     "example.fst",
#'     package = "GAnnoViz")
#'
#' res <- anno_fst_dmr(
#'   gff_file = gff_file,
#'   format = "auto",
#'   genomic_ranges = fst_table,
#'   chrom_col = "CHROM",
#'   start_col = "BIN_START",
#'   end_col = "BIN_END",
#'   upstream = 2000,
#'   downstream = 200,
#'   ignore_strand = TRUE,
#'   features = c("promoter", "UTR5", "gene", "exon", "intron", "CDS", "UTR3", "intergenic")
#' )
#' head(res)
#'
#' # Annotate DMR
#' dmr_table <- system.file(
#'     "extdata",
#'     "example.dmr",
#'     package = "GAnnoViz")
#'
#' res <- anno_fst_dmr(
#'   gff_file = gff_file,
#'   format = "auto",
#'   genomic_ranges = dmr_table,
#'   chrom_col = "chr",
#'   start_col = "start",
#'   end_col = "end",
#'   upstream = 2000,
#'   downstream = 200,
#'   ignore_strand = TRUE,
#'   features = c("promoter", "UTR5", "gene", "exon", "intron", "CDS", "UTR3", "intergenic")
#' )
#' head(res)
#'
anno_fst_dmr <- function(gff_file,
						 format = "auto",
						 genomic_ranges,
						 chrom_col = "CHROM",
						 start_col = "BIN_START",
						 end_col = "BIN_END",
						 upstream = 2000,
						 downstream = 200,
						 ignore_strand = TRUE,
						 features = c("promoter", "UTR5", "gene", "exon", "intron", "CDS", "UTR3", "intergenic")) {
	# TXDB
	fmt <- resolve_gff_format(gff_file, format)
	txdb <- suppressWarnings(txdbmaker::makeTxDbFromGFF(file = gff_file, format = fmt))

	# Features
	promoter <- if ("promoter" %in% features)
		suppressWarnings(
			GenomicFeatures::promoters(
				txdb,
				upstream = upstream,
				downstream = downstream,
				use.names = TRUE
			)
		)
	else
		NULL
	gene <- if ("gene" %in% features)
		suppressWarnings(GenomicFeatures::genes(txdb))
	else
		NULL
	utr5 <- if ("UTR5" %in% features)
		suppressWarnings(unlist(
			GenomicFeatures::fiveUTRsByTranscript(txdb, use.names = FALSE)
		))
	else
		NULL
	cds <- if ("CDS" %in% features)
		suppressWarnings(GenomicFeatures::cds(txdb))
	else
		NULL
	exon <- if ("exon" %in% features)
		suppressWarnings(GenomicFeatures::exons(txdb))
	else
		NULL
	intron <- if ("intron" %in% features)
		suppressWarnings(unlist(
			GenomicFeatures::intronsByTranscript(txdb, use.names = FALSE)
		))
	else
		NULL
	utr3 <- if ("UTR3" %in% features)
		suppressWarnings(unlist(
			GenomicFeatures::threeUTRsByTranscript(txdb, use.names = FALSE)
		))
	else
		NULL

	# Map promoter transcript name/id to gene_id
	if (!is.null(promoter)) {
		tx_by_gene <- suppressWarnings(GenomicFeatures::transcriptsBy(txdb, by = "gene"))
		if (!is.null(tx_by_gene) && length(tx_by_gene) > 0) {
			map_df <- do.call(rbind, lapply(names(tx_by_gene), function(gid) {
				gr <- tx_by_gene[[gid]]
				tx_nm <- as.character(S4Vectors::mcols(gr)$tx_name)
				if (is.null(tx_nm))
					tx_nm <- as.character(S4Vectors::mcols(gr)$tx_id)
				data.frame(
					tx_name = tx_nm,
					gene_id = gid,
					stringsAsFactors = FALSE
				)
			}))
			if (!is.null(S4Vectors::mcols(promoter)$tx_name)) {
				idx <- match(as.character(S4Vectors::mcols(promoter)$tx_name),
							 map_df$tx_name)
				S4Vectors::mcols(promoter)$gene_id <- map_df$gene_id[idx]
			} else if (!is.null(S4Vectors::mcols(promoter)$tx_id)) {
				idx <- match(as.character(S4Vectors::mcols(promoter)$tx_id),
							 map_df$tx_name)
				S4Vectors::mcols(promoter)$gene_id <- map_df$gene_id[idx]
			}
		}
	}

  # Genomic ranges file
  range_df <- utils::read.table(
    genomic_ranges,
    header = TRUE,
    sep = "\t",
    fill = TRUE,
    stringsAsFactors = FALSE
  )
	if (!all(c(chrom_col, start_col, end_col) %in% colnames(range_df)))
		stop("genomic_ranges must contain columns: ", paste(c(chrom_col, start_col, end_col), collapse = ", "))

	# Genomic Ranges
	gene_range <- GenomicRanges::GRanges(
		seqnames = range_df[[chrom_col]],
		ranges = IRanges::IRanges(start = as.integer(range_df[[start_col]]), end = as.integer(range_df[[end_col]]))
	)

	# Chroms
	chroms <- unique(as.character(range_df[[chrom_col]]))
	filt <- function(gr) {
		if (is.null(gr))
			return(NULL)
		gr[as.character(GenomicRanges::seqnames(gr)) %in% chroms]
	}
	promoter <- filt(promoter)
	gene <- filt(gene)
	utr5 <- filt(utr5)
	cds <- filt(cds)
	exon <- filt(exon)
	intron <- filt(intron)
	utr3 <- filt(utr3)

	# Overlap
	get_overlap_info <- function(query,
								 subject,
								 type_label,
								 id_col = NULL,
								 id_suffix = NULL) {
		hits <- GenomicRanges::findOverlaps(query, subject, ignore.strand = ignore_strand)
		if (length(hits) == 0)
			return(NULL)
		ids <- if (!is.null(id_col) &&
				   id_col %in% names(S4Vectors::mcols(subject))) {
			paste0(as.character(S4Vectors::mcols(subject)[[id_col]][S4Vectors::subjectHits(hits)]), id_suffix)
		} else
			rep(NA_character_, length(S4Vectors::subjectHits(hits)))
		data.frame(
			idx = S4Vectors::queryHits(hits),
			anno_type = type_label,
			gene_id = ids,
			stringsAsFactors = FALSE
		)
	}

	anno_promoter <- if (!is.null(promoter))
		get_overlap_info(gene_range,
						 promoter,
						 "promoter",
						 id_col = "gene_id",
						 id_suffix = "(p)")
	else
		NULL
	anno_gene <- if (!is.null(gene))
		get_overlap_info(gene_range,
						 gene,
						 "gene",
						 id_col = "gene_id",
						 id_suffix = "(g)")
	else
		NULL
	anno_utr5 <- if (!is.null(utr5))
		get_overlap_info(gene_range, utr5, "UTR5")
	else
		NULL
	anno_cds <- if (!is.null(cds))
		get_overlap_info(gene_range, cds, "CDS")
	else
		NULL
	anno_exon <- if (!is.null(exon))
		get_overlap_info(gene_range, exon, "exon")
	else
		NULL
	anno_intron <- if (!is.null(intron))
		get_overlap_info(gene_range, intron, "intron")
	else
		NULL
	anno_utr3 <- if (!is.null(utr3))
		get_overlap_info(gene_range, utr3, "UTR3")
	else
		NULL

	# Annotation
	anno_all <- dplyr::bind_rows(anno_promoter,
								 anno_gene,
								 anno_utr5,
								 anno_cds,
								 anno_exon,
								 anno_intron,
								 anno_utr3)
	if (is.null(anno_all) || nrow(anno_all) == 0) {
		anno_full <- data.frame(
			idx = seq_len(length(gene_range)),
			anno_type = "intergenic",
			gene_id = "",
			stringsAsFactors = FALSE
		)
	} else {
		anno_summary <- dplyr::summarise(
			dplyr::group_by(anno_all, idx),
			anno_type = paste(unique(anno_type), collapse = ","),
			gene_id = paste(stats::na.omit(unique(gene_id)), collapse = ","),
			.groups = "drop"
		)
		anno_full <- data.frame(idx = seq_len(length(gene_range)),
								stringsAsFactors = FALSE)
		anno_full <- dplyr::left_join(anno_full, anno_summary, by = "idx")
		anno_full <- dplyr::mutate(
			anno_full,
			anno_type = ifelse(is.na(anno_type), "intergenic", anno_type),
			gene_id = ifelse(is.na(gene_id), "", gene_id)
		)
	}

	# Bind
	res <- cbind(range_df, anno_full[, c("anno_type", "gene_id")])
	res
}
