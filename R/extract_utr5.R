#' @title Extract 5'UTR ranges from GFF or GTF
#' @description Extract \bold{\emph{5'UTR}} ranges from GFF or GTF.
#' @author benben-miao
#'
#' @return A \bold{\emph{GRanges}} object.
#' @param gff_file Genomic structural annotation \bold{\file{GFF3/GTF}} file path.
#' @param format Format of GFF3/GTF file. (\bold{\emph{"auto"}}, "gff3", "gtf").
#' @param utr5_info mRNA information. (\bold{\emph{"all"}}, "chrom_id", "utr5_id", "utr5_range").
#'
#' @export
#'
#' @examples
#' # Example GFF3 file in GAnnoViz
#' gff_file <- system.file(
#'   "extdata",
#'   "example.gff3.gz",
#'   package = "GAnnoViz")
#'
#' # Extract 5'UTR
#' utr5 <- extract_utr5(
#'   gff_file = gff_file,
#'   format = "auto",
#'   utr5_info = "all")
#' utr5
#'
#' # 5'UTR info: utr5_range
#' utr5_range <- extract_utr5(
#'   gff_file = gff_file,
#'   format = "auto",
#'   utr5_info = "utr5_range")
#' head(utr5_range)
#'
extract_utr5 <- function(
		gff_file,
		format = "auto",
		utr5_info = "all") {

  # GFF3/GTF -> TXDB
  fmt <- resolve_gff_format(gff_file, format)
  txdb <- suppressWarnings(txdbmaker::makeTxDbFromGFF(
  	file = gff_file,
  	format = fmt))

  # Extract UTR5
  utr5 <- suppressWarnings(unlist(
  	GenomicFeatures::fiveUTRsByTranscript(txdb, use.names = FALSE)))

  # UTR5 info
  if (utr5_info == "all") {
  	utr5
  } else if (utr5_info == "chrom_id") {
  	levels(utr5@seqnames@values)
  } else if (utr5_info == "utr5_id") {
  	utr5$exon_name
  } else if (utr5_info == "utr5_range") {
  	utr5@ranges
  }
}
