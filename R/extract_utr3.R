#' @title Extract 3'UTR ranges from GFF or GTF
#' @description Extract 3'UTR ranges from GFF or GTF.
#' @author benben-miao
#'
#' @return A \emph{GRanges} object.
#' @param gff_file Genomic structural annotation GFF3/GTF file path.
#' @param format Format of GFF3/GTF file. (\emph{"auto"}, "gff3", "gtf").
#' @param utr3_info mRNA information. (\emph{"all"}, "chrom_id", "utr3_id", "utr3_range").
#'
#' @export
#'
#' @examples
#' # Example GFF3 file in SlideAnno
#' gff_file <- system.file(
#'   "extdata",
#'   "example.gff",
#'   package = "SlideAnno")
#'
#' # Extract 3'UTR
#' utr3 <- extract_utr3(
#'   gff_file = gff_file,
#'   format = "auto",
#'   utr3_info = "all")
#' utr3
#'
#' # 3'UTR info: utr3_id
#' utr3_id <- extract_utr3(
#'   gff_file = gff_file,
#'   format = "auto",
#'   utr3_info = "utr3_id")
#' head(utr3_id)
#'
#' # 3'UTR info: utr3_range
#' utr3_range <- extract_utr3(
#'   gff_file = gff_file,
#'   format = "auto",
#'   utr3_info = "utr3_range")
#' head(utr3_range)
#'
extract_utr3 <- function(
		gff_file,
		format = "auto",
		utr3_info = "all") {
  # GFF3/GTF -> TXDB
  txdb <- suppressWarnings(txdbmaker::makeTxDbFromGFF(
  	file = gff_file,
  	format = format))

  # Extract UTR3
  utr3 <- suppressWarnings(unlist(
  	GenomicFeatures::threeUTRsByTranscript(txdb, use.names = FALSE)))

  # UTR3 info
  if (utr3_info == "all") {
  	utr3
  } else if (utr3_info == "chrom_id") {
  	levels(utr3@seqnames@values)
  } else if (utr3_info == "utr3_id") {
  	utr3$exon_name
  } else if (utr3_info == "utr3_range") {
  	utr3@ranges
  }
}