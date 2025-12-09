#' @title Extract promoter ranges from GFF or GTF
#' @description Extract promoter ranges from GFF or GTF.
#' @author benben-miao
#'
#' @return A \emph{GRanges} object.
#' @param gff_file Genomic structural annotation GFF3/GTF file path.
#' @param format Format of GFF3/GTF file. (\emph{"auto"}, "gff3", "gtf").
#' @param upstream Promoter upstream (bp). (\emph{2000}).
#' @param downstream Promoter downstream (bp). (\emph{200}).
#' @param promoter_info Promoter information. (\emph{"all"}, "chrom_id", "promoter_id", "promoter_range").
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
#' # Extract Promoters
#' promoters <- extract_promoters(
#'   gff_file = gff_file,
#'   format = "auto",
#'   upstream = 2000,
#'   downstream = 200,
#'   promoter_info = "all")
#' promoters
#'
#' # Promoter info: promoter_id
#' promoter_id <- extract_promoters(
#'   gff_file = gff_file,
#'   format = "auto",
#'   upstream = 2000,
#'   downstream = 200,
#'   promoter_info = "promoter_id")
#' head(promoter_id)
#'
#' # Promoter info: promoter_range
#' promoter_range <- extract_promoters(
#'   gff_file = gff_file,
#'   format = "auto",
#'   upstream = 2000,
#'   downstream = 200,
#'   promoter_info = "promoter_range")
#' head(promoter_range)
#'
extract_promoters <- function(
		gff_file,
		format = "auto",
		upstream = 2000,
		downstream = 200,
		promoter_info = "all") {
  # GFF3/GTF -> TXDB
  txdb <- suppressWarnings(txdbmaker::makeTxDbFromGFF(
  	file = gff_file,
  	format = format))

  # Extract Promoters
  promoters <- suppressWarnings(GenomicFeatures::promoters(
  	txdb,
  	upstream = 2000,
  	downstream = 200,
  	use.names = TRUE))

  # Promoter info
  if (promoter_info == "all") {
  	promoters
  } else if (promoter_info == "chrom_id") {
  	levels(promoters@seqnames@values)
  } else if (promoter_info == "promoter_id") {
  	promoters$tx_name
  } else if (promoter_info == "promoter_range") {
  	promoters@ranges
  }
}