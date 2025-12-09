#' @title Extract CDS ranges from GFF or GTF
#' @description Extract CDS ranges from GFF or GTF.
#' @author benben-miao
#'
#' @return A \emph{GRanges} object.
#' @param gff_file Genomic structural annotation GFF3/GTF file path.
#' @param format Format of GFF3/GTF file. (\emph{"auto"}, "gff3", "gtf").
#' @param cds_info CDS information. (\emph{"all"}, "chrom_id", "cds_range").
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
#' # Extract CDS
#' cds <- extract_cds(
#'   gff_file = gff_file,
#'   format = "auto",
#'   cds_info = "all")
#' cds
#'
#' # CDS info: cds_range
#' cds_range <- extract_cds(
#'   gff_file = gff_file,
#'   format = "auto",
#'   cds_info = "cds_range")
#' head(cds_range)
#'
extract_cds <- function(
		gff_file,
		format = "auto",
		cds_info = "all") {
  # GFF3/GTF -> TXDB
  txdb <- suppressWarnings(txdbmaker::makeTxDbFromGFF(
  	file = gff_file,
  	format = format))

  # Extract CDS
  cds <- suppressWarnings(GenomicFeatures::cds(txdb))

  # CDS info
  if (cds_info == "all") {
  	cds
  } else if (cds_info == "chrom_id") {
  	levels(cds@seqnames@values)
  } else if (cds_info == "cds_range") {
  	cds@ranges
  }
}