#' @title Extract mRNA transcripts from GFF or GTF
#' @description Extract mRNA transcripts from GFF or GTF.
#' @author benben-miao
#'
#' @return A \emph{GRanges} object.
#' @param gff_file Genomic structural annotation GFF3/GTF file path.
#' @param format Format of GFF3/GTF file. (\emph{"auto"}, "gff3", "gtf").
#' @param mrna_info mRNA information. (\emph{"all"}, "chrom_id", "mrna_id", "mrna_range").
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
#' # Extract mRNAs
#' mrnas <- extract_mrnas(
#'   gff_file = gff_file,
#'   format = "auto",
#'   mrna_info = "all")
#' mrnas
#'
#' # mRNA info: mrna_id
#' mrna_id <- extract_mrnas(
#'   gff_file = gff_file,
#'   format = "auto",
#'   mrna_info = "mrna_id")
#' head(mrna_id)
#'
#' # mRNA info: mrna_range
#' mrna_range <- extract_mrnas(
#'   gff_file = gff_file,
#'   format = "auto",
#'   mrna_info = "mrna_range")
#' head(mrna_range)
#'
extract_mrnas <- function(
		gff_file,
		format = "auto",
		mrna_info = "all") {
  # GFF3/GTF -> TXDB
  txdb <- suppressWarnings(txdbmaker::makeTxDbFromGFF(
  	file = gff_file,
  	format = format))

  # Extract Genes
  mrnas <- suppressWarnings(GenomicFeatures::transcripts(txdb))

  # Genes info
  if (mrna_info == "all") {
  	mrnas
  } else if (mrna_info == "chrom_id") {
  	levels(mrnas@seqnames@values)
  } else if (mrna_info == "mrna_id") {
  	mrnas$tx_name
  } else if (mrna_info == "mrna_range") {
  	mrnas@ranges
  }
}