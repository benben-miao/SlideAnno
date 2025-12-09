#' @title Extract Exons ranges from GFF or GTF
#' @description Extract Exons ranges from GFF or GTF.
#' @author benben-miao
#'
#' @return A \emph{GRanges} object.
#' @param gff_file Genomic structural annotation GFF3/GTF file path.
#' @param format Format of GFF3/GTF file. (\emph{"auto"}, "gff3", "gtf").
#' @param exon_info Exon information. (\emph{"all"}, "chrom_id", "exon_range").
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
#' # Extract Exons
#' exons <- extract_exons(
#'   gff_file = gff_file,
#'   format = "auto",
#'   exon_info = "all")
#' exons
#'
#' # Exon info: exon_range
#' exon_range <- extract_exons(
#'   gff_file = gff_file,
#'   format = "auto",
#'   exon_info = "exon_range")
#' head(exon_range)
#'
extract_exons <- function(
		gff_file,
		format = "auto",
		exon_info = "all") {
  # GFF3/GTF -> TXDB
  txdb <- suppressWarnings(txdbmaker::makeTxDbFromGFF(
  	file = gff_file,
  	format = format))

  # Extract Exons
  exons <- suppressWarnings(GenomicFeatures::exons(txdb))

  # Exon info
  if (exon_info == "all") {
  	exons
  } else if (exon_info == "chrom_id") {
  	levels(exons@seqnames@values)
  } else if (exon_info == "exon_range") {
  	exons@ranges
  }
}