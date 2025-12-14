#' @title Extract CDS ranges from GFF or GTF
#' @description Extract \bold{\emph{CDS}} ranges from GFF or GTF.
#' @author benben-miao
#'
#' @return A \bold{\emph{GRanges}} object.
#' @param gff_file Genomic structural annotation \bold{\file{GFF3/GTF}} file path.
#' @param format Format of GFF3/GTF file. (\bold{\emph{"auto"}}, "gff3", "gtf").
#' @param cds_info CDS information. (\bold{\emph{"all"}}, "chrom_id", "cds_range").
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
  fmt <- resolve_gff_format(gff_file, format)
  txdb <- suppressWarnings(txdbmaker::makeTxDbFromGFF(
  	file = gff_file,
  	format = fmt))

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
