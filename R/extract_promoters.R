#' @title Extract promoter ranges from GFF or GTF
#' @description Extract \bold{\emph{promoter}} ranges from GFF or GTF.
#' @author benben-miao
#'
#' @return A \bold{\emph{GRanges}} object.
#' @param gff_file Genomic structural annotation \bold{\file{GFF3/GTF}} file path.
#' @param format Format of GFF3/GTF file. (\bold{\emph{"auto"}}, "gff3", "gtf").
#' @param upstream Promoter upstream (bp). (\bold{\emph{2000}}).
#' @param downstream Promoter downstream (bp). (\bold{\emph{200}}).
#' @param promoter_info Promoter information. (\bold{\emph{"all"}}, "chrom_id", "promoter_id", "promoter_range").
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
#' # Extract Promoters
#' promoters <- extract_promoters(
#'   gff_file = gff_file,
#'   format = "auto",
#'   upstream = 2000,
#'   downstream = 200,
#'   promoter_info = "all")
#' promoters
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
  fmt <- resolve_gff_format(gff_file, format)
  txdb <- suppressWarnings(txdbmaker::makeTxDbFromGFF(
  	file = gff_file,
  	format = fmt))

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
