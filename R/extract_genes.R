#' @title Extract genes information from GFF or GTF
#' @description Extract \bold{\emph{genes}} information from GFF or GTF.
#' @author benben-miao
#'
#' @return A \bold{\emph{GRanges}} object.
#' @param gff_file Genomic structural annotation \bold{\file{GFF3/GTF}} file path.
#' @param format Format of GFF3/GTF file. (\bold{\emph{"auto"}}, "gff3", "gtf").
#' @param gene_info Gene information. (\bold{\emph{"all"}}, "chrom_id", "gene_id", "gene_range").
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
#' # Extract Genes
#' genes <- extract_genes(
#'   gff_file = gff_file,
#'   format = "auto",
#'   gene_info = "all")
#' genes
#'
#' # Gene info: gene_range
#' gene_range <- extract_genes(
#'   gff_file = gff_file,
#'   format = "auto",
#'   gene_info = "gene_range")
#' head(gene_range)
#'
extract_genes <- function(
		gff_file,
		format = "auto",
		gene_info = "all") {

  # GFF3/GTF -> TXDB
  fmt <- resolve_gff_format(gff_file, format)
  txdb <- suppressWarnings(txdbmaker::makeTxDbFromGFF(
  	file = gff_file,
  	format = fmt))

  # Extract Genes
  genes <- suppressWarnings(GenomicFeatures::genes(txdb))

  # Genes info
  if (gene_info == "all") {
  	genes
  } else if (gene_info == "chrom_id") {
  	levels(genes@seqnames@values)
  } else if (gene_info == "gene_id") {
  	genes$gene_id
  } else if (gene_info == "gene_range") {
  	genes@ranges
  }
}
