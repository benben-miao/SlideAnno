#' @title Annotate differentially expressed genes (DEGs) with chromosome positions
#' @description Annotate \bold{\emph{differentially expressed genes (DEGs)}} with chromosome positions.
#'
#' @return A data.frame with columns: (\bold{\emph{chrom, start, end, GeneID, log2FoldChange, strand}}).
#' @param deg_file DEG table from \bold{\emph{DESeq2}} analysis.
#' @param gff_file Genomic structural annotation \bold{\file{GFF3/GTF}} file path.
#' @param format Format of GFF3/GTF file. (\bold{\emph{"auto"}}, "gff3", "gtf").
#' @param id_col Gene IDs column name. (\bold{\emph{"GeneID"}}).
#' @param fc_col Log2(fold change) column name. (\bold{\emph{"log2FoldChange"}}).
#' @param use_strand Whether to resput actual strand instead of \code{"*"}. (\bold{\emph{FALSE}}).
#' @param drop_unmapped Whether to drop DEG entries not found in annotation. (\bold{\emph{TRUE}}).
#'
#' @export
#'
#' @examples
#' # Example DEGs and GFF in SlideAnno
#' deg_file <- system.file(
#'   "extdata",
#'   "example.deg",
#'   package = "SlideAnno")
#'
#' gff_file <- system.file(
#'   "extdata",
#'   "example.gff",
#'   package = "SlideAnno")
#'
#' # Annotate DEGs with chromosome positions
#' res <- anno_deg_chrom(
#'   deg_file = deg_file,
#'   gff_file = gff_file,
#'   format = "auto",
#'   id_col = "GeneID",
#'   fc_col = "log2FoldChange",
#'   use_strand = FALSE,
#'   drop_unmapped = TRUE
#' )
#' head(res)
#'
anno_deg_chrom <- function(deg_file,
                           gff_file,
                           format = "auto",
                           id_col = "GeneID",
                           fc_col = "log2FoldChange",
                           use_strand = FALSE,
                           drop_unmapped = TRUE) {
  # DEG results
  deg <- utils::read.table(
    deg_file,
    header = TRUE,
    sep = "\t",
    fill = TRUE,
    na.strings = "NA",
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  # Genes
  genes <- extract_genes(gff_file = gff_file,
                         format = format,
                         gene_info = "all")
  df_genes <- data.frame(
    chrom = as.character(GenomicRanges::seqnames(genes)),
    start = BiocGenerics::start(genes),
    end = BiocGenerics::end(genes),
    gene_id = S4Vectors::mcols(genes)$gene_id,
    strand = as.character(BiocGenerics::strand(genes)),
    stringsAsFactors = FALSE
  )
  deg_cols <- deg[, c(id_col, fc_col), drop = FALSE]
  colnames(deg_cols) <- c("gene_id", "score")

  # Results
  res <- dplyr::inner_join(df_genes, deg_cols, by = "gene_id")
  if (!use_strand)
    res$strand <- "*"
  res <- res[, c("chrom", "start", "end", "gene_id", "score", "strand")]

  if (!drop_unmapped) {
    all_ids <- unique(deg[[id_col]])
    missing_ids <- setdiff(all_ids, res$gene_id)
    if (length(missing_ids) > 0) {
      add <- data.frame(
        chrom = NA,
        start = NA,
        end = NA,
        gene_id = missing_ids,
        score = deg_cols$score[match(missing_ids, deg_cols$gene_id)],
        strand = "*",
        stringsAsFactors = FALSE
      )
      res <- dplyr::bind_rows(res, add)
    }
  }
  res
}
