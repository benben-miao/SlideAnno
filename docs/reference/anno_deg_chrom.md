# Annotate differentially expressed genes (DEGs) with chromosome positions

Annotate ***differentially expressed genes (DEGs)*** with chromosome
positions.

## Usage

``` r
anno_deg_chrom(
  deg_file,
  gff_file,
  format = "auto",
  id_col = "GeneID",
  fc_col = "log2FoldChange",
  use_strand = FALSE,
  drop_unmapped = TRUE
)
```

## Arguments

- deg_file:

  DEG table from ***DESeq2*** analysis.

- gff_file:

  Genomic structural annotation **`GFF3/GTF`** file path.

- format:

  Format of GFF3/GTF file. (***"auto"***, "gff3", "gtf").

- id_col:

  Gene IDs column name. (***"GeneID"***).

- fc_col:

  Log2(fold change) column name. (***"log2FoldChange"***).

- use_strand:

  Whether to resput actual strand instead of `"*"`. (***FALSE***).

- drop_unmapped:

  Whether to drop DEG entries not found in annotation. (***TRUE***).

## Value

A data.frame with columns: (***chrom, start, end, GeneID,
log2FoldChange, strand***).

## Examples

``` r
# Example DEGs and GFF in GAnnoViz
deg_file <- system.file(
  "extdata",
  "example.deg",
  package = "GAnnoViz")

gff_file <- system.file(
  "extdata",
  "example.gff",
  package = "GAnnoViz")

# Annotate DEGs with chromosome positions
res <- anno_deg_chrom(
  deg_file = deg_file,
  gff_file = gff_file,
  format = "auto",
  id_col = "GeneID",
  fc_col = "log2FoldChange",
  use_strand = FALSE,
  drop_unmapped = TRUE
)
#> Import genomic features from the file as a GRanges object ... 
#> OK
#> Prepare the 'metadata' data frame ... 
#> OK
#> Make the TxDb object ... 
#> OK
head(res)
#>   chrom    start      end   gene_id     score strand
#> 1 chr14 12802589 12807112 HdF000031  1.408356      *
#> 2 chr14 12979855 12984114 HdF000041 -1.114777      *
#> 3 chr14 13092225 13094734 HdF000049  1.271638      *
#> 4 chr14 13915421 13919063 HdF000066 -6.180426      *
#> 5 chr14 15038730 15050966 HdF000091 -1.250090      *
#> 6 chr14 15455652 15456778 HdF000106 -6.085891      *
```
