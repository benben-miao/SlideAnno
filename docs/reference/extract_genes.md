# Extract genes information from GFF or GTF

Extract ***genes*** information from GFF or GTF.

## Usage

``` r
extract_genes(gff_file, format = "auto", gene_info = "all")
```

## Arguments

- gff_file:

  Genomic structural annotation **`GFF3/GTF`** file path.

- format:

  Format of GFF3/GTF file. (***"auto"***, "gff3", "gtf").

- gene_info:

  Gene information. (***"all"***, "chrom_id", "gene_id", "gene_range").

## Value

A ***GRanges*** object.

## Author

benben-miao

## Examples

``` r
# Example GFF3 file in GAnnoViz
gff_file <- system.file(
  "extdata",
  "example.gff",
  package = "GAnnoViz")

# Extract Genes
genes <- extract_genes(
  gff_file = gff_file,
  format = "auto",
  gene_info = "all")
#> Import genomic features from the file as a GRanges object ... 
#> OK
#> Prepare the 'metadata' data frame ... 
#> OK
#> Make the TxDb object ... 
#> OK
genes
#> GRanges object with 28263 ranges and 1 metadata column:
#>                    seqnames            ranges strand |     gene_id
#>                       <Rle>         <IRanges>  <Rle> | <character>
#>    HdF000001          chr14 10724003-10725460      + |   HdF000001
#>    HdF000002          chr14 10808174-10823308      - |   HdF000002
#>    HdF000003          chr14 10911391-10930352      + |   HdF000003
#>    HdF000004          chr14 10939872-10941326      - |   HdF000004
#>    HdF000005          chr14 10982117-10999259      - |   HdF000005
#>          ...            ...               ...    ... .         ...
#>    ggene-ND3 mith1tg000589c         1578-1931      + |   ggene-ND3
#>    ggene-ND4 mith1tg000589c         8450-9742      - |   ggene-ND4
#>   ggene-ND4L mith1tg000589c        9838-10137      - |  ggene-ND4L
#>    ggene-ND5 mith1tg000589c         6599-8341      - |   ggene-ND5
#>    ggene-ND6 mith1tg000589c       11509-12015      - |   ggene-ND6
#>   -------
#>   seqinfo: 89 sequences from an unspecified genome; no seqlengths

# Gene info: gene_id
gene_id <- extract_genes(
  gff_file = gff_file,
  format = "auto",
  gene_info = "gene_id")
#> Import genomic features from the file as a GRanges object ... 
#> OK
#> Prepare the 'metadata' data frame ... 
#> OK
#> Make the TxDb object ... 
#> OK
head(gene_id)
#> [1] "HdF000001" "HdF000002" "HdF000003" "HdF000004" "HdF000005" "HdF000006"

# Gene info: gene_range
gene_range <- extract_genes(
  gff_file = gff_file,
  format = "auto",
  gene_info = "gene_range")
#> Import genomic features from the file as a GRanges object ... 
#> OK
#> Prepare the 'metadata' data frame ... 
#> OK
#> Make the TxDb object ... 
#> OK
head(gene_range)
#> IRanges object with 6 ranges and 0 metadata columns:
#>                 start       end     width
#>             <integer> <integer> <integer>
#>   HdF000001  10724003  10725460      1458
#>   HdF000002  10808174  10823308     15135
#>   HdF000003  10911391  10930352     18962
#>   HdF000004  10939872  10941326      1455
#>   HdF000005  10982117  10999259     17143
#>   HdF000006  11012358  11024009     11652
```
