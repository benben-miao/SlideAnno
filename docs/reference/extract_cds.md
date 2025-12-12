# Extract CDS ranges from GFF or GTF

Extract ***CDS*** ranges from GFF or GTF.

## Usage

``` r
extract_cds(gff_file, format = "auto", cds_info = "all")
```

## Arguments

- gff_file:

  Genomic structural annotation **`GFF3/GTF`** file path.

- format:

  Format of GFF3/GTF file. (***"auto"***, "gff3", "gtf").

- cds_info:

  CDS information. (***"all"***, "chrom_id", "cds_range").

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

# Extract CDS
cds <- extract_cds(
  gff_file = gff_file,
  format = "auto",
  cds_info = "all")
#> Import genomic features from the file as a GRanges object ... 
#> OK
#> Prepare the 'metadata' data frame ... 
#> OK
#> Make the TxDb object ... 
#> OK
cds
#> GRanges object with 174657 ranges and 1 metadata column:
#>              seqnames        ranges strand |    cds_id
#>                 <Rle>     <IRanges>  <Rle> | <integer>
#>        [1]       chr1 810081-810090      + |         1
#>        [2]       chr1 816466-816734      + |         2
#>        [3]       chr1 818765-818971      + |         3
#>        [4]       chr1 819680-819740      + |         4
#>        [5]       chr1 821042-821214      + |         5
#>        ...        ...           ...    ... .       ...
#>   [174653] superscaf8 374132-374321      + |    174653
#>   [174654] superscaf8 375441-375559      + |    174654
#>   [174655] superscaf8 377150-377320      + |    174655
#>   [174656] superscaf8 246044-246595      - |    174656
#>   [174657] superscaf8 356310-356849      - |    174657
#>   -------
#>   seqinfo: 89 sequences from an unspecified genome; no seqlengths

# CDS info: cds_range
cds_range <- extract_cds(
  gff_file = gff_file,
  format = "auto",
  cds_info = "cds_range")
#> Import genomic features from the file as a GRanges object ... 
#> OK
#> Prepare the 'metadata' data frame ... 
#> OK
#> Make the TxDb object ... 
#> OK
head(cds_range)
#> IRanges object with 6 ranges and 0 metadata columns:
#>           start       end     width
#>       <integer> <integer> <integer>
#>   [1]    810081    810090        10
#>   [2]    816466    816734       269
#>   [3]    818765    818971       207
#>   [4]    819680    819740        61
#>   [5]    821042    821214       173
#>   [6]    825267    825329        63
```
