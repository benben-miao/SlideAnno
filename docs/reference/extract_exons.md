# Extract Exons ranges from GFF or GTF

Extract ***Exons*** ranges from GFF or GTF.

## Usage

``` r
extract_exons(gff_file, format = "auto", exon_info = "all")
```

## Arguments

- gff_file:

  Genomic structural annotation **`GFF3/GTF`** file path.

- format:

  Format of GFF3/GTF file. (***"auto"***, "gff3", "gtf").

- exon_info:

  Exon information. (***"all"***, "chrom_id", "exon_range").

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

# Extract Exons
exons <- extract_exons(
  gff_file = gff_file,
  format = "auto",
  exon_info = "all")
#> Import genomic features from the file as a GRanges object ... 
#> OK
#> Prepare the 'metadata' data frame ... 
#> OK
#> Make the TxDb object ... 
#> OK
exons
#> GRanges object with 200859 ranges and 1 metadata column:
#>              seqnames        ranges strand |   exon_id
#>                 <Rle>     <IRanges>  <Rle> | <integer>
#>        [1]       chr1 809987-810090      + |         1
#>        [2]       chr1 816466-816734      + |         2
#>        [3]       chr1 818765-818971      + |         3
#>        [4]       chr1 819680-819740      + |         4
#>        [5]       chr1 821042-821214      + |         5
#>        ...        ...           ...    ... .       ...
#>   [200855] superscaf8 240788-242074      - |    200855
#>   [200856] superscaf8 243394-243450      - |    200856
#>   [200857] superscaf8 245265-246722      - |    200857
#>   [200858] superscaf8 351316-352212      - |    200858
#>   [200859] superscaf8 355498-356949      - |    200859
#>   -------
#>   seqinfo: 89 sequences from an unspecified genome; no seqlengths

# Exon info: exon_range
exon_range <- extract_exons(
  gff_file = gff_file,
  format = "auto",
  exon_info = "exon_range")
#> Import genomic features from the file as a GRanges object ... 
#> OK
#> Prepare the 'metadata' data frame ... 
#> OK
#> Make the TxDb object ... 
#> OK
head(exon_range)
#> IRanges object with 6 ranges and 0 metadata columns:
#>           start       end     width
#>       <integer> <integer> <integer>
#>   [1]    809987    810090       104
#>   [2]    816466    816734       269
#>   [3]    818765    818971       207
#>   [4]    819680    819740        61
#>   [5]    821042    821214       173
#>   [6]    825267    825329        63
```
