# Extract 3'UTR ranges from GFF or GTF

Extract ***3'UTR*** ranges from GFF or GTF.

## Usage

``` r
extract_utr3(gff_file, format = "auto", utr3_info = "all")
```

## Arguments

- gff_file:

  Genomic structural annotation **`GFF3/GTF`** file path.

- format:

  Format of GFF3/GTF file. (***"auto"***, "gff3", "gtf").

- utr3_info:

  mRNA information. (***"all"***, "chrom_id", "utr3_id", "utr3_range").

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

# Extract 3'UTR
utr3 <- extract_utr3(
  gff_file = gff_file,
  format = "auto",
  utr3_info = "all")
#> Import genomic features from the file as a GRanges object ... 
#> OK
#> Prepare the 'metadata' data frame ... 
#> OK
#> Make the TxDb object ... 
#> OK
utr3
#> GRanges object with 25405 ranges and 3 metadata columns:
#>           seqnames          ranges strand |   exon_id              exon_name
#>              <Rle>       <IRanges>  <Rle> | <integer>            <character>
#>       1       chr1   850987-851590      + |         9 HdF051577.8.cds.exon.9
#>       2       chr1   897489-897666      + |        15 HdF051578.5.cds.exon.6
#>       3       chr1 1006434-1007073      + |        47 HdF051580.30.cds.exo..
#>       4       chr1 1115499-1117817      + |        79 HdF027208.31.cds.exo..
#>       6       chr1 1317155-1318472      + |       106 HdF027212.25.cds.exo..
#>     ...        ...             ...    ... .       ...                    ...
#>   28262 superscaf8   245265-246043      - |    200857 HdF055315.1.3utr.exo..
#>   28262 superscaf8   243394-243450      - |    200856 HdF055315.2.3utr.exo..
#>   28262 superscaf8   240788-242074      - |    200855 HdF055315.3.3utr.exo..
#>   28263 superscaf8   355498-356309      - |    200859 HdF055189.1.3utr.exo..
#>   28263 superscaf8   351316-352212      - |    200858 HdF055189.2.3utr.exo..
#>         exon_rank
#>         <integer>
#>       1         9
#>       2         6
#>       3        32
#>       4        32
#>       6        26
#>     ...       ...
#>   28262         1
#>   28262         2
#>   28262         3
#>   28263         1
#>   28263         2
#>   -------
#>   seqinfo: 89 sequences from an unspecified genome; no seqlengths

# 3'UTR info: utr3_id
utr3_id <- extract_utr3(
  gff_file = gff_file,
  format = "auto",
  utr3_info = "utr3_id")
#> Import genomic features from the file as a GRanges object ... 
#> OK
#> Prepare the 'metadata' data frame ... 
#> OK
#> Make the TxDb object ... 
#> OK
head(utr3_id)
#> [1] "HdF051577.8.cds.exon.9"   "HdF051578.5.cds.exon.6"  
#> [3] "HdF051580.30.cds.exon.32" "HdF027208.31.cds.exon.32"
#> [5] "HdF027212.25.cds.exon.26" "HdF027218.4.cds.exon.4"  

# 3'UTR info: utr3_range
utr3_range <- extract_utr3(
  gff_file = gff_file,
  format = "auto",
  utr3_info = "utr3_range")
#> Import genomic features from the file as a GRanges object ... 
#> OK
#> Prepare the 'metadata' data frame ... 
#> OK
#> Make the TxDb object ... 
#> OK
head(utr3_range)
#> IRanges object with 6 ranges and 0 metadata columns:
#>         start       end     width
#>     <integer> <integer> <integer>
#>   1    850987    851590       604
#>   2    897489    897666       178
#>   3   1006434   1007073       640
#>   4   1115499   1117817      2319
#>   6   1317155   1318472      1318
#>   8   1875146   1875151         6
```
