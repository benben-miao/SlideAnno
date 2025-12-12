# Extract 5'UTR ranges from GFF or GTF

Extract ***5'UTR*** ranges from GFF or GTF.

## Usage

``` r
extract_utr5(gff_file, format = "auto", utr5_info = "all")
```

## Arguments

- gff_file:

  Genomic structural annotation **`GFF3/GTF`** file path.

- format:

  Format of GFF3/GTF file. (***"auto"***, "gff3", "gtf").

- utr5_info:

  mRNA information. (***"all"***, "chrom_id", "utr5_id", "utr5_range").

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

# Extract 5'UTR
utr5 <- extract_utr5(
  gff_file = gff_file,
  format = "auto",
  utr5_info = "all")
#> Import genomic features from the file as a GRanges object ... 
#> OK
#> Prepare the 'metadata' data frame ... 
#> OK
#> Make the TxDb object ... 
#> OK
utr5
#> GRanges object with 40374 ranges and 3 metadata columns:
#>           seqnames          ranges strand |   exon_id              exon_name
#>              <Rle>       <IRanges>  <Rle> | <integer>            <character>
#>       1       chr1   809987-810080      + |         1 HdF051577.1.5utr.exo..
#>       2       chr1   884426-884645      + |        10 HdF051578.1.5utr.exo..
#>       3       chr1   951197-951558      + |        16 HdF051580.1.5utr.exo..
#>       3       chr1   952741-952756      + |        17 HdF051580.2.5utr.exo..
#>       4       chr1 1060746-1060848      + |        48 HdF027208.1.5utr.exo..
#>     ...        ...             ...    ... .       ...                    ...
#>   28259 superscaf8   227656-227795      + |    200841 HdF055316.2.5utr.exo..
#>   28260 superscaf8   310498-310703      + |    200843 HdF055188.1.5utr.exo..
#>   28260 superscaf8   339024-339144      + |    200844 HdF055188.2.5utr.exo..
#>   28262 superscaf8   246596-246722      - |    200857 HdF055315.1.3utr.exo..
#>   28263 superscaf8   356850-356949      - |    200859 HdF055189.1.3utr.exo..
#>         exon_rank
#>         <integer>
#>       1         1
#>       2         1
#>       3         1
#>       3         2
#>       4         1
#>     ...       ...
#>   28259         2
#>   28260         1
#>   28260         2
#>   28262         1
#>   28263         1
#>   -------
#>   seqinfo: 89 sequences from an unspecified genome; no seqlengths

# 5'UTR info: utr5_id
utr5_id <- extract_utr5(
  gff_file = gff_file,
  format = "auto",
  utr5_info = "utr5_id")
#> Import genomic features from the file as a GRanges object ... 
#> OK
#> Prepare the 'metadata' data frame ... 
#> OK
#> Make the TxDb object ... 
#> OK
head(utr5_id)
#> [1] "HdF051577.1.5utr.exon.1" "HdF051578.1.5utr.exon.1"
#> [3] "HdF051580.1.5utr.exon.1" "HdF051580.2.5utr.exon.2"
#> [5] "HdF027208.1.5utr.exon.1" "HdF027212.1.5utr.exon.1"

# 5'UTR info: utr5_range
utr5_range <- extract_utr5(
  gff_file = gff_file,
  format = "auto",
  utr5_info = "utr5_range")
#> Import genomic features from the file as a GRanges object ... 
#> OK
#> Prepare the 'metadata' data frame ... 
#> OK
#> Make the TxDb object ... 
#> OK
head(utr5_range)
#> IRanges object with 6 ranges and 0 metadata columns:
#>         start       end     width
#>     <integer> <integer> <integer>
#>   1    809987    810080        94
#>   2    884426    884645       220
#>   3    951197    951558       362
#>   3    952741    952756        16
#>   4   1060746   1060848       103
#>   6   1261096   1261128        33
```
