# Extract mRNA transcripts from GFF or GTF

Extract ***mRNA*** transcripts from GFF or GTF.

## Usage

``` r
extract_mrnas(gff_file, format = "auto", mrna_info = "all")
```

## Arguments

- gff_file:

  Genomic structural annotation **`GFF3/GTF`** file path.

- format:

  Format of GFF3/GTF file. (***"auto"***, "gff3", "gtf").

- mrna_info:

  mRNA information. (***"all"***, "chrom_id", "mrna_id", "mrna_range").

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

# Extract mRNAs
mrnas <- extract_mrnas(
  gff_file = gff_file,
  format = "auto",
  mrna_info = "all")
#> Import genomic features from the file as a GRanges object ... 
#> OK
#> Prepare the 'metadata' data frame ... 
#> OK
#> Make the TxDb object ... 
#> OK
mrnas
#> GRanges object with 28263 ranges and 2 metadata columns:
#>             seqnames          ranges strand |     tx_id         tx_name
#>                <Rle>       <IRanges>  <Rle> | <integer>     <character>
#>       [1]       chr1   809987-851590      + |         1   nbis-mrna-126
#>       [2]       chr1   884426-897666      + |         2   nbis-mrna-127
#>       [3]       chr1  951197-1007073      + |         3   nbis-mrna-129
#>       [4]       chr1 1060746-1117817      + |         4   nbis-mrna-130
#>       [5]       chr1 1186302-1187372      + |         5   nbis-mrna-132
#>       ...        ...             ...    ... .       ...             ...
#>   [28259] superscaf8   226875-235847      + |     28259 nbis-mrna-11014
#>   [28260] superscaf8   310498-342926      + |     28260 nbis-mrna-11016
#>   [28261] superscaf8   367010-377320      + |     28261 nbis-mrna-11018
#>   [28262] superscaf8   240788-246722      - |     28262 nbis-mrna-11015
#>   [28263] superscaf8   351316-356949      - |     28263 nbis-mrna-11017
#>   -------
#>   seqinfo: 89 sequences from an unspecified genome; no seqlengths

# mRNA info: mrna_id
mrna_id <- extract_mrnas(
  gff_file = gff_file,
  format = "auto",
  mrna_info = "mrna_id")
#> Import genomic features from the file as a GRanges object ... 
#> OK
#> Prepare the 'metadata' data frame ... 
#> OK
#> Make the TxDb object ... 
#> OK
head(mrna_id)
#> [1] "nbis-mrna-126" "nbis-mrna-127" "nbis-mrna-129" "nbis-mrna-130"
#> [5] "nbis-mrna-132" "nbis-mrna-134"

# mRNA info: mrna_range
mrna_range <- extract_mrnas(
  gff_file = gff_file,
  format = "auto",
  mrna_info = "mrna_range")
#> Import genomic features from the file as a GRanges object ... 
#> OK
#> Prepare the 'metadata' data frame ... 
#> OK
#> Make the TxDb object ... 
#> OK
head(mrna_range)
#> IRanges object with 6 ranges and 0 metadata columns:
#>           start       end     width
#>       <integer> <integer> <integer>
#>   [1]    809987    851590     41604
#>   [2]    884426    897666     13241
#>   [3]    951197   1007073     55877
#>   [4]   1060746   1117817     57072
#>   [5]   1186302   1187372      1071
#>   [6]   1261096   1318472     57377
```
