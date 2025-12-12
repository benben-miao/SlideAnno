# Extract promoter ranges from GFF or GTF

Extract ***promoter*** ranges from GFF or GTF.

## Usage

``` r
extract_promoters(
  gff_file,
  format = "auto",
  upstream = 2000,
  downstream = 200,
  promoter_info = "all"
)
```

## Arguments

- gff_file:

  Genomic structural annotation **`GFF3/GTF`** file path.

- format:

  Format of GFF3/GTF file. (***"auto"***, "gff3", "gtf").

- upstream:

  Promoter upstream (bp). (***2000***).

- downstream:

  Promoter downstream (bp). (***200***).

- promoter_info:

  Promoter information. (***"all"***, "chrom_id", "promoter_id",
  "promoter_range").

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

# Extract Promoters
promoters <- extract_promoters(
  gff_file = gff_file,
  format = "auto",
  upstream = 2000,
  downstream = 200,
  promoter_info = "all")
#> Import genomic features from the file as a GRanges object ... 
#> OK
#> Prepare the 'metadata' data frame ... 
#> OK
#> Make the TxDb object ... 
#> OK
promoters
#> GRanges object with 28263 ranges and 2 metadata columns:
#>                     seqnames          ranges strand |     tx_id         tx_name
#>                        <Rle>       <IRanges>  <Rle> | <integer>     <character>
#>     nbis-mrna-126       chr1   807987-810186      + |         1   nbis-mrna-126
#>     nbis-mrna-127       chr1   882426-884625      + |         2   nbis-mrna-127
#>     nbis-mrna-129       chr1   949197-951396      + |         3   nbis-mrna-129
#>     nbis-mrna-130       chr1 1058746-1060945      + |         4   nbis-mrna-130
#>     nbis-mrna-132       chr1 1184302-1186501      + |         5   nbis-mrna-132
#>               ...        ...             ...    ... .       ...             ...
#>   nbis-mrna-11014 superscaf8   224875-227074      + |     28259 nbis-mrna-11014
#>   nbis-mrna-11016 superscaf8   308498-310697      + |     28260 nbis-mrna-11016
#>   nbis-mrna-11018 superscaf8   365010-367209      + |     28261 nbis-mrna-11018
#>   nbis-mrna-11015 superscaf8   246523-248722      - |     28262 nbis-mrna-11015
#>   nbis-mrna-11017 superscaf8   356750-358949      - |     28263 nbis-mrna-11017
#>   -------
#>   seqinfo: 89 sequences from an unspecified genome; no seqlengths

# Promoter info: promoter_id
promoter_id <- extract_promoters(
  gff_file = gff_file,
  format = "auto",
  upstream = 2000,
  downstream = 200,
  promoter_info = "promoter_id")
#> Import genomic features from the file as a GRanges object ... 
#> OK
#> Prepare the 'metadata' data frame ... 
#> OK
#> Make the TxDb object ... 
#> OK
head(promoter_id)
#> [1] "nbis-mrna-126" "nbis-mrna-127" "nbis-mrna-129" "nbis-mrna-130"
#> [5] "nbis-mrna-132" "nbis-mrna-134"

# Promoter info: promoter_range
promoter_range <- extract_promoters(
  gff_file = gff_file,
  format = "auto",
  upstream = 2000,
  downstream = 200,
  promoter_info = "promoter_range")
#> Import genomic features from the file as a GRanges object ... 
#> OK
#> Prepare the 'metadata' data frame ... 
#> OK
#> Make the TxDb object ... 
#> OK
head(promoter_range)
#> IRanges object with 6 ranges and 0 metadata columns:
#>                     start       end     width
#>                 <integer> <integer> <integer>
#>   nbis-mrna-126    807987    810186      2200
#>   nbis-mrna-127    882426    884625      2200
#>   nbis-mrna-129    949197    951396      2200
#>   nbis-mrna-130   1058746   1060945      2200
#>   nbis-mrna-132   1184302   1186501      2200
#>   nbis-mrna-134   1259096   1261295      2200
```
