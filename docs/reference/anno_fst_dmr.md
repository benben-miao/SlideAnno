# Annotate FST/DMR slide windows with genomic features

Annotate ***FST/DMR*** slide windows with ***genomic features***.

## Usage

``` r
anno_fst_dmr(
  gff_file,
  format = "auto",
  genomic_ranges,
  chrom_col = "CHROM",
  start_col = "BIN_START",
  end_col = "BIN_END",
  upstream = 2000,
  downstream = 200,
  ignore_strand = TRUE,
  features = c("promoter", "UTR5", "gene", "exon", "intron", "CDS", "UTR3", "intergenic")
)
```

## Arguments

- gff_file:

  Genomic structural annotation **`GFF3/GTF`** file path.

- format:

  Format of GFF3/GTF file. (***"auto"***, "gff3", "gtf").

- genomic_ranges:

  Genomic ranges file (e.g., FST, DMR). (***chromosome (prefix: chr),
  start, end***).

- chrom_col:

  Chromosome column name. (***FST: "CHROM", DMR: "chr"***).

- start_col:

  Start column name. (***FST: "BIN_START", DMR: "start"***).

- end_col:

  End column name. (***FST: "BIN_END", DMR: "end"***).

- upstream:

  Promoter upstream (bp). (***2000***).

- downstream:

  Promoter downstream (bp). (***200***).

- ignore_strand:

  Whether to ignore strand when computing overlaps. (***TRUE***).

- features:

  Which feature types to annotate. c(***"promoter", "UTR5", "gene",
  "exon", "intron", "CDS", "UTR3", "intergenic"***).

## Value

Data.frame with input intervals and ***annotation results***.

## Author

benben-miao

## Examples

``` r
# Example GFF3 file in GAnnoViz
gff_file <- system.file(
  "extdata",
  "example.gff",
  package = "GAnnoViz")

# Annotate FST
fst_table <- system.file(
    "extdata",
    "example.fst",
    package = "GAnnoViz")

fst <- read.table(
  file = fst_table,
  header = TRUE,
  sep = "\t",
  na.strings = NA,
  stringsAsFactors = FALSE
)
head(fst)
#>   CHROM BIN_START BIN_END N_VARIANTS WEIGHTED_FST    MEAN_FST
#> 1 chr11     10001   20000          7   0.04366440  0.02391140
#> 2 chr11     15001   25000          7   0.04366440  0.02391140
#> 3 chr11     25001   35000          4   0.01503910  0.01619730
#> 4 chr11     30001   40000          4   0.01503910  0.01619730
#> 5 chr11     55001   65000          1  -0.00444874 -0.00444874
#> 6 chr11     60001   70000          1  -0.00444874 -0.00444874

res <- anno_fst_dmr(
  gff_file = gff_file,
  format = "auto",
  genomic_ranges = fst_table,
  chrom_col = "CHROM",
  start_col = "BIN_START",
  end_col = "BIN_END",
  upstream = 2000,
  downstream = 200,
  ignore_strand = TRUE,
  features = c("promoter", "UTR5", "gene", "exon", "intron", "CDS", "UTR3", "intergenic")
)
#> Import genomic features from the file as a GRanges object ... 
#> OK
#> Prepare the 'metadata' data frame ... 
#> OK
#> Make the TxDb object ... 
#> OK
head(res)
#>   CHROM BIN_START BIN_END N_VARIANTS WEIGHTED_FST    MEAN_FST
#> 1 chr11     10001   20000          7   0.04366440  0.02391140
#> 2 chr11     15001   25000          7   0.04366440  0.02391140
#> 3 chr11     25001   35000          4   0.01503910  0.01619730
#> 4 chr11     30001   40000          4   0.01503910  0.01619730
#> 5 chr11     55001   65000          1  -0.00444874 -0.00444874
#> 6 chr11     60001   70000          1  -0.00444874 -0.00444874
#>                                 anno_type                   gene_id
#> 1 promoter,gene,UTR5,CDS,exon,intron,UTR3 HdF041849(p),HdF041849(g)
#> 2          promoter,gene,UTR5,exon,intron HdF041849(p),HdF041849(g)
#> 3                              intergenic                          
#> 4                              intergenic                          
#> 5                              intergenic                          
#> 6                              intergenic                          

# Annotate DMR
dmr_table <- system.file(
    "extdata",
    "example.dmr",
    package = "GAnnoViz")

dmr <- read.table(
  file = dmr_table,
  header = TRUE,
  sep = "\t",
  na.strings = NA,
  stringsAsFactors = FALSE
)
head(dmr)
#>    chr   start     end strand        pvalue        qvalue meth.diff
#> 1 chr1  466001  468000      *  4.157118e-03  1.206712e-02 -27.31707
#> 2 chr1  660001  662000      *  4.041133e-09  3.622605e-08 -33.37925
#> 3 chr1 1454001 1456000      *  8.939178e-12  1.083862e-10 -26.55631
#> 4 chr1 2750001 2752000      * 2.716109e-129 5.710428e-127  25.04745
#> 5 chr1 3428001 3430000      *  2.049072e-10  2.144328e-09 -46.33124
#> 6 chr1 3604001 3606000      *  2.137854e-09  1.984104e-08  37.85743

res <- anno_fst_dmr(
  gff_file = gff_file,
  format = "auto",
  genomic_ranges = dmr_table,
  chrom_col = "chr",
  start_col = "start",
  end_col = "end",
  upstream = 2000,
  downstream = 200,
  ignore_strand = TRUE,
  features = c("promoter", "UTR5", "gene", "exon", "intron", "CDS", "UTR3", "intergenic")
)
#> Import genomic features from the file as a GRanges object ... 
#> OK
#> Prepare the 'metadata' data frame ... 
#> OK
#> Make the TxDb object ... 
#> OK
#> Warning: Each of the 2 combined objects has sequence levels not in the other:
#>   - in 'x': h1tg000177l, h1tg000356l, h1tg000377l, h1tg000386l, h1tg000411l, h1tg000632l, h1tg000725l, h1tg000791l, h1tg000855l, h1tg000859l, superscaf9, h1tg000571c, h1tg000573l, h1tg000597l, h1tg000610l, h1tg000665l, h1tg000675l, h1tg000822l
#>   - in 'y': h1tg000310c, h1tg000374l, h1tg000441l, h1tg000520l, h1tg000541l, h1tg000549l, h1tg000554l, h1tg000609l, h1tg000618l, h1tg000620l, h1tg000625l, h1tg000633l, h1tg000686l, h1tg000739l, h1tg000742l, h1tg000750l, h1tg000756l, h1tg000761l, h1tg000769l, h1tg000774l, h1tg000782l, h1tg000790l, h1tg000797l, h1tg000798l, h1tg000799l, h1tg000831l, h1tg000840l, h1tg000847l, h1tg000852l, h1tg000856l, h1tg000864l, h1tg000873l, h1tg000875l, h1tg000878l, h1tg000890l, h1tg000906l, h1tg000912l, h1tg000913l, mith1tg000589c, superscaf10, superscaf13, superscaf27, superscaf28, superscaf29, superscaf3, superscaf31, superscaf32, superscaf33, superscaf35, superscaf37, superscaf38, superscaf40
#>   Make sure to always combine/compare objects based on the same reference
#>   genome (use suppressWarnings() to suppress this warning).
#> Warning: Each of the 2 combined objects has sequence levels not in the other:
#>   - in 'x': h1tg000177l, h1tg000356l, h1tg000377l, h1tg000386l, h1tg000411l, h1tg000632l, h1tg000725l, h1tg000791l, h1tg000855l, h1tg000859l, superscaf9, h1tg000571c, h1tg000573l, h1tg000597l, h1tg000610l, h1tg000665l, h1tg000675l, h1tg000822l
#>   - in 'y': h1tg000310c, h1tg000374l, h1tg000441l, h1tg000520l, h1tg000541l, h1tg000549l, h1tg000554l, h1tg000609l, h1tg000618l, h1tg000620l, h1tg000625l, h1tg000633l, h1tg000686l, h1tg000739l, h1tg000742l, h1tg000750l, h1tg000756l, h1tg000761l, h1tg000769l, h1tg000774l, h1tg000782l, h1tg000790l, h1tg000797l, h1tg000798l, h1tg000799l, h1tg000831l, h1tg000840l, h1tg000847l, h1tg000852l, h1tg000856l, h1tg000864l, h1tg000873l, h1tg000875l, h1tg000878l, h1tg000890l, h1tg000906l, h1tg000912l, h1tg000913l, mith1tg000589c, superscaf10, superscaf13, superscaf27, superscaf28, superscaf29, superscaf3, superscaf31, superscaf32, superscaf33, superscaf35, superscaf37, superscaf38, superscaf40
#>   Make sure to always combine/compare objects based on the same reference
#>   genome (use suppressWarnings() to suppress this warning).
#> Warning: Each of the 2 combined objects has sequence levels not in the other:
#>   - in 'x': h1tg000177l, h1tg000356l, h1tg000377l, h1tg000386l, h1tg000411l, h1tg000632l, h1tg000725l, h1tg000791l, h1tg000855l, h1tg000859l, superscaf9, h1tg000571c, h1tg000573l, h1tg000597l, h1tg000610l, h1tg000665l, h1tg000675l, h1tg000822l
#>   - in 'y': h1tg000310c, h1tg000374l, h1tg000441l, h1tg000520l, h1tg000541l, h1tg000549l, h1tg000554l, h1tg000609l, h1tg000618l, h1tg000620l, h1tg000625l, h1tg000633l, h1tg000686l, h1tg000739l, h1tg000742l, h1tg000750l, h1tg000756l, h1tg000761l, h1tg000769l, h1tg000774l, h1tg000782l, h1tg000790l, h1tg000797l, h1tg000798l, h1tg000799l, h1tg000831l, h1tg000840l, h1tg000847l, h1tg000852l, h1tg000856l, h1tg000864l, h1tg000873l, h1tg000875l, h1tg000878l, h1tg000890l, h1tg000906l, h1tg000912l, h1tg000913l, mith1tg000589c, superscaf10, superscaf13, superscaf27, superscaf28, superscaf29, superscaf3, superscaf31, superscaf32, superscaf33, superscaf35, superscaf37, superscaf38, superscaf40
#>   Make sure to always combine/compare objects based on the same reference
#>   genome (use suppressWarnings() to suppress this warning).
#> Warning: Each of the 2 combined objects has sequence levels not in the other:
#>   - in 'x': h1tg000177l, h1tg000356l, h1tg000377l, h1tg000386l, h1tg000411l, h1tg000632l, h1tg000725l, h1tg000791l, h1tg000855l, h1tg000859l, superscaf9, h1tg000571c, h1tg000573l, h1tg000597l, h1tg000610l, h1tg000665l, h1tg000675l, h1tg000822l
#>   - in 'y': h1tg000310c, h1tg000374l, h1tg000441l, h1tg000520l, h1tg000541l, h1tg000549l, h1tg000554l, h1tg000609l, h1tg000618l, h1tg000620l, h1tg000625l, h1tg000633l, h1tg000686l, h1tg000739l, h1tg000742l, h1tg000750l, h1tg000756l, h1tg000761l, h1tg000769l, h1tg000774l, h1tg000782l, h1tg000790l, h1tg000797l, h1tg000798l, h1tg000799l, h1tg000831l, h1tg000840l, h1tg000847l, h1tg000852l, h1tg000856l, h1tg000864l, h1tg000873l, h1tg000875l, h1tg000878l, h1tg000890l, h1tg000906l, h1tg000912l, h1tg000913l, mith1tg000589c, superscaf10, superscaf13, superscaf27, superscaf28, superscaf29, superscaf3, superscaf31, superscaf32, superscaf33, superscaf35, superscaf37, superscaf38, superscaf40
#>   Make sure to always combine/compare objects based on the same reference
#>   genome (use suppressWarnings() to suppress this warning).
#> Warning: Each of the 2 combined objects has sequence levels not in the other:
#>   - in 'x': h1tg000177l, h1tg000356l, h1tg000377l, h1tg000386l, h1tg000411l, h1tg000632l, h1tg000725l, h1tg000791l, h1tg000855l, h1tg000859l, superscaf9, h1tg000571c, h1tg000573l, h1tg000597l, h1tg000610l, h1tg000665l, h1tg000675l, h1tg000822l
#>   - in 'y': h1tg000310c, h1tg000374l, h1tg000441l, h1tg000520l, h1tg000541l, h1tg000549l, h1tg000554l, h1tg000609l, h1tg000618l, h1tg000620l, h1tg000625l, h1tg000633l, h1tg000686l, h1tg000739l, h1tg000742l, h1tg000750l, h1tg000756l, h1tg000761l, h1tg000769l, h1tg000774l, h1tg000782l, h1tg000790l, h1tg000797l, h1tg000798l, h1tg000799l, h1tg000831l, h1tg000840l, h1tg000847l, h1tg000852l, h1tg000856l, h1tg000864l, h1tg000873l, h1tg000875l, h1tg000878l, h1tg000890l, h1tg000906l, h1tg000912l, h1tg000913l, mith1tg000589c, superscaf10, superscaf13, superscaf27, superscaf28, superscaf29, superscaf3, superscaf31, superscaf32, superscaf33, superscaf35, superscaf37, superscaf38, superscaf40
#>   Make sure to always combine/compare objects based on the same reference
#>   genome (use suppressWarnings() to suppress this warning).
#> Warning: Each of the 2 combined objects has sequence levels not in the other:
#>   - in 'x': h1tg000177l, h1tg000356l, h1tg000377l, h1tg000386l, h1tg000411l, h1tg000632l, h1tg000725l, h1tg000791l, h1tg000855l, h1tg000859l, superscaf9, h1tg000571c, h1tg000573l, h1tg000597l, h1tg000610l, h1tg000665l, h1tg000675l, h1tg000822l
#>   - in 'y': h1tg000310c, h1tg000374l, h1tg000441l, h1tg000520l, h1tg000541l, h1tg000549l, h1tg000554l, h1tg000609l, h1tg000618l, h1tg000620l, h1tg000625l, h1tg000633l, h1tg000686l, h1tg000739l, h1tg000742l, h1tg000750l, h1tg000756l, h1tg000761l, h1tg000769l, h1tg000774l, h1tg000782l, h1tg000790l, h1tg000797l, h1tg000798l, h1tg000799l, h1tg000831l, h1tg000840l, h1tg000847l, h1tg000852l, h1tg000856l, h1tg000864l, h1tg000873l, h1tg000875l, h1tg000878l, h1tg000890l, h1tg000906l, h1tg000912l, h1tg000913l, mith1tg000589c, superscaf10, superscaf13, superscaf27, superscaf28, superscaf29, superscaf3, superscaf31, superscaf32, superscaf33, superscaf35, superscaf37, superscaf38, superscaf40
#>   Make sure to always combine/compare objects based on the same reference
#>   genome (use suppressWarnings() to suppress this warning).
#> Warning: Each of the 2 combined objects has sequence levels not in the other:
#>   - in 'x': h1tg000177l, h1tg000356l, h1tg000377l, h1tg000386l, h1tg000411l, h1tg000632l, h1tg000725l, h1tg000791l, h1tg000855l, h1tg000859l, superscaf9, h1tg000571c, h1tg000573l, h1tg000597l, h1tg000610l, h1tg000665l, h1tg000675l, h1tg000822l
#>   - in 'y': h1tg000310c, h1tg000374l, h1tg000441l, h1tg000520l, h1tg000541l, h1tg000549l, h1tg000554l, h1tg000609l, h1tg000618l, h1tg000620l, h1tg000625l, h1tg000633l, h1tg000686l, h1tg000739l, h1tg000742l, h1tg000750l, h1tg000756l, h1tg000761l, h1tg000769l, h1tg000774l, h1tg000782l, h1tg000790l, h1tg000797l, h1tg000798l, h1tg000799l, h1tg000831l, h1tg000840l, h1tg000847l, h1tg000852l, h1tg000856l, h1tg000864l, h1tg000873l, h1tg000875l, h1tg000878l, h1tg000890l, h1tg000906l, h1tg000912l, h1tg000913l, mith1tg000589c, superscaf10, superscaf13, superscaf27, superscaf28, superscaf29, superscaf3, superscaf31, superscaf32, superscaf33, superscaf35, superscaf37, superscaf38, superscaf40
#>   Make sure to always combine/compare objects based on the same reference
#>   genome (use suppressWarnings() to suppress this warning).
head(res)
#>    chr   start     end strand        pvalue        qvalue meth.diff
#> 1 chr1  466001  468000      *  4.157118e-03  1.206712e-02 -27.31707
#> 2 chr1  660001  662000      *  4.041133e-09  3.622605e-08 -33.37925
#> 3 chr1 1454001 1456000      *  8.939178e-12  1.083862e-10 -26.55631
#> 4 chr1 2750001 2752000      * 2.716109e-129 5.710428e-127  25.04745
#> 5 chr1 3428001 3430000      *  2.049072e-10  2.144328e-09 -46.33124
#> 6 chr1 3604001 3606000      *  2.137854e-09  1.984104e-08  37.85743
#>        anno_type      gene_id
#> 1    gene,intron HdF051575(g)
#> 2    gene,intron HdF051575(g)
#> 3 gene,exon,UTR3 HdF027216(g)
#> 4     intergenic             
#> 5    gene,intron HdF027234(g)
#> 6    gene,intron HdF027242(g)
```
