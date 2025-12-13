# SeuratVisPro Tutorial

## GAnnoViz

### 1. Introduction

**SourceCode:** <https://github.com/benben-miao/GAnnoViz/>

**Website API**: <https://benben-miao.github.io/GAnnoViz/>

GAnnoViz: a R package for genomic annotation and visualization.

#### Key features

- Comprehensive GFF/GTF parsing into tidy GRanges workflows
- Publication-grade genomic visualizations: gene structure, interval
  structure, and chromosome layouts
- Signal overlays for population and epigenetic analyses (FST, DMR) with
  robust window aggregation
- Innovative genome-wide density maps and circos-style rings for
  exploratory analysis
- Modular API with consistent theming for reproducible, high-quality
  figures

#### Typical inputs

- Structural annotations: GFF3/GTF
- Population genetics: sliding-window FST tables
- Epigenetics: methylKit-derived DMR tables
- Differential expression: DESeq2-style DEG tables

#### Design principles

- Consistent chromosome ordering and coordinate handling across
  functions
- Overlap-length weighted aggregation for window-based statistics
- Publication-ready defaults, minimal external dependencies,
  reproducible outputs

### 2. Installation

### 3. Shiny App

``` r

# GAnnoVizApp()
```

### 4. Extract Features

#### Extract Genes

``` r

# Example GFF3 file in GAnnoViz
gff_file <- system.file(
  "extdata",
  "example.gff",
  package = "GAnnoViz")

# Extract Genes
# genes <- extract_genes(
#   gff_file = gff_file,
#   format = "auto",
#   gene_info = "all")
# genes

# Gene info: gene_range
gene_range <- extract_genes(
  gff_file = gff_file,
  format = "auto",
  gene_info = "gene_range")
# head(gene_range)

head(gene_range) %>%
  kbl(
    caption = "Table. Gene range results", 
    escape = FALSE, 
    align = "c") %>%
  kable_classic(
    "striped", 
    full_width = TRUE, 
    html_font = "Times") %>%
  row_spec(0, bold = TRUE) %>%
  row_spec(1:5, extra_css = "line-height: 2em;")
```

|  start   |   end    | width |   names   |
|:--------:|:--------:|:-----:|:---------:|
| 10724003 | 10725460 | 1458  | HdF000001 |
| 10808174 | 10823308 | 15135 | HdF000002 |
| 10911391 | 10930352 | 18962 | HdF000003 |
| 10939872 | 10941326 | 1455  | HdF000004 |
| 10982117 | 10999259 | 17143 | HdF000005 |
| 11012358 | 11024009 | 11652 | HdF000006 |

Table. Gene range results {.table .lightable-classic .lightable-striped
style="font-family: Times; margin-left: auto; margin-right: auto;"}

#### Extract CDS

``` r

# Extract CDS
# cds <- extract_cds(
#   gff_file = gff_file,
#   format = "auto",
#   cds_info = "all")
# cds

# CDS info: cds_range
cds_range <- extract_cds(
  gff_file = gff_file,
  format = "auto",
  cds_info = "cds_range")
# head(cds_range)

head(cds_range) %>%
  kbl(
    caption = "Table. CDS range results", 
    escape = FALSE, 
    align = "c") %>%
  kable_classic(
    "striped", 
    full_width = TRUE, 
    html_font = "Times") %>%
  row_spec(0, bold = TRUE) %>%
  row_spec(1:5, extra_css = "line-height: 2em;")
```

| start  |  end   | width |
|:------:|:------:|:-----:|
| 810081 | 810090 |  10   |
| 816466 | 816734 |  269  |
| 818765 | 818971 |  207  |
| 819680 | 819740 |  61   |
| 821042 | 821214 |  173  |
| 825267 | 825329 |  63   |

Table. CDS range results {.table .lightable-classic .lightable-striped
style="font-family: Times; margin-left: auto; margin-right: auto;"}

#### Extract Promoters

``` r

# Extract Promoters
# promoters <- extract_promoters(
#   gff_file = gff_file,
#   format = "auto",
#   upstream = 2000,
#   downstream = 200,
#   promoter_info = "all")
# promoters

# Promoter info: promoter_range
promoter_range <- extract_promoters(
  gff_file = gff_file,
  format = "auto",
  upstream = 2000,
  downstream = 200,
  promoter_info = "promoter_range")
# head(promoter_range)

head(promoter_range) %>%
  kbl(
    caption = "Table. Promoter range results", 
    escape = FALSE, 
    align = "c") %>%
  kable_classic(
    "striped", 
    full_width = TRUE, 
    html_font = "Times") %>%
  row_spec(0, bold = TRUE) %>%
  row_spec(1:5, extra_css = "line-height: 2em;")
```

|  start  |   end   | width |     names     |
|:-------:|:-------:|:-----:|:-------------:|
| 807987  | 810186  | 2200  | nbis-mrna-126 |
| 882426  | 884625  | 2200  | nbis-mrna-127 |
| 949197  | 951396  | 2200  | nbis-mrna-129 |
| 1058746 | 1060945 | 2200  | nbis-mrna-130 |
| 1184302 | 1186501 | 2200  | nbis-mrna-132 |
| 1259096 | 1261295 | 2200  | nbis-mrna-134 |

Table. Promoter range results {.table .lightable-classic
.lightable-striped
style="font-family: Times; margin-left: auto; margin-right: auto;"}

#### Extract 5’UTR

``` r

# Extract 5'UTR
# utr5 <- extract_utr5(
#   gff_file = gff_file,
#   format = "auto",
#   utr5_info = "all")
# utr5

# 5'UTR info: utr5_range
utr5_range <- extract_utr5(
  gff_file = gff_file,
  format = "auto",
  utr5_info = "utr5_range")
# head(utr5_range)

head(utr5_range) %>%
  kbl(
    caption = "Table. 5'UTR range results", 
    escape = FALSE, 
    align = "c") %>%
  kable_classic(
    "striped", 
    full_width = TRUE, 
    html_font = "Times") %>%
  row_spec(0, bold = TRUE) %>%
  row_spec(1:5, extra_css = "line-height: 2em;")
```

|  start  |   end   | width | names |
|:-------:|:-------:|:-----:|:-----:|
| 809987  | 810080  |  94   |   1   |
| 884426  | 884645  |  220  |   2   |
| 951197  | 951558  |  362  |   3   |
| 952741  | 952756  |  16   |   3   |
| 1060746 | 1060848 |  103  |   4   |
| 1261096 | 1261128 |  33   |   6   |

Table. 5’UTR range results {.table .lightable-classic .lightable-striped
style="font-family: Times; margin-left: auto; margin-right: auto;"}

### 5. Plot Structure

#### Plot gene stats for chromosomes

``` r

# Plot gene stats
plot_gene_stats(
    gff_file = gff_file,
    format = "auto",
    bar_width = 0.7,
    bar_color = "#0055ff55",
    lable_size = 3)
```

![](Tutorial_files/figure-html/plot_gene_stats-1.png)

#### Plot gene structure (Promoter, 3’UTR, Exon, Intron, 5’UTR)

``` r

# Plot gene structure
plot_gene_structure(
  gff_file = gff_file,
  format = "auto",
  gene_id = "HdF029609",
  upstream = 2000,
  downstream = 200,
  feature_alpha = 0.8,
  intron_width = 1,
  x_breaks = 10,
  arrow_length = 5,
  arrow_count = 1,
  arrow_unit = "pt",
  promoter_color = "#ff8800",
  utr5_color = "#008833",
  utr3_color = "#ff0033",
  exon_color = "#0033ff",
  intron_color = "#333333"
)
```

![](Tutorial_files/figure-html/plot_gene_structure-1.png)

#### Plot gene structures for a genomic interval

``` r

# Plot interval structure
plot_interval_structure(
  gff_file = gff_file,
  format = "auto",
  chrom_id = "chr1",
  start = 950000,
  end = 1180000,
  x_breaks = 10,
  upstream = 2000,
  downstream = 200,
  feature_alpha = 0.8,
  intron_width = 1,
  arrow_count = 1,
  arrow_length = 5,
  arrow_unit = "pt",
  promoter_color = "#ff8800",
  utr5_color = "#008833",
  utr3_color = "#ff0033",
  exon_color = "#0033ff",
  intron_color = "#333333"
)
```

![](Tutorial_files/figure-html/plot_interval_structure-1.png)

#### Plot interval flank structure around a focal gene

``` r

# Neighborhood around a focal gene on its chromosome
plot_interval_flank(
  gff_file = gff_file,
  format = "auto",
  gene_id = "HdF029609",
  flank_upstream = 200000,
  flank_downstream = 200000,
  show_promoters = TRUE,
  upstream = 2000,
  downstream = 200,
  arrow_length = 5,
  arrow_unit = "pt",
  gene_color = "#0088ff",
  promoter_color = "#ff8800",
  label_size = 3
)
```

![](Tutorial_files/figure-html/plot_interval_flank-1.png)

#### Plot chromosome structures and gene stats

``` r

# Plot chrom structure
plot_chrom_structure(
  gff_file = gff_file,
  format = "auto",
  orientation = "vertical",
  bar_width = 0.6,
  chrom_alpha = 0.1,
  gene_width = 0.5,
  chrom_color = "#008888",
  gene_color = "#0088ff",
  telomere_color = "#ff0000",
  label_size = 3
)
```

![](Tutorial_files/figure-html/plot_chrom_structure-1.png)

#### Plot chromosome structures and gene annotation

``` r

genes <- data.frame(
  gene_id = c("HdF029609", "HdF029610"),
  gene_name = c("GeneA", "GeneB"))

# Vertical, annotate by name
plot_chrom_genes(
  gff_file = gff_file,
  gene_table = genes,
  annotate = "name",
  orientation = "vertical")
```

![](Tutorial_files/figure-html/plot_chrom_genes-1.png)

### 6. DEG Anno & Viz

#### Annotate differentially expressed genes (DEGs) with chromosome positions

``` r

# Example DEGs GAnnoViz
deg_file <- system.file(
  "extdata",
  "example.deg",
  package = "GAnnoViz")

deg <- read.table(
  file = deg_file,
  header = TRUE,
  sep = "\t",
  na.strings = NA,
  stringsAsFactors = FALSE
)
# head(deg)

head(deg) %>%
  kbl(
    caption = "Table. DEGs results from DESeq2", 
    escape = FALSE, 
    align = "c") %>%
  kable_classic(
    "striped", 
    full_width = TRUE, 
    html_font = "Times") %>%
  row_spec(0, bold = TRUE) %>%
  row_spec(1:5, extra_css = "line-height: 2em;")
```

|  GeneID   |  baseMean   | log2FoldChange |   lfcSE   |   stat    |  pvalue   |   padj    |
|:---------:|:-----------:|:--------------:|:---------:|:---------:|:---------:|:---------:|
| HdF054777 |  16.999581  |   -1.480620    | 0.7189789 | -2.059337 | 0.0394619 | 0.5479174 |
| HdF055254 |  10.012632  |    3.039624    | 1.3467639 | 2.256984  | 0.0240091 | 0.4341369 |
| HdF055463 |  3.926022   |    5.319215    | 2.2806921 | 2.332281  | 0.0196859 |    NA     |
| HdF027237 | 1433.976916 |   -1.802171    | 0.3336465 | -5.401438 | 0.0000001 | 0.0000283 |
| HdF027261 | 178.124271  |   -1.049145    | 0.4681808 | -2.240897 | 0.0250327 | 0.4406169 |
| HdF002748 | 185.015989  |    2.051539    | 0.9494449 | 2.160778  | 0.0307125 | 0.4865774 |

Table. DEGs results from DESeq2 {.table .lightable-classic
.lightable-striped
style="font-family: Times; margin-left: auto; margin-right: auto;"}

``` r

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
# head(res)

head(res) %>%
  kbl(
    caption = "Table. DEGs annotation results", 
    escape = FALSE, 
    align = "c") %>%
  kable_classic(
    "striped", 
    full_width = TRUE, 
    html_font = "Times") %>%
  row_spec(0, bold = TRUE) %>%
  row_spec(1:5, extra_css = "line-height: 2em;")
```

[TABLE]

Table. DEGs annotation results {.table .lightable-classic
.lightable-striped
style="font-family: Times; margin-left: auto; margin-right: auto;"}

#### Plot genomic feature density heatmap

``` r

# Gene density heatmap
plot_chrom_heatmap(
  gff_file = gff_file,
  format = "auto",
  feature = "gene",
  bin_size = 1e6,
  orientation = "horizontal",
  palette = c("#ffffff", "#0055aa"),
  alpha = 0.9
)
```

![](Tutorial_files/figure-html/plot_chrom_heatmap-1.png)

#### Plot differentially expressed genes (DEGs) hyper/hypo distributions by chromosome

``` r

# Plot chrom DEGs
plot_chrom_deg(
  deg_file = deg_file,
  gff_file = gff_file,
  format = "auto",
  id_col = "GeneID",
  fc_col = "log2FoldChange",
  violin_scale = "count",
  violin_border = 0.5,
  point_shape = 16,
  point_size = 2,
  jitter_width = 0.2,
  hyper_color = "#ff000088",
  hypo_color = "#00880088"
)
```

![](Tutorial_files/figure-html/plot_chrom_deg-1.png)

### 7. SNP Anno & Plot

#### Annotate FST slide windows with genomic features

``` r

# Annotate FST
fst_table <- system.file(
    "extdata",
    "example.fst",
    package = "GAnnoViz")

fst <- read.table(
  file = fst_table,
  header = TRUE,
  sep = "\t",
  na.strings = "NA",
  stringsAsFactors = FALSE
)
# head(fst)

head(fst) %>%
  kbl(
    caption = "Table. FST results from VCFtools", 
    escape = FALSE, 
    align = "c") %>%
  kable_classic(
    "striped", 
    full_width = TRUE, 
    html_font = "Times") %>%
  row_spec(0, bold = TRUE) %>%
  row_spec(1:5, extra_css = "line-height: 2em;")
```

| CHROM | BIN_START | BIN_END | N_VARIANTS | WEIGHTED_FST |  MEAN_FST  |
|:-----:|:---------:|:-------:|:----------:|:------------:|:----------:|
| chr11 |   10001   |  20000  |     7      |  0.0436644   | 0.0239114  |
| chr11 |   15001   |  25000  |     7      |  0.0436644   | 0.0239114  |
| chr11 |   25001   |  35000  |     4      |  0.0150391   | 0.0161973  |
| chr11 |   30001   |  40000  |     4      |  0.0150391   | 0.0161973  |
| chr11 |   55001   |  65000  |     1      |  -0.0044487  | -0.0044487 |
| chr11 |   60001   |  70000  |     1      |  -0.0044487  | -0.0044487 |

Table. FST results from VCFtools {.table .lightable-classic
.lightable-striped
style="font-family: Times; margin-left: auto; margin-right: auto;"}

``` r

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
# head(res)

head(res) %>%
  kbl(
    caption = "Table. FST annotation results", 
    escape = FALSE, 
    align = "c") %>%
  kable_classic(
    "striped", 
    full_width = TRUE, 
    html_font = "Times") %>%
  row_spec(0, bold = TRUE) %>%
  row_spec(1:5, extra_css = "line-height: 2em;")
```

| CHROM | BIN_START | BIN_END | N_VARIANTS | WEIGHTED_FST | MEAN_FST | anno_type | gene_id |
|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|
| chr11 | 10001 | 20000 | 7 | 0.0436644 | 0.0239114 | promoter,gene,UTR5,CDS,exon,intron,UTR3 | HdF041849(p),HdF041849(g) |
| chr11 | 15001 | 25000 | 7 | 0.0436644 | 0.0239114 | promoter,gene,UTR5,exon,intron | HdF041849(p),HdF041849(g) |
| chr11 | 25001 | 35000 | 4 | 0.0150391 | 0.0161973 | intergenic |  |
| chr11 | 30001 | 40000 | 4 | 0.0150391 | 0.0161973 | intergenic |  |
| chr11 | 55001 | 65000 | 1 | -0.0044487 | -0.0044487 | intergenic |  |
| chr11 | 60001 | 70000 | 1 | -0.0044487 | -0.0044487 | intergenic |  |

Table. FST annotation results {.table .lightable-classic
.lightable-striped
style="font-family: Times; margin-left: auto; margin-right: auto;"}

#### Annotate DMR slide windows with genomic features

``` r

# Annotate DMR
dmr_table <- system.file(
    "extdata",
    "example.dmr",
    package = "GAnnoViz")

dmr <- read.table(
  file = dmr_table,
  header = TRUE,
  sep = "\t",
  na.strings = "NA",
  stringsAsFactors = FALSE
)
# head(dmr)

head(dmr) %>%
  kbl(
    caption = "Table. DMR results from MethylKit", 
    escape = FALSE, 
    align = "c") %>%
  kable_classic(
    "striped", 
    full_width = TRUE, 
    html_font = "Times") %>%
  row_spec(0, bold = TRUE) %>%
  row_spec(1:5, extra_css = "line-height: 2em;")
```

[TABLE]

Table. DMR results from MethylKit {.table .lightable-classic
.lightable-striped
style="font-family: Times; margin-left: auto; margin-right: auto;"}

``` r

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
# head(res)

head(res) %>%
  kbl(
    caption = "Table. DMR annotation results", 
    escape = FALSE, 
    align = "c") %>%
  kable_classic(
    "striped", 
    full_width = TRUE, 
    html_font = "Times") %>%
  row_spec(0, bold = TRUE) %>%
  row_spec(1:5, extra_css = "line-height: 2em;")
```

[TABLE]

Table. DMR annotation results {.table .lightable-classic
.lightable-striped
style="font-family: Times; margin-left: auto; margin-right: auto;"}

#### Plot SNP density at chromosome level

``` r

# Plot SNP density
plot_snp_density(
  fst_file = fst_table,
  LOG10 = FALSE,
  bin_size = 1e6,
  density_color = c("#0088ff", "#ff8800", "#ff0000")
)
```

![](Tutorial_files/figure-html/plot_snp_density-1.png)

#### Plot genomic weighted FST heatmap

``` r

# Plot weighted FST
plot_snp_fst(
  fst_file = fst_table,
  bin_size = 1e6,
  metric = "fst_mean",
  orientation = "horizontal",
  palette = c("#ffffff", "#aa00aa"),
  alpha = 0.9
)
```

![](Tutorial_files/figure-html/plot_snp_fst-1.png)

#### Plot genomic FST with Top-N gene annotations

``` r

# Chromosome FST with Top-20 gene annotations on chr11
plot_snp_anno(
  fst_file = fst_table,
  gff_file = gff_file,
  format = "auto",
  chrom_id = "chr2",
  top_n = 20,
  orientation = "vertical",
  smooth_span = 0.5,
  fst_color = "#0088ff",
  point_size = 1,
  point_alpha = 0.3,
  label_size = 3,
  connector_dx1 = 2e4,
  connector_dx2 = 4e4,
  gap_frac = 0.05
)
```

![](Tutorial_files/figure-html/plot_snp_anno-1.png)

### 8. DMG Anno & Plot

#### Plot differentially methylated regions (DMRs) hyper/hypo distributions by chromosome

``` r

# Plot chrom DMRs
plot_dmg_chrom(
  dmr_file = dmr_table,
  violin_scale = "count",
  violin_border = 0.5,
  point_shape = 8,
  point_size = 2,
  jitter_width = 0.2,
  hyper_color = "#ff880088",
  hypo_color = "#0088ff88"
)
```

![](Tutorial_files/figure-html/plot_chrom_dmr-1.png)

#### Plot chromosomal DMGs trend

``` r

# Plot DMR trend
plot_dmg_trend(
  chrom_id = "chr1",
  dmr_file = dmr_table,
  smooth_span = 0.1,
  hyper_color = "#ff0000",
  hypo_color = "#008800",
  point_size = 3,
  point_alpha = 0.5
)
```

![](Tutorial_files/figure-html/plot_dmg_trend-1.png)
