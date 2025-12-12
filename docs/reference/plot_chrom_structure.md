# Plot chromosome structures and gene stats

Plot ***chromosome structures and gene stats***.

## Usage

``` r
plot_chrom_structure(
  gff_file,
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

## Arguments

- gff_file:

  Genomic structural annotation **`GFF3/GTF`** file path.

- format:

  Format of GFF3/GTF file. (***"auto"***, "gff3", "gtf").

- orientation:

  Coordinate orientation: (***"vertical"***, "horizontal").

- bar_width:

  Chromosome bars relative width. (***0.6***).

- chrom_alpha:

  Chromosome bars alpha. (***0.2***).

- gene_width:

  Gene bar relative width. (***0.5***).

- chrom_color:

  Chromosome bar color. (***"#008888"***).

- gene_color:

  Gene rectangle color. (***"#0088ff"***).

- telomere_color:

  Telomere color. (***"#ff0000"***).

- label_size:

  Label text size. (***3***).

## Value

Plot of chromosome structure.

## Examples

``` r
# Example GFF3 file in GAnnoViz
gff_file <- system.file(
  "extdata",
  "example.gff",
  package = "GAnnoViz")

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
#> Import genomic features from the file as a GRanges object ... 
#> OK
#> Prepare the 'metadata' data frame ... 
#> OK
#> Make the TxDb object ... 
#> OK

```
