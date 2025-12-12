# Plot gene stats for chromosomes

***Plot*** gene stats for chromosomes.

## Usage

``` r
plot_gene_stats(
  gff_file,
  format = "auto",
  bar_width = 0.7,
  bar_color = "#0055ff55",
  lable_size = 3
)
```

## Arguments

- gff_file:

  Path to **`GFF3/GTF`** file as input.

- format:

  Format of GFF3/GTF file. (***"auto"***, "gff3", "gtf").

- bar_width:

  Bar width percent. (***0.7***).

- bar_color:

  Bar color with name or hex code. (***"#0055ff55"***).

- lable_size:

  Lable text size. (***3***).

## Value

A plot object of gene stats for chromosomes.

## Author

benben-miao

## Examples

``` r
# Example GFF3 file in GAnnoViz
gff_file <- system.file(
    "extdata",
    "example.gff",
    package = "GAnnoViz")

# Plot gene stats
plot_gene_stats(
    gff_file = gff_file,
    format = "auto",
    bar_width = 0.7,
    bar_color = "#0055ff55",
    lable_size = 3)
#> Import genomic features from the file as a GRanges object ... 
#> OK
#> Prepare the 'metadata' data frame ... 
#> OK
#> Make the TxDb object ... 
#> OK

```
