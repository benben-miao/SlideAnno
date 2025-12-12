# Plot differentially methylated regions (DMRs) hyper/hypo distributions by chromosome

Plot ***differentially methylated regions (DMRs)*** hyper/hypo
distributions by chromosome.

## Usage

``` r
plot_dmg_chrom(
  dmr_file,
  violin_scale = "count",
  violin_border = 0.5,
  point_shape = 8,
  point_size = 2,
  jitter_width = 0.2,
  hyper_color = "#ff880088",
  hypo_color = "#0088ff88"
)
```

## Arguments

- dmr_file:

  DEG table from ***MethylKit*** analysis.

- violin_scale:

  Violin scale mode. (***"count"***, "area", "width").

- violin_border:

  Violin border width. (***0.5***).

- point_shape:

  Points shape (0-25). (***8***).

- point_size:

  Point size. (***2***).

- jitter_width:

  Horizontal jitter width. (***0.2***).

- hyper_color:

  Color for hyper-methylated points. (***"#ff000088"***).

- hypo_color:

  Color for hypo-methylated points. (***"#00880088"***).

## Value

A ***ggplot object*** of chromosome-wise DMR distributions.

## Examples

``` r
# DMR results from Methylkit
dmr_file <- system.file(
    "extdata",
    "example.dmr",
    package = "GAnnoViz")

# Plot
plot_dmg_chrom(
  dmr_file = dmr_file,
  violin_scale = "count",
  violin_border = 0.5,
  point_shape = 8,
  point_size = 2,
  jitter_width = 0.2,
  hyper_color = "#ff880088",
  hypo_color = "#0088ff88"
)

```
