library(shiny)
library(bs4Dash)
library(GAnnoViz)
library(ggplot2)
library(patchwork)
library(DT)
library(dplyr)

ui <- bs4DashPage(
  title = "GAnnoViz",
  fullscreen = TRUE,
  help = TRUE,
  dark = FALSE,
  scrollToTop = TRUE,
  header = bs4DashNavbar(skin = "light"),
  sidebar = bs4DashSidebar(
    disable = FALSE,
    skin = "dark",
    status = "warning",
    elevation = 3,
    collapsed = FALSE,
    minified = TRUE,
    expandOnHover = TRUE,
    fixed = TRUE,
    width = "300px",
    bs4SidebarMenu(
      id = "main_menu",
      bs4SidebarUserPanel(name = "GAnnoViz", image = "logo.png"),
      bs4SidebarHeader(title = "GAnnoViz"),
      bs4SidebarMenuItem(
        text = "Extract Features",
        icon = icon("database"),
        startExpanded = FALSE,
        bs4SidebarMenuSubItem(
          text = "extract_promoters",
          icon = icon("r-project"),
          tabName = "extract_promoters"
        ),
        bs4SidebarMenuSubItem(
          text = "extract_utr5",
          icon = icon("r-project"),
          tabName = "extract_utr5"
        ),
        bs4SidebarMenuSubItem(
          text = "extract_genes",
          icon = icon("r-project"),
          tabName = "extract_genes"
        ),
        bs4SidebarMenuSubItem(
          text = "extract_mrnas",
          icon = icon("r-project"),
          tabName = "extract_mrnas"
        ),
        bs4SidebarMenuSubItem(
          text = "extract_cds",
          icon = icon("r-project"),
          tabName = "extract_cds"
        ),
        bs4SidebarMenuSubItem(
          text = "extract_exons",
          icon = icon("r-project"),
          tabName = "extract_exons"
        ),
        bs4SidebarMenuSubItem(
          text = "extract_utr3",
          icon = icon("r-project"),
          tabName = "extract_utr3"
        )
      ),
      bs4SidebarMenuItem(
        text = "Plot Structure",
        icon = icon("circle-nodes"),
        startExpanded = TRUE,
        bs4SidebarMenuSubItem(
          text = "plot_gene_domains",
          icon = icon("r-project"),
          tabName = "plot_gene_domains"
        ),
        bs4SidebarMenuSubItem(
          text = "plot_gene_stats",
          icon = icon("r-project"),
          tabName = "plot_gene_stats"
        ),
        bs4SidebarMenuSubItem(
          text = "plot_gene_structure",
          icon = icon("r-project"),
          tabName = "plot_gene_structure"
        ),
        bs4SidebarMenuSubItem(
          text = "plot_interval_structure",
          icon = icon("r-project"),
          tabName = "plot_interval_structure"
        ),
        bs4SidebarMenuSubItem(
          text = "plot_interval_flank",
          icon = icon("r-project"),
          tabName = "plot_interval_flank"
        ),
        bs4SidebarMenuSubItem(
          text = "plot_chrom_structure",
          icon = icon("r-project"),
          tabName = "plot_chrom_structure"
        ),
        bs4SidebarMenuSubItem(
          text = "plot_chrom_genes",
          icon = icon("r-project"),
          tabName = "plot_chrom_genes"
        ),
        bs4SidebarMenuSubItem(
          text = "plot_chrom_heatmap",
          icon = icon("r-project"),
          tabName = "plot_chrom_heatmap"
        )
      ),
      bs4SidebarMenuItem(
        text = "DEG Anno & Viz",
        icon = icon("chart-simple"),
        startExpanded = FALSE,
        bs4SidebarMenuSubItem(
          text = "anno_deg_chrom",
          icon = icon("r-project"),
          tabName = "anno_deg_chrom"
        ),
        bs4SidebarMenuSubItem(
          text = "plot_deg_chrom",
          icon = icon("r-project"),
          tabName = "plot_deg_chrom"
        ),
        bs4SidebarMenuSubItem(
          text = "plot_deg_exp",
          icon = icon("r-project"),
          tabName = "plot_deg_exp"
        ),
        bs4SidebarMenuSubItem(
          text = "plot_deg_volcano",
          icon = icon("r-project"),
          tabName = "plot_deg_volcano"
        )
      ),
      bs4SidebarMenuItem(
        text = "SNP Anno & Plot",
        icon = icon("arrows-to-dot"),
        startExpanded = FALSE,
        bs4SidebarMenuSubItem(
          text = "plot_snp_density",
          icon = icon("r-project"),
          tabName = "plot_snp_density"
        ),
        bs4SidebarMenuSubItem(
          text = "plot_snp_fst",
          icon = icon("r-project"),
          tabName = "plot_snp_fst"
        ),
        bs4SidebarMenuSubItem(
          text = "plot_snp_anno",
          icon = icon("r-project"),
          tabName = "plot_snp_anno"
        )
      ),
      bs4SidebarMenuItem(
        text = "DMG Anno & Plot",
        icon = icon("atom"),
        startExpanded = FALSE,
        bs4SidebarMenuSubItem(
          text = "anno_fst_dmr",
          icon = icon("r-project"),
          tabName = "anno_fst_dmr"
        ),
        bs4SidebarMenuSubItem(
          text = "plot_dmg_chrom",
          icon = icon("r-project"),
          tabName = "plot_dmg_chrom"
        ),
        bs4SidebarMenuSubItem(
          text = "plot_dmg_exp",
          icon = icon("r-project"),
          tabName = "plot_dmg_exp"
        ),
        bs4SidebarMenuSubItem(
          text = "plot_dmg_trend",
          icon = icon("r-project"),
          tabName = "plot_dmg_trend"
        )
      )
    )
  ),
  body = bs4DashBody(
    tags$head(
      tags$link(rel = "stylesheet", href = "style.css"),
      tags$style(
        HTML(
          "
            .card-body {
              max-height: 800px;
              overflow-y: auto;
              scrollbar-width: thin;
            }
            .action-button {
            	width: 100%;
            	background-color: #ff880088;
            	border-radius: 10px;
            	color: #333333;
            	font-weight: bold;
            }
            .form-control, .selectize-input, .shiny-input-number, .shiny-input-select, .custom-file-input, .custom-file-label {
              border-radius: 10px !important;
              background-color: #f7f7f7 !important;
            }
            .irs-line, .irs-bar, .irs-handle {
              border-radius: 10px !important;
            }
          "
        )
      )
    ),
    bs4TabItems(
      bs4TabItem(tabName = "plot_gene_structure", fluidRow(
        column(
          width = 3,
          bs4Card(
            title = "Parameters",
            status = "info",
            solidHeader = TRUE,
            width = 12,
            actionButton(
              inputId = "run_plot_gene_structure",
              label = "Run",
              icon = icon("circle-play"),
              style = "width: 100%",
              class = "btn-block"
            ),
            br(),
            fileInput(inputId = "gff_file_gene", label = "GFF/GTF file"),
            selectInput(
              inputId = "gff_format_gene",
              label = "Format",
              choices = c("auto", "gff3", "gtf"),
              selected = "auto"
            ),
            textInput(
              inputId = "gene_id",
              label = "Gene ID",
              value = "HdF029609"
            ),
            numericInput(
              inputId = "upstream",
              label = "Promoter upstream",
              value = 2000,
              min = 0,
              step = 100
            ),
            numericInput(
              inputId = "downstream",
              label = "Promoter downstream",
              value = 200,
              min = 0,
              step = 50
            ),
            sliderInput(
              inputId = "feature_alpha",
              label = "Feature alpha",
              min = 0,
              max = 1,
              value = 0.8,
              step = 0.05
            ),
            numericInput(
              inputId = "intron_width",
              label = "Intron width",
              value = 1,
              min = 0,
              step = 0.5
            ),
            numericInput(
              inputId = "x_breaks",
              label = "X breaks",
              value = 10,
              min = 2,
              step = 1
            ),
            numericInput(
              inputId = "arrow_length",
              label = "Arrow length",
              value = 5,
              min = 1,
              step = 1
            ),
            numericInput(
              inputId = "arrow_count",
              label = "Arrow count",
              value = 1,
              min = 0,
              step = 1
            ),
            selectInput(
              inputId = "arrow_unit",
              label = "Arrow unit",
              choices = c("pt", "mm", "cm", "inches"),
              selected = "pt"
            ),
            textInput(
              inputId = "promoter_color",
              label = "Promoter color",
              value = "#ff8800"
            ),
            textInput(
              inputId = "utr5_color",
              label = "5'UTR color",
              value = "#008833"
            ),
            textInput(
              inputId = "utr3_color",
              label = "3'UTR color",
              value = "#ff0033"
            ),
            textInput(
              inputId = "exon_color",
              label = "Exon color",
              value = "#0033ff"
            ),
            textInput(
              inputId = "intron_color",
              label = "Intron color",
              value = "#333333"
            )
          )
        ), column(
          width = 9,
          bs4Card(
            title = "Plot",
            status = "danger",
            solidHeader = TRUE,
            width = 12,
            plotOutput("plot_gene_structure", height = "600px")
          ),

        )
      )),

      bs4TabItem(tabName = "anno_deg_chrom", fluidRow(
        column(
          width = 3,
          bs4Card(
            title = "Parameters",
            status = "info",
            solidHeader = TRUE,
            width = 12,
            actionButton(
              inputId = "run_anno_deg_chrom",
              label = "Run",
              icon = icon("circle-play"),
              style = "width: 100%",
              class = "btn-block"
            ),
            downloadButton(
              outputId = "download_table_anno_deg_chrom",
              label = "Download",
              icon = icon("download"),
              style = "width: 100%",
              class = "btn btn-success btn-block"
            ),
            br(),
            fileInput(inputId = "deg_file_anno", label = "DEG table"),
            fileInput(inputId = "gff_file_deg_anno", label = "GFF/GTF file"),
            selectInput(
              inputId = "gff_format_deg_anno",
              label = "Format",
              choices = c("auto", "gff3", "gtf"),
              selected = "auto"
            ),
            textInput(
              inputId = "id_col_anno",
              label = "ID column",
              value = "GeneID"
            ),
            textInput(
              inputId = "fc_col_anno",
              label = "FC column",
              value = "log2FoldChange"
            ),
            checkboxInput(
              inputId = "use_strand_anno",
              label = "Use strand",
              value = FALSE
            ),
            checkboxInput(
              inputId = "drop_unmapped_anno",
              label = "Drop unmapped",
              value = TRUE
            )
          )
        ), column(
          width = 9,
          bs4Card(
            title = "Data",
            status = "primary",
            solidHeader = TRUE,
            width = 12,
            DT::dataTableOutput("table_anno_deg_chrom")
          )
        )
      )),

      bs4TabItem(tabName = "plot_gene_domains", fluidRow(
        column(
          width = 3,
          bs4Card(
            title = "Parameters",
            status = "info",
            solidHeader = TRUE,
            width = 12,
            actionButton(
              inputId = "run_plot_gene_domains",
              label = "Run",
              icon = icon("circle-play"),
              style = "width: 100%",
              class = "btn-block"
            ),
            br(),
            textInput(
              inputId = "gene_name_domains",
              label = "Gene name",
              value = "TP53"
            ),
            textInput(
              inputId = "species_domains",
              label = "Species",
              value = "hsapiens"
            ),
            textInput(
              inputId = "transcript_id_domains",
              label = "Transcript ID",
              value = ""
            ),
            selectInput(
              inputId = "transcript_choice_domains",
              label = "Transcript choice",
              choices = c("longest", "canonical"),
              selected = "longest"
            ),
            selectInput(
              inputId = "palette_domains",
              label = "Palette",
              choices = c("Set 2", "Set 3", "Warm", "Cold", "Dynamic", "Viridis", "Plasma", "Inferno", "Rocket", "Mako"),
              selected = "Set 2"
            ),
            numericInput(
              inputId = "legend_ncol_domains",
              label = "Legend columns",
              value = 2,
              min = 1,
              step = 1
            )
          )
        ), column(
          width = 9,
          bs4Card(
            title = "Plot",
            status = "danger",
            solidHeader = TRUE,
            width = 12,
            plotOutput("plot_gene_domains", height = "600px")
          )
        )
      )),

      bs4TabItem(tabName = "plot_dmg_trend", fluidRow(
        column(
          width = 3,
          bs4Card(
            title = "Parameters",
            status = "info",
            solidHeader = TRUE,
            width = 12,
            actionButton(
              inputId = "run_plot_dmg_trend",
              label = "Run",
              icon = icon("circle-play"),
              class = "btn-block"
            ),
            br(),
            fileInput(inputId = "dmr_file_trend", label = "DMR table"),
            textInput(
              inputId = "chrom_id_trend",
              label = "Chrom ID",
              value = "chr1"
            ),
            sliderInput(
              inputId = "smooth_span_trend",
              label = "Smooth span",
              min = 0.05,
              max = 1,
              value = 0.1,
              step = 0.05
            ),
            textInput(
              inputId = "hyper_color_trend",
              label = "Hyper color",
              value = "#ff000055"
            ),
            textInput(
              inputId = "hypo_color_trend",
              label = "Hypo color",
              value = "#00880055"
            ),
            numericInput(
              inputId = "point_size_trend",
              label = "Point size",
              value = 3,
              min = 0.5,
              step = 0.5
            ),
            sliderInput(
              inputId = "point_alpha_trend",
              label = "Point alpha",
              min = 0,
              max = 1,
              value = 0.5,
              step = 0.05
            )
          )
        ), column(
          width = 9,
          bs4Card(
            title = "Plot",
            status = "danger",
            solidHeader = TRUE,
            width = 12,
            plotOutput("plot_dmg_trend", height = "600px")
          ),

        )
      )),

      bs4TabItem(tabName = "plot_interval_structure", fluidRow(
        column(
          width = 3,
          bs4Card(
            title = "Parameters",
            status = "info",
            solidHeader = TRUE,
            width = 12,
            actionButton(
              inputId = "run_plot_interval_structure",
              label = "Run",
              icon = icon("circle-play"),
              style = "width: 100%",
              class = "btn-block"
            ),
            br(),
            fileInput(inputId = "gff_file_interval", label = "GFF/GTF file"),
            selectInput(
              inputId = "gff_format_interval",
              label = "Format",
              choices = c("auto", "gff3", "gtf"),
              selected = "auto"
            ),
            textInput(
              inputId = "chrom_id",
              label = "Chrom ID",
              value = "chr1"
            ),
            numericInput(
              inputId = "win_start",
              label = "Start",
              value = 950000,
              min = 0,
              step = 1000
            ),
            numericInput(
              inputId = "win_end",
              label = "End",
              value = 1180000,
              min = 0,
              step = 1000
            ),
            numericInput(
              inputId = "x_breaks_interval",
              label = "X breaks",
              value = 10,
              min = 2,
              step = 1
            ),
            numericInput(
              inputId = "upstream_interval",
              label = "Promoter upstream",
              value = 2000,
              min = 0,
              step = 100
            ),
            numericInput(
              inputId = "downstream_interval",
              label = "Promoter downstream",
              value = 200,
              min = 0,
              step = 50
            ),
            sliderInput(
              inputId = "feature_alpha_interval",
              label = "Feature alpha",
              min = 0,
              max = 1,
              value = 0.8,
              step = 0.05
            ),
            numericInput(
              inputId = "intron_width_interval",
              label = "Intron width",
              value = 1,
              min = 0,
              step = 0.5
            ),
            numericInput(
              inputId = "arrow_length_interval",
              label = "Arrow length",
              value = 5,
              min = 1,
              step = 1
            ),
            numericInput(
              inputId = "arrow_count_interval",
              label = "Arrow count",
              value = 1,
              min = 0,
              step = 1
            ),
            selectInput(
              inputId = "arrow_unit_interval",
              label = "Arrow unit",
              choices = c("pt", "mm", "cm", "inches"),
              selected = "pt"
            ),
            textInput(
              inputId = "promoter_color_interval",
              label = "Promoter color",
              value = "#ff8800"
            ),
            textInput(
              inputId = "utr5_color_interval",
              label = "5'UTR color",
              value = "#008833"
            ),
            textInput(
              inputId = "utr3_color_interval",
              label = "3'UTR color",
              value = "#ff0033"
            ),
            textInput(
              inputId = "exon_color_interval",
              label = "Exon color",
              value = "#0033ff"
            ),
            textInput(
              inputId = "intron_color_interval",
              label = "Intron color",
              value = "#333333"
            )
          )
        ), column(
          width = 9,
          bs4Card(
            title = "Plot",
            status = "danger",
            solidHeader = TRUE,
            width = 12,
            plotOutput("plot_interval_structure", height = "600px")
          ),

        )
      )),

      bs4TabItem(tabName = "plot_interval_flank", fluidRow(
        column(
          width = 3,
          bs4Card(
            title = "Parameters",
            status = "info",
            solidHeader = TRUE,
            width = 12,
            actionButton(
              inputId = "run_plot_interval_flank",
              label = "Run",
              icon = icon("circle-play"),
              style = "width: 100%",
              class = "btn-block"
            ),
            br(),
            fileInput(inputId = "gff_file_flank", label = "GFF/GTF file"),
            selectInput(
              inputId = "gff_format_flank",
              label = "Format",
              choices = c("auto", "gff3", "gtf"),
              selected = "auto"
            ),
            textInput(
              inputId = "flank_gene_id",
              label = "Gene ID",
              value = "HdF029609"
            ),
            numericInput(
              inputId = "flank_upstream",
              label = "Flank upstream",
              value = 200000,
              min = 0,
              step = 1000
            ),
            numericInput(
              inputId = "flank_downstream",
              label = "Flank downstream",
              value = 200000,
              min = 0,
              step = 1000
            ),
            checkboxInput(
              inputId = "show_promoters",
              label = "Show promoters",
              value = TRUE
            ),
            numericInput(
              inputId = "upstream_flank",
              label = "Promoter upstream",
              value = 2000,
              min = 0,
              step = 100
            ),
            numericInput(
              inputId = "downstream_flank",
              label = "Promoter downstream",
              value = 200,
              min = 0,
              step = 50
            ),
            numericInput(
              inputId = "arrow_length_flank",
              label = "Arrow length",
              value = 5,
              min = 1,
              step = 1
            ),
            selectInput(
              inputId = "arrow_unit_flank",
              label = "Arrow unit",
              choices = c("pt", "mm", "cm", "inches"),
              selected = "pt"
            ),
            textInput(
              inputId = "gene_color_flank",
              label = "Gene color",
              value = "#0088ff"
            ),
            textInput(
              inputId = "promoter_color_flank",
              label = "Promoter color",
              value = "#ff8800"
            ),
            numericInput(
              inputId = "label_size_flank",
              label = "Label size",
              value = 3,
              min = 1,
              step = 1
            )
          )
        ), column(
          width = 9,
          bs4Card(
            title = "Plot",
            status = "danger",
            solidHeader = TRUE,
            width = 12,
            plotOutput("plot_interval_flank", height = "600px")
          ),

        )
      )),

      bs4TabItem(tabName = "plot_chrom_structure", fluidRow(
        column(
          width = 3,
          bs4Card(
            title = "Parameters",
            status = "info",
            solidHeader = TRUE,
            width = 12,
            actionButton(
              inputId = "run_plot_chrom_structure",
              label = "Run",
              icon = icon("circle-play"),
              style = "width: 100%",
              class = "btn-block"
            ),
            br(),
            fileInput(inputId = "gff_file_chrom", label = "GFF/GTF file"),
            selectInput(
              inputId = "gff_format_chrom",
              label = "Format",
              choices = c("auto", "gff3", "gtf"),
              selected = "auto"
            ),
            selectInput(
              inputId = "chrom_orientation",
              label = "Orientation",
              choices = c("vertical", "horizontal"),
              selected = "vertical"
            ),
            sliderInput(
              inputId = "bar_width",
              label = "Bar width",
              min = 0.1,
              max = 1,
              value = 0.6,
              step = 0.05
            ),
            sliderInput(
              inputId = "chrom_alpha",
              label = "Chrom alpha",
              min = 0,
              max = 1,
              value = 0.1,
              step = 0.05
            ),
            sliderInput(
              inputId = "gene_width",
              label = "Gene width",
              min = 0.1,
              max = 1,
              value = 0.5,
              step = 0.05
            ),
            textInput(
              inputId = "chrom_color",
              label = "Chrom color",
              value = "#008888"
            ),
            textInput(
              inputId = "gene_color_chrom",
              label = "Gene color",
              value = "#0088ff"
            ),
            textInput(
              inputId = "telomere_color",
              label = "Telomere color",
              value = "#ff0000"
            ),
            numericInput(
              inputId = "label_size_chrom",
              label = "Label size",
              value = 3,
              min = 1,
              step = 1
            )
          )
        ), column(
          width = 9,
          bs4Card(
            title = "Plot",
            status = "danger",
            solidHeader = TRUE,
            width = 12,
            plotOutput("plot_chrom_structure", height = "600px")
          ),

        )
      )),

      bs4TabItem(tabName = "plot_chrom_genes", fluidRow(
        column(
          width = 3,
          bs4Card(
            title = "Parameters",
            status = "info",
            solidHeader = TRUE,
            width = 12,
            actionButton(
              inputId = "run_plot_chrom_genes",
              label = "Run",
              icon = icon("circle-play"),
              style = "width: 100%",
              class = "btn-block"
            ),
            br(),
            fileInput(inputId = "gff_file_chrom_genes", label = "GFF/GTF file"),
            selectInput(
              inputId = "gff_format_chrom_genes",
              label = "Format",
              choices = c("auto", "gff3", "gtf"),
              selected = "auto"
            ),
            fileInput(inputId = "gene_table_file", label = "Gene table (2 cols: id,name)"),
            selectInput(
              inputId = "annotate_mode",
              label = "Annotate",
              choices = c("id", "name"),
              selected = "id"
            ),
            selectInput(
              inputId = "chrom_genes_orientation",
              label = "Orientation",
              choices = c("vertical", "horizontal"),
              selected = "vertical"
            ),
            sliderInput(
              inputId = "min_gap_frac",
              label = "Min gap frac",
              min = 0.005,
              max = 0.1,
              value = 0.02,
              step = 0.005
            ),
            sliderInput(
              inputId = "bar_width_genes",
              label = "Bar width",
              min = 0.1,
              max = 1,
              value = 0.6,
              step = 0.05
            ),
            sliderInput(
              inputId = "chrom_alpha_genes",
              label = "Chrom alpha",
              min = 0,
              max = 1,
              value = 0.1,
              step = 0.05
            ),
            sliderInput(
              inputId = "gene_width_genes",
              label = "Gene width",
              min = 0.1,
              max = 1,
              value = 0.5,
              step = 0.05
            ),
            textInput(
              inputId = "chrom_color_genes",
              label = "Chrom color",
              value = "#008888"
            ),
            textInput(
              inputId = "gene_color_genes",
              label = "Gene color",
              value = "#0088ff"
            ),
            textInput(
              inputId = "telomere_color_genes",
              label = "Telomere color",
              value = "#ff0000"
            ),
            numericInput(
              inputId = "label_size_genes",
              label = "Label size",
              value = 3,
              min = 1,
              step = 1
            ),
            numericInput(
              inputId = "connector_dx1_genes",
              label = "Connector dx1",
              value = 0.2,
              min = 0,
              step = 0.05
            ),
            numericInput(
              inputId = "connector_dx2_genes",
              label = "Connector dx2",
              value = 0.2,
              min = 0,
              step = 0.05
            )
          )
        ), column(
          width = 9,
          bs4Card(
            title = "Plot",
            status = "danger",
            solidHeader = TRUE,
            width = 12,
            plotOutput("plot_chrom_genes", height = "600px")
          ),

        )
      )),

      bs4TabItem(tabName = "plot_chrom_heatmap", fluidRow(
        column(
          width = 3,
          bs4Card(
            title = "Parameters",
            status = "info",
            solidHeader = TRUE,
            width = 12,
            actionButton(
              inputId = "run_plot_chrom_heatmap",
              label = "Run",
              icon = icon("circle-play"),
              style = "width: 100%",
              class = "btn-block"
            ),
            br(),
            fileInput(inputId = "gff_file_heatmap", label = "GFF/GTF file"),
            selectInput(
              inputId = "gff_format_heatmap",
              label = "Format",
              choices = c("auto", "gff3", "gtf"),
              selected = "auto"
            ),
            selectInput(
              inputId = "feature_type",
              label = "Feature",
              choices = c("gene", "exon", "CDS", "promoter"),
              selected = "gene"
            ),
            numericInput(
              inputId = "bin_size_heatmap",
              label = "Bin size",
              value = 1e6,
              min = 1e4,
              step = 1e5
            ),
            selectInput(
              inputId = "orientation_heatmap",
              label = "Orientation",
              choices = c("horizontal", "vertical"),
              selected = "horizontal"
            ),
            textInput(
              inputId = "palette_start",
              label = "Palette start",
              value = "#ffffff"
            ),
            textInput(
              inputId = "palette_end",
              label = "Palette end",
              value = "#0055aa"
            ),
            sliderInput(
              inputId = "alpha_heatmap",
              label = "Alpha",
              min = 0,
              max = 1,
              value = 0.9,
              step = 0.05
            )
          )
        ), column(
          width = 9,
          bs4Card(
            title = "Plot",
            status = "danger",
            solidHeader = TRUE,
            width = 12,
            plotOutput("plot_chrom_heatmap", height = "600px")
          ),

        )
      )),

      bs4TabItem(tabName = "plot_deg_chrom", fluidRow(
        column(
          width = 3,
          bs4Card(
            title = "Parameters",
            status = "info",
            solidHeader = TRUE,
            width = 12,
            actionButton(
              inputId = "run_plot_deg_chrom",
              label = "Run",
              icon = icon("circle-play"),
              style = "width: 100%",
              class = "btn-block"
            ),
            downloadButton(
              outputId = "download_plot_deg_chrom",
              label = "Download",
              icon = icon("download"),
              style = "width: 100%",
              class = "btn btn-success btn-block"
            ),
            br(),
            fileInput(inputId = "deg_file", label = "DEG table"),
            fileInput(inputId = "gff_file_deg", label = "GFF/GTF file"),
            selectInput(
              inputId = "gff_format_deg",
              label = "Format",
              choices = c("auto", "gff3", "gtf"),
              selected = "auto"
            ),
            textInput(
              inputId = "id_col",
              label = "ID column",
              value = "GeneID"
            ),
            textInput(
              inputId = "fc_col",
              label = "FC column",
              value = "log2FoldChange"
            ),
            selectInput(
              inputId = "violin_scale",
              label = "Violin scale",
              choices = c("count", "area", "width"),
              selected = "count"
            ),
            sliderInput(
              inputId = "violin_border",
              label = "Violin border",
              min = 0,
              max = 2,
              value = 0.5,
              step = 0.1
            ),
            numericInput(
              inputId = "point_shape_deg",
              label = "Point shape",
              value = 16,
              min = 0,
              step = 1
            ),
            numericInput(
              inputId = "point_size_deg",
              label = "Point size",
              value = 2,
              min = 0.5,
              step = 0.5
            ),
            sliderInput(
              inputId = "jitter_width_deg",
              label = "Jitter width",
              min = 0,
              max = 1,
              value = 0.2,
              step = 0.05
            ),
            textInput(
              inputId = "hyper_color_deg",
              label = "Hyper color",
              value = "#ff000088"
            ),
            textInput(
              inputId = "hypo_color_deg",
              label = "Hypo color",
              value = "#00880088"
            ),
            numericInput(
              inputId = "plot_width_deg_chrom",
              label = "Width (in)",
              value = 10,
              min = 1,
              step = 0.5
            ),
            numericInput(
              inputId = "plot_height_deg_chrom",
              label = "Height (in)",
              value = 6,
              min = 1,
              step = 0.5
            ),
            numericInput(
              inputId = "plot_dpi_deg_chrom",
              label = "DPI",
              value = 300,
              min = 72,
              step = 10
            ),
            selectInput(
              inputId = "plot_format_deg_chrom",
              label = "Format",
              choices = c("pdf", "jpeg"),
              selected = "pdf"
            )
          )
        ), column(
          width = 9,
          bs4Card(
            title = "Plot",
            status = "danger",
            solidHeader = TRUE,
            width = 12,
            plotOutput("plot_deg_chrom", height = "600px")
          ),

        )
      )),

      bs4TabItem(tabName = "plot_deg_exp", fluidRow(
        column(
          width = 3,
          bs4Card(
            title = "Parameters",
            status = "info",
            solidHeader = TRUE,
            width = 12,
            actionButton(
              inputId = "run_plot_deg_exp",
              label = "Run",
              icon = icon("circle-play"),
              style = "width: 100%",
              class = "btn-block"
            ),
            downloadButton(
              outputId = "download_plot_deg_exp",
              label = "Download",
              icon = icon("download"),
              style = "width: 100%",
              class = "btn btn-success btn-block"
            ),
            br(),
            fileInput(inputId = "deg_file_exp", label = "DEG table"),
            fileInput(inputId = "gff_file_deg_exp", label = "GFF/GTF file"),
            selectInput(
              inputId = "gff_format_deg_exp",
              label = "Format",
              choices = c("auto", "gff3", "gtf"),
              selected = "auto"
            ),
            textInput(
              inputId = "id_col_exp",
              label = "ID column",
              value = "GeneID"
            ),
            textInput(
              inputId = "fc_col_exp",
              label = "FC column",
              value = "log2FoldChange"
            ),
            selectInput(
              inputId = "orientation_deg_exp",
              label = "Orientation",
              choices = c("horizontal", "vertical"),
              selected = "horizontal"
            ),
            sliderInput(
              inputId = "chrom_alpha_deg_exp",
              label = "Chrom alpha",
              min = 0,
              max = 1,
              value = 0.1,
              step = 0.05
            ),
            numericInput(
              inputId = "bar_height_deg_exp",
              label = "Bar height",
              value = 0.8,
              min = 0.1,
              step = 0.1
            ),
            textInput(
              inputId = "chrom_color_deg_exp",
              label = "Chrom color",
              value = "#008888"
            ),
            numericInput(
              inputId = "point_size_deg_exp",
              label = "Point size",
              value = 1,
              min = 0.5,
              step = 0.5
            ),
            sliderInput(
              inputId = "point_alpha_deg_exp",
              label = "Point alpha",
              min = 0,
              max = 1,
              value = 0.3,
              step = 0.05
            ),
            textInput(
              inputId = "up_color_deg_exp",
              label = "Up color",
              value = "#ff0000"
            ),
            textInput(
              inputId = "down_color_deg_exp",
              label = "Down color",
              value = "#008800"
            ),
            selectInput(
              inputId = "mark_style_deg_exp",
              label = "Mark style",
              choices = c("point", "line"),
              selected = "point"
            ),
            numericInput(
              inputId = "line_width_deg_exp",
              label = "Line width",
              value = 0.6,
              min = 0.1,
              step = 0.1
            ),
            numericInput(
              inputId = "line_height_deg_exp",
              label = "Line height",
              value = 0.8,
              min = 0.1,
              step = 0.1
            ),
            numericInput(
              inputId = "plot_width_deg_exp",
              label = "Width (in)",
              value = 10,
              min = 1,
              step = 0.5
            ),
            numericInput(
              inputId = "plot_height_deg_exp",
              label = "Height (in)",
              value = 6,
              min = 1,
              step = 0.5
            ),
            numericInput(
              inputId = "plot_dpi_deg_exp",
              label = "DPI",
              value = 300,
              min = 72,
              step = 10
            ),
            selectInput(
              inputId = "plot_format_deg_exp",
              label = "Format",
              choices = c("pdf", "jpeg"),
              selected = "pdf"
            )
          )
        ), column(
          width = 9,
          bs4Card(
            title = "Plot",
            status = "danger",
            solidHeader = TRUE,
            width = 12,
            plotOutput("plot_deg_exp", height = "600px")
          )
        )
      )),

      bs4TabItem(tabName = "plot_deg_volcano", fluidRow(
        column(
          width = 3,
          bs4Card(
            title = "Parameters",
            status = "info",
            solidHeader = TRUE,
            width = 12,
            actionButton(
              inputId = "run_plot_deg_volcano",
              label = "Run",
              icon = icon("circle-play"),
              style = "width: 100%",
              class = "btn-block"
            ),
            downloadButton(
              outputId = "download_plot_deg_volcano",
              label = "Download",
              icon = icon("download"),
              style = "width: 100%",
              class = "btn btn-success btn-block"
            ),
            br(),
            fileInput(inputId = "deg_file_volcano", label = "DEG table"),
            textInput(
              inputId = "id_col_volcano",
              label = "ID column",
              value = "GeneID"
            ),
            textInput(
              inputId = "fc_col_volcano",
              label = "FC column",
              value = "log2FoldChange"
            ),
            textInput(
              inputId = "sig_col_volcano",
              label = "Sig column",
              value = "padj"
            ),
            numericInput(
              inputId = "fc_threshold_volcano",
              label = "FC threshold",
              value = 1,
              min = 0,
              step = 0.1
            ),
            numericInput(
              inputId = "sig_threshold_volcano",
              label = "Sig threshold",
              value = 0.05,
              min = 0,
              max = 1,
              step = 0.01
            ),
            numericInput(
              inputId = "point_size_volcano",
              label = "Point size",
              value = 1.5,
              min = 0.5,
              step = 0.5
            ),
            sliderInput(
              inputId = "point_alpha_volcano",
              label = "Point alpha",
              min = 0,
              max = 1,
              value = 0.6,
              step = 0.05
            ),
            textInput(
              inputId = "up_color_volcano",
              label = "Up color",
              value = "#ff0000"
            ),
            textInput(
              inputId = "down_color_volcano",
              label = "Down color",
              value = "#008800"
            ),
            textInput(
              inputId = "ns_color_volcano",
              label = "NS color",
              value = "#999999"
            ),
            numericInput(
              inputId = "plot_width_deg_volcano",
              label = "Width (in)",
              value = 10,
              min = 1,
              step = 0.5
            ),
            numericInput(
              inputId = "plot_height_deg_volcano",
              label = "Height (in)",
              value = 6,
              min = 1,
              step = 0.5
            ),
            numericInput(
              inputId = "plot_dpi_deg_volcano",
              label = "DPI",
              value = 300,
              min = 72,
              step = 10
            ),
            selectInput(
              inputId = "plot_format_deg_volcano",
              label = "Format",
              choices = c("pdf", "jpeg"),
              selected = "pdf"
            )
          )
        ), column(
          width = 9,
          bs4Card(
            title = "Plot",
            status = "danger",
            solidHeader = TRUE,
            width = 12,
            plotOutput("plot_deg_volcano", height = "600px")
          )
        )
      )),

      bs4TabItem(tabName = "plot_snp_fst", fluidRow(
        column(
          width = 3,
          bs4Card(
            title = "Parameters",
            status = "info",
            solidHeader = TRUE,
            width = 12,
            actionButton(
              inputId = "run_plot_snp_fst",
              label = "Run",
              icon = icon("circle-play"),
              style = "width: 100%",
              class = "btn-block"
            ),
            br(),
            fileInput(inputId = "fst_file_heat", label = "FST table"),
            numericInput(
              inputId = "bin_size_fst",
              label = "Bin size",
              value = 1e6,
              min = 1e4,
              step = 1e5
            ),
            selectInput(
              inputId = "metric",
              label = "Metric",
              choices = c("fst_mean", "variant_count"),
              selected = "fst_mean"
            ),
            selectInput(
              inputId = "orientation_fst",
              label = "Orientation",
              choices = c("horizontal", "vertical"),
              selected = "horizontal"
            ),
            textInput(
              inputId = "palette_start_fst",
              label = "Palette start",
              value = "#ffffff"
            ),
            textInput(
              inputId = "palette_end_fst",
              label = "Palette end",
              value = "#aa00aa"
            ),
            sliderInput(
              inputId = "alpha_fst",
              label = "Alpha",
              min = 0,
              max = 1,
              value = 0.9,
              step = 0.05
            )
          )
        ), column(
          width = 9,
          bs4Card(
            title = "Plot",
            status = "danger",
            solidHeader = TRUE,
            width = 12,
            plotOutput("plot_snp_fst", height = "600px")
          ),

        )
      )),

      bs4TabItem(tabName = "plot_snp_anno", fluidRow(
        column(
          width = 3,
          bs4Card(
            title = "Parameters",
            status = "info",
            solidHeader = TRUE,
            width = 12,
            actionButton(
              inputId = "run_plot_snp_anno",
              label = "Run",
              icon = icon("circle-play"),
              style = "width: 100%",
              class = "btn-block"
            ),
            br(),
            fileInput(inputId = "fst_file_anno", label = "FST table"),
            fileInput(inputId = "gff_file_fst", label = "GFF/GTF file"),
            selectInput(
              inputId = "gff_format_fst",
              label = "Format",
              choices = c("auto", "gff3", "gtf"),
              selected = "auto"
            ),
            textInput(
              inputId = "chrom_id_fst",
              label = "Chrom ID",
              value = "chr2"
            ),
            numericInput(
              inputId = "top_n",
              label = "Top N",
              value = 20,
              min = 1,
              step = 1
            ),
            selectInput(
              inputId = "orientation_fst_anno",
              label = "Orientation",
              choices = c("vertical", "horizontal"),
              selected = "vertical"
            ),
            sliderInput(
              inputId = "smooth_span",
              label = "Smooth span",
              min = 0.05,
              max = 1,
              value = 0.5,
              step = 0.05
            ),
            numericInput(
              inputId = "point_size_fst",
              label = "Point size",
              value = 1,
              min = 0.5,
              step = 0.5
            ),
            sliderInput(
              inputId = "point_alpha_fst",
              label = "Point alpha",
              min = 0,
              max = 1,
              value = 0.3,
              step = 0.05
            ),
            numericInput(
              inputId = "label_size_fst",
              label = "Label size",
              value = 3,
              min = 1,
              step = 1
            ),
            numericInput(
              inputId = "connector_dx1",
              label = "Connector dx1",
              value = 2e4,
              min = 0,
              step = 1e3
            ),
            numericInput(
              inputId = "connector_dx2",
              label = "Connector dx2",
              value = 4e4,
              min = 0,
              step = 1e3
            ),
            sliderInput(
              inputId = "gap_frac",
              label = "Gap frac",
              min = 0.01,
              max = 0.2,
              value = 0.05,
              step = 0.01
            ),
            textInput(
              inputId = "fst_color",
              label = "FST color",
              value = "#0088ff"
            )
          )
        ), column(
          width = 9,
          bs4Card(
            title = "Plot",
            status = "danger",
            solidHeader = TRUE,
            width = 12,
            plotOutput("plot_snp_anno", height = "600px")
          ),

        )
      )),

      bs4TabItem(tabName = "plot_dmg_chrom", fluidRow(
        column(
          width = 3,
          bs4Card(
            title = "Parameters",
            status = "info",
            solidHeader = TRUE,
            width = 12,
            actionButton(
              inputId = "run_plot_dmg_chrom",
              label = "Run",
              icon = icon("circle-play"),
              style = "width: 100%",
              class = "btn-block"
            ),
            br(),
            fileInput(inputId = "dmr_file", label = "DMR table"),
            selectInput(
              inputId = "violin_scale_dmr",
              label = "Violin scale",
              choices = c("count", "area", "width"),
              selected = "count"
            ),
            sliderInput(
              inputId = "violin_border_dmr",
              label = "Violin border",
              min = 0,
              max = 2,
              value = 0.5,
              step = 0.1
            ),
            numericInput(
              inputId = "point_shape_dmr",
              label = "Point shape",
              value = 8,
              min = 0,
              step = 1
            ),
            numericInput(
              inputId = "point_size_dmr",
              label = "Point size",
              value = 2,
              min = 0.5,
              step = 0.5
            ),
            sliderInput(
              inputId = "jitter_width_dmr",
              label = "Jitter width",
              min = 0,
              max = 1,
              value = 0.2,
              step = 0.05
            ),
            textInput(
              inputId = "hyper_color_dmr",
              label = "Hyper color",
              value = "#ff880088"
            ),
            textInput(
              inputId = "hypo_color_dmr",
              label = "Hypo color",
              value = "#0088ff88"
            )
          )
        ), column(
          width = 9,
          bs4Card(
            title = "Plot",
            status = "danger",
            solidHeader = TRUE,
            width = 12,
            plotOutput("plot_dmg_chrom", height = "600px")
          ),

        )
      )),
      bs4TabItem(tabName = "plot_dmg_exp", fluidRow(
        column(
          width = 3,
          bs4Card(
            title = "Parameters",
            status = "info",
            solidHeader = TRUE,
            width = 12,
            actionButton(
              inputId = "run_plot_dmg_exp",
              label = "Run",
              icon = icon("circle-play"),
              style = "width: 100%",
              class = "btn-block"
            ),
            br(),
            fileInput(inputId = "dmr_file_exp", label = "DMR table"),
            selectInput(
              inputId = "orientation_dmg_exp",
              label = "Orientation",
              choices = c("horizontal", "vertical"),
              selected = "horizontal"
            ),
            sliderInput(
              inputId = "chrom_alpha_dmg_exp",
              label = "Chrom alpha",
              min = 0,
              max = 1,
              value = 0.1,
              step = 0.05
            ),
            numericInput(
              inputId = "bar_height_dmg_exp",
              label = "Bar height",
              value = 0.8,
              min = 0.1,
              step = 0.1
            ),
            textInput(
              inputId = "chrom_color_dmg_exp",
              label = "Chrom color",
              value = "#008888"
            ),
            numericInput(
              inputId = "point_size_dmg_exp",
              label = "Point size",
              value = 1,
              min = 0.5,
              step = 0.5
            ),
            sliderInput(
              inputId = "point_alpha_dmg_exp",
              label = "Point alpha",
              min = 0,
              max = 1,
              value = 0.3,
              step = 0.05
            ),
            textInput(
              inputId = "hyper_color_dmg_exp",
              label = "Hyper color",
              value = "#ff0000"
            ),
            textInput(
              inputId = "hypo_color_dmg_exp",
              label = "Hypo color",
              value = "#008800"
            ),
            selectInput(
              inputId = "mark_style_dmg_exp",
              label = "Mark style",
              choices = c("point", "line"),
              selected = "point"
            ),
            numericInput(
              inputId = "line_width_dmg_exp",
              label = "Line width",
              value = 0.6,
              min = 0.1,
              step = 0.1
            ),
            numericInput(
              inputId = "line_height_dmg_exp",
              label = "Line height",
              value = 0.8,
              min = 0.1,
              step = 0.1
            )
          )
        ), column(
          width = 9,
          bs4Card(
            title = "Plot",
            status = "danger",
            solidHeader = TRUE,
            width = 12,
            plotOutput("plot_dmg_exp", height = "600px")
          )
        )
      ))
      ,
      bs4TabItem(tabName = "plot_gene_stats", fluidRow(
        column(
          width = 3,
          bs4Card(
            title = "Parameters",
            status = "info",
            solidHeader = TRUE,
            width = 12,
            actionButton(
              inputId = "run_plot_gene_stats",
              label = "Run",
              icon = icon("circle-play"),
              style = "width: 100%",
              class = "btn-block"
            ),
            br(),
            fileInput(inputId = "gff_file_gene_stats", label = "GFF/GTF file"),
            selectInput(
              inputId = "gff_format_gene_stats",
              label = "Format",
              choices = c("auto", "gff3", "gtf"),
              selected = "auto"
            ),
            sliderInput(
              inputId = "bar_width_gene_stats",
              label = "Bar width",
              min = 0.1,
              max = 1,
              value = 0.7,
              step = 0.05
            ),
            textInput(
              inputId = "bar_color_gene_stats",
              label = "Bar color",
              value = "#0055ff55"
            ),
            numericInput(
              inputId = "label_size_gene_stats",
              label = "Label size",
              value = 3,
              min = 1,
              step = 1
            )
          )
        ), column(
          width = 9,
          bs4Card(
            title = "Plot",
            status = "danger",
            solidHeader = TRUE,
            width = 12,
            plotOutput("plot_gene_stats", height = "600px")
          )
        )
      )),
      bs4TabItem(tabName = "plot_snp_density", fluidRow(
        column(
          width = 3,
          bs4Card(
            title = "Parameters",
            status = "info",
            solidHeader = TRUE,
            width = 12,
            actionButton(
              inputId = "run_plot_snp_density",
              label = "Run",
              icon = icon("circle-play"),
              style = "width: 100%",
              class = "btn-block"
            ),
            br(),
            fileInput(inputId = "fst_file_density", label = "FST table"),
            checkboxInput(
              inputId = "log10_density",
              label = "LOG10",
              value = FALSE
            ),
            numericInput(
              inputId = "bin_size_density",
              label = "Bin size",
              value = 1e6,
              min = 1e4,
              step = 1e5
            ),
            textInput(
              inputId = "density_color1",
              label = "Density color 1",
              value = "#0088ff"
            ),
            textInput(
              inputId = "density_color2",
              label = "Density color 2",
              value = "#ff8800"
            ),
            textInput(
              inputId = "density_color3",
              label = "Density color 3",
              value = "#ff0000"
            )
          )
        ), column(
          width = 9,
          bs4Card(
            title = "Plot",
            status = "danger",
            solidHeader = TRUE,
            width = 12,
            plotOutput("plot_snp_density", height = "600px")
          )
        )
      )),
      bs4TabItem(tabName = "anno_fst_dmr", fluidRow(
        column(
          width = 3,
          bs4Card(
            title = "Parameters",
            status = "info",
            solidHeader = TRUE,
            width = 12,
            actionButton(
              inputId = "run_anno_fst_dmr",
              label = "Run",
              icon = icon("circle-play"),
              style = "width: 100%",
              class = "btn-block"
            ),
            br(),
            fileInput(inputId = "gff_file_anno_ranges", label = "GFF/GTF file"),
            selectInput(
              inputId = "gff_format_anno_ranges",
              label = "Format",
              choices = c("auto", "gff3", "gtf"),
              selected = "auto"
            ),
            fileInput(inputId = "genomic_ranges_file", label = "Genomic ranges (FST/DMR)"),
            textInput(
              inputId = "chrom_col_ranges",
              label = "Chrom column",
              value = "CHROM"
            ),
            textInput(
              inputId = "start_col_ranges",
              label = "Start column",
              value = "BIN_START"
            ),
            textInput(
              inputId = "end_col_ranges",
              label = "End column",
              value = "BIN_END"
            ),
            numericInput(
              inputId = "upstream_ranges",
              label = "Promoter upstream",
              value = 2000,
              min = 0,
              step = 100
            ),
            numericInput(
              inputId = "downstream_ranges",
              label = "Promoter downstream",
              value = 200,
              min = 0,
              step = 50
            ),
            checkboxInput(
              inputId = "ignore_strand_ranges",
              label = "Ignore strand",
              value = TRUE
            ),
            selectInput(
              inputId = "features_ranges",
              label = "Features",
              multiple = TRUE,
              choices = c(
                "promoter",
                "UTR5",
                "gene",
                "exon",
                "intron",
                "CDS",
                "UTR3",
                "intergenic"
              ),
              selected = c("promoter", "UTR5", "gene", "exon", "intron", "CDS", "UTR3")
            )
          )
        ), column(
          width = 9,
          bs4Card(
            title = "Data",
            status = "primary",
            solidHeader = TRUE,
            width = 12,
            DT::dataTableOutput("table_anno_fst_dmr")
          )
        )
      )),
      bs4TabItem(tabName = "extract_promoters", fluidRow(
        column(
          width = 3,
          bs4Card(
            title = "Parameters",
            status = "info",
            solidHeader = TRUE,
            width = 12,
            actionButton(
              inputId = "run_extract_promoters",
              label = "Run",
              icon = icon("circle-play"),
              style = "width: 100%",
              class = "btn-block"
            ),
            br(),
            fileInput(inputId = "gff_file_extract_promoters", label = "GFF/GTF file"),
            selectInput(
              inputId = "gff_format_extract_promoters",
              label = "Format",
              choices = c("auto", "gff3", "gtf"),
              selected = "auto"
            ),
            numericInput(
              inputId = "upstream_extract_promoters",
              label = "Promoter upstream",
              value = 2000,
              min = 0,
              step = 100
            ),
            numericInput(
              inputId = "downstream_extract_promoters",
              label = "Promoter downstream",
              value = 200,
              min = 0,
              step = 50
            ),
            selectInput(
              inputId = "promoter_info",
              label = "Info",
              choices = c("all", "chrom_id", "promoter_id", "promoter_range"),
              selected = "all"
            )
          )
        ), column(
          width = 9,
          bs4Card(
            title = "Data",
            status = "primary",
            solidHeader = TRUE,
            width = 12,
            DT::dataTableOutput("table_extract_promoters")
          )
        )
      )),
      bs4TabItem(tabName = "extract_utr5", fluidRow(
        column(
          width = 3,
          bs4Card(
            title = "Parameters",
            status = "info",
            solidHeader = TRUE,
            width = 12,
            actionButton(
              inputId = "run_extract_utr5",
              label = "Run",
              icon = icon("circle-play"),
              style = "width: 100%",
              class = "btn-block"
            ),
            br(),
            fileInput(inputId = "gff_file_extract_utr5", label = "GFF/GTF file"),
            selectInput(
              inputId = "gff_format_extract_utr5",
              label = "Format",
              choices = c("auto", "gff3", "gtf"),
              selected = "auto"
            ),
            selectInput(
              inputId = "utr5_info",
              label = "Info",
              choices = c("all", "chrom_id", "utr5_range"),
              selected = "all"
            )
          )
        ), column(
          width = 9,
          bs4Card(
            title = "Data",
            status = "primary",
            solidHeader = TRUE,
            width = 12,
            DT::dataTableOutput("table_extract_utr5")
          )
        )
      )),
      bs4TabItem(tabName = "extract_genes", fluidRow(
        column(
          width = 3,
          bs4Card(
            title = "Parameters",
            status = "info",
            solidHeader = TRUE,
            width = 12,
            actionButton(
              inputId = "run_extract_genes",
              label = "Run",
              icon = icon("circle-play"),
              style = "width: 100%",
              class = "btn-block"
            ),
            br(),
            fileInput(inputId = "gff_file_extract_genes", label = "GFF/GTF file"),
            selectInput(
              inputId = "gff_format_extract_genes",
              label = "Format",
              choices = c("auto", "gff3", "gtf"),
              selected = "auto"
            ),
            selectInput(
              inputId = "gene_info_opt",
              label = "Info",
              choices = c("all", "chrom_id", "gene_id", "gene_range"),
              selected = "all"
            )
          )
        ), column(
          width = 9,
          bs4Card(
            title = "Data",
            status = "primary",
            solidHeader = TRUE,
            width = 12,
            DT::dataTableOutput("table_extract_genes")
          )
        )
      )),
      bs4TabItem(tabName = "extract_mrnas", fluidRow(
        column(
          width = 3,
          bs4Card(
            title = "Parameters",
            status = "info",
            solidHeader = TRUE,
            width = 12,
            actionButton(
              inputId = "run_extract_mrnas",
              label = "Run",
              icon = icon("circle-play"),
              style = "width: 100%",
              class = "btn-block"
            ),
            br(),
            fileInput(inputId = "gff_file_extract_mrnas", label = "GFF/GTF file"),
            selectInput(
              inputId = "gff_format_extract_mrnas",
              label = "Format",
              choices = c("auto", "gff3", "gtf"),
              selected = "auto"
            ),
            selectInput(
              inputId = "mrna_info_opt",
              label = "Info",
              choices = c("all", "chrom_id", "mrna_id", "mrna_range"),
              selected = "all"
            )
          )
        ), column(
          width = 9,
          bs4Card(
            title = "Data",
            status = "primary",
            solidHeader = TRUE,
            width = 12,
            DT::dataTableOutput("table_extract_mrnas")
          )
        )
      )),
      bs4TabItem(tabName = "extract_cds", fluidRow(
        column(
          width = 3,
          bs4Card(
            title = "Parameters",
            status = "info",
            solidHeader = TRUE,
            width = 12,
            actionButton(
              inputId = "run_extract_cds",
              label = "Run",
              icon = icon("circle-play"),
              style = "width: 100%",
              class = "btn-block"
            ),
            br(),
            fileInput(inputId = "gff_file_extract_cds", label = "GFF/GTF file"),
            selectInput(
              inputId = "gff_format_extract_cds",
              label = "Format",
              choices = c("auto", "gff3", "gtf"),
              selected = "auto"
            ),
            selectInput(
              inputId = "cds_info_opt",
              label = "Info",
              choices = c("all", "chrom_id", "cds_id", "cds_range"),
              selected = "all"
            )
          )
        ), column(
          width = 9,
          bs4Card(
            title = "Data",
            status = "primary",
            solidHeader = TRUE,
            width = 12,
            DT::dataTableOutput("table_extract_cds")
          )
        )
      )),
      bs4TabItem(tabName = "extract_exons", fluidRow(
        column(
          width = 3,
          bs4Card(
            title = "Parameters",
            status = "info",
            solidHeader = TRUE,
            width = 12,
            actionButton(
              inputId = "run_extract_exons",
              label = "Run",
              icon = icon("circle-play"),
              style = "width: 100%",
              class = "btn-block"
            ),
            br(),
            fileInput(inputId = "gff_file_extract_exons", label = "GFF/GTF file"),
            selectInput(
              inputId = "gff_format_extract_exons",
              label = "Format",
              choices = c("auto", "gff3", "gtf"),
              selected = "auto"
            ),
            selectInput(
              inputId = "exon_info_opt",
              label = "Info",
              choices = c("all", "chrom_id", "exon_id", "exon_range"),
              selected = "all"
            )
          )
        ), column(
          width = 9,
          bs4Card(
            title = "Data",
            status = "primary",
            solidHeader = TRUE,
            width = 12,
            DT::dataTableOutput("table_extract_exons")
          )
        )
      )),
      bs4TabItem(tabName = "extract_utr3", fluidRow(
        column(
          width = 3,
          bs4Card(
            title = "Parameters",
            status = "info",
            solidHeader = TRUE,
            width = 12,
            actionButton(
              inputId = "run_extract_utr3",
              label = "Run",
              icon = icon("circle-play"),
              style = "width: 100%",
              class = "btn-block"
            ),
            br(),
            fileInput(inputId = "gff_file_extract_utr3", label = "GFF/GTF file"),
            selectInput(
              inputId = "gff_format_extract_utr3",
              label = "Format",
              choices = c("auto", "gff3", "gtf"),
              selected = "auto"
            ),
            selectInput(
              inputId = "utr3_info_opt",
              label = "Info",
              choices = c("all", "chrom_id", "utr3_range"),
              selected = "all"
            )
          )
        ), column(
          width = 9,
          bs4Card(
            title = "Data",
            status = "primary",
            solidHeader = TRUE,
            width = 12,
            DT::dataTableOutput("table_extract_utr3")
          )
        )
      ))
    )
  )
)

server <- function(input, output, session) {
  observeEvent(TRUE, {
    updateTabItems(session, "main_menu", "plot_interval_structure")
  }, ignoreInit = FALSE, once = TRUE)

  getGff <- function(infile) {
    if (!is.null(infile))
      infile$datapath
    else
      system.file("extdata", "example.gff3.gz", package = "GAnnoViz")
  }
  getDeg <- function(infile) {
    if (!is.null(infile))
      infile$datapath
    else
      system.file("extdata", "example.deg", package = "GAnnoViz")
  }
  getFst <- function(infile) {
    if (!is.null(infile))
      infile$datapath
    else
      system.file("extdata", "example.fst", package = "GAnnoViz")
  }
  getDmr <- function(infile) {
    if (!is.null(infile))
      infile$datapath
    else
      system.file("extdata", "example.dmr", package = "GAnnoViz")
  }
  readGeneTable <- function(infile) {
    if (is.null(infile)) {
      data.frame(
        gene_id = c("ENSMUSG00000042414", "ENSMUSG00000025935", "ENSMUSG00000048701", "ENSMUSG00000035385"),
        gene_name = c("Prdm14", "Tram1", "Ccdc6", "Ccl2"),
        stringsAsFactors = FALSE
      )
    } else {
      utils::read.table(
        infile$datapath,
        header = TRUE,
        sep = "\t",
        stringsAsFactors = FALSE,
        check.names = FALSE
      )
    }
  }
  grToDf <- function(gr) {
    if (is.null(gr))
      return(data.frame())
    data.frame(
      chrom = as.character(GenomicRanges::seqnames(gr)),
      start = BiocGenerics::start(gr),
      end = BiocGenerics::end(gr),
      as.data.frame(S4Vectors::mcols(gr)),
      stringsAsFactors = FALSE
    )
  }

  dfGeneFeatures <- function(gff_file,
                             format,
                             gene_id,
                             upstream,
                             downstream) {
    txdb <- suppressWarnings(txdbmaker::makeTxDbFromGFF(file = gff_file, format = format))
    genes <- suppressWarnings(GenomicFeatures::genes(txdb))
    if (!(gene_id %in% genes$gene_id))
      return(data.frame())
    exons_by_tx <- suppressWarnings(GenomicFeatures::exonsBy(txdb, by = "tx"))
    introns_by_tx <- suppressWarnings(GenomicFeatures::intronsByTranscript(txdb, use.names = TRUE))
    utr5_by_tx <- suppressWarnings(GenomicFeatures::fiveUTRsByTranscript(txdb, use.names = TRUE))
    utr3_by_tx <- suppressWarnings(GenomicFeatures::threeUTRsByTranscript(txdb, use.names = TRUE))
    promoters_tx <- suppressWarnings(
      GenomicFeatures::promoters(
        txdb,
        upstream = upstream,
        downstream = downstream,
        use.names = TRUE
      )
    )
    tx_by_gene <- suppressWarnings(GenomicFeatures::transcriptsBy(txdb, by = "gene"))
    tx_names <- if (gene_id %in% names(tx_by_gene))
      tx_by_gene[[gene_id]]$tx_name
    else
      character(0)
    build_df <- function(gr_list, feature) {
      if (is.null(gr_list) || length(gr_list) == 0)
        return(data.frame())
      lst <- split(gr_list, gr_list$tx_name)
      lst <- lst[names(lst) %in% tx_names]
      if (length(lst) == 0)
        return(data.frame())
      dfs <- lapply(names(lst), function(nm) {
        gr <- lst[[nm]]
        data.frame(
          tx_name = nm,
          feature = feature,
          start = BiocGenerics::start(gr),
          end = BiocGenerics::end(gr),
          strand = as.character(GenomicRanges::strand(gr)),
          stringsAsFactors = FALSE
        )
      })
      do.call(rbind, dfs)
    }
    df <- do.call(rbind,
                  list(
                    build_df(exons_by_tx, "exon"),
                    build_df(introns_by_tx, "intron"),
                    build_df(utr5_by_tx, "utr5"),
                    build_df(utr3_by_tx, "utr3"),
                    build_df(promoters_tx, "promoter")
                  ))
    df
  }

  plot_gene_structure_ev <- eventReactive(input$run_plot_gene_structure, {
    plot_gene_structure(
      gff_file = getGff(input$gff_file_gene),
      format = input$gff_format_gene,
      gene_id = input$gene_id,
      upstream = input$upstream,
      downstream = input$downstream,
      feature_alpha = input$feature_alpha,
      intron_width = input$intron_width,
      x_breaks = input$x_breaks,
      arrow_length = input$arrow_length,
      arrow_count = input$arrow_count,
      arrow_unit = input$arrow_unit,
      promoter_color = input$promoter_color,
      utr5_color = input$utr5_color,
      utr3_color = input$utr3_color,
      exon_color = input$exon_color,
      intron_color = input$intron_color
    )
  })
  output$plot_gene_structure <- renderPlot({
    req(plot_gene_structure_ev())
    print(plot_gene_structure_ev())
  })

  plot_interval_structure_ev <- eventReactive(input$run_plot_interval_structure, {
    plot_interval_structure(
      gff_file = getGff(input$gff_file_interval),
      format = input$gff_format_interval,
      chrom_id = input$chrom_id,
      start = input$win_start,
      end = input$win_end,
      x_breaks = input$x_breaks_interval,
      upstream = input$upstream_interval,
      downstream = input$downstream_interval,
      feature_alpha = input$feature_alpha_interval,
      intron_width = input$intron_width_interval,
      arrow_length = input$arrow_length_interval,
      arrow_count = input$arrow_count_interval,
      arrow_unit = input$arrow_unit_interval,
      promoter_color = input$promoter_color_interval,
      utr5_color = input$utr5_color_interval,
      utr3_color = input$utr3_color_interval,
      exon_color = input$exon_color_interval,
      intron_color = input$intron_color_interval
    )
  })
  output$plot_interval_structure <- renderPlot({
    req(plot_interval_structure_ev())
    print(plot_interval_structure_ev())
  })

  plot_interval_flank_ev <- eventReactive(input$run_plot_interval_flank, {
    plot_interval_flank(
      gff_file = getGff(input$gff_file_flank),
      format = input$gff_format_flank,
      gene_id = input$flank_gene_id,
      flank_upstream = input$flank_upstream,
      flank_downstream = input$flank_downstream,
      show_promoters = input$show_promoters,
      upstream = input$upstream_flank,
      downstream = input$downstream_flank,
      arrow_length = input$arrow_length_flank,
      arrow_unit = input$arrow_unit_flank,
      gene_color = input$gene_color_flank,
      promoter_color = input$promoter_color_flank,
      label_size = input$label_size_flank
    )
  })
  output$plot_interval_flank <- renderPlot({
    req(plot_interval_flank_ev())
    print(plot_interval_flank_ev())
  })

  plot_chrom_structure_ev <- eventReactive(input$run_plot_chrom_structure, {
    plot_chrom_structure(
      gff_file = getGff(input$gff_file_chrom),
      format = input$gff_format_chrom,
      orientation = input$chrom_orientation,
      bar_width = input$bar_width,
      chrom_alpha = input$chrom_alpha,
      gene_width = input$gene_width,
      chrom_color = input$chrom_color,
      gene_color = input$gene_color_chrom,
      telomere_color = input$telomere_color,
      label_size = input$label_size_chrom
    )
  })
  output$plot_chrom_structure <- renderPlot({
    req(plot_chrom_structure_ev())
    print(plot_chrom_structure_ev())
  })

  plot_chrom_genes_ev <- eventReactive(input$run_plot_chrom_genes, {
    plot_chrom_genes(
      gff_file = getGff(input$gff_file_chrom_genes),
      gene_table = readGeneTable(input$gene_table_file),
      format = input$gff_format_chrom_genes,
      annotate = input$annotate_mode,
      orientation = input$chrom_genes_orientation,
      bar_width = input$bar_width_genes,
      chrom_alpha = input$chrom_alpha_genes,
      gene_width = input$gene_width_genes,
      chrom_color = input$chrom_color_genes,
      gene_color = input$gene_color_genes,
      telomere_color = input$telomere_color_genes,
      label_size = input$label_size_genes,
      connector_dx1 = input$connector_dx1_genes,
      connector_dx2 = input$connector_dx2_genes,
      min_gap_frac = input$min_gap_frac
    )
  })
  output$plot_chrom_genes <- renderPlot({
    req(plot_chrom_genes_ev())
    print(plot_chrom_genes_ev())
  })

  plot_chrom_heatmap_ev <- eventReactive(input$run_plot_chrom_heatmap, {
    plot_chrom_heatmap(
      gff_file = getGff(input$gff_file_heatmap),
      format = input$gff_format_heatmap,
      feature = input$feature_type,
      bin_size = input$bin_size_heatmap,
      orientation = input$orientation_heatmap,
      palette = c(input$palette_start, input$palette_end),
      alpha = input$alpha_heatmap
    )
  })
  output$plot_chrom_heatmap <- renderPlot({
    req(plot_chrom_heatmap_ev())
    print(plot_chrom_heatmap_ev())
  })

  plot_deg_chrom_ev <- eventReactive(input$run_plot_deg_chrom, {
    plot_deg_chrom(
      deg_file = getDeg(input$deg_file),
      gff_file = getGff(input$gff_file_deg),
      format = input$gff_format_deg,
      id_col = input$id_col,
      fc_col = input$fc_col,
      violin_scale = input$violin_scale,
      violin_border = input$violin_border,
      point_shape = input$point_shape_deg,
      point_size = input$point_size_deg,
      jitter_width = input$jitter_width_deg,
      hyper_color = input$hyper_color_deg,
      hypo_color = input$hypo_color_deg
    )
  })
  output$plot_deg_chrom <- renderPlot({
    req(plot_deg_chrom_ev())
    print(plot_deg_chrom_ev())
  })
  output$download_plot_deg_chrom <- downloadHandler(
    filename = function() {
      fmt <- input$plot_format_deg_chrom
      sprintf("plot_deg_chrom.%s", fmt)
    },
    content = function(file) {
      p <- plot_deg_chrom_ev()
      req(p)
      ggplot2::ggsave(
        filename = file,
        plot = p,
        width = input$plot_width_deg_chrom,
        height = input$plot_height_deg_chrom,
        dpi = input$plot_dpi_deg_chrom,
        device = input$plot_format_deg_chrom
      )
    }
  )

  plot_deg_exp_ev <- eventReactive(input$run_plot_deg_exp, {
    plot_deg_exp(
      deg_file = getDeg(input$deg_file_exp),
      gff_file = getGff(input$gff_file_deg_exp),
      format = input$gff_format_deg_exp,
      id_col = input$id_col_exp,
      fc_col = input$fc_col_exp,
      orientation = input$orientation_deg_exp,
      chrom_alpha = input$chrom_alpha_deg_exp,
      chrom_color = input$chrom_color_deg_exp,
      bar_height = input$bar_height_deg_exp,
      point_size = input$point_size_deg_exp,
      point_alpha = input$point_alpha_deg_exp,
      up_color = input$up_color_deg_exp,
      down_color = input$down_color_deg_exp,
      mark_style = input$mark_style_deg_exp,
      line_width = input$line_width_deg_exp,
      line_height = input$line_height_deg_exp
    )
  })
  output$plot_deg_exp <- renderPlot({
    req(plot_deg_exp_ev())
    print(plot_deg_exp_ev())
  })
  output$download_plot_deg_exp <- downloadHandler(
    filename = function() {
      fmt <- input$plot_format_deg_exp
      sprintf("plot_deg_exp.%s", fmt)
    },
    content = function(file) {
      p <- plot_deg_exp_ev()
      req(p)
      ggplot2::ggsave(
        filename = file,
        plot = p,
        width = input$plot_width_deg_exp,
        height = input$plot_height_deg_exp,
        dpi = input$plot_dpi_deg_exp,
        device = input$plot_format_deg_exp
      )
    }
  )

  plot_deg_volcano_ev <- eventReactive(input$run_plot_deg_volcano, {
    plot_deg_volcano(
      deg_file = getDeg(input$deg_file_volcano),
      id_col = input$id_col_volcano,
      fc_col = input$fc_col_volcano,
      sig_col = input$sig_col_volcano,
      fc_threshold = input$fc_threshold_volcano,
      sig_threshold = input$sig_threshold_volcano,
      point_size = input$point_size_volcano,
      point_alpha = input$point_alpha_volcano,
      up_color = input$up_color_volcano,
      down_color = input$down_color_volcano,
      ns_color = input$ns_color_volcano
    )
  })
  output$plot_deg_volcano <- renderPlot({
    req(plot_deg_volcano_ev())
    print(plot_deg_volcano_ev())
  })
  output$download_plot_deg_volcano <- downloadHandler(
    filename = function() {
      fmt <- input$plot_format_deg_volcano
      sprintf("plot_deg_volcano.%s", fmt)
    },
    content = function(file) {
      p <- plot_deg_volcano_ev()
      req(p)
      ggplot2::ggsave(
        filename = file,
        plot = p,
        width = input$plot_width_deg_volcano,
        height = input$plot_height_deg_volcano,
        dpi = input$plot_dpi_deg_volcano,
        device = input$plot_format_deg_volcano
      )
    }
  )

  plot_snp_fst_ev <- eventReactive(input$run_plot_snp_fst, {
    plot_snp_fst(
      fst_file = getFst(input$fst_file_heat),
      bin_size = input$bin_size_fst,
      metric = input$metric,
      orientation = input$orientation_fst,
      palette = c(input$palette_start_fst, input$palette_end_fst),
      alpha = input$alpha_fst
    )
  })
  output$plot_snp_fst <- renderPlot({
    req(plot_snp_fst_ev())
    print(plot_snp_fst_ev())
  })

  plot_snp_anno_ev <- eventReactive(input$run_plot_snp_anno, {
    plot_snp_anno(
      fst_file = getFst(input$fst_file_anno),
      gff_file = getGff(input$gff_file_fst),
      format = input$gff_format_fst,
      chrom_id = input$chrom_id_fst,
      top_n = input$top_n,
      orientation = input$orientation_fst_anno,
      smooth_span = input$smooth_span,
      fst_color = input$fst_color,
      point_size = input$point_size_fst,
      point_alpha = input$point_alpha_fst,
      label_size = input$label_size_fst,
      connector_dx1 = input$connector_dx1,
      connector_dx2 = input$connector_dx2,
      gap_frac = input$gap_frac
    )
  })
  output$plot_snp_anno <- renderPlot({
    req(plot_snp_anno_ev())
    print(plot_snp_anno_ev())
  })

  plot_dmg_chrom_ev <- eventReactive(input$run_plot_dmg_chrom, {
    plot_dmg_chrom(
      dmr_file = getDmr(input$dmr_file),
      violin_scale = input$violin_scale_dmr,
      violin_border = input$violin_border_dmr,
      point_shape = input$point_shape_dmr,
      point_size = input$point_size_dmr,
      jitter_width = input$jitter_width_dmr,
      hyper_color = input$hyper_color_dmr,
      hypo_color = input$hypo_color_dmr
    )
  })
  output$plot_dmg_chrom <- renderPlot({
    req(plot_dmg_chrom_ev())
    print(plot_dmg_chrom_ev())
  })

  plot_dmg_trend_ev <- eventReactive(input$run_plot_dmg_trend, {
    plot_dmg_trend(
      chrom_id = input$chrom_id_trend,
      dmr_file = getDmr(input$dmr_file_trend),
      smooth_span = input$smooth_span_trend,
      hyper_color = input$hyper_color_trend,
      hypo_color = input$hypo_color_trend,
      point_size = input$point_size_trend,
      point_alpha = input$point_alpha_trend
    )
  })
  output$plot_dmg_trend <- renderPlot({
    req(plot_dmg_trend_ev())
    print(plot_dmg_trend_ev())
  })

  plot_dmg_exp_ev <- eventReactive(input$run_plot_dmg_exp, {
    plot_dmg_exp(
      dmr_file = getDmr(input$dmr_file_exp),
      orientation = input$orientation_dmg_exp,
      chrom_alpha = input$chrom_alpha_dmg_exp,
      chrom_color = input$chrom_color_dmg_exp,
      bar_height = input$bar_height_dmg_exp,
      point_size = input$point_size_dmg_exp,
      point_alpha = input$point_alpha_dmg_exp,
      hyper_color = input$hyper_color_dmg_exp,
      hypo_color = input$hypo_color_dmg_exp,
      mark_style = input$mark_style_dmg_exp,
      line_width = input$line_width_dmg_exp,
      line_height = input$line_height_dmg_exp
    )
  })
  output$plot_dmg_exp <- renderPlot({
    req(plot_dmg_exp_ev())
    print(plot_dmg_exp_ev())
  })

  plot_gene_stats_ev <- eventReactive(input$run_plot_gene_stats, {
    plot_gene_stats(
      gff_file = getGff(input$gff_file_gene_stats),
      format = input$gff_format_gene_stats,
      bar_width = input$bar_width_gene_stats,
      bar_color = input$bar_color_gene_stats,
      lable_size = input$label_size_gene_stats
    )
  })
  output$plot_gene_stats <- renderPlot({
    req(plot_gene_stats_ev())
    print(plot_gene_stats_ev())
  })

  plot_snp_density_ev <- eventReactive(input$run_plot_snp_density, {
    plot_snp_density(
      fst_file = getFst(input$fst_file_density),
      LOG10 = input$log10_density,
      bin_size = input$bin_size_density,
      density_color = c(
        input$density_color1,
        input$density_color2,
        input$density_color3
      )
    )
  })
  output$plot_snp_density <- renderPlot({
    req(plot_snp_density_ev())
    print(plot_snp_density_ev())
  })

  anno_fst_dmr_ev <- eventReactive(input$run_anno_fst_dmr, {
    anno_fst_dmr(
      gff_file = getGff(input$gff_file_anno_ranges),
      format = input$gff_format_anno_ranges,
      genomic_ranges = if (!is.null(input$genomic_ranges_file))
        input$genomic_ranges_file$datapath
      else
        getFst(NULL),
      chrom_col = input$chrom_col_ranges,
      start_col = input$start_col_ranges,
      end_col = input$end_col_ranges,
      upstream = input$upstream_ranges,
      downstream = input$downstream_ranges,
      ignore_strand = input$ignore_strand_ranges,
      features = input$features_ranges
    )
  })
  output$table_anno_fst_dmr <- DT::renderDataTable({
    req(anno_fst_dmr_ev())
    DT::datatable(
      anno_fst_dmr_ev(),
      options = list(pageLength = 10, scrollX = TRUE),
      rownames = FALSE
    )
  })

  anno_deg_chrom_ev <- eventReactive(input$run_anno_deg_chrom, {
    anno_deg_chrom(
      deg_file = getDeg(input$deg_file_anno),
      gff_file = getGff(input$gff_file_deg_anno),
      format = input$gff_format_deg_anno,
      id_col = input$id_col_anno,
      fc_col = input$fc_col_anno,
      use_strand = input$use_strand_anno,
      drop_unmapped = input$drop_unmapped_anno
    )
  })
  output$table_anno_deg_chrom <- DT::renderDataTable({
    req(anno_deg_chrom_ev())
    DT::datatable(
      anno_deg_chrom_ev(),
      options = list(pageLength = 10, scrollX = TRUE),
      rownames = FALSE
    )
  })
  output$download_table_anno_deg_chrom <- downloadHandler(
    filename = function() sprintf("anno_deg_chrom.txt"),
    content = function(file) {
      df <- anno_deg_chrom_ev()
      if (is.null(df)) df <- data.frame()
      utils::write.table(df, file, sep = "\t", row.names = FALSE, quote = FALSE)
    }
  )

  plot_gene_domains_ev <- eventReactive(input$run_plot_gene_domains, {
    plot_gene_domains(
      gene_name = if (nzchar(input$transcript_id_domains)) NULL else input$gene_name_domains,
      species = input$species_domains,
      transcript_id = if (nzchar(input$transcript_id_domains)) input$transcript_id_domains else NULL,
      transcript_choice = input$transcript_choice_domains,
      palette = input$palette_domains,
      legend_ncol = input$legend_ncol_domains,
      return_data = FALSE
    )
  })
  output$plot_gene_domains <- renderPlot({
    req(plot_gene_domains_ev())
    print(plot_gene_domains_ev())
  })

  extract_promoters_ev <- eventReactive(input$run_extract_promoters, {
    extract_promoters(
      gff_file = getGff(input$gff_file_extract_promoters),
      format = input$gff_format_extract_promoters,
      upstream = input$upstream_extract_promoters,
      downstream = input$downstream_extract_promoters,
      promoter_info = input$promoter_info
    )
  })
  output$table_extract_promoters <- DT::renderDataTable({
    res <- extract_promoters_ev()
    if (is.null(res))
      return(DT::datatable(data.frame()))
    if (is(res, "GRanges")) {
      DT::datatable(
        grToDf(res),
        options = list(pageLength = 10, scrollX = TRUE),
        rownames = FALSE
      )
    } else if (is(res, "IRanges")) {
      DT::datatable(
        as.data.frame(res),
        options = list(pageLength = 10, scrollX = TRUE),
        rownames = FALSE
      )
    } else {
      DT::datatable(
        data.frame(value = res),
        options = list(pageLength = 10, scrollX = TRUE),
        rownames = FALSE
      )
    }
  })

  extract_utr5_ev <- eventReactive(input$run_extract_utr5, {
    extract_utr5(
      gff_file = getGff(input$gff_file_extract_utr5),
      format = input$gff_format_extract_utr5,
      utr5_info = input$utr5_info
    )
  })
  output$table_extract_utr5 <- DT::renderDataTable({
    res <- extract_utr5_ev()
    if (is.null(res))
      return(DT::datatable(data.frame()))
    if (is(res, "GRanges")) {
      DT::datatable(
        grToDf(res),
        options = list(pageLength = 10, scrollX = TRUE),
        rownames = FALSE
      )
    } else if (is(res, "IRanges")) {
      DT::datatable(
        as.data.frame(res),
        options = list(pageLength = 10, scrollX = TRUE),
        rownames = FALSE
      )
    } else {
      DT::datatable(
        data.frame(value = res),
        options = list(pageLength = 10, scrollX = TRUE),
        rownames = FALSE
      )
    }
  })

  extract_genes_ev <- eventReactive(input$run_extract_genes, {
    extract_genes(
      gff_file = getGff(input$gff_file_extract_genes),
      format = input$gff_format_extract_genes,
      gene_info = input$gene_info_opt
    )
  })
  output$table_extract_genes <- DT::renderDataTable({
    res <- extract_genes_ev()
    if (is.null(res))
      return(DT::datatable(data.frame()))
    if (is(res, "GRanges")) {
      DT::datatable(
        grToDf(res),
        options = list(pageLength = 10, scrollX = TRUE),
        rownames = FALSE
      )
    } else if (is(res, "IRanges")) {
      DT::datatable(
        as.data.frame(res),
        options = list(pageLength = 10, scrollX = TRUE),
        rownames = FALSE
      )
    } else {
      DT::datatable(
        data.frame(value = res),
        options = list(pageLength = 10, scrollX = TRUE),
        rownames = FALSE
      )
    }
  })

  extract_mrnas_ev <- eventReactive(input$run_extract_mrnas, {
    extract_mrnas(
      gff_file = getGff(input$gff_file_extract_mrnas),
      format = input$gff_format_extract_mrnas,
      mrna_info = input$mrna_info_opt
    )
  })
  output$table_extract_mrnas <- DT::renderDataTable({
    res <- extract_mrnas_ev()
    if (is.null(res))
      return(DT::datatable(data.frame()))
    if (is(res, "GRanges")) {
      DT::datatable(
        grToDf(res),
        options = list(pageLength = 10, scrollX = TRUE),
        rownames = FALSE
      )
    } else if (is(res, "IRanges")) {
      DT::datatable(
        as.data.frame(res),
        options = list(pageLength = 10, scrollX = TRUE),
        rownames = FALSE
      )
    } else {
      DT::datatable(
        data.frame(value = res),
        options = list(pageLength = 10, scrollX = TRUE),
        rownames = FALSE
      )
    }
  })

  extract_cds_ev <- eventReactive(input$run_extract_cds, {
    extract_cds(
      gff_file = getGff(input$gff_file_extract_cds),
      format = input$gff_format_extract_cds,
      cds_info = input$cds_info_opt
    )
  })
  output$table_extract_cds <- DT::renderDataTable({
    res <- extract_cds_ev()
    if (is.null(res))
      return(DT::datatable(data.frame()))
    if (is(res, "GRanges")) {
      DT::datatable(
        grToDf(res),
        options = list(pageLength = 10, scrollX = TRUE),
        rownames = FALSE
      )
    } else if (is(res, "IRanges")) {
      DT::datatable(
        as.data.frame(res),
        options = list(pageLength = 10, scrollX = TRUE),
        rownames = FALSE
      )
    } else {
      DT::datatable(
        data.frame(value = res),
        options = list(pageLength = 10, scrollX = TRUE),
        rownames = FALSE
      )
    }
  })

  extract_exons_ev <- eventReactive(input$run_extract_exons, {
    extract_exons(
      gff_file = getGff(input$gff_file_extract_exons),
      format = input$gff_format_extract_exons,
      exon_info = input$exon_info_opt
    )
  })
  output$table_extract_exons <- DT::renderDataTable({
    res <- extract_exons_ev()
    if (is.null(res))
      return(DT::datatable(data.frame()))
    if (is(res, "GRanges")) {
      DT::datatable(
        grToDf(res),
        options = list(pageLength = 10, scrollX = TRUE),
        rownames = FALSE
      )
    } else if (is(res, "IRanges")) {
      DT::datatable(
        as.data.frame(res),
        options = list(pageLength = 10, scrollX = TRUE),
        rownames = FALSE
      )
    } else {
      DT::datatable(
        data.frame(value = res),
        options = list(pageLength = 10, scrollX = TRUE),
        rownames = FALSE
      )
    }
  })

  extract_utr3_ev <- eventReactive(input$run_extract_utr3, {
    extract_utr3(
      gff_file = getGff(input$gff_file_extract_utr3),
      format = input$gff_format_extract_utr3,
      utr3_info = input$utr3_info_opt
    )
  })
  output$table_extract_utr3 <- DT::renderDataTable({
    res <- extract_utr3_ev()
    if (is.null(res))
      return(DT::datatable(data.frame()))
    if (is(res, "GRanges")) {
      DT::datatable(
        grToDf(res),
        options = list(pageLength = 10, scrollX = TRUE),
        rownames = FALSE
      )
    } else if (is(res, "IRanges")) {
      DT::datatable(
        as.data.frame(res),
        options = list(pageLength = 10, scrollX = TRUE),
        rownames = FALSE
      )
    } else {
      DT::datatable(
        data.frame(value = res),
        options = list(pageLength = 10, scrollX = TRUE),
        rownames = FALSE
      )
    }
  })
}

shinyApp(ui, server)
