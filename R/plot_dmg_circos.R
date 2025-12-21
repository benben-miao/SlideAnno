#' @title Plot DMGs circos across chromosomes
#' @description Plot genome-wide DMGs \bold{\emph{circos}} with separate tracks for hyper/hypo and top-gene annotations.
#' @author benben-miao
#'
#' @return Plot of DMGs circos across chromosomes.
#' @param dmr_file DEG table from \bold{\emph{MethylKit}} analysis.
#' @param gff_file Genomic structural annotation \bold{\file{GFF3/GTF}} file path.
#' @param format Format of GFF3/GTF file. (\bold{\emph{"auto"}}, "gff3", "gtf").
#' @param label_type Label by gene name or gene id. (\bold{\emph{"name"}}, "id").
#' @param gene_table Optional gene ID/name mapping table (first two columns: id, name).
#' @param y_transform Y-axis transformation for \code{meth.diff}. (\bold{\emph{"none"}}, "log2", "log10").
#'
#' @param chrom_height Chrom track height. (\bold{\emph{0.08}}).
#' @param chrom_color Chromosome track fill color. (\bold{\emph{"#eeeeee"}}).
#' @param chrom_border Chromosome track border color. (\bold{\emph{"#333333"}}).
#' @param chrom_cex Chrom name cex. (\bold{\emph{0.8}}).
#' @param gap_degree Gap degree between chromosomes. (\bold{\emph{1}}).
#' @param x_tick_by Chromosome axis tick step (bp). (\bold{\emph{2e7}}).
#' @param axis_cex Axis text cex. (\bold{\emph{0.5}}).
#' @param last_gap_degree Gap degree between last and first chromosome. (\bold{\emph{3}}).
#'
#' @param scatter_height Scatter track height. (\bold{\emph{0.15}}).
#' @param top_up Number of top hyper-methylated DMGs to annotate. (\bold{\emph{30}}).
#' @param top_down Number of top hypo-methylated DMGs to annotate. (\bold{\emph{30}}).
#' @param up_color Point color for hyper-methylated DMGs. (\bold{\emph{"#ff0000"}}).
#' @param down_color Point color for hypo-methylated DMGs. (\bold{\emph{"#008800"}}).
#' @param point_cex Point size (cex). (\bold{\emph{0.5}}).
#'
#' @param annotation_height Track height for annotation rings. (\bold{\emph{0.15}}).
#' @param label_cex Label size (cex) for top-gene annotations. (\bold{\emph{0.5}}).
#' @param label_rotate Label rotation degree (counterclockwise). (\bold{\emph{0}}).
#' @param label_font Font for chromosome and gene labels. (\bold{\emph{1}}).
#' @param connector_lwd Line width for point-label connectors. (\bold{\emph{0.5}}).
#' @param connector_col Color for point-label connectors. (\bold{\emph{"#333333"}}).
#' @param connector_len Connector length. (\bold{\emph{0.2}}).
#' @param connector_elbow Connector elbow. (\bold{\emph{0.8}}).
#'
#'
#' @importFrom magrittr %>%
#' @export
#'
#' @examples
#' dmr_file <- system.file(
#'   "extdata",
#'   "example.dmr",
#'   package = "GAnnoViz")
#'
#' gff_file <- system.file(
#'   "extdata",
#'   "example.gff3.gz",
#'   package = "GAnnoViz")
#'
#' plot_dmg_circos(
#'   dmr_file = dmr_file,
#'   gff_file = gff_file,
#'   format = "auto",
#'   label_type = "name",
#'   gene_table = NULL,
#'   y_transform = "none",
#'   chrom_height = 0.08,
#'   chrom_color = "#eeeeee",
#'   chrom_border = "#cccccc",
#'   chrom_cex = 0.8,
#'   gap_degree = 1,
#'   x_tick_by = 2e7,
#'   axis_cex = 0.5,
#'   last_gap_degree = 3,
#'   scatter_height = 0.15,
#'   top_up = 30,
#'   top_down = 30,
#'   up_color = "#ff0000",
#'   down_color = "#008800",
#'   point_cex = 0.5,
#'   annotation_height = 0.18,
#'   label_cex = 0.5,
#'   label_rotate = 0,
#'   label_font = 1,
#'   connector_lwd = 0.5,
#'   connector_col = "#333333",
#'   connector_len = 0.2,
#'   connector_elbow = 0.8
#' )
#'
plot_dmg_circos <- function(dmr_file,
                            gff_file = NULL,
                            format = "auto",
                            label_type = "name",
                            gene_table = NULL,
                            y_transform = "none",
                            chrom_height = 0.08,
                            chrom_color = "#eeeeee",
                            chrom_border = "#cccccc",
                            chrom_cex = 0.8,
                            gap_degree = 1,
                            x_tick_by = 2e7,
                            axis_cex = 0.5,
                            last_gap_degree = 3,
                            scatter_height = 0.15,
                            top_up = 30,
                            top_down = 30,
                            up_color = "#ff0000",
                            down_color = "#008800",
                            point_cex = 0.5,
                            annotation_height = 0.18,
                            label_cex = 0.5,
                            label_rotate = 0,
                            label_font = 1,
                            connector_lwd = 0.5,
                            connector_col = "#333333",
                            connector_len = 0.2,
                            connector_elbow = 0.8) {
  y_transform <- match.arg(y_transform)
  label_type <- match.arg(label_type)
  if (!requireNamespace("circlize", quietly = TRUE))
    stop("Package 'circlize' is required for plot_dmg_circos(). Please install it first.")
  point_alpha <- 0.5
  highlight_alpha <- 0.4
  start_degree <- 90
  op <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(op), add = TRUE)
  graphics::par(family = "sans", xpd = NA)
  connector_lwd <- suppressWarnings(as.numeric(connector_lwd))
  if (!is.finite(connector_lwd) || connector_lwd <= 0)
    connector_lwd <- 0.6
  label_rotate <- suppressWarnings(as.numeric(label_rotate))
  if (!is.finite(label_rotate))
    label_rotate <- 0
  connector_len <- suppressWarnings(as.numeric(connector_len))
  if (!is.finite(connector_len))
    connector_len <- 0.2
  connector_len <- max(0.05, min(0.95, connector_len))
  connector_elbow <- suppressWarnings(as.numeric(connector_elbow))
  if (!is.finite(connector_elbow))
    connector_elbow <- 0.8
  connector_elbow <- max(0.05, min(0.95, connector_elbow))
  axis_cex <- suppressWarnings(as.numeric(axis_cex))
  if (!is.finite(axis_cex) || axis_cex <= 0)
    axis_cex <- 0.45
  chrom_cex <- suppressWarnings(as.numeric(chrom_cex))
  if (!is.finite(chrom_cex) || chrom_cex <= 0)
    chrom_cex <- 0.65
  annotation_height <- suppressWarnings(as.numeric(annotation_height))
  if (!is.finite(annotation_height) || annotation_height <= 0)
    annotation_height <- 0.18
  scatter_height <- suppressWarnings(as.numeric(scatter_height))
  if (!is.finite(scatter_height) || scatter_height <= 0)
    scatter_height <- 0.13
  chrom_height <- suppressWarnings(as.numeric(chrom_height))
  if (!is.finite(chrom_height) || chrom_height <= 0)
    chrom_height <- 0.075
  point_cex <- suppressWarnings(as.numeric(point_cex))
  if (!is.finite(point_cex) || point_cex <= 0)
    point_cex <- 0.8
  label_cex <- suppressWarnings(as.numeric(label_cex))
  if (!is.finite(label_cex) || label_cex <= 0)
    label_cex <- 0.8

  dmr <- utils::read.table(
    dmr_file,
    header = TRUE,
    sep = "\t",
    fill = TRUE,
    na.strings = "NA",
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  cols <- colnames(dmr)
  needed <- c("chr", "start", "meth.diff")
  missing <- setdiff(needed, cols)
  if (length(missing) > 0)
    stop(paste0("Missing required columns: ", paste(missing, collapse = ", ")))
  dmr <- dmr[stats::complete.cases(dmr$chr, dmr$start, dmr$meth.diff), , drop = FALSE]
  if (nrow(dmr) == 0)
    stop("No DMGs with valid positions found")

  chrom_raw <- as.character(dmr$chr)
  chr_prefix <- any(stringr::str_detect(chrom_raw, "^(?i)chr"))
  keep_mask <- if (chr_prefix)
    stringr::str_detect(chrom_raw, "^(?i)chr")
  else
    rep(TRUE, length(chrom_raw))
  dmr <- dmr[keep_mask, , drop = FALSE]
  if (nrow(dmr) == 0)
    stop("No DMGs found on selected chromosomes")

  dmr$end <- if ("end" %in% colnames(dmr))
    dmr$end
  else
    dmr$start
  dmr$pos <- (dmr$start + dmr$end) / 2

  chrom_chr <- as.character(dmr$chr)
  chrom_num <- stringr::str_extract(chrom_chr, "\\d+")
  suppressWarnings({
    chrom_num <- as.numeric(chrom_num)
  })
  order_index <- order(is.na(chrom_num), chrom_num, chrom_chr)
  chrom_levels <- unique(chrom_chr[order_index])

  length_map <- stats::aggregate(end ~ chr,
                                 data = data.frame(chr = dmr$chr, end = dmr$end),
                                 FUN = max)
  length_map <- length_map[match(chrom_levels, length_map$chr), , drop = FALSE]

  y_fun <- switch(
    y_transform,
    none = function(x)
      abs(x),
    log2 = function(x)
      log(1 + abs(x), base = 2),
    log10 = function(x)
      log10(1 + abs(x))
  )
  dmr$abs_y <- y_fun(dmr$meth.diff)
  dmr$state <- ifelse(dmr$meth.diff >= 0, "up", "down")

  gene_id_mapped <- NULL
  label_vec <- NULL
  if (label_type == "name") {
    auto_cols <- c("gene_name", "GeneName", "SYMBOL", "symbol", "gene")
    hit <- auto_cols[auto_cols %in% colnames(dmr)][1]
    if (!is.null(hit) && !is.na(hit))
      label_vec <- as.character(dmr[[hit]])
    if (is.null(label_vec)) {
      auto_cols <- c("gene_id", "GeneID", "ID", "gene")
      hit <- auto_cols[auto_cols %in% colnames(dmr)][1]
      if (!is.null(hit) && !is.na(hit)) {
        gene_id_mapped <- as.character(dmr[[hit]])
        label_vec <- gene_id_mapped
      }
    }
  }
  if (is.null(label_vec) && label_type == "id") {
    auto_cols <- c("gene_id", "GeneID", "ID", "gene")
    hit <- auto_cols[auto_cols %in% colnames(dmr)][1]
    if (!is.null(hit) && !is.na(hit))
      label_vec <- as.character(dmr[[hit]])
  }

  if (is.null(label_vec) &&
      !missing(gff_file) && !is.null(gff_file)) {
    genes <- extract_genes(gff_file = gff_file,
                           format = format,
                           gene_info = "all")
    genes <- genes[as.character(GenomicRanges::seqnames(genes)) %in% chrom_levels]
    if (length(genes) > 0) {
      gr_dmr <- GenomicRanges::GRanges(
        seqnames = dmr$chr,
        ranges = IRanges::IRanges(start = as.integer(dmr$start), end = as.integer(dmr$end))
      )
      hits <- GenomicRanges::findOverlaps(gr_dmr, genes, ignore.strand = TRUE)
      gene_id <- rep(NA_character_, length(gr_dmr))
      if (length(hits) > 0) {
        q <- S4Vectors::queryHits(hits)
        s <- S4Vectors::subjectHits(hits)
        keep_first <- !duplicated(q)
        q <- q[keep_first]
        s <- s[keep_first]
        gene_id[q] <- as.character(genes$gene_id[s])
      }
      na_idx <- which(is.na(gene_id) | gene_id == "")
      if (length(na_idx) > 0) {
        near <- GenomicRanges::nearest(gr_dmr[na_idx], genes, ignore.strand = TRUE)
        ok <- !is.na(near)
        if (any(ok))
          gene_id[na_idx[ok]] <- as.character(genes$gene_id[near[ok]])
      }

      gene_id_mapped <- gene_id
      label_vec <- gene_id
    }
  }

  if (label_type == "name" && !is.null(gene_id_mapped)) {
    gene_name <- rep(NA_character_, length(gene_id_mapped))
    if (!is.null(gene_table)) {
      gene_table <- as.data.frame(gene_table, stringsAsFactors = FALSE)
      if (ncol(gene_table) >= 2) {
        names(gene_table)[1:2] <- c("gene_id", "gene_name")
        idx <- match(gene_id_mapped, gene_table$gene_id)
        ok <- !is.na(idx)
        if (any(ok))
          gene_name[ok] <- as.character(gene_table$gene_name[idx[ok]])
      }
    }
    if (all(is.na(gene_name) |
            gene_name == "") &&
        !missing(gff_file) && !is.null(gff_file)) {
      infer_gene_name_map <- function(gff_file, ids, fmt) {
        ids <- unique(as.character(ids))
        ids <- ids[!is.na(ids) & ids != ""]
        if (length(ids) == 0)
          return(stats::setNames(character(0), character(0)))
        ids_norm <- sub("^.+?:", "", ids)
        wanted <- unique(ids_norm)
        map <- stats::setNames(rep(NA_character_, length(wanted)), wanted)
        con <- if (grepl("\\.gz$", gff_file, ignore.case = TRUE))
          gzfile(gff_file, open = "rt")
        else
          file(gff_file, open = "rt")
        on.exit(try(close(con), silent = TRUE)
                , add = TRUE)
        repeat {
          lines <- readLines(con, n = 50000, warn = FALSE)
          if (length(lines) == 0)
            break
          lines <- lines[!startsWith(lines, "#")]
          if (length(lines) == 0)
            next
          parts <- strsplit(lines, "\t", fixed = TRUE)
          attrs <- vapply(parts, function(z)
            if (length(z) >= 9)
              z[[9]]
            else
              NA_character_, character(1))
          attrs <- attrs[!is.na(attrs)]
          if (length(attrs) == 0)
            next
          if (identical(fmt, "gtf")) {
            gid <- stringr::str_match(attrs, 'gene_id\\s+"([^"]+)"')[, 2]
            gnm <- stringr::str_match(attrs, 'gene_name\\s+"([^"]+)"')[, 2]
            ok <- !is.na(gid) & gid != "" & !is.na(gnm) & gnm != ""
            if (any(ok)) {
              gid_norm <- sub("^.+?:", "", gid[ok])
              hit <- gid_norm %in% wanted
              if (any(hit)) {
                idx <- match(gid_norm[hit], names(map))
                fill <- is.na(map[idx]) | map[idx] == ""
                if (any(fill))
                  map[idx[fill]] <- gnm[ok][hit][fill]
              }
            }
          } else {
            gid <- stringr::str_match(attrs, "(?:^|;)ID=([^;]+)")[, 2]
            gnm <- stringr::str_match(attrs, "(?:^|;)Name=([^;]+)")[, 2]
            if (all(is.na(gnm) | gnm == "")) {
              gnm <- stringr::str_match(attrs, "(?:^|;)gene_name=([^;]+)")[, 2]
            }
            ok <- !is.na(gid) & gid != "" & !is.na(gnm) & gnm != ""
            if (any(ok)) {
              gid_norm <- sub("^.+?:", "", gid[ok])
              hit <- gid_norm %in% wanted
              if (any(hit)) {
                idx <- match(gid_norm[hit], names(map))
                fill <- is.na(map[idx]) | map[idx] == ""
                if (any(fill))
                  map[idx[fill]] <- gnm[ok][hit][fill]
              }
            }
          }
          if (all(!is.na(map) & map != ""))
            break
        }
        map[!is.na(map) & map != ""]
      }

      fmt <- resolve_gff_format(gff_file, format)
      gmap <- infer_gene_name_map(gff_file, gene_id_mapped, fmt)
      if (length(gmap) > 0) {
        gid_norm <- sub("^.+?:", "", as.character(gene_id_mapped))
        idx <- match(gid_norm, names(gmap))
        ok <- !is.na(idx)
        if (any(ok))
          gene_name[ok] <- unname(gmap[idx[ok]])
      }
    }
    label_vec <- ifelse(!is.na(gene_name) &
                          gene_name != "",
                        gene_name,
                        gene_id_mapped)
  }

  if (is.null(label_vec)) {
    label_vec <- paste0(dmr$chr, ":", dmr$start, "-", dmr$end)
  } else {
    interval_label <- paste0(dmr$chr, ":", dmr$start, "-", dmr$end)
    empty <- is.na(label_vec) | label_vec == ""
    label_vec[empty] <- interval_label[empty]
  }
  dmr$label <- label_vec

  top_up_n <- max(0, as.integer(top_up))
  top_down_n <- max(0, as.integer(top_down))
  df_up <- dmr %>%
    dplyr::filter(is.finite(meth.diff) & meth.diff > 0) %>%
    dplyr::arrange(dplyr::desc(meth.diff))
  df_down <- dmr %>%
    dplyr::filter(is.finite(meth.diff) & meth.diff < 0) %>%
    dplyr::arrange(meth.diff)
  df_top_up <- utils::head(df_up, top_up_n)
  df_top_down <- utils::head(df_down, top_down_n)

  dmr$chr <- factor(dmr$chr, levels = chrom_levels)
  df_up$chr <- factor(df_up$chr, levels = chrom_levels)
  df_down$chr <- factor(df_down$chr, levels = chrom_levels)
  df_top_up$chr <- factor(df_top_up$chr, levels = chrom_levels)
  df_top_down$chr <- factor(df_top_down$chr, levels = chrom_levels)

  last_gap <- suppressWarnings(as.numeric(last_gap_degree))
  if (!is.finite(last_gap) || last_gap < 0)
    last_gap <- gap_degree * 3
  scale_sector <- ".scale"
  sectors <- c(chrom_levels, scale_sector)
  gap_after <- rep(gap_degree, length(sectors))
  if (length(sectors) >= 2) {
    half <- last_gap / 2
    if (!is.finite(half) || half < 0)
      half <- gap_degree
    gap_after[length(sectors) - 1] <- half
    gap_after[length(sectors)] <- half
  }

  circlize::circos.clear()
  on.exit(circlize::circos.clear(), add = TRUE)
  circlize::circos.par(
    start.degree = start_degree,
    gap.after = gap_after,
    cell.padding = c(0, 0, 0, 0),
    track.margin = c(0.005, 0.005),
    points.overflow.warning = FALSE
  )

  circlize::circos.initialize(factors = sectors, xlim = rbind(cbind(rep(
    0, nrow(length_map)
  ), length_map$end), c(0, 1)))

  circlize::circos.trackPlotRegion(
    factors = sectors,
    ylim = c(0, 1),
    track.height = chrom_height,
    bg.border = NA,
    panel.fun = function(x, y) {
      sector.index <- circlize::get.cell.meta.data("sector.index")
      xlim <- circlize::get.cell.meta.data("xlim")
      if (identical(sector.index, scale_sector))
        return()
      circlize::circos.rect(xlim[1],
                            0,
                            xlim[2],
                            1,
                            col = chrom_color,
                            border = chrom_border,
                            lwd = 0.6)
      tick_by <- suppressWarnings(as.numeric(x_tick_by))
      if (!is.finite(tick_by) || tick_by <= 0)
        tick_by <- 2e7
      tick_by <- max(2e7, tick_by)
      if (xlim[2] / tick_by > 15) {
        tick_by <- ceiling((xlim[2] / 15) / 1e7) * 1e7
        tick_by <- max(2e7, tick_by)
      }
      major.at <- seq(0, xlim[2], by = tick_by)
      if (length(major.at) == 0)
        major.at <- c(0, xlim[2])
      if (major.at[length(major.at)] < xlim[2])
        major.at <- c(major.at, xlim[2])
      mb_lab <- as.character(as.integer(round(major.at / 1e6)))
      circlize::circos.axis(
        h = "top",
        major.at = major.at,
        labels = mb_lab,
        labels.cex = axis_cex,
        major.tick.length = circlize::convert_y(0.3, "mm"),
        minor.ticks = 0
      )
      circlize::circos.text(
        x = mean(xlim),
        y = 0.5,
        labels = sector.index,
        cex = chrom_cex,
        facing = "inside",
        niceFacing = TRUE,
        font = label_font
      )
    }
  )
  track_chr <- circlize::get.current.track.index()

  ymax_all <- max(dmr$abs_y, na.rm = TRUE)
  if (!is.finite(ymax_all) || ymax_all == 0)
    ymax_all <- 1
  max_abs_raw <- max(abs(dmr$meth.diff), na.rm = TRUE)
  if (!is.finite(max_abs_raw) || max_abs_raw == 0)
    max_abs_raw <- 1
  tick_raw <- pretty(c(0, max_abs_raw), n = 4)
  tick_raw <- tick_raw[tick_raw >= 0 & tick_raw <= max_abs_raw]
  if (length(tick_raw) == 0)
    tick_raw <- c(0, max_abs_raw)
  tick_y <- y_fun(tick_raw)
  tick_ok <- is.finite(tick_y) & tick_y >= 0 & tick_y <= ymax_all
  tick_raw <- tick_raw[tick_ok]
  tick_y <- tick_y[tick_ok]
  label_y_up <- connector_len
  label_y_down <- 1 - connector_len

  if (nrow(df_top_up) > 0) {
    col_up_hl <- grDevices::adjustcolor(up_color, alpha.f = highlight_alpha)
    circlize::circos.trackPlotRegion(
      factors = sectors,
      ylim = c(0, 1),
      track.height = annotation_height,
      bg.border = NA,
      panel.fun = function(x, y) {
        sector.index <- circlize::get.cell.meta.data("sector.index")
        if (identical(sector.index, scale_sector))
          return()
        sub <- df_top_up[df_top_up$chr == sector.index, , drop = FALSE]
        if (nrow(sub) == 0)
          return()
        for (i in seq_len(nrow(sub))) {
          s <- sub$start[i]
          e <- sub$end[i]
          if (!is.finite(s) || !is.finite(e))
            next
          if (e < s) {
            tmp <- s
            s <- e
            e <- tmp
          }
          circlize::circos.rect(s, 0, e, 1, col = col_up_hl, border = NA)
        }
      }
    )
    track_up_anno <- circlize::get.current.track.index()
  }

  if (nrow(df_up) > 0) {
    col_up <- grDevices::adjustcolor(up_color, alpha.f = point_alpha)
    circlize::circos.trackPlotRegion(
      factors = sectors,
      ylim = c(0, ymax_all),
      track.height = scatter_height,
      bg.border = NA,
      panel.fun = function(x, y) {
        sector.index <- circlize::get.cell.meta.data("sector.index")
        if (identical(sector.index, scale_sector)) {
          if (length(tick_y) > 0) {
            circlize::circos.yaxis(
              side = "right",
              at = tick_y,
              labels = tick_raw,
              labels.cex = axis_cex,
              tick.length = circlize::convert_x(0.5, "mm")
            )
          }
          return()
        }
        sub <- df_up[df_up$chr == sector.index, , drop = FALSE]
        if (nrow(sub) == 0)
          return()
        circlize::circos.points(sub$pos,
                                sub$abs_y,
                                pch = 16,
                                cex = point_cex,
                                col = col_up)
      }
    )
    track_up_scatter <- circlize::get.current.track.index()
  }

  if (nrow(df_down) > 0) {
    col_down <- grDevices::adjustcolor(down_color, alpha.f = point_alpha)
    circlize::circos.trackPlotRegion(
      factors = sectors,
      ylim = c(0, ymax_all),
      track.height = scatter_height,
      bg.border = NA,
      panel.fun = function(x, y) {
        sector.index <- circlize::get.cell.meta.data("sector.index")
        if (identical(sector.index, scale_sector)) {
          if (length(tick_y) > 0) {
            circlize::circos.yaxis(
              side = "right",
              at = ymax_all - tick_y,
              labels = -tick_raw,
              labels.cex = axis_cex,
              tick.length = circlize::convert_x(0.5, "mm")
            )
          }
          return()
        }
        sub <- df_down[df_down$chr == sector.index, , drop = FALSE]
        if (nrow(sub) == 0)
          return()
        y_plot <- ymax_all - sub$abs_y
        circlize::circos.points(sub$pos,
                                y_plot,
                                pch = 16,
                                cex = point_cex,
                                col = col_down)
      }
    )
    track_down_scatter <- circlize::get.current.track.index()
  }

  if (nrow(df_top_down) > 0) {
    col_down_hl <- grDevices::adjustcolor(down_color, alpha.f = highlight_alpha)
    circlize::circos.trackPlotRegion(
      factors = sectors,
      ylim = c(0, 1),
      track.height = annotation_height,
      bg.border = NA,
      panel.fun = function(x, y) {
        sector.index <- circlize::get.cell.meta.data("sector.index")
        if (identical(sector.index, scale_sector))
          return()
        sub <- df_top_down[df_top_down$chr == sector.index, , drop = FALSE]
        if (nrow(sub) == 0)
          return()
        for (i in seq_len(nrow(sub))) {
          s <- sub$start[i]
          e <- sub$end[i]
          if (!is.finite(s) || !is.finite(e))
            next
          if (e < s) {
            tmp <- s
            s <- e
            e <- tmp
          }
          circlize::circos.rect(s, 0, e, 1, col = col_down_hl, border = NA)
        }
      }
    )
    track_down_anno <- circlize::get.current.track.index()
  }

  track_up_anno <- if (exists("track_up_anno", inherits = FALSE))
    track_up_anno
  else
    NA_integer_
  track_up_scatter <- if (exists("track_up_scatter", inherits = FALSE))
    track_up_scatter
  else
    NA_integer_
  track_down_scatter <- if (exists("track_down_scatter", inherits = FALSE))
    track_down_scatter
  else
    NA_integer_
  track_down_anno <- if (exists("track_down_anno", inherits = FALSE))
    track_down_anno
  else
    NA_integer_

  adjust_label_x <- function(x, x_max, min_gap) {
    x <- suppressWarnings(as.numeric(x))
    if (length(x) == 0)
      return(x)
    x[!is.finite(x)] <- NA_real_
    x <- pmin(pmax(x, 0), x_max)
    ord <- order(x, na.last = TRUE)
    x_sorted <- x[ord]
    if (sum(is.finite(x_sorted)) <= 1)
      return(x)
    idx <- which(is.finite(x_sorted))
    for (k in idx[-1]) {
      prev <- x_sorted[k - 1]
      cur <- x_sorted[k]
      if (is.finite(prev) &&
          is.finite(cur) && (cur - prev) < min_gap)
        x_sorted[k] <- prev + min_gap
    }
    if (max(x_sorted[idx], na.rm = TRUE) > x_max) {
      x_sorted[idx] <- pmin(x_sorted[idx], x_max)
      for (k in rev(idx[-length(idx)])) {
        nxt <- x_sorted[k + 1]
        cur <- x_sorted[k]
        if (is.finite(nxt) &&
            is.finite(cur) && (nxt - cur) < min_gap)
          x_sorted[k] <- nxt - min_gap
      }
      x_sorted[idx] <- pmax(x_sorted[idx], 0)
    }
    x_adj <- x
    x_adj[ord] <- x_sorted
    x_adj
  }

  pos2coord <- function(sector, x, y, track, ylim_top) {
    xlim <- circlize::get.cell.meta.data("xlim", sector.index = sector, track.index = track)
    if (!is.numeric(xlim) ||
        length(xlim) < 2 ||
        !all(is.finite(xlim)) || xlim[2] <= xlim[1])
      return(c(NA_real_, NA_real_))
    start_deg <- circlize::get.cell.meta.data("cell.start.degree",
                                              sector.index = sector,
                                              track.index = track)
    end_deg <- circlize::get.cell.meta.data("cell.end.degree",
                                            sector.index = sector,
                                            track.index = track)
    if (!is.finite(start_deg) || !is.finite(end_deg))
      return(c(NA_real_, NA_real_))
    bot <- circlize::get.cell.meta.data("cell.bottom.radius",
                                        sector.index = sector,
                                        track.index = track)
    top <- circlize::get.cell.meta.data("cell.top.radius",
                                        sector.index = sector,
                                        track.index = track)
    if (!is.finite(bot) || !is.finite(top))
      return(c(NA_real_, NA_real_))
    frac_x <- (x - xlim[1]) / (xlim[2] - xlim[1])
    frac_x <- max(0, min(1, frac_x))
    theta <- start_deg + frac_x * (end_deg - start_deg)
    frac_y <- y / ylim_top
    frac_y <- max(0, min(1, frac_y))
    r <- bot + frac_y * (top - bot)
    th <- theta * pi / 180
    c(r * cos(th), r * sin(th))
  }

  draw_gene_labels <- function(df_top,
                               is_down,
                               track_scatter,
                               track_anno,
                               label_y) {
    if (nrow(df_top) == 0 ||
        !is.finite(track_scatter) || !is.finite(track_anno))
      return()
    sectors_unique <- unique(as.character(df_top$chr))
    for (sector in sectors_unique) {
      if (is.na(sector) ||
          sector == "" || identical(sector, scale_sector))
        next
      sub <- df_top[as.character(df_top$chr) == sector, , drop = FALSE]
      if (nrow(sub) == 0)
        next
      sector_len <- length_map$end[match(sector, length_map$chr)]
      if (!is.finite(sector_len) || sector_len <= 0)
        next
      x_point <- if ("pos" %in% colnames(sub))
        sub$pos
      else
        (sub$start + sub$end) / 2
      x_point <- suppressWarnings(as.numeric(x_point))
      x_point <- pmin(pmax(x_point, 0), sector_len)
      min_gap <- if (is.null(NULL))
        max(1e6, sector_len * 0.02)
      else
        NULL
      ord <- order(x_point, na.last = TRUE)
      x_sorted <- x_point[ord]
      dx_step <- if (!is.null(NULL))
        NULL
      else
        max(5e5, min(min_gap * 0.3, sector_len * 0.005))
      dx_max <- if (!is.null(NULL))
        NULL
      else
        max(dx_step, min_gap * 1.2, sector_len * 0.03)
      x_label <- x_point
      y_label <- rep(label_y, length(x_point))
      if (sum(is.finite(x_sorted)) > 1) {
        groups <- integer(length(x_sorted))
        g <- 1L
        groups[1] <- g
        for (k in 2:length(x_sorted)) {
          if (!is.finite(x_sorted[k]) || !is.finite(x_sorted[k - 1])) {
            g <- g + 1L
          } else if ((x_sorted[k] - x_sorted[k - 1]) >= min_gap) {
            g <- g + 1L
          }
          groups[k] <- g
        }
        for (gg in unique(groups)) {
          idx_sorted <- which(groups == gg & is.finite(x_sorted))
          if (length(idx_sorted) == 0)
            next
          m <- length(idx_sorted)
          center <- (m + 1) / 2
          for (jj in seq_len(m)) {
            ii <- ord[idx_sorted[jj]]
            dx <- (jj - center) * dx_step
            dx <- max(-dx_max, min(dx_max, dx))
            x_label[ii] <- pmin(pmax(x_point[ii] + dx, 0), sector_len)
            lev <- (jj - 1) %% 1
            if (is_down) {
              y_label[ii] <- max(0, label_y - lev * 0.08)
            } else {
              y_label[ii] <- min(1, label_y + lev * 0.08)
            }
          }
        }
      }
      for (i in seq_len(nrow(sub))) {
        x_p <- x_point[i]
        x_l <- x_label[i]
        if (!is.finite(x_p) || !is.finite(x_l))
          next
        y_p <- if (is_down)
          (ymax_all - sub$abs_y[i])
        else
          sub$abs_y[i]
        if (!is.finite(y_p))
          next
        p <- pos2coord(
          sector = sector,
          x = x_p,
          y = y_p,
          track = track_scatter,
          ylim_top = ymax_all
        )
        x_e <- x_p + (x_l - x_p) * connector_elbow
        x_e <- pmin(pmax(x_e, 0), sector_len)
        e <- pos2coord(
          sector = sector,
          x = x_e,
          y = label_y,
          track = track_anno,
          ylim_top = 1
        )
        l <- pos2coord(
          sector = sector,
          x = x_l,
          y = y_label[i],
          track = track_anno,
          ylim_top = 1
        )
        if (length(p) < 2 ||
            length(l) < 2 ||
            any(!is.finite(p)) || any(!is.finite(l)))
          next
        if (length(e) < 2 || any(!is.finite(e)))
          next
        graphics::segments(p[1], p[2], e[1], e[2], col = connector_col, lwd = connector_lwd)
        graphics::segments(e[1], e[2], l[1], l[2], col = connector_col, lwd = connector_lwd)
        theta_deg <- (atan2(l[2], l[1]) * 180 / pi) %% 360
        base_srt <- if (is_down)
          (theta_deg + 180 + label_rotate) %% 360
        else
          (theta_deg + label_rotate) %% 360
        flip <- base_srt > 90 && base_srt < 270
        srt <- if (flip)
          (base_srt + 180) %% 360
        else
          base_srt
        adj_x <- if (flip)
          1
        else
          0
        graphics::text(
          l[1],
          l[2],
          labels = sub$label[i],
          cex = label_cex,
          font = label_font,
          srt = srt,
          adj = c(adj_x, 0.5),
          col = connector_col
        )
      }
    }
  }

  if (nrow(df_top_up) > 0)
    draw_gene_labels(
      df_top_up,
      is_down = FALSE,
      track_scatter = track_up_scatter,
      track_anno = track_up_anno,
      label_y = label_y_up
    )
  if (nrow(df_top_down) > 0)
    draw_gene_labels(
      df_top_down,
      is_down = TRUE,
      track_scatter = track_down_scatter,
      track_anno = track_down_anno,
      label_y = label_y_down
    )

  invisible(list(
    dmr = dmr,
    chrom = length_map,
    top_up = df_top_up,
    top_down = df_top_down
  ))
}
