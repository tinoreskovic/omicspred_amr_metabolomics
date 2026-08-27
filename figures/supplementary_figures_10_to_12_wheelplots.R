#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(circlize)
})

results_dir <- Sys.getenv("OMICSPRED_RESULTS_DIR")
out_dir <- Sys.getenv("OMICSPRED_FIGURE_DIR")
if (!nzchar(results_dir) || !nzchar(out_dir)) stop("Set OMICSPRED_RESULTS_DIR and OMICSPRED_FIGURE_DIR")
input_csv <- file.path(results_dir, "aou_disease_associations.csv")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

dat <- read.csv(input_csv, check.names = FALSE, stringsAsFactors = FALSE)

subclasses <- c(
  "chylomicrons and extremely large VLDL",
  "very large VLDL",
  "large VLDL",
  "medium VLDL",
  "small VLDL",
  "very small VLDL",
  "IDL",
  "large LDL",
  "medium LDL",
  "small LDL",
  "very large HDL",
  "large HDL",
  "medium HDL",
  "small HDL"
)

particle_names <- c(
  "Concentration of chylomicrons and extremely large VLDL particles",
  "Concentration of very large VLDL particles",
  "Concentration of large VLDL particles",
  "Concentration of medium VLDL particles",
  "Concentration of small VLDL particles",
  "Concentration of very small VLDL particles",
  "Concentration of IDL particles",
  "Concentration of large LDL particles",
  "Concentration of medium LDL particles",
  "Concentration of small LDL particles",
  "Concentration of very large HDL particles",
  "Concentration of large HDL particles",
  "Concentration of medium HDL particles",
  "Concentration of small HDL particles"
)

component_names <- function(component) {
  paste(component, "in", subclasses)
}

subclass_labels <- c("XXL", "XL", "L", "M", "S", "XS", "IDL", "L", "M", "S", "XL", "L", "M", "S")

cholesterol_names <- c(
  "Total cholesterol",
  "VLDL cholesterol",
  component_names("Cholesterol")[1:7],
  "LDL cholesterol",
  component_names("Cholesterol")[8:10],
  "HDL cholesterol",
  component_names("Cholesterol")[11:14]
)

cholesterol_labels <- c(
  "Tot", "Tot", subclass_labels[1:7],
  "Tot", subclass_labels[8:10],
  "Tot", subclass_labels[11:14]
)

free_cholesterol_names <- c(
  "Total free cholesterol",
  component_names("Free cholesterol")
)

free_cholesterol_labels <- c("Tot", subclass_labels)

cholesteryl_ester_names <- c(
  "Total esterified cholesterol",
  component_names("Cholesteryl esters")
)

cholesteryl_ester_labels <- c("Tot", subclass_labels)

triglyceride_names <- c(
  "Total triglycerides",
  "Triglycerides in VLDL",
  component_names("Triglycerides")[1:7],
  "Triglycerides in LDL",
  component_names("Triglycerides")[8:10],
  "Triglycerides in HDL",
  component_names("Triglycerides")[11:14]
)

triglyceride_labels <- c(
  "Tot", "Tot", subclass_labels[1:7],
  "Tot", subclass_labels[8:10],
  "Tot", subclass_labels[11:14]
)

trait_order <- c(
  particle_names,
  cholesterol_names,
  free_cholesterol_names,
  cholesteryl_ester_names,
  triglyceride_names,
  component_names("Phospholipids"),
  component_names("Total lipids"),
  "Average diameter for VLDL particles",
  "Average diameter for LDL particles",
  "Average diameter for HDL particles",
  "Apolipoprotein A1",
  "Apolipoprotein B",
  "Polyunsaturated fatty acids",
  "Monounsaturated fatty acids",
  "Saturated fatty acids",
  "Docosahexaenoic acid",
  "Linoleic acid",
  "Omega-3 fatty acids",
  "Omega-6 fatty acids",
  "Total fatty acids",
  "Degree of unsaturation",
  "Total cholines",
  "Phosphatidylcholines",
  "Sphingomyelins",
  "Phosphoglycerides",
  "Glucose",
  "Alanine",
  "Glutamine",
  "Glycine",
  "Histidine",
  "Isoleucine",
  "Leucine",
  "Valine",
  "Phenylalanine",
  "Tyrosine",
  "Acetoacetate",
  "Creatinine",
  "Glycoprotein acetyls"
)

if (length(trait_order) != 139L) {
  stop("Trait order has ", length(trait_order), " entries, expected 139.")
}

missing_traits <- setdiff(trait_order, unique(dat$Biomarker.Name))
extra_traits <- setdiff(unique(dat$Biomarker.Name), trait_order)
if (length(missing_traits) > 0L || length(extra_traits) > 0L) {
  stop("Trait mapping does not match the AoU source data: ", length(missing_traits), " missing and ", length(extra_traits), " extra traits")
}

labs_extra <- c(
  "VLDL-D", "LDL-D", "HDL-D", "Apo-AI", "Apo-B",
  "PUFA", "MUFA", "SFA", "DHA", "LA", "Omega-3", "Omega-6", "TotFA", "Unsat.",
  "TotCho", "PC", "SM", "PG",
  "Glc",
  "Ala", "Gln", "Gly", "His", "Ile", "Leu", "Val", "Phe", "Tyr",
  "AcAce", "Crea", "Glyc-A"
)

trait_key <- data.frame(
  id_name_s = seq_along(trait_order),
  Biomarker.Name = trait_order,
  label = c(
    subclass_labels,
    cholesterol_labels,
    free_cholesterol_labels,
    cholesteryl_ester_labels,
    triglyceride_labels,
    subclass_labels,
    subclass_labels,
    labs_extra
  ),
  stringsAsFactors = FALSE
)

blocks <- data.frame(
  start = c(1, 15, 33, 48, 63, 81, 95, 109, 114, 123, 127, 137),
  end = c(14, 32, 47, 62, 80, 94, 108, 113, 122, 126, 136, 139),
  label = c(
    "Lipoprotein particles",
    "Cholesterol",
    "Free cholesterol",
    "Cholesteryl esters",
    "Triglycerides",
    "Phospholipids",
    "Total lipids",
    "lipoprotein particle sizes",
    "Fatty acids",
    "Other lipids",
    "Glycolysis & amino acids",
    "Ketone bodies, fluid balance,"
  ),
  label2 = c(
    rep(NA_character_, 7),
    "Apolipoproteins and",
    rep(NA_character_, 3),
    "and inflammation"
  ),
  label_x_shift = c(rep(0, 7), 2.0, rep(0, 3), -1.0),
  label2_x_shift = c(rep(0, 7), 0, rep(0, 3), 0),
  label_adj = c(rep(0.5, 7), 1, rep(0.5, 3), 0),
  stringsAsFactors = FALSE
)

make_subgroups <- function(block, labels, starts, ends) {
  fill_map <- c(
    All = "#5F5F5F80",
    VLDL = "#BDBDBD66",
    IDL = "#6E6E6E73",
    LDL = "#D9D9D966",
    HDL = "#EFEFEF99"
  )
  data.frame(
    block = block,
    start = starts,
    end = ends,
    label = labels,
    fill = unname(fill_map[labels]),
    stringsAsFactors = FALSE
  )
}

subgroup_segments <- do.call(rbind, list(
  make_subgroups("Lipoprotein particles", c("VLDL", "IDL", "LDL", "HDL"), c(1, 7, 8, 11), c(6, 7, 10, 14)),
  make_subgroups("Cholesterol", c("All", "VLDL", "IDL", "LDL", "HDL"), c(15, 16, 23, 24, 28), c(15, 22, 23, 27, 32)),
  make_subgroups("Free cholesterol", c("All", "VLDL", "IDL", "LDL", "HDL"), c(33, 34, 40, 41, 44), c(33, 39, 40, 43, 47)),
  make_subgroups("Cholesteryl esters", c("All", "VLDL", "IDL", "LDL", "HDL"), c(48, 49, 55, 56, 59), c(48, 54, 55, 58, 62)),
  make_subgroups("Triglycerides", c("All", "VLDL", "IDL", "LDL", "HDL"), c(63, 64, 71, 72, 76), c(63, 70, 71, 75, 80)),
  make_subgroups("Phospholipids", c("VLDL", "IDL", "LDL", "HDL"), c(81, 87, 88, 91), c(86, 87, 90, 94)),
  make_subgroups("Total lipids", c("VLDL", "IDL", "LDL", "HDL"), c(95, 101, 102, 105), c(100, 101, 104, 108))
))
subgroup_segments$x <- (subgroup_segments$start + subgroup_segments$end) / 2

subgroup_radial_boundaries <- sort(unique(c(
  subgroup_segments$start - 0.5,
  subgroup_segments$end + 0.5
)))

disease_specs <- data.frame(
  Disease_Name = c("Type 2 Diabetes", "Ischemic Heart Disease", "Chronic Kidney Disease"),
  disease_label = c("Type 2 diabetes", "Ischaemic heart disease", "Chronic kidney disease"),
  disease_file = c("t2d", "ihd", "ckd"),
  stringsAsFactors = FALSE
)

score_specs <- data.frame(
  score_set = c("MCPS-trained", "INTERVAL-trained"),
  score_label = c("MCPS-trained", "INTERVAL-trained"),
  score_file = c("mcps", "interval"),
  color = c("#0072B2", "#D55E00"),
  x_offset = c(-0.19, 0.19),
  stringsAsFactors = FALSE
)

axis_specs <- list(
  "Type 2 Diabetes" = list(
    ylim_or = c(0.85, 1.40),
    ycuts = c(0.85, 0.90, 1.00, 1.15, 1.30, 1.40)
  ),
  "Ischemic Heart Disease" = list(
    ylim_or = c(0.90, 1.15),
    ycuts = c(0.90, 0.95, 1.00, 1.05, 1.10, 1.15)
  ),
  "Chronic Kidney Disease" = list(
    ylim_or = c(0.87, 1.22),
    ycuts = c(0.90, 0.95, 1.00, 1.05, 1.10, 1.20)
  )
)

format_axis_label <- function(x) {
  sub("\\.?0+$", "", sprintf("%.2f", x))
}

prepare_pair_data <- function(disease_name) {
  sub <- dat[dat$Disease_Name == disease_name & dat$score_set %in% score_specs$score_set, ]
  sub <- merge(sub, trait_key, by = "Biomarker.Name", all.x = TRUE, sort = FALSE)
  sub <- merge(sub, score_specs, by = "score_set", all.x = TRUE, sort = FALSE)
  sub <- sub[order(sub$id_name_s, sub$score_file), ]
  if (any(is.na(sub$id_name_s)) || any(is.na(sub$color))) {
    stop("Missing trait/model mapping after merge for ", disease_name)
  }
  sub$x_plot <- sub$id_name_s + sub$x_offset
  sub$log_or <- log(sub$estimate)
  sub$log_l95 <- log(sub$L95)
  sub$log_u95 <- log(sub$U95)
  sub$sig <- sub$fdr < 0.05
  sub
}

add_center_legend <- function(title_text) {
  text(0, 0.13, title_text, cex = 0.75, font = 2)

  cex_legend <- 0.56
  point_cex <- 0.9
  x_left <- -0.34
  x_right <- 0.08
  y1 <- -0.10
  y2 <- -0.18

  points(x_left, y1, pch = 21, col = score_specs$color[score_specs$score_set == "MCPS-trained"],
         bg = score_specs$color[score_specs$score_set == "MCPS-trained"], cex = point_cex, lwd = 1.1)
  text(x_left + 0.045, y1, "MCPS-trained", adj = 0, cex = cex_legend)

  points(x_left, y2, pch = 21, col = score_specs$color[score_specs$score_set == "INTERVAL-trained"],
         bg = score_specs$color[score_specs$score_set == "INTERVAL-trained"], cex = point_cex, lwd = 1.1)
  text(x_left + 0.045, y2, "INTERVAL-trained", adj = 0, cex = cex_legend)

  points(x_right, y1, pch = 21, col = "gray20", bg = "gray20", cex = point_cex, lwd = 1.1)
  text(x_right + 0.045, y1, expression(FDR~italic(P)~"<"~0.05), adj = 0, cex = cex_legend)

  points(x_right, y2, pch = 21, col = "gray20", bg = "white", cex = point_cex, lwd = 1.1)
  text(x_right + 0.045, y2, expression(FDR~italic(P) >= 0.05), adj = 0, cex = cex_legend)
}

draw_lipoprotein_subgroup_ring <- function(fact, ylim, yspan) {
  ring_inner <- max(ylim) + yspan * 0.060
  ring_outer <- max(ylim) + yspan * 0.125
  divider_inner <- ring_inner
  divider_outer <- ring_outer

  for (i in seq_len(nrow(subgroup_segments))) {
    pos <- circlize(
      c(subgroup_segments$start[i] - 0.5, subgroup_segments$end[i] + 0.5),
      c(ring_inner, ring_outer),
      sector.index = fact,
      track.index = 1
    )
    draw.sector(
      pos[1, "theta"], pos[2, "theta"], pos[1, "rou"], pos[2, "rou"],
      clock.wise = TRUE,
      col = subgroup_segments$fill[i],
      border = "#FFFFFFB3"
    )
  }

  for (x_boundary in unique(c(subgroup_segments$start - 0.5, subgroup_segments$end + 0.5))) {
    circos.segments(
      x0 = x_boundary,
      y0 = divider_inner,
      x1 = x_boundary,
      y1 = divider_outer,
      col = "#55555580",
      lwd = 0.75
    )
  }

  for (x_boundary in subgroup_radial_boundaries) {
    circos.segments(
      x0 = x_boundary,
      y0 = min(ylim),
      x1 = x_boundary,
      y1 = max(ylim),
      col = "#66666680",
      lwd = 0.70
    )
  }
}

draw_paired_wheel <- function(plot_dat, title_text, output_file, axis_spec, show_subgroup_ring = TRUE) {
  pdf(output_file, width = 11.5, height = 11.5, useDingbats = FALSE, bg = "white")
  on.exit({
    circos.clear()
    dev.off()
  }, add = TRUE)

  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par), add = TRUE)
  par(xpd = NA, mar = rep(0, 4), oma = rep(2.1, 4), family = "serif", bg = "white")

  fact <- "mets"
  len_data <- length(trait_order)
  xlim <- c(1, len_data)
  ylim_or <- axis_spec$ylim_or
  ylim <- log(ylim_or)
  ycuts <- log(axis_spec$ycuts)
  ycut_labs <- format_axis_label(exp(ycuts))
  yspan <- diff(ylim)
  ylab_trait <- max(ylim) + yspan * 0.27
  ylab_subgroup <- max(ylim) + yspan * 0.39
  ylab_group <- max(ylim) + yspan * 0.58
  ylab_group_inner <- max(ylim) + yspan * 0.49

  adj <- 0.004
  clamp <- function(x) pmin(pmax(x, min(ylim) + abs(min(ylim)) * adj), max(ylim) - abs(max(ylim)) * adj)
  plot_dat$log_or_plot <- clamp(plot_dat$log_or)
  plot_dat$log_l95_plot <- clamp(plot_dat$log_l95)
  plot_dat$log_u95_plot <- clamp(plot_dat$log_u95)

  circos.clear()
  circos.par(
    track.height = 0.58,
    cell.padding = c(0, 0, 0, 0),
    gap.degree = 50,
    start.degree = 90,
    unit.circle.segments = 5000,
    canvas.xlim = c(-1.25, 1.25),
    canvas.ylim = c(-1.25, 1.25),
    points.overflow.warning = FALSE
  )
  circos.initialize(factors = fact, xlim = c(0, len_data))
  circos.track(factors = fact, ylim = ylim, bg.border = NA)

  for (i in seq_len(nrow(blocks))) {
    if (i %% 2 == 1) {
      pos <- circlize(
        c(blocks$start[i] - 0.5, blocks$end[i] + 0.5),
        c(min(ylim), max(ylim)),
        sector.index = fact,
        track.index = 1
      )
      draw.sector(pos[1, "theta"], pos[2, "theta"], pos[1, "rou"], pos[2, "rou"],
                  clock.wise = TRUE, col = "#D9D9D966", border = NA)
    }
  }

  circos.track(track.index = 1, bg.border = "white", factors = fact, ylim = ylim, panel.fun = function(x, y) {
    circos.segments(x0 = min(xlim) - 0.5, y0 = max(ylim), x1 = min(xlim) - 0.5, y1 = min(ylim), col = "black", lwd = 1.2)
    circos.segments(x0 = max(xlim) + 0.5, y0 = max(ylim), x1 = max(xlim) + 0.5, y1 = min(ylim), col = "black", lwd = 1.2)
    circos.segments(x0 = min(xlim) - 0.75, y0 = max(ylim), x1 = max(xlim) + 0.5, y1 = max(ylim), col = "black", lwd = 1.2)
    circos.segments(x0 = min(xlim) - 0.75, y0 = min(ylim), x1 = max(xlim) + 0.5, y1 = min(ylim), col = "black", lwd = 1.2)

    for (yc in ycuts) {
      circos.segments(
        x0 = min(xlim) - 0.75,
        y0 = yc,
        x1 = max(xlim) + 0.5,
        y1 = yc,
        col = ifelse(abs(yc) < 1e-8, "black", "gray75"),
        lwd = ifelse(abs(yc) < 1e-8, 1.3, 0.8)
      )
    }

    if (show_subgroup_ring) {
      draw_lipoprotein_subgroup_ring(fact, ylim, yspan)
    }

    for (score_name in score_specs$score_set) {
      z <- plot_dat[plot_dat$score_set == score_name, ]
      ci_col <- adjustcolor(z$color, alpha.f = 0.78)
      circos.segments(z$x_plot, z$log_l95_plot, z$x_plot, z$log_u95_plot, col = ci_col, lwd = 0.85)
      circos.segments(z$x_plot - 0.055, z$log_l95_plot, z$x_plot + 0.055, z$log_l95_plot, col = ci_col, lwd = 0.85)
      circos.segments(z$x_plot - 0.055, z$log_u95_plot, z$x_plot + 0.055, z$log_u95_plot, col = ci_col, lwd = 0.85)
    }

    circos.points(
      plot_dat$x_plot,
      plot_dat$log_or_plot,
      pch = 21,
      col = plot_dat$color,
      bg = ifelse(plot_dat$sig, plot_dat$color, "white"),
      cex = 0.50
    )

    circos.yaxis(at = ycuts, labels = ycut_labs, labels.cex = 0.58, tick = FALSE, col = "white")

    for (i in seq_len(nrow(blocks))) {
      mid <- (blocks$start[i] + blocks$end[i]) / 2 + blocks$label_x_shift[i]
      circos.text(mid, ylab_group, labels = blocks$label[i], facing = "bending.inside",
                  niceFacing = TRUE, cex = 0.66, adj = c(blocks$label_adj[i], 0.5), font = 2)
      if (!is.na(blocks$label2[i])) {
        circos.text(mid + blocks$label2_x_shift[i], ylab_group_inner, labels = blocks$label2[i], facing = "bending.inside",
                    niceFacing = TRUE, cex = 0.66, adj = c(blocks$label_adj[i], 0.5), font = 2)
      }
    }

    circos.text(subgroup_segments$x, ylab_subgroup, labels = subgroup_segments$label,
                facing = "bending.inside", niceFacing = TRUE, cex = 0.52, adj = c(0.5, 0.5), font = 2)

    circos.text(trait_key$id_name_s + 0.25, ylab_trait, labels = trait_key$label,
                facing = "clockwise", niceFacing = TRUE, cex = 0.46, adj = 0, font = 1)
  })

  add_center_legend(title_text)
}

for (i in seq_len(nrow(disease_specs))) {
  plot_dat <- prepare_pair_data(disease_specs$Disease_Name[i])
  title_text <- paste(
    paste0(disease_specs$disease_label[i], " and"),
    "metabolomic traits predicted by",
    "MCPS- and INTERVAL-trained models",
    sep = "\n"
  )
  draw_paired_wheel(
    plot_dat,
    title_text,
    file.path(out_dir, paste0("supplementary_figure_", i + 9L, ".pdf")),
    axis_specs[[disease_specs$Disease_Name[i]]],
    show_subgroup_ring = TRUE
  )
}

message("Wrote paired model wheel plots to: ", normalizePath(out_dir))
