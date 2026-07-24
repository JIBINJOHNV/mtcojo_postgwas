#!/usr/bin/env Rscript
# plot_gwas.R - Multi-trait Manhattan & Q-Q Plotting via rMVP
# Called from Python wrapper in report_generator.py

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Usage: Rscript plot_gwas.R <out_dir> <merged_gwas_summary_tsv> <rg_csv>")
}

out_dir      <- args[1]
merged_tsv   <- args[2]
rg_csv       <- args[3]

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

cat("[R] Loading merged GWAS summary data...\n")
if (!file.exists(merged_tsv)) {
  stop(paste("Merged summary TSV not found:", merged_tsv))
}

if (!requireNamespace("data.table", quietly = TRUE)) {
  stop("[R] data.table is required for fast large-file plotting. Install with: conda install -n mtcojo_postgwas r-data.table")
}
merged_df <- data.table::fread(
  merged_tsv,
  sep = "\t",
  header = TRUE,
  data.table = FALSE,
  showProgress = FALSE
)

# Ensure CHR, BP, and SNP are correct types
merged_df$CHR <- as.numeric(merged_df$CHR)
merged_df$BP  <- as.numeric(merged_df$BP)
merged_df$SNP <- as.character(merged_df$SNP)

# Filter out missing positions
merged_df <- merged_df[!is.na(merged_df$CHR) & !is.na(merged_df$BP), ]

p_cols <- colnames(merged_df)[!(colnames(merged_df) %in% c("SNP", "CHR", "BP"))]

cat(paste("[R] Found", length(p_cols), "trait tracks for plotting:", paste(p_cols, collapse=", "), "\n"))

# Replace NAs with 1.0 and clamp plot p-values to the valid range.
for (col in p_cols) {
  merged_df[is.na(merged_df[[col]]), col] <- 1.0
  p <- suppressWarnings(as.numeric(merged_df[[col]]))
  p[is.finite(p) & p < 1e-300] <- 1e-300
  p[is.finite(p) & p > 1.0] <- 1.0
  merged_df[[col]] <- p
}

plot_df  <- merged_df
n_traits <- length(p_cols)

# ── Unified publication palette ──────────────────────────────────────────────
# Manhattan chromosome bands: neutral charcoal / cool grey (reads cleanly in print & greyscale)
manhattan_band_cols <- c("#242C37", "#98A2B3")
# Trait colors: Okabe-Ito colorblind-safe set, extended with ColorBrewer Dark2 purple
# instead of black, since black collides visually with the charcoal Manhattan band
trait_cols <- c("#0072B2", "#D55E00", "#E69F00", "#009E73", "#CC79A7", "#56B4E9", "#7570B3")
# Single crimson accent used consistently for genome-wide significance across all plots
sig_col       <- "#B2182B"
suggestive_col <- "#6B7280"
# Diverging palette for the rg heatmap (ColorBrewer RdBu), sharing the crimson accent hue
heatmap_low  <- "#2166AC"
heatmap_high <- "#B2182B"

keep_for_manhattan <- rep(FALSE, nrow(plot_df))
if (length(p_cols) > 0) {
  for (col in p_cols) {
    p <- suppressWarnings(as.numeric(plot_df[[col]]))
    keep_for_manhattan <- keep_for_manhattan | (is.finite(p) & p < 0.1)
  }
}
keep_for_manhattan <- keep_for_manhattan | ((seq_len(nrow(plot_df)) %% 20) == 1)
plot_df_manhattan <- plot_df[keep_for_manhattan, ]
cat(paste("[R] Manhattan plotting retained", nrow(plot_df_manhattan), "of", nrow(plot_df), "variants after thinning 95% of P >= 0.1 background points.\n"))

choose_plot_file <- function(files) {
  files <- files[file.exists(files)]
  if (length(files) == 0) {
    return(NA_character_)
  }
  info <- file.info(files)
  files[order(info$mtime, info$size, decreasing = TRUE, na.last = TRUE)[1]]
}

make_overlaid_manhattan_with_legend <- function(df, p_cols, out_file) {
  plot_df <- df[!is.na(df$CHR) & !is.na(df$BP), ]
  plot_df <- plot_df[order(plot_df$CHR, plot_df$BP), ]
  if (nrow(plot_df) == 0 || length(p_cols) == 0) {
    return(FALSE)
  }

  chr_levels <- sort(unique(plot_df$CHR))
  offset <- 0
  axis_at <- numeric(0)
  axis_labels <- character(0)
  plot_df$plot_pos <- NA_real_

  for (chr in chr_levels) {
    idx <- which(plot_df$CHR == chr)
    chr_bp <- plot_df$BP[idx]
    plot_df$plot_pos[idx] <- chr_bp + offset
    axis_at <- c(axis_at, mean(range(plot_df$plot_pos[idx], na.rm = TRUE)))
    axis_labels <- c(axis_labels, as.character(chr))
    offset <- max(plot_df$plot_pos[idx], na.rm = TRUE) + 1e6
  }

  max_y <- 0
  for (col in p_cols) {
    p <- suppressWarnings(as.numeric(plot_df[[col]]))
    p <- p[is.finite(p) & p > 0]
    if (length(p) > 0) {
      max_y <- max(max_y, -log10(min(p)))
    }
  }
  y_lim <- c(0, max(5.5, min(max_y * 1.15 + 0.5, 30)))

  jpeg(out_file, width = 15, height = 6.2, units = "in", res = 300, quality = 95)
  par(mar = c(5.2, 5.2, 3.8, 8.8), xpd = NA)
  plot(
    NA,
    xlim = range(plot_df$plot_pos, na.rm = TRUE),
    ylim = y_lim,
    xaxt = "n",
    xlab = "Chromosome",
    ylab = expression(-log[10](italic(P))),
    main = "Multi-Trait Overlaid Manhattan Plot",
    bty = "l"
  )
  grid(nx = NA, ny = NULL, col = "#D0D5DD", lty = "dotted")
  abline(h = -log10(1e-5), col = suggestive_col, lty = "dotted", lwd = 1.2)
  abline(h = -log10(5e-8), col = sig_col, lty = "dashed", lwd = 1.2)
  axis(1, at = axis_at, labels = axis_labels, tick = FALSE, cex.axis = 0.85)

  for (i in seq_along(p_cols)) {
    col_name <- p_cols[i]
    p <- suppressWarnings(as.numeric(plot_df[[col_name]]))
    keep <- is.finite(p) & p > 0 & is.finite(plot_df$plot_pos)
    if (any(keep)) {
      points(
        plot_df$plot_pos[keep],
        -log10(p[keep]),
        pch = 20,
        cex = 0.28,
        col = adjustcolor(trait_cols[(i - 1) %% length(trait_cols) + 1], alpha.f = 0.34)
      )
    }
  }

  legend(
    "topright",
    inset = c(-0.17, 0),
    legend = gsub("_", " ", p_cols),
    col = trait_cols[seq_along(p_cols)],
    pch = 20,
    pt.cex = 1.2,
    cex = 0.9,
    bty = "n",
    title = "Phenotype"
  )
  dev.off()
  TRUE
}

# ─────────────────────────────────────────────────────────────────────────────
# rMVP Multi-Track Manhattan & 3-Track Q-Q Plots
# ─────────────────────────────────────────────────────────────────────────────
has_rmvp   <- requireNamespace("rMVP", quietly = TRUE)
has_cmplot <- requireNamespace("CMplot", quietly = TRUE)

if (has_rmvp) {
  cat("[R] rMVP package found. Generating plots with rMVP::MVP.Report...\n")
  setwd(out_dir)
  
  # 1. Multi-Track Manhattan Plot (crisp points with vivid chromosome bands)
  tryCatch({
    rMVP::MVP.Report(
      plot_df_manhattan,
      plot.type = "m",
      col = manhattan_band_cols,
      multracks = TRUE,
      threshold = c(1e-6, 1e-4),
      threshold.lty = c(1, 2),
      threshold.lwd = c(1, 1),
      threshold.col = c(suggestive_col, sig_col),
      amplify = TRUE,
      bin.size = 1e6,
      chr.den.col = c("#56B4E9", "#E69F00", sig_col),
      signal.col = c("#0072B2", "#D55E00"),
      signal.cex = c(0.3, 0.3),
      width = 16,
      height = max(7.0, 3.4 * n_traits),
      file.type = "jpg",
      memo = "mt",
      dpi = 300
    )
  }, error = function(e) {
    cat(paste("[R] Warning in rMVP Manhattan plot:", e$message, "\n"))
  })
  
  # Stacked multracks Manhattan plot (pick file with longest filename containing all traits)
  man_tracks <- list.files(pattern = ".*Multracks.*Manhattan.*\\.jpg$")
  man_tracks <- man_tracks[man_tracks != "Rectangular-Manhattan.multi_trait.jpg"]
  if (length(man_tracks) > 0) {
    best_file <- man_tracks[order(nchar(man_tracks), decreasing = TRUE)[1]]
    file.copy(best_file, "Rectangular-Manhattan.multi_trait.jpg", overwrite = TRUE)
    cat(paste("[R] Saved Stacked Manhattan plot → Rectangular-Manhattan.multi_trait.jpg from", best_file, "\n"))
  }
  
  # Overlaid multraits Manhattan plot
  man_overlaid <- list.files(pattern = ".*Multraits.*Manhattan.*\\.jpg$")
  man_overlaid <- man_overlaid[man_overlaid != "Overlaid-Manhattan.multi_trait.jpg"]
  if (length(man_overlaid) > 0) {
    best_file <- man_overlaid[order(nchar(man_overlaid), decreasing = TRUE)[1]]
    file.copy(best_file, "Overlaid-Manhattan.multi_trait.jpg", overwrite = TRUE)
    cat(paste("[R] Saved Overlaid Manhattan plot → Overlaid-Manhattan.multi_trait.jpg from", best_file, "\n"))
  }
  
  # 2. 3-Track Separate Q-Q Plot (multracks = TRUE for 3 distinct panel tracks)
  tryCatch({
    rMVP::MVP.Report(
      plot_df,
      plot.type = "q",
      col = trait_cols,
      threshold = 1e6,
      signal.pch = 19,
      signal.cex = 0.25,
      signal.col = sig_col,
      conf.int = TRUE,
      box = FALSE,
      multracks = TRUE,
      file.type = "jpg",
      memo = "mt",
      dpi = 300
    )
  }, error = function(e) {
    cat(paste("[R] Warning in rMVP Q-Q plot (multracks=TRUE):", e$message, "\n"))
  })
  
  qq_multracks <- list.files(pattern = ".*Multracks.*QQ.*\\.jpg$")
  qq_multracks <- qq_multracks[qq_multracks != "QQ-Plot.multi_trait.jpg"]
  if (length(qq_multracks) > 0) {
    best_file <- qq_multracks[order(nchar(qq_multracks), decreasing = TRUE)[1]]
    file.copy(best_file, "QQ-Plot.multi_trait.jpg", overwrite = TRUE)
    cat(paste("[R] Saved 3-Panel Separate Q-Q plot → QQ-Plot.multi_trait.jpg from", best_file, "\n"))
  }

  # 3. Combined Single-Panel Overlaid Q-Q Plot (multracks = FALSE for all traits together)
  tryCatch({
    rMVP::MVP.Report(
      plot_df,
      plot.type = "q",
      col = trait_cols,
      threshold = 1e6,
      signal.pch = 19,
      signal.cex = 0.25,
      signal.col = sig_col,
      conf.int = TRUE,
      box = FALSE,
      multracks = FALSE,
      file.type = "jpg",
      memo = "combined",
      dpi = 300
    )
  }, error = function(e) {
    cat(paste("[R] Warning in rMVP Q-Q plot (multracks=FALSE):", e$message, "\n"))
  })

  qq_combined <- list.files(pattern = ".*Multraits.*QQ.*\\.jpg$|.*combined.*QQ.*\\.jpg$")
  qq_combined <- qq_combined[qq_combined != "QQ-Plot.combined.jpg" & qq_combined != "QQ-Plot.multi_trait.jpg"]
  if (length(qq_combined) > 0) {
    best_file <- qq_combined[order(nchar(qq_combined), decreasing = TRUE)[1]]
    file.copy(best_file, "QQ-Plot.combined.jpg", overwrite = TRUE)
    cat(paste("[R] Saved Combined Single-Panel Q-Q plot → QQ-Plot.combined.jpg from", best_file, "\n"))
  }
  
} else if (has_cmplot) {
  cat("[R] CMplot package found. Generating plots with CMplot...\n")
  setwd(out_dir)
  CMplot::CMplot(
    plot_df_manhattan,
    type = "p",
    plot.type = "m",
    LOG10 = TRUE,
    threshold = c(1e-5, 5e-8),
    threshold.col = c(suggestive_col, sig_col),
    col = manhattan_band_cols,
    cex = 0.3,
    file = "jpg",
    dpi = 300,
    file.output = TRUE,
    verbose = TRUE,
    multracks = TRUE,
    width = 12,
    height = max(7.0, 3.8 * n_traits)
  )
  CMplot::CMplot(
    plot_df,
    plot.type = "q",
    col = trait_cols,
    file = "jpg",
    dpi = 300,
    file.output = TRUE,
    verbose = TRUE,
    multracks = TRUE,
    width = 7.5,
    height = 7
  )
  man_files <- list.files(pattern = "^Multi-tracks_Manhtn.*\\.jpg$|^Rectangular-Manhattan.*\\.jpg$")
  man_files <- man_files[man_files != "Rectangular-Manhattan.multi_trait.jpg"]
  if (length(man_files) > 0) {
    best_file <- choose_plot_file(man_files)
    file.copy(best_file, "Rectangular-Manhattan.multi_trait.jpg", overwrite = TRUE)
    cat(paste("[R] Saved Stacked Manhattan plot → Rectangular-Manhattan.multi_trait.jpg from", best_file, "\n"))
  }
  qq_files <- list.files(pattern = "^Multi-tracks_QQplot.*\\.jpg$|^QQ-Plot.*\\.jpg$")
  qq_files <- qq_files[qq_files != "QQ-Plot.multi_trait.jpg"]
  if (length(qq_files) > 0) {
    best_file <- choose_plot_file(qq_files)
    file.copy(best_file, "QQ-Plot.multi_trait.jpg", overwrite = TRUE)
    cat(paste("[R] Saved 3-Panel Separate Q-Q plot → QQ-Plot.multi_trait.jpg from", best_file, "\n"))
  }
} else {
  cat("[R] Neither rMVP nor CMplot is installed.\n")
}

tryCatch({
  ok <- make_overlaid_manhattan_with_legend(plot_df_manhattan, p_cols, file.path(out_dir, "Overlaid-Manhattan.multi_trait.jpg"))
  if (ok) {
    cat("[R] Saved Overlaid Manhattan plot with phenotype legend → Overlaid-Manhattan.multi_trait.jpg\n")
  }
}, error = function(e) {
  cat(paste("[R] Warning in custom overlaid Manhattan legend plot:", e$message, "\n"))
})

# ─────────────────────────────────────────────────────────────────────────────
# 3. Genetic Correlation Matrix Heatmap
# ─────────────────────────────────────────────────────────────────────────────
if (file.exists(rg_csv) && file.info(rg_csv)$size > 10) {
  cat("[R] Generating genetic correlation heatmap...\n")
  rg <- read.csv(rg_csv, stringsAsFactors = FALSE)
  if (nrow(rg) > 0) {
    all_traits <- sort(unique(c(rg$p1, rg$p2)))
    n <- length(all_traits)
    mat <- matrix(1.0, nrow = n, ncol = n, dimnames = list(all_traits, all_traits))
    
    for (row_idx in seq_len(nrow(rg))) {
      t1 <- rg$p1[row_idx]
      t2 <- rg$p2[row_idx]
      mat[t1, t2] <- rg$rg[row_idx]
      mat[t2, t1] <- rg$rg[row_idx]
    }
    
    png(file.path(out_dir, "ldsc_rg_heatmap.png"), width = 600, height = 500, res = 120)
    par(mar = c(6, 6, 4, 4))
    image(1:n, 1:n, mat, col = colorRampPalette(c(heatmap_low, "white", heatmap_high))(100), 
          zlim = c(-1, 1), axes = FALSE, xlab = "", ylab = "")
    axis(1, at = 1:n, labels = all_traits, las = 2)
    axis(2, at = 1:n, labels = all_traits, las = 2)
    title("LDSC Genetic Correlation (rg)")
    
    for (i in 1:n) {
      for (j in 1:n) {
        text(j, i, sprintf("%.3f", mat[i, j]), cex = 0.8, font = 2)
      }
    }
    dev.off()
  }
}

cat("[R] Plotting script complete.\n")