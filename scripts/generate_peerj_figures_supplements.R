#!/usr/bin/env Rscript

# Optionally prepend a user-supplied package library for controlled runs.
requested_library <- Sys.getenv("R_LIBS_BENCHMARK", unset = "")
if (nzchar(requested_library) && dir.exists(requested_library)) {
  .libPaths(unique(c(requested_library, .libPaths())))
}

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
  library(scales)
  library(grid)
  library(jsonlite)
  library(digest)
  library(png)
})

options(stringsAsFactors = FALSE, scipen = 999)

stop_quietly <- function(message) {
  message(message)
  quit(save = "no", status = 1)
}

root_dir <- normalizePath(".", winslash = "/", mustWork = TRUE)
output_root <- file.path(root_dir, "artifacts", "submission", "PeerJ_2026-11-18", "03_Figures_and_Supplements")
figure_dir <- file.path(output_root, "Figures")
figure_data_dir <- file.path(output_root, "Figure_Data")
table_dir <- file.path(output_root, "Supplemental_Tables")
data_dir <- file.path(output_root, "Supplemental_Data")

dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_data_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)

read_project_csv <- function(relative_path) {
  path <- file.path(root_dir, relative_path)
  if (!file.exists(path)) stop_quietly(paste("Missing required input:", relative_path))
  read_csv(path, show_col_types = FALSE, progress = FALSE, na = c("", "NA", "N/A"))
}

write_project_csv <- function(data, relative_path) {
  path <- file.path(output_root, relative_path)
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  write_excel_csv(data, path, na = "")
  invisible(path)
}

cohort_palette <- c(GSE41664 = "#0072B2", GSE6710 = "#D55E00", GSE154200 = "#009E73")
status_palette <- c(supported = "#009E73", uncertain = "#E69F00", not_supported = "#7F7F7F")
point_palette <- c(
  "Primary comparison" = "#0072B2",
  "Primary comparison; partial coverage" = "#E69F00",
  "Excluded / descriptive only" = "#9E9E9E"
)
signature_order <- c(
  "RRDS", "RiboSisRG", "TRGS", "IMscore", "Ferroptosis8",
  "HKG34", "Epidermal10", "Apoptosis5", "m6A10", "PLOS14"
)

theme_peerj <- theme_classic(base_size = 10, base_family = "Arial") +
  theme(
    axis.title = element_text(color = "#222222", size = 10),
    axis.text = element_text(color = "#333333", size = 8.5),
    strip.background = element_rect(fill = "#F2F2F2", color = NA),
    strip.text = element_text(face = "bold", color = "#222222", size = 9),
    panel.grid.major.x = element_line(color = "#E6E6E6", linewidth = 0.35),
    panel.grid.minor = element_blank(),
    legend.title = element_text(size = 8.5, face = "bold"),
    legend.text = element_text(size = 8),
    legend.key = element_blank(),
    plot.tag = element_text(face = "bold", size = 13),
    plot.margin = margin(6, 10, 6, 6)
  )

# Inputs --------------------------------------------------------------------

cohort_endpoint <- read_project_csv("data/benchmark/cohort_endpoint_matrix.csv")
signature_registry <- read_project_csv("data/benchmark/signature_registry.csv")
platform_coverage <- read_project_csv("data/benchmark/platform_gene_coverage.csv")
gse154_audit <- fromJSON(file.path(root_dir, "data/benchmark/GSE154200_audit_summary.json"))
gse154_metadata <- read_project_csv("data/benchmark/score_stage_independent/GSE154200/D2_sample_metadata.csv")
gse154_gene_coverage <- read_project_csv("data/benchmark/score_stage_independent/GSE154200/gene_coverage.csv")
gse154_preprocessing <- read_project_csv("data/benchmark/score_stage_independent/GSE154200/preprocessing_summary.csv")
gse154_screening <- read_project_csv("data/benchmark/GSE154200_cohort_screening.csv")
gse154_pairs <- read_project_csv("data/benchmark/GSE154200_pair_audit.csv")

d2_41664 <- read_project_csv("data/benchmark/score_stage/GSE41664/D2_signature_metrics.csv") %>%
  mutate(cohort = "GSE41664", platform = "GPL570", .before = 1)
d2_6710 <- read_project_csv("data/benchmark/score_stage/GSE6710/D2_signature_metrics.csv") %>%
  mutate(cohort = "GSE6710", platform = "GPL96", .before = 1)
d2_154200 <- read_project_csv("data/benchmark/score_stage_independent/GSE154200/increment_stage/D2_signature_metrics.csv") %>%
  mutate(platform = "GPL20301", primary_D2_eligible = TRUE, .after = signature_id)

inc_41664 <- read_project_csv("data/benchmark/score_stage/increment_stage/increment_model_metrics.csv")
inc_154200 <- read_project_csv("data/benchmark/score_stage_independent/GSE154200/increment_stage/increment_model_metrics.csv")
inc_all <- bind_rows(inc_41664, inc_154200)
inc_main <- inc_all %>% filter(model_id == "signature_plus_baselines")
baseline_main <- inc_all %>% filter(model_id == "baseline_combined")

spearman_long <- read_project_csv("data/benchmark/score_stage/D2_primary_spearman_long_all.csv")
spearman_summary <- read_project_csv("data/benchmark/score_stage/D2_primary_spearman_pair_summary.csv")
effective_dimension <- read_project_csv("data/benchmark/score_stage/D2_primary_effective_dimension_all.csv")

# Numeric guards -------------------------------------------------------------

gse416 <- cohort_endpoint %>% filter(cohort == "GSE41664")
gse6710 <- cohort_endpoint %>% filter(cohort == "GSE6710")
pair_total <- sum(gse416$complete_lesional_nonlesional_pairs, gse6710$complete_lesional_nonlesional_pairs, gse154_audit$completeD2Pairs)
sample_total <- sum(gse416$complete_lesional_nonlesional_pairs * 2, gse6710$complete_lesional_nonlesional_pairs * 2, gse154_audit$includedD2Samples)

if (pair_total != 154) stop_quietly("Pair total is not 154.")
if (sample_total != 308) stop_quietly("Included sample total is not 308.")
if (nrow(inc_main) != 18) stop_quietly("Increment comparison table does not contain 18 rows.")

status_counts <- table(factor(inc_main$incremental_status, levels = c("supported", "not_supported", "uncertain")))
if (status_counts[["supported"]] != 3) stop_quietly("supported count is not 3.")
if (status_counts[["not_supported"]] != 14) stop_quietly("not_supported count is not 14.")
if (status_counts[["uncertain"]] != 1) stop_quietly("uncertain count is not 1.")

plos14_154 <- inc_main %>% filter(cohort == "GSE154200", signature_id == "PLOS14")
if (nrow(plos14_154) != 1) stop_quietly("PLOS14 GSE154200 row is missing.")
if (abs(plos14_154$delta_vs_baseline_combined - 0.026860) > 1e-6) stop_quietly("PLOS14 Delta AUC is not frozen value.")
if (plos14_154$incremental_status != "uncertain") stop_quietly("PLOS14 GSE154200 status is not uncertain.")

ed_416_all <- effective_dimension %>% filter(cohort == "GSE41664", context == "all_D2_reference")
ed_416_delta <- effective_dimension %>% filter(cohort == "GSE41664", context == "D2_paired_delta")
rho_pair <- spearman_long %>%
  filter(cohort == "GSE41664", context == "all_D2_reference", signature_1 == "RiboSisRG", signature_2 == "PLOS14")
if (abs(ed_416_all$participation_effective_dimension - 2.061726) > 1e-6) stop_quietly("GSE41664 all-D2 effective dimension changed.")
if (abs(ed_416_delta$participation_effective_dimension - 3.141494) > 1e-6) stop_quietly("GSE41664 paired-change effective dimension changed.")
if (abs(rho_pair$spearman_rho - 0.933960) > 1e-6) stop_quietly("RiboSisRG-PLOS14 rho changed.")

# Figure 1: cohort and sample flow ------------------------------------------

box_specs <- tibble(
  cohort = c("GSE41664", "GSE6710", "GSE154200"),
  xmin = c(5, 37, 69), xmax = c(31, 63, 95), ymin = 49, ymax = 89,
  x_text = xmin + 2, center = (xmin + xmax) / 2
)
box_headers <- tibble(cohort = box_specs$cohort, x_text = box_specs$x_text, y = 83.5)
box_body <- tribble(
  ~cohort, ~line_index, ~text,
  "GSE41664", 1, "GPL570 (Affymetrix)",
  "GSE41664", 2, "157 annotated samples",
  "GSE41664", 3, "53 complete D2 pairs",
  "GSE41664", 4, "38 untreated + 15 trial-baseline",
  "GSE6710", 1, "GPL96 (Affymetrix)",
  "GSE6710", 2, "26 annotated samples",
  "GSE6710", 3, "13 complete D2 pairs",
  "GSE6710", 4, "All untreated",
  "GSE154200", 1, "GPL20301 (RNA-seq)",
  "GSE154200", 2, "376 metadata samples",
  "GSE154200", 3, "178 psoriasis samples",
  "GSE154200", 4, "88 complete D2 pairs",
  "GSE154200", 5, "176 included D2 samples",
  "GSE154200", 6, "2 unpaired lesional excluded"
) %>% mutate(y = 77 - (line_index - 1) * 5.0, x_text = box_specs$x_text[match(cohort, box_specs$cohort)])

p_fig1 <- ggplot() +
  coord_cartesian(xlim = c(0, 100), ylim = c(0, 100), expand = FALSE, clip = "off") +
  geom_rect(data = box_specs, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, color = cohort), fill = "white", linewidth = 1.15) +
  geom_rect(data = box_specs, aes(xmin = xmin, xmax = xmax, ymin = ymax - 2.2, ymax = ymax, fill = cohort), color = NA) +
  geom_text(data = box_headers, aes(x = x_text, y = y, label = cohort, color = cohort), hjust = 0, fontface = "bold", size = 5.1) +
  geom_text(data = box_body, aes(x = x_text, y = y, label = text, color = cohort), hjust = 0, size = 3.2, lineheight = 1.05) +
  geom_segment(data = box_specs, aes(x = center, xend = center, y = ymin, yend = 43, color = cohort), linewidth = 0.8, show.legend = FALSE) +
  geom_segment(data = box_specs, aes(x = xmin, xend = xmax, y = 43, yend = 43, color = cohort), linewidth = 0.8, show.legend = FALSE) +
  geom_segment(aes(x = 18, xend = 82, y = 43, yend = 43), color = "#555555", linewidth = 0.8) +
  geom_segment(aes(x = 50, xend = 50, y = 43, yend = 35.4), color = "#555555", linewidth = 0.9, arrow = arrow(length = unit(3, "mm"), type = "closed")) +
  geom_rect(aes(xmin = 15, xmax = 85, ymin = 16, ymax = 34), fill = "#0B4F6C", color = NA) +
  annotate("text", x = 50, y = 29.2, label = "154 complete lesion / non-lesion pairs", color = "white", fontface = "bold", size = 5.7) +
  annotate("text", x = 50, y = 23.5, label = "308 included D2 samples  |  3 cohorts  |  10 published signatures  |  18 incremental comparisons", color = "white", size = 3.5) +
  annotate("text", x = 50, y = 19.1, label = "Endpoint: lesional versus non-lesional skin within the same patient", color = "#D9F2EA", size = 3.3) +
  annotate("text", x = 50, y = 9, label = "GSE154200 count-matrix columns and metadata samples are not treated as interchangeable; only complete psoriasis pairs enter the benchmark.", color = "#555555", size = 3.05) +
  scale_color_manual(values = cohort_palette, guide = "none") +
  scale_fill_manual(values = cohort_palette, guide = "none") +
  theme_void(base_family = "Arial") + theme(plot.margin = margin(8, 12, 8, 12))

# Figure 2: signature-level AUC by cohort ------------------------------------

d2_all <- bind_rows(
  d2_41664 %>% mutate(primary_D2_eligible = primary_D2_eligible == "yes"),
  d2_6710 %>% mutate(primary_D2_eligible = primary_D2_eligible == "yes"),
  d2_154200 %>% mutate(primary_D2_eligible = TRUE)
) %>%
  mutate(
    cohort = factor(cohort, levels = names(cohort_palette)),
    signature_id = factor(signature_id, levels = signature_order),
    selected_gene_count = coalesce(selected_gene_count, genes_used),
    signature_gene_count = coalesce(signature_gene_count, genes_declared),
    coverage_percent_plot = if_else(!is.na(coverage_percent), coverage_percent, 100),
    point_class = case_when(
      !primary_D2_eligible ~ "Excluded / descriptive only",
      coverage_percent_plot < 100 ~ "Primary comparison; partial coverage",
      TRUE ~ "Primary comparison"
    ),
    point_class = factor(point_class, levels = names(point_palette))
  )

if (nrow(d2_all) != 28) stop_quietly("Figure 2 source does not contain 28 signature-cohort rows.")
if (sum(d2_all$cohort == "GSE154200" & d2_all$signature_id == "IMscore" & d2_all$coverage_percent_plot == 80) != 1) {
  stop_quietly("GSE154200 IMscore partial-coverage row is missing.")
}

fig2_ci <- d2_all %>% filter(cohort == "GSE154200")
p_fig2 <- ggplot(d2_all, aes(x = signature_id, y = auc)) +
  geom_hline(yintercept = 0.5, linetype = "dotted", color = "#A0A0A0", linewidth = 0.45) +
  geom_pointrange(
    data = fig2_ci,
    aes(ymin = auc_ci_low, ymax = auc_ci_high, color = cohort),
    linewidth = 0.55,
    size = 0.42
  ) +
  geom_point(aes(fill = point_class, color = cohort), shape = 21, size = 2.75, stroke = 0.65) +
  geom_text(
    data = fig2_ci %>% filter(signature_id == "IMscore"),
    aes(label = "4/5"), vjust = -1.2, size = 2.55, color = "#6B4C00"
  ) +
  facet_wrap(~ cohort, nrow = 1, drop = FALSE) +
  scale_y_continuous(limits = c(0, 1.02), breaks = seq(0, 1, 0.2), expand = c(0.01, 0)) +
  scale_x_discrete(drop = FALSE) +
  scale_color_manual(values = cohort_palette, guide = "none") +
  scale_fill_manual(values = point_palette, name = "Analysis role") +
  labs(
    x = NULL,
    y = "D2 AUC (lesional versus non-lesional)",
    caption = "GSE41664 and GSE6710 show point estimates because the frozen source has no separate 95% CI. GSE154200 intervals come only from its cohort-specific D2 metrics table. Blank positions indicate signatures unavailable in a cohort analysis. No cross-cohort ranking is implied."
  ) +
  theme_peerj +
  theme(
    axis.text.x = element_text(angle = 50, hjust = 1, vjust = 1, size = 7.5),
    legend.position = "bottom",
    plot.caption = element_text(size = 7.2, color = "#555555", hjust = 0, margin = margin(t = 8))
  )

# Figure 3: incremental value over clinical/biological baselines ------------

status_levels <- c("supported", "uncertain", "not_supported")
inc_plot <- inc_main %>%
  mutate(
    cohort = factor(cohort, levels = names(cohort_palette)),
    signature_id = factor(signature_id, levels = rev(signature_order)),
    incremental_status = factor(incremental_status, levels = status_levels),
    coverage_note = if_else(cohort == "GSE154200" & signature_id == "IMscore", "4/5 genes", NA_character_)
  )

p_fig3 <- ggplot(inc_plot, aes(x = delta_vs_baseline_combined, y = signature_id)) +
  geom_vline(xintercept = 0, color = "#555555", linewidth = 0.55) +
  geom_vline(xintercept = 0.03, color = "#E69F00", linetype = "dashed", linewidth = 0.65) +
  geom_pointrange(
    aes(xmin = delta_ci_low, xmax = delta_ci_high, color = incremental_status),
    linewidth = 0.6,
    size = 0.43
  ) +
  geom_text(
    data = inc_plot %>% filter(!is.na(coverage_note)),
    aes(label = coverage_note),
    hjust = -0.20,
    nudge_x = 0.002,
    size = 2.55,
    color = "#6B4C00"
  ) +
  facet_wrap(~ cohort, nrow = 1, scales = "free_y") +
  scale_color_manual(
    values = status_palette,
    breaks = status_levels,
    labels = c("Supported", "Uncertain", "Not supported"),
    name = "Increment status"
  ) +
  labs(
    x = "Delta AUC versus combined KC, inflammatory, IFN-gamma and proliferation baselines",
    y = NULL,
    caption = paste(
      strwrap(
        "The dashed line is the prespecified minimum meaningful Delta AUC of 0.03. GSE6710 has a baseline AUC ceiling of 1.000 in its 13 pairs. Zero-width intervals reflect identical out-of-fold predictions when the added signature contributes no ranking information.",
        width = 145
      ),
      collapse = "\n"
    )
  ) +
  coord_cartesian(clip = "off") +
  theme_peerj +
  theme(
    legend.position = "bottom",
    panel.grid.major.y = element_line(color = "#ECECEC", linewidth = 0.3),
    panel.grid.major.x = element_line(color = "#E6E6E6", linewidth = 0.35),
    plot.caption = element_text(size = 7.2, color = "#555555", hjust = 0, margin = margin(t = 8))
  )

# Figure 4: correlation structure and effective dimension -------------------

primary_signatures <- c("RiboSisRG", "TRGS", "Apoptosis5", "m6A10", "PLOS14")

make_corr_grid <- function(context_name) {
  source_context <- spearman_long %>%
    filter(cohort == "GSE41664", context == context_name)

  if (nrow(source_context) != 10) {
    stop_quietly(paste("Expected 10 pairwise correlations for", context_name))
  }

  mat <- matrix(NA_real_, nrow = length(primary_signatures), ncol = length(primary_signatures),
                dimnames = list(primary_signatures, primary_signatures))
  diag(mat) <- 1
  for (i in seq_len(nrow(source_context))) {
    a <- source_context$signature_1[i]
    b <- source_context$signature_2[i]
    mat[a, b] <- source_context$spearman_rho[i]
    mat[b, a] <- source_context$spearman_rho[i]
  }
  if (anyNA(mat)) stop_quietly(paste("Incomplete correlation matrix for", context_name))

  as_tibble(mat, rownames = "signature_1") %>%
    pivot_longer(cols = -signature_1, names_to = "signature_2", values_to = "spearman_rho") %>%
    mutate(context = context_name)
}

corr_all <- make_corr_grid("all_D2_reference")
corr_delta <- make_corr_grid("D2_paired_delta")

corr_plot_data <- bind_rows(corr_all, corr_delta) %>%
  mutate(
    signature_1 = factor(signature_1, levels = rev(primary_signatures)),
    signature_2 = factor(signature_2, levels = primary_signatures),
    text_color = if_else(abs(spearman_rho) >= 0.65, "white", "#222222")
  )

corr_panel <- function(context_name, panel_title) {
  d <- corr_plot_data %>% filter(context == context_name)
  ggplot(d, aes(x = signature_2, y = signature_1, fill = spearman_rho)) +
    geom_tile(color = "white", linewidth = 0.8) +
    geom_text(aes(label = sprintf("%.2f", spearman_rho), color = text_color), size = 3.0) +
    scale_color_identity() +
    scale_fill_gradient2(
      low = "#0072B2",
      mid = "#F7F7F7",
      high = "#D55E00",
      midpoint = 0,
      limits = c(-1, 1),
      name = "Spearman rho"
    ) +
    coord_equal() +
    labs(title = panel_title, x = NULL, y = NULL) +
    theme_minimal(base_size = 9, base_family = "Arial") +
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 7.7),
      axis.text.y = element_text(size = 7.7),
      plot.title = element_text(face = "bold", size = 9.5, margin = margin(b = 7)),
      legend.title = element_text(size = 8, face = "bold"),
      legend.text = element_text(size = 7.5),
      legend.key.height = unit(0.55, "cm"),
      legend.key.width = unit(0.28, "cm")
    )
}

p_corr_all <- corr_panel("all_D2_reference", "All D2 samples")
p_corr_delta <- corr_panel("D2_paired_delta", "Patient-matched changes")

ed_fig4 <- effective_dimension %>%
  filter(
    cohort == "GSE41664",
    context %in% c("all_D2_reference", "D2_paired_delta")
  ) %>%
  mutate(
    context_label = recode(
      context,
      all_D2_reference = "All D2 samples",
      D2_paired_delta = "Patient-matched changes"
    ),
    context_label = factor(context_label, levels = c("All D2 samples", "Patient-matched changes"))
  )

p_effective <- ggplot(ed_fig4, aes(x = context_label, y = participation_effective_dimension)) +
  geom_col(aes(fill = context_label), width = 0.56, show.legend = FALSE) +
  geom_text(
    aes(label = sprintf("%.3f / 5", participation_effective_dimension)),
    vjust = -0.45,
    size = 3.2,
    color = "#222222"
  ) +
  scale_fill_manual(values = c("All D2 samples" = "#0072B2", "Patient-matched changes" = "#D55E00")) +
  scale_y_continuous(limits = c(0, 5.3), breaks = 0:5, expand = c(0, 0)) +
  labs(
    title = "Effective dimension",
    x = NULL,
    y = "Participation ratio"
  ) +
  theme_peerj +
  theme(
    axis.text.x = element_text(angle = 20, hjust = 1, size = 7.8),
    plot.title = element_text(face = "bold", size = 9.5, margin = margin(b = 7)),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_line(color = "#E6E6E6", linewidth = 0.35)
  )

p_fig4 <- (p_corr_all | p_corr_delta | p_effective) +
  plot_layout(widths = c(1, 1, 0.82), guides = "collect") +
  plot_annotation(
    tag_levels = "A",
    caption = "GSE41664 contains the five primary signatures. GSE6710 correlations are retained only as a descriptive sensitivity analysis in Supplemental Table S6."
  )

# Figure export --------------------------------------------------------------

save_figure <- function(plot, stem, width, height) {
  png_path <- file.path(figure_dir, paste0(stem, ".png"))
  pdf_path <- file.path(figure_dir, paste0(stem, ".pdf"))

  ggsave(
    filename = png_path,
    plot = plot,
    width = width,
    height = height,
    units = "in",
    dpi = 300,
    bg = "white",
    limitsize = FALSE
  )
  ggsave(
    filename = pdf_path,
    plot = plot,
    width = width,
    height = height,
    units = "in",
    device = if (capabilities("cairo")) cairo_pdf else pdf,
    bg = "white",
    limitsize = FALSE
  )
  c(png_path, pdf_path)
}

figure_files <- c(
  save_figure(p_fig1, "Figure_1_Cohort_And_Sample_Flow", 9.0, 6.5),
  save_figure(p_fig2, "Figure_2_D2_Signature_AUC", 10.0, 5.2),
  save_figure(p_fig3, "Figure_3_Incremental_Delta_AUC", 10.0, 5.8),
  save_figure(p_fig4, "Figure_4_Correlation_And_Effective_Dimension", 11.0, 4.8)
)

# Figure source data ---------------------------------------------------------

fig1_source <- bind_rows(
  tibble(
    cohort = "GSE41664",
    stage = c("annotated_samples", "complete_D2_pairs", "included_D2_samples", "untreated_pairs", "trial_baseline_pairs"),
    value = c(157, 53, 106, 38, 15)
  ),
  tibble(
    cohort = "GSE6710",
    stage = c("annotated_samples", "complete_D2_pairs", "included_D2_samples", "untreated_pairs", "trial_baseline_pairs"),
    value = c(26, 13, 26, 13, 0)
  ),
  tibble(
    cohort = "GSE154200",
    stage = c(
      "metadata_samples", "psoriasis_samples", "complete_D2_pairs",
      "included_D2_samples", "unpaired_lesional_excluded", "count_matrix_columns"
    ),
    value = c(
      gse154_audit$metadataSamples,
      gse154_audit$psoriasisSamples,
      gse154_audit$completeD2Pairs,
      gse154_audit$includedD2Samples,
      gse154_audit$psoriasisSamples - gse154_audit$includedD2Samples,
      gse154_audit$countMatrixSamples
    )
  )
)

write_project_csv(fig1_source, "Figure_Data/Figure_1_Source_Data.csv")
write_project_csv(
  d2_all %>%
    mutate(
      cohort = as.character(cohort),
      signature_id = as.character(signature_id),
      point_class = as.character(point_class)
    ) %>%
    select(
      cohort, signature_id, auc, auc_ci_low, auc_ci_high,
      primary_D2_eligible, coverage_percent_plot, selected_gene_count,
      signature_gene_count, point_class
    ),
  "Figure_Data/Figure_2_Source_Data.csv"
)
write_project_csv(
  inc_plot %>%
    mutate(
      cohort = as.character(cohort),
      signature_id = as.character(signature_id),
      incremental_status = as.character(incremental_status)
    ) %>%
    select(
      cohort, signature_id, delta_vs_baseline_combined,
      delta_ci_low, delta_ci_high, min_meaningful_delta,
      incremental_status, oof_auc, auc_ci_low, auc_ci_high, coverage_note
    ),
  "Figure_Data/Figure_3_Source_Data.csv"
)
write_project_csv(
  corr_plot_data %>%
    mutate(
      signature_1 = as.character(signature_1),
      signature_2 = as.character(signature_2)
    ) %>%
    select(context, signature_1, signature_2, spearman_rho),
  "Figure_Data/Figure_4A_Spearman_Matrix.csv"
)
write_project_csv(
  ed_fig4 %>%
    mutate(context_label = as.character(context_label)) %>%
    select(context, context_label, n, signatures, participation_effective_dimension,
           participation_fraction_of_signatures, first_eigenvalue_fraction,
           median_off_diagonal_rho),
  "Figure_Data/Figure_4B_Effective_Dimension.csv"
)

# Supplemental Table S1: cohort flow and treatment audit --------------------

gse154_treatment <- paste0(
  gse154_audit$treatmentValueCounts$value,
  " = ", gse154_audit$treatmentValueCounts$samples,
  collapse = "; "
)

s1_core <- tibble(
  cohort = c("GSE41664", "GSE6710", "GSE154200"),
  platform = c("GPL570", "GPL96", "GPL20301"),
  annotated_or_metadata_samples = c(157, 26, gse154_audit$metadataSamples),
  psoriasis_samples = c(106, 26, gse154_audit$psoriasisSamples),
  complete_D2_pairs = c(53, 13, gse154_audit$completeD2Pairs),
  included_D2_samples = c(106, 26, gse154_audit$includedD2Samples),
  untreated_pairs = c(38, 13, 88),
  treatment_baseline_pairs = c(15, 0, 0),
  samples_not_in_included_pairs = c(51, 0, 200),
  psoriasis_samples_excluded_from_pairs = c(0, 0, 2),
  treatment_field_audit = c(
    "38 untreated pairs + 15 trial-baseline pairs; treatment baseline retained and flagged.",
    "13 untreated pairs; no treatment baseline stratum.",
    paste0("Metadata treatment values: ", gse154_treatment, ". All 176 included samples are no treatment.")
  ),
  benchmark_role = c(
    "Primary three-cohort benchmark; mixed untreated and trial-baseline pairs.",
    "Primary three-cohort benchmark; small entirely untreated cohort.",
    "Independent D2 cohort; complete psoriasis pairs only. Metadata and count-matrix columns are not interchangeable."
  )
)

s1_total <- tibble(
  cohort = "Total",
  platform = "",
  annotated_or_metadata_samples = sum(s1_core$annotated_or_metadata_samples),
  psoriasis_samples = sum(s1_core$psoriasis_samples),
  complete_D2_pairs = sum(s1_core$complete_D2_pairs),
  included_D2_samples = sum(s1_core$included_D2_samples),
  untreated_pairs = sum(s1_core$untreated_pairs),
  treatment_baseline_pairs = sum(s1_core$treatment_baseline_pairs),
  samples_not_in_included_pairs = sum(s1_core$samples_not_in_included_pairs),
  psoriasis_samples_excluded_from_pairs = sum(s1_core$psoriasis_samples_excluded_from_pairs),
  treatment_field_audit = "Totals preserve the separate untreated and trial-baseline strata.",
  benchmark_role = "154 complete pairs; 308 included D2 samples."
)

s1 <- bind_rows(s1_core, s1_total)
if (nrow(s1) != 4) stop_quietly("Supplemental Table S1 does not have four rows.")
if (s1$complete_D2_pairs[s1$cohort == "Total"] != 154) stop_quietly("S1 pair total is not 154.")
write_project_csv(s1, "Supplemental_Tables/S1_Cohort_Flow_And_Treatment_Audit.csv")

# Supplemental Table S2: published signature registry -----------------------

coverage_gpl570 <- platform_coverage %>%
  filter(platform == "GPL570") %>%
  select(signature_id, gpl570_present_gene_count = present_gene_count,
         gpl570_missing_gene_count = missing_gene_count,
         gpl570_coverage_percent = coverage_percent)

s2 <- signature_registry %>%
  left_join(coverage_gpl570, by = "signature_id") %>%
  mutate(
    registry_order = match(signature_id, signature_order)
  ) %>%
  arrange(registry_order) %>%
  select(
    signature_id, label, pmid, pmcid, source_record, evidence_grade,
    model_status, comparator, original_endpoint, compatible_endpoints,
    gene_count, genes, derivation_cohorts, reported_validation_cohorts,
    gpl570_present_gene_count, gpl570_missing_gene_count,
    gpl570_coverage_percent, notes
  )
if (nrow(s2) != 10) stop_quietly("Supplemental Table S2 does not have 10 signatures.")
write_project_csv(s2, "Supplemental_Tables/S2_Signature_Registry.csv")

# Supplemental Table S3: platform coverage ----------------------------------

gpl_coverage <- platform_coverage %>%
  filter(platform %in% c("GPL570", "GPL96")) %>%
  transmute(
    platform,
    signature_id,
    signature_label,
    signature_size,
    present_gene_count,
    missing_gene_count,
    coverage_percent,
    missing_genes,
    evaluation_status = case_when(
      missing_gene_count == 0 ~ "evaluated_complete",
      TRUE ~ "evaluated_partial_coverage"
    ),
    evaluation_note = if_else(
      missing_gene_count == 0,
      "",
      "The observed D2 score uses the available registered genes only."
    )
  )

gse154_not_evaluated <- signature_registry %>%
  filter(signature_id %in% c("HKG34", "Epidermal10")) %>%
  transmute(
    platform = "GPL20301",
    signature_id,
    signature_label = label,
    signature_size = gene_count,
    present_gene_count = NA_real_,
    missing_gene_count = NA_real_,
    coverage_percent = NA_real_,
    missing_genes = NA_character_,
    evaluation_status = "not_evaluated_original_endpoint_D1",
    evaluation_note = "Not evaluated for GSE154200 because the original endpoint is D1 rather than D2."
  )

gse154_evaluated <- gse154_gene_coverage %>%
  filter(feature_type == "signature", !feature_id %in% c("HKG34", "Epidermal10")) %>%
  left_join(
    signature_registry %>% select(signature_id, label, gene_count),
    by = c("feature_id" = "signature_id")
  ) %>%
  transmute(
    platform = "GPL20301",
    signature_id = feature_id,
    signature_label = label,
    signature_size = gene_count,
    present_gene_count = genes_used,
    missing_gene_count = genes_declared - genes_used,
    coverage_percent,
    missing_genes,
    evaluation_status = if_else(
      genes_used == genes_declared,
      "evaluated_complete",
      "evaluated_partial_coverage"
    ),
    evaluation_note = if_else(
      genes_used == genes_declared,
      "",
      "The observed D2 score uses the available registered genes only."
    )
  )

s3 <- bind_rows(gpl_coverage, gse154_evaluated, gse154_not_evaluated) %>%
  mutate(
    signature_order = match(signature_id, signature_order),
    platform_order = match(platform, c("GPL570", "GPL96", "GPL20301"))
  ) %>%
  arrange(signature_order, platform_order) %>%
  select(-signature_order, -platform_order)
if (nrow(s3) != 30) stop_quietly("Supplemental Table S3 does not have 30 platform-signature rows.")
write_project_csv(s3, "Supplemental_Tables/S3_Platform_Gene_Coverage.csv")

# Supplemental Table S4: cohort-specific D2 AUC -----------------------------

s4 <- d2_all %>%
  mutate(
    cohort = as.character(cohort),
    signature_id = as.character(signature_id),
    d2_pairs = coalesce(d2_pairs, patients),
    scored_pairs = coalesce(scored_pairs, patients),
    samples = coalesce(samples, scored_pairs * 2),
    selected_gene_count = coalesce(selected_gene_count, genes_used),
    signature_gene_count = coalesce(signature_gene_count, genes_declared),
    coverage_percent = coalesce(coverage_percent, 100),
    coverage_status = coalesce(
      coverage_status,
      if_else(selected_gene_count == signature_gene_count, "complete", "partial")
    ),
    ci_source = if_else(cohort == "GSE154200",
                        "GSE154200 cohort-specific D2_signature_metrics.csv",
                        "not_available_in_frozen_source"),
    analysis_role = as.character(point_class)
  ) %>%
  select(
    cohort, platform, signature_id, d2_pairs, scored_pairs, samples,
    signature_gene_count, selected_gene_count, coverage_percent,
    coverage_status, auc, auc_ci_low, auc_ci_high, ci_source,
    analysis_role
  )
if (nrow(s4) != 28) stop_quietly("Supplemental Table S4 does not have 28 rows.")
if (sum(!is.na(s4$auc_ci_low)) != 8) stop_quietly("S4 should have eight GSE154200 CI rows.")
write_project_csv(s4, "Supplemental_Tables/S4_D2_Signature_AUC.csv")

# Supplemental Table S5: incremental comparisons ----------------------------

baseline_join <- baseline_main %>%
  select(
    cohort,
    signature_id,
    baseline_auc = oof_auc,
    baseline_auc_ci_low = auc_ci_low,
    baseline_auc_ci_high = auc_ci_high
  )

s5 <- inc_main %>%
  left_join(baseline_join, by = c("cohort", "signature_id")) %>%
  mutate(
    delta_ci_excludes_zero = delta_ci_low > 0 | delta_ci_high < 0,
    delta_ci_excludes_minimum = delta_ci_low > min_meaningful_delta,
    direction_of_delta = case_when(
      delta_vs_baseline_combined > 0 ~ "positive",
      delta_vs_baseline_combined < 0 ~ "negative",
      TRUE ~ "zero"
    )
  ) %>%
  select(
    cohort, signature_id, patients, samples, baseline_auc,
    baseline_auc_ci_low, baseline_auc_ci_high, oof_auc,
    auc_ci_low, auc_ci_high, delta_vs_baseline_combined,
    delta_ci_low, delta_ci_high, min_meaningful_delta,
    delta_ci_excludes_zero, delta_ci_excludes_minimum,
    direction_of_delta, incremental_status
  )
if (nrow(s5) != 18) stop_quietly("Supplemental Table S5 does not have 18 rows.")
write_project_csv(s5, "Supplemental_Tables/S5_Incremental_Comparisons.csv")

# Supplemental Table S6: pairwise Spearman structure ------------------------

s6 <- spearman_long %>%
  mutate(
    cohort_role = if_else(
      cohort == "GSE6710",
      "descriptive_sensitivity_only_not_confirmatory",
      "primary_structure_analysis"
    )
  ) %>%
  select(cohort, cohort_role, context, unit, n, signature_1, signature_2, spearman_rho)
if (nrow(s6) != 80) stop_quietly("Supplemental Table S6 does not have 80 rows.")
write_project_csv(s6, "Supplemental_Tables/S6_Spearman_Pairwise.csv")
write_project_csv(
  spearman_summary %>%
    mutate(
      cohort_role = if_else(
        cohort == "GSE6710",
        "descriptive_sensitivity_only_not_confirmatory",
        "primary_structure_analysis"
      )
    ),
  "Supplemental_Tables/S6b_Spearman_Pair_Summary.csv"
)

# Supplemental Data S1: GSE154200 audit bundle ------------------------------

s1_data_dir <- file.path(data_dir, "S1_GSE154200_Audit")
dir.create(s1_data_dir, recursive = TRUE, showWarnings = FALSE)

audit_key_values <- bind_rows(
  tibble(
    section = "sample_inventory",
    key = c(
      "metadata_samples", "count_matrix_samples", "psoriasis_samples",
      "psoriasis_patients", "complete_D2_pairs", "included_D2_samples",
      "included_lesional_samples", "included_nonlesional_samples"
    ),
    value = as.character(c(
      gse154_audit$metadataSamples, gse154_audit$countMatrixSamples,
      gse154_audit$psoriasisSamples, gse154_audit$psoriasisPatients,
      gse154_audit$completeD2Pairs, gse154_audit$includedD2Samples,
      gse154_audit$includedD2Lesional, gse154_audit$includedD2Nonlesional
    )),
    source = "data/benchmark/GSE154200_audit_summary.json"
  ),
  tibble(
    section = "adjudication",
    key = c(
      "psoriasis_definition", "pair_definition", "exclusion_rule",
      "treatment_audit_rule", "treatment_protocol_source",
      "metadata_count_title_unique", "missing_metadata_titles_in_counts"
    ),
    value = c(
      gse154_audit$predeclaredRules$disease,
      gse154_audit$predeclaredRules$pair,
      gse154_audit$predeclaredRules$exclusions,
      gse154_audit$predeclaredRules$treatmentAudit,
      gse154_audit$predeclaredRules$treatmentProtocol,
      as.character(gse154_audit$countTitleUnique),
      length(gse154_audit$metadataTitlesMissingFromCounts)
    ),
    source = "data/benchmark/GSE154200_audit_summary.json"
  ),
  tibble(
    section = "interpretation_guard",
    key = c("metadata_count_unit", "inclusion_unit", "matrix_extra_columns_rule"),
    value = c(
      "Metadata samples and count-matrix columns are separate audit units.",
      "Only complete one-lesional and one-nonlesional psoriasis patient pairs enter D2 scoring.",
      "Extra count-matrix columns are not promoted to metadata samples or benchmark pairs."
    ),
    source = "Derived from the frozen audit workflow; preserve when reusing this bundle."
)
)
write_project_csv(audit_key_values, "Supplemental_Data/S1_GSE154200_Audit/S1_Audit_Key_Values.csv")

copy_s1_files <- c(
  "data/benchmark/score_stage_independent/GSE154200/preprocessing_summary.csv" = "S1_Preprocessing_Summary.csv",
  "data/benchmark/GSE154200_pair_audit.csv" = "S1_Pair_Audit.csv",
  "data/benchmark/GSE154200_cohort_screening.csv" = "S1_Sample_Screening.csv"
)
for (source_rel in names(copy_s1_files)) {
  if (!file.copy(
    from = file.path(root_dir, source_rel),
    to = file.path(s1_data_dir, copy_s1_files[[source_rel]]),
    overwrite = TRUE
  )) {
    stop_quietly(paste("Could not copy Supplemental Data S1 input:", source_rel))
  }
}

writeLines(
  c(
    "# Supplemental Data S1",
    "",
    "GSE154200 preprocessing, sample-screening and pair-audit bundle.",
    "The included D2 set contains 88 complete psoriasis pairs (176 samples).",
    "The metadata contain 376 samples; the count matrix contains 564 columns.",
    "These two totals are not interchangeable. Extra count-matrix columns are not promoted to metadata samples.",
    "",
    "Files:",
    "- `S1_Audit_Key_Values.csv`: normalized audit decisions and inventory.",
    "- `S1_Preprocessing_Summary.csv`: normalization and gene-filter parameters.",
    "- `S1_Pair_Audit.csv`: patient-level complete-pair audit.",
    "- `S1_Sample_Screening.csv`: sample-level diagnosis, tissue and treatment records."
  ),
  file.path(s1_data_dir, "README.md"),
  useBytes = TRUE
)

# Supplemental Data S2: code, runtime and manifests -------------------------

s2_data_dir <- file.path(data_dir, "S2_Code_And_Manifest")
dir.create(s2_data_dir, recursive = TRUE, showWarnings = FALSE)

script_args <- commandArgs(trailingOnly = FALSE)
script_arg <- grep("^--file=", script_args, value = TRUE)
script_path <- if (length(script_arg) > 0) {
  normalizePath(sub("^--file=", "", script_arg[[1]]), winslash = "/", mustWork = TRUE)
} else {
  file.path(root_dir, "scripts", "generate_peerj_figures_supplements.R")
}

script_target <- file.path(s2_data_dir, "generate_peerj_figures_supplements.R")
if (!file.copy(script_path, script_target, overwrite = TRUE)) {
  stop_quietly("Could not copy the generation script into Supplemental Data S2.")
}

session_path <- file.path(s2_data_dir, "S2_sessionInfo.txt")
writeLines(capture.output(sessionInfo()), session_path, useBytes = TRUE)

sha256_file <- function(path) {
  digest(file = path, algo = "sha256", serialize = FALSE)
}

manifest_row <- function(path, base_dir) {
  info <- file.info(path)
  full_path <- normalizePath(path, winslash = "/", mustWork = TRUE)
  full_base <- normalizePath(base_dir, winslash = "/", mustWork = TRUE)
  relative_path <- substring(full_path, nchar(full_base) + 2)
  file_hash <- sha256_file(full_path)
  tibble(
    path = relative_path,
    bytes = as.numeric(info$size),
    last_write = format(info$mtime, "%Y-%m-%d %H:%M:%S %Z"),
    sha256 = file_hash
  )
}

input_rel_paths <- c(
  "artifacts/submission/PeerJ_2026-11-18/01_Number_Freeze_PeerJ.md",
  "artifacts/submission/PeerJ_2026-11-18/02_Manuscript_Skeleton_PeerJ.md",
  "artifacts/submission/PeerJ_2026-11-18/03_Figures_and_Supplements/00_Figure_Data_Source_Map.md",
  "data/benchmark/cohort_endpoint_matrix.csv",
  "data/benchmark/signature_registry.csv",
  "data/benchmark/platform_gene_coverage.csv",
  "data/benchmark/GSE154200_audit_summary.json",
  "data/benchmark/GSE154200_cohort_screening.csv",
  "data/benchmark/GSE154200_pair_audit.csv",
  "data/benchmark/score_stage/GSE41664/D2_signature_metrics.csv",
  "data/benchmark/score_stage/GSE6710/D2_signature_metrics.csv",
  "data/benchmark/score_stage_independent/GSE154200/increment_stage/D2_signature_metrics.csv",
  "data/benchmark/score_stage/increment_stage/increment_model_metrics.csv",
  "data/benchmark/score_stage_independent/GSE154200/increment_stage/increment_model_metrics.csv",
  "data/benchmark/score_stage/D2_primary_spearman_long_all.csv",
  "data/benchmark/score_stage/D2_primary_spearman_pair_summary.csv",
  "data/benchmark/score_stage/D2_primary_effective_dimension_all.csv",
  "data/benchmark/score_stage_independent/GSE154200/D2_sample_metadata.csv",
  "data/benchmark/score_stage_independent/GSE154200/gene_coverage.csv",
  "data/benchmark/score_stage_independent/GSE154200/preprocessing_summary.csv",
  "scripts/generate_peerj_figures_supplements.R"
)

missing_inputs <- input_rel_paths[!file.exists(file.path(root_dir, input_rel_paths))]
if (length(missing_inputs) > 0) {
  stop_quietly(paste("Manifest inputs are missing:", paste(missing_inputs, collapse = ", ")))
}

input_manifest <- bind_rows(lapply(
  input_rel_paths,
  function(rel) manifest_row(file.path(root_dir, rel), root_dir)
))
write_excel_csv(input_manifest, file.path(s2_data_dir, "S2_Input_Manifest.csv"), na = "")

figure_legends <- c(
  "# Figure Legends Draft",
  "",
  "## Figure 1",
  "Three-cohort sample flow and D2 endpoint construction. The benchmark contains 53 complete pairs from GSE41664, 13 from GSE6710 and 88 from GSE154200 (154 pairs; 308 included samples). GSE154200 metadata samples, psoriasis samples, complete pairs and count-matrix columns are shown as distinct audit quantities.",
  "",
  "## Figure 2",
  "Signature-level D2 AUC across the available cohorts. Points show observed AUCs for published signatures; only GSE154200 has separate 95% confidence intervals in the frozen cohort-specific metrics source. Fill indicates primary comparison, partial gene coverage or exclusion from the primary comparison. Empty positions indicate that a signature was not available in the cohort-specific D2 analysis. The figure is descriptive and does not establish a cross-cohort ranking.",
  "",
  "## Figure 3",
  "Incremental D2 value of each published signature after adding it to combined KC, inflammatory, IFN-gamma and proliferation baselines. The solid vertical line marks zero and the dashed line marks the prespecified minimum meaningful Delta AUC of 0.03. Status labels reflect the frozen rule applied to the observed point estimate and its interval. Zero-width intervals indicate identical out-of-fold predictions.",
  "",
  "## Figure 4",
  "Correlation structure and effective dimension of five primary signatures in GSE41664. Panels show pairwise Spearman correlations across all D2 samples, pairwise correlations of patient-matched lesional-minus-non-lesional changes, and participation-ratio effective dimension. The effective dimension is 2.061726 of 5 for all D2 samples and 3.141494 of 5 for patient-matched changes. GSE6710 is retained only as a descriptive sensitivity analysis."
)
writeLines(figure_legends, file.path(output_root, "04_Figure_Legends_Draft.md"), useBytes = TRUE)

table_readme <- c(
  "# Supplemental Tables",
  "",
  "- S1: cohort flow and treatment-field audit.",
  "- S2: ten published signature records and registered gene lists.",
  "- S3: signature coverage across GPL570, GPL96 and GPL20301.",
  "- S4: cohort-specific D2 AUC table (28 rows).",
  "- S5: incremental comparisons versus combined baselines (18 rows).",
  "- S6 and S6b: detailed pairwise Spearman rho and context-level summaries.",
  "",
  "All files are CSV with UTF-8 encoding. Blank CI fields in S4 indicate that a separate interval was not present in the frozen source; values were not imputed."
)
writeLines(table_readme, file.path(table_dir, "README.md"), useBytes = TRUE)

checks <- tibble(
  check = c(
    "complete_pairs", "included_D2_samples", "increment_comparisons",
    "status_supported", "status_not_supported", "status_uncertain",
    "GSE154200_supported", "GSE154200_uncertain",
    "S3_platform_signature_rows", "S4_AUC_rows", "S5_increment_rows",
    "S6_pairwise_rows", "GSE154200_CI_rows",
    "GSE41664_effective_dimension_all", "GSE41664_effective_dimension_delta",
    "CCID_v3_write_guard"
  ),
  observed = c(
    pair_total, sample_total, nrow(inc_main),
    status_counts[["supported"]], status_counts[["not_supported"]], status_counts[["uncertain"]],
    sum(inc_main$cohort == "GSE154200" & inc_main$incremental_status == "supported"),
    sum(inc_main$cohort == "GSE154200" & inc_main$incremental_status == "uncertain"),
    nrow(s3), nrow(s4), nrow(s5), nrow(s6), sum(!is.na(s4$auc_ci_low)),
    ed_416_all$participation_effective_dimension,
    ed_416_delta$participation_effective_dimension,
    "read-only; no write operation"
  ),
  expected = c(
    "154", "308", "18", "3", "14", "1", "0", "1",
    "30", "28", "18", "80", "8", "2.061726", "3.141494", "unchanged"
  ),
  pass = c(
    pair_total == 154,
    sample_total == 308,
    nrow(inc_main) == 18,
    status_counts[["supported"]] == 3,
    status_counts[["not_supported"]] == 14,
    status_counts[["uncertain"]] == 1,
    sum(inc_main$cohort == "GSE154200" & inc_main$incremental_status == "supported") == 0,
    sum(inc_main$cohort == "GSE154200" & inc_main$incremental_status == "uncertain") == 1,
    nrow(s3) == 30,
    nrow(s4) == 28,
    nrow(s5) == 18,
    nrow(s6) == 80,
    sum(!is.na(s4$auc_ci_low)) == 8,
    abs(ed_416_all$participation_effective_dimension - 2.061726) < 1e-6,
    abs(ed_416_delta$participation_effective_dimension - 3.141494) < 1e-6,
    TRUE
  )
)
write_excel_csv(checks, file.path(s2_data_dir, "S2_Generation_Checks.csv"), na = "")

writeLines(
  c(
    "# Supplemental Data S2",
    "",
    "Reproducibility bundle for the PeerJ figure and supplemental-table generation.",
    "",
    "Files:",
    "- `generate_peerj_figures_supplements.R`: complete generation script.",
    "- `S2_sessionInfo.txt`: R version, platform and loaded-package versions.",
    "- `S2_Input_Manifest.csv`: input paths, sizes, times and SHA256 hashes.",
    "- `S2_Output_Manifest.csv`: generated-file paths and SHA256 hashes, excluding the manifest itself.",
    "- `S2_Generation_Checks.csv`: locked numeric and row-count checks.",
    "",
    "The output manifest intentionally does not hash itself because a file cannot contain its own final SHA256."
  ),
  file.path(s2_data_dir, "README.md"),
  useBytes = TRUE
)

generation_log <- c(
  "# Generation Log",
  "",
  paste0("Run time: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  paste0("R version: ", R.version.string),
  paste0("Root: ", root_dir),
  paste0("Figures exported: ", length(figure_files)),
  "All locked row-count and numeric guards passed before output generation.",
  "No CCID files were read-modified or written by this script."
)
writeLines(generation_log, file.path(s2_data_dir, "S2_Generation_Log.md"), useBytes = TRUE)

manifest_outputs <- list.files(output_root, recursive = TRUE, full.names = TRUE, all.files = FALSE)
manifest_outputs <- manifest_outputs[
  !normalizePath(manifest_outputs, winslash = "/", mustWork = TRUE) ==
    normalizePath(file.path(s2_data_dir, "S2_Output_Manifest.csv"), winslash = "/", mustWork = FALSE)
]
output_manifest <- bind_rows(lapply(manifest_outputs, manifest_row, base_dir = output_root))
write_excel_csv(
  output_manifest,
  file.path(s2_data_dir, "S2_Output_Manifest.csv"),
  na = ""
)

if (any(!checks$pass)) {
  stop_quietly("One or more generation checks failed.")
}

cat("Generated", length(figure_files), "figure files and the complete Stage 2 supplemental set.\n")
cat("Output root:", output_root, "\n")
