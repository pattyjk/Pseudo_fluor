#!/usr/bin/env Rscript
# =============================================================================
# GO Term Figures — Strong vs Weak enriched genes
# =============================================================================

library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)
library(readr)
library(forcats)

SIG_FILE <- "~/Documents/GitHub/Pseudo_fluor/significant_annotated.tsv"
OUT_DIR  <- "~/Documents/GitHub/Pseudo_fluor/figures/"
dir.create(OUT_DIR, showWarnings = FALSE)

dat <- read_tsv(SIG_FILE, show_col_types = FALSE)

# ── Shared theme ──────────────────────────────────────────────────────────────
theme_clean <- function() {
  theme_bw(base_size = 12) +
    theme(
      panel.grid.minor   = element_blank(),
      panel.grid.major.y = element_blank(),
      strip.background   = element_rect(fill = "grey92", color = NA),
      strip.text         = element_text(face = "bold", size = 11),
      legend.position    = "bottom",
      plot.title         = element_text(face = "bold", size = 13),
      plot.subtitle      = element_text(color = "grey40", size = 10),
      axis.text          = element_text(size = 10)
    )
}

direction_colors <- c("Enriched in Strong" = "#3B4EA6", "Enriched in Weak" = "#B5532A")

# ── 1. Annotation rate bar chart ─────────────────────────────────────────────
# One row per gene — collapse multi-GO rows
gene_summary <- dat %>%
  group_by(Gene, Direction) %>%
  summarise(annotated = any(!is.na(go_id)), .groups = "drop")

annotation_summary <- gene_summary %>%
  count(Direction, annotated) %>%
  mutate(status = ifelse(annotated, "Annotated", "Unannotated"))

p1 <- ggplot(annotation_summary, aes(x = Direction, y = n, fill = status)) +
  geom_col(position = "stack", width = 0.6) +
  geom_text(
    aes(label = n),
    position = position_stack(vjust = 0.5),
    color = "white", fontface = "bold", size = 4
  ) +
  scale_fill_manual(values = c("Annotated" = "#3B6D11", "Unannotated" = "grey65"),
                    name = NULL) +
  scale_x_discrete(labels = c("Enriched in Strong" = "Strong", "Enriched in Weak" = "Weak")) +
  labs(
    title    = "GO annotation coverage by enrichment direction",
    subtitle = "Number of significant genes with and without GO annotation",
    x = NULL, y = "Number of genes"
  ) +
  theme_clean()

ggsave(file.path(OUT_DIR, "fig1_annotation_rate.pdf"), p1, width = 5, height = 5)

# ── 2. Top GO terms per aspect — faceted dot plot (Strong vs Weak) ────────────
# Parse aspect from go_name prefix (C: / F: / P:)
go_dat <- dat %>%
  filter(!is.na(go_name), go_name != "NA") %>%
  mutate(
    aspect     = str_extract(go_name, "^[CFP](?=:)"),
    term_clean = str_remove(go_name, "^[CFP]:") %>% str_trim(),
    aspect_label = case_when(
      aspect == "P" ~ "Biological process",
      aspect == "F" ~ "Molecular function",
      aspect == "C" ~ "Cellular component",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(aspect))

# Pick top N terms per aspect (union across both directions)
TOP_N <- 15

top_terms <- go_dat %>%
  count(aspect_label, term_clean) %>%
  group_by(aspect_label) %>%
  slice_max(n, n = TOP_N, with_ties = FALSE) %>%
  pull(term_clean) %>%
  unique()

go_counts <- go_dat %>%
  filter(term_clean %in% top_terms) %>%
  count(Direction, aspect_label, term_clean)

# Compute total per direction for percentage
totals <- gene_summary %>% count(Direction, name = "total")

go_counts <- go_counts %>%
  left_join(totals, by = "Direction") %>%
  mutate(pct = n / total * 100)

# Order terms by total count
term_order <- go_counts %>%
  group_by(term_clean) %>%
  summarise(total = sum(n)) %>%
  arrange(total) %>%
  pull(term_clean)

go_counts <- go_counts %>%
  mutate(term_clean = factor(term_clean, levels = term_order))

p2 <- ggplot(go_counts, aes(x = pct, y = term_clean, color = Direction)) +
  geom_line(aes(group = term_clean), color = "grey75", linewidth = 0.5) +
  geom_point(size = 3, alpha = 0.9) +
  scale_color_manual(values = direction_colors, name = NULL) +
  scale_x_continuous(labels = function(x) paste0(round(x, 1), "%")) +
  facet_wrap(~ aspect_label, scales = "free_y", ncol = 1) +
  labs(
    title    = "Top GO terms: strong vs weak enriched genes",
    subtitle = "Point position = % of genes in that direction with the term",
    x = "% of genes in direction", y = NULL
  ) +
  theme_clean() +
  theme(axis.text.y = element_text(size = 9))

ggsave(file.path(OUT_DIR, "fig2_go_dotplot.pdf"), p2, width = 8, height = 14)

# ── 3. Diverging bar chart — terms enriched in Strong vs Weak ─────────────────
# For each term, compute log2 ratio of % strong / % weak
go_wide <- go_counts %>%
  dplyr::select(Direction, aspect_label, term_clean, pct) %>%
  pivot_wider(names_from = Direction, values_from = pct, values_fill = 0) %>%
  rename(strong = `Enriched in Strong`, weak = `Enriched in Weak`) %>%
  mutate(
    log2fc = log2((strong + 0.5) / (weak + 0.5)),  # pseudocount avoids Inf
    dominant = ifelse(log2fc > 0, "Enriched in Strong", "Enriched in Weak")
  )

# Take top 12 most divergent terms per aspect
top_divergent <- go_wide %>%
  group_by(aspect_label) %>%
  slice_max(abs(log2fc), n = 12, with_ties = FALSE) %>%
  ungroup()

term_order2 <- top_divergent %>%
  arrange(aspect_label, log2fc) %>%
  pull(term_clean)

top_divergent <- top_divergent %>%
  mutate(term_clean = factor(term_clean, levels = unique(term_order2)))

p3 <- ggplot(top_divergent, aes(x = log2fc, y = term_clean, fill = dominant)) +
  geom_col(width = 0.7) +
  geom_vline(xintercept = 0, linewidth = 0.4, color = "grey30") +
  scale_fill_manual(values = direction_colors, name = NULL) +
  facet_wrap(~ aspect_label, scales = "free_y", ncol = 1) +
  labs(
    title    = "Divergence of GO terms between strong and weak enriched genes",
    subtitle = "log2(% strong / % weak) — positive = more common in strong, negative = more common in weak",
    x = "log2 fold-difference", y = NULL
  ) +
  theme_clean() +
  theme(axis.text.y = element_text(size = 9))

ggsave(file.path(OUT_DIR, "fig3_go_diverging.pdf"), p3, width = 8, height = 14)

# ── 4. Stacked bar — aspect composition per direction (annotated genes only) ──
aspect_summary <- go_dat %>%
  # one row per gene × aspect (avoid double-counting a gene for same aspect)
  distinct(Gene, Direction, aspect_label) %>%
  count(Direction, aspect_label)

aspect_totals <- gene_summary %>%
  filter(annotated) %>%
  count(Direction, name = "total")

aspect_summary <- aspect_summary %>%
  left_join(aspect_totals, by = "Direction") %>%
  mutate(pct = n / total * 100)

p4 <- ggplot(aspect_summary,
             aes(x = Direction, y = pct,
                 fill = aspect_label)) +
  geom_col(position = "stack", width = 0.6) +
  geom_text(aes(label = paste0(round(pct), "%")),
            position = position_stack(vjust = 0.5),
            color = "white", fontface = "bold", size = 3.5) +
  scale_fill_manual(
    values = c("Biological process"  = "#3B4EA6",
               "Molecular function"  = "#0F6E56",
               "Cellular component"  = "#185FA5"),
    name = NULL
  ) +
  scale_x_discrete(labels = c("Enriched in Strong" = "Strong",
                              "Enriched in Weak"   = "Weak")) +
  labs(
    title    = "GO aspect composition of annotated genes",
    subtitle = "% of annotated genes assigned to each GO domain (gene may appear in multiple)",
    x = NULL, y = "% of annotated genes"
  ) +
  theme_clean()

ggsave(file.path(OUT_DIR, "fig4_aspect_composition.pdf"), p4, width = 5, height = 5)

