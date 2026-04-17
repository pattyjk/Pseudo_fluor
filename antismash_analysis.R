library(tidyverse)
library(dunn.test)
library(ggpubr)
library(ggrepel)

# ── 0. Colour palettes ────────────────────────────────────────────────────────
pond_colors       <- c("Sixty Lake"   = "#0072B2",
                       "Conness Pond" = "#E69F00")

inhibition_colors <- c("Weakly Inhibitory"   = "#56B4E9",
                       "Strongly Inhibitory" = "#E69F00")

# 7-color Okabe-Ito palette for phage families
family_colors <- c("#E69F00", "#56B4E9", "#009E73",
                   "#CC79A7", "#0072B2", "#D55E00", "#F0E442")

# ── 1. Load data ──────────────────────────────────────────────────────────────
blast    <- read_tsv("~/Documents/GitHub/Pseudo_fluor/viral_blast_hits.txt")
antismash <- read_tsv("~/Documents/GitHub/Pseudo_fluor/antismash.txt")
inhibition <- read_tsv("~/Documents/GitHub/Pseudo_fluor/isolate_inhibi_category.txt") %>%
  dplyr::rename(Isolate = Group)

# ── 2. Parse BGC types ────────────────────────────────────────────────────────
# Some BGC regions have compound types (e.g. "NRPS,azole-containing-RiPP")
# Split these so each BGC type gets its own row for counting

bgc <- antismash %>%
  mutate(Type = str_replace_all(Type, "\\s", "")) %>%
  separate_rows(Type, sep = ",") %>%           # one row per BGC type
  left_join(inhibition, by = "Isolate")        # attach inhibition metadata

# ── 3. Parse prophage taxonomy ────────────────────────────────────────────────
safe_nth <- function(lst, i) {
  map_chr(lst, ~ if (length(.x) >= i) .x[[i]] else NA_character_)
}

blast_tax <- blast %>%
  mutate(
    FullTaxonomy = str_replace_all(FullTaxonomy, "\\s+", " ") %>% str_trim(),
    tax_levels   = str_split(FullTaxonomy, ";\\s*")
  ) %>%
  mutate(
    family       = safe_nth(tax_levels, 7),
    genus        = safe_nth(tax_levels, 8),
    lowest_taxon = map_chr(tax_levels, ~ str_remove(tail(.x, 1), "\\.$"))
  ) %>%
  select(-tax_levels)

# Deduplicate: best hit per isolate × subject phage
blast_dedup <- blast_tax %>%
  group_by(Pond, Isolate, subject_id, Hit, family, genus, lowest_taxon) %>%
  dplyr::slice_max(bit_score, n = 1, with_ties = FALSE) %>%
  ungroup()

# ── 4. Build unified per-isolate summary table ────────────────────────────────
# Counts: unique phage hits, total BGCs, BGCs by broad category

phage_counts <- blast_dedup %>%
  dplyr::distinct(Pond, Isolate, subject_id) %>%
  dplyr::count(Pond, Isolate, name = "n_phages")

bgc_counts <- bgc %>%
  group_by(Isolate) %>%
  dplyr::summarise(
    n_bgc_total   = dplyr::n(),
    n_bgc_unique  = dplyr::n_distinct(Type),
    .groups = "drop"
  )

isolate_data <- inhibition %>%
  left_join(phage_counts, by = "Isolate") %>%
  left_join(bgc_counts,   by = "Isolate") %>%
  mutate(
    Inhibition_Flag = factor(Inhibition_Flag,
                             levels = c("Weakly Inhibitory",
                                        "Strongly Inhibitory"))
  )

print(isolate_data)

# ── 5. BGC type counts per isolate ───────────────────────────────────────────
bgc_type_counts <- bgc %>%
  dplyr::count(Isolate, Type) %>%
  left_join(inhibition %>% dplyr::select(Isolate, Pond, Inhibition_Flag),
            by = "Isolate")

# ── 6. Plot 1: Stacked bar — prophage family per isolate ─────────────────────
blast_dedup %>%
  dplyr::distinct(Pond, Isolate, subject_id, family) %>%
  left_join(inhibition %>% dplyr::select(Isolate, Inhibition_Flag),
            by = "Isolate") %>%
  ggplot(aes(x = Isolate, fill = family)) +
  geom_bar() +
  facet_wrap(~ Pond, scales = "free_x") +
  scale_fill_manual(values = family_colors) +
  labs(title = "Unique Prophage Hits per Isolate",
       x = "Isolate", y = "Number of Unique Phage Hits", fill = "Family") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom")

# ── 7. Plot 2: Stacked bar — BGC types per isolate ───────────────────────────
pal<-c("#771155", "#CC99BB", "#114477", "#4477AA", "#117777", "#44AAAA", "#77CCCC", "#117744", "#44AA77", "#88CCAA", "#777711", "#AAAA44", "#DDDD77", "#774411", "#AA7744", "#DDAA77", "#771122", "#AA4455", "#DD7788","#41AB5D", "#252525", "#525252", "#737373", "#969696")
bgc_type_counts %>%
  mutate(Isolate = fct_reorder(Isolate, n, .fun = sum, .desc = TRUE)) %>%
  ggplot(aes(x = Isolate, y = n, fill = Type)) +
  geom_col() +
  scale_fill_manual(values = pal)+
  facet_wrap(~ Pond, scales = "free_x") +
  labs(title = "Biosynthetic Gene Clusters per Isolate",
       x = "Isolate", y = "Number of BGC Regions", fill = "BGC Type") +
  theme_bw() +
  theme(axis.text.x  = element_text(angle = 45, hjust = 1),
        legend.position = "bottom",
        legend.text  = element_text(size = 7)) +
  guides(fill = guide_legend(ncol = 3))

# ── 8. Plot 3: Box/violin — BGC count by inhibition category ─────────────────
comparisons_inhib <- list(c("Weakly Inhibitory", "Strongly Inhibitory"))

ggplot(isolate_data, aes(x = Inhibition_Flag, y = n_bgc_total,
                         fill = Inhibition_Flag)) +
  #geom_violin(alpha = 0.4, trim = FALSE) +
  geom_boxplot(width = 0.2, outlier.shape = NA, alpha = 0.7) +
  #geom_jitter(aes(color = Pond), width = 0.1, size = 2.5, alpha = 0.9) +
  #geom_text_repel(aes(label = Isolate), size = 2.8, max.overlaps = 15) +
  stat_compare_means(comparisons = comparisons_inhib,
                    method = "wilcox.test", label = "p.format") +
  scale_fill_manual(values  = inhibition_colors) +
  scale_color_manual(values = pond_colors) +
  facet_wrap(~ Pond.x) +
  labs(title = "Total BGC Count by Inhibition Category",
       x = "Inhibition Category", y = "Total BGC Regions",
       fill = "Inhibition Category", color = "Pond") +
  theme_bw() +
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 25, hjust = 1))

# ── 9. Plot 4: Scatter — BGC count vs. mean inhibition ───────────────────────
ggplot(isolate_data, aes(x = mean, y = n_bgc_total,
                         color = Pond.x, shape = Inhibition_Flag)) +
  geom_errorbarh(aes(xmin = mean - se, xmax = mean + se),
                 height = 0.2, alpha = 0.4) +
  geom_point(size = 3.5, alpha = 0.9) +
  geom_smooth(aes(group = Pond.x), method = "lm", se = TRUE,
              linetype = "dashed", alpha = 0.15) +
  geom_text_repel(aes(label = Isolate), size = 2.8, max.overlaps = 15) +
  scale_color_manual(values = pond_colors) +
  labs(title  = "BGC Count vs. Mean Chytrid Inhibition",
       x      = "Mean Inhibition (%)", y = "Total BGC Regions",
       color  = "Pond.x", shape = "Inhibition Category") +
  theme_bw() +
  theme(legend.position = "bottom")

# ── 10. Plot 5: Scatter — # phages vs. mean inhibition ───────────────────────
ggplot(isolate_data, aes(x = mean, y = n_phages,
                         color = Pond.x, shape = Inhibition_Flag)) +
  geom_errorbarh(aes(xmin = mean - se, xmax = mean + se),
                 height = 0.2, alpha = 0.4) +
  geom_point(size = 3.5, alpha = 0.9) +
  geom_smooth(aes(group = Pond.x), method = "lm", se = TRUE,
              linetype = "dashed", alpha = 0.15) +
  geom_text_repel(aes(label = Isolate), size = 2.8, max.overlaps = 15) +
  scale_color_manual(values = pond_colors) +
  labs(title  = "Prophage Count vs. Mean Chytrid Inhibition",
       x      = "Mean Inhibition (%)", y = "Number of Unique Phage Hits",
       color  = "Pond", shape = "Inhibition Category") +
  theme_bw() +
  theme(legend.position = "bottom")

# ── 11. Plot 6: Scatter — # phages vs. total BGCs ────────────────────────────
ggplot(isolate_data, aes(x = n_bgc_total, y = n_phages,
                         color = Pond.x, shape = Inhibition_Flag)) +
  geom_point(size = 3.5, alpha = 0.9) +
  geom_smooth(aes(group = Pond.x), method = "lm", se = TRUE,
              linetype = "dashed", alpha = 0.15) +
  geom_text_repel(aes(label = Isolate), size = 2.8, max.overlaps = 15) +
  scale_color_manual(values = pond_colors) +
  labs(title  = "Prophage Count vs. Total BGC Regions",
       x      = "Total BGC Regions", y = "Number of Unique Phage Hits",
       color  = "Pond", shape = "Inhibition Category") +
  theme_bw() +
  theme(legend.position = "bottom")

# ── 12. Plot 7: Heatmap — BGC type presence per isolate ──────────────────────
# Useful to spot shared/unique biosynthetic capacity across isolates

bgc_matrix <- bgc %>%
  dplyr::distinct(Isolate, Type) %>%
  mutate(present = 1) %>%
  pivot_wider(names_from = Type, values_from = present, values_fill = 0) %>%
  pivot_longer(-Isolate, names_to = "BGC_Type", values_to = "present") %>%
  left_join(inhibition %>% dplyr::select(Isolate, Pond, Inhibition_Flag,  mean),
            by = "Isolate") %>%
  mutate(Isolate = fct_reorder(Isolate, mean))   # order by inhibition strength

ggplot(bgc_matrix, aes(x = BGC_Type, y = Isolate, fill = factor(present))) +
  geom_tile(color = "white", linewidth = 0.4) +
  facet_grid(Pond ~ ., scales = "free_y", space = "free_y") +
  scale_fill_manual(values = c("0" = "grey92", "1" = "#009E73"),
                    labels = c("Absent", "Present")) +
  labs(title = "BGC Type Presence/Absence per Isolate\n(ordered by mean inhibition)",
       x = "BGC Type", y = "Isolate", fill = "") +
  theme_bw() +
  theme(axis.text.x  = element_text(angle = 45, hjust = 1, size = 8),
        legend.position = "bottom")

# ── 13. Statistical tests ──────────────────────────────────────────────────────

# 13a. Kruskal-Wallis: phage count ~ pond
cat("\n── Kruskal-Wallis: Phage Count ~ Pond ──\n")
kruskal.test(n_phages ~ Pond, data = isolate_data)

# 13b. Kruskal-Wallis: phage count ~ inhibition category
cat("\n── Kruskal-Wallis: Phage Count ~ Inhibition Category ──\n")
kruskal.test(n_phages ~ Inhibition_Flag, data = isolate_data)

# 13c. Kruskal-Wallis: BGC count ~ inhibition category
cat("\n── Kruskal-Wallis: BGC Count ~ Inhibition Category ──\n")
kruskal.test(n_bgc_total ~ Inhibition_Flag, data = isolate_data)

# 13d. Kruskal-Wallis: BGC count ~ pond
cat("\n── Kruskal-Wallis: BGC Count ~ Pond ──\n")
kruskal.test(n_bgc_total ~ Pond, data = isolate_data)

# 13e. Dunn post-hoc tests
cat("\n── Dunn Test (Bonferroni): Phage Count ~ Pond ──\n")
dunn.test(isolate_data$n_phages, isolate_data$Pond,
          method = "bonferroni", kw = TRUE, label = TRUE)

cat("\n── Dunn Test (Bonferroni): Phage Count ~ Inhibition Category ──\n")
dunn.test(isolate_data$n_phages, isolate_data$Inhibition_Flag,
          method = "bonferroni", kw = TRUE, label = TRUE)

cat("\n── Dunn Test (Bonferroni): BGC Count ~ Inhibition Category ──\n")
dunn.test(isolate_data$n_bgc_total, isolate_data$Inhibition_Flag,
          method = "bonferroni", kw = TRUE, label = TRUE)

# 13f. Spearman correlations
cat("\n── Spearman: n_phages vs. mean inhibition ──\n")
cor.test(isolate_data$n_phages, isolate_data$mean, method = "spearman")

cat("\n── Spearman: n_bgc_total vs. mean inhibition ──\n")
cor.test(isolate_data$n_bgc_total, isolate_data$mean, method = "spearman")

cat("\n── Spearman: n_phages vs. n_bgc_total ──\n")
cor.test(isolate_data$n_phages, isolate_data$n_bgc_total, method = "spearman")