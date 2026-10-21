library(tidyverse)
library(dunn.test)
library(ggpubr)
library(ggrepel)

#Colour palettes
pond_colors       <- c("Sixty Lake"   = "#0072B2",
                       "Conness Pond" = "#E69F00")

inhibition_colors <- c("Weakly Inhibitory"   = "#56B4E9",
                       "Strongly Inhibitory" = "#E69F00")

# 7-color palette for phage families
family_colors <- c("#E69F00", "#56B4E9", "#009E73",
                   "#CC79A7", "#0072B2", "#D55E00", "#F0E442", "#771155", "#CC99BB", "#114477")

#Load data
blast    <- read_tsv("~/Documents/GitHub/Pseudo_fluor/viral_blast_hits.txt")
antismash <- read_tsv("~/Documents/GitHub/Pseudo_fluor/antismash.txt")
inhibition <- read_tsv("~/Documents/GitHub/Pseudo_fluor/isolate_inhibi_category.txt") %>%
  dplyr::rename(Isolate = Group)

#Parse BGC types
# Some BGC regions have compound types (e.g. "NRPS,azole-containing-RiPP")
# Split these so each BGC type gets its own row for counting

bgc <- antismash %>%
  mutate(Type = str_replace_all(Type, "\\s", "")) %>%
  separate_rows(Type, sep = ",") %>%           # one row per BGC type
  left_join(inhibition, by = "Isolate")        # attach inhibition metadata

#Parse prophage taxonomy
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
  ) 

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

#BGC type counts per isolate
bgc_type_counts <- bgc %>%
  dplyr::count(Isolate, Type) %>%
  left_join(inhibition %>% dplyr::select(Isolate, Pond, Inhibition_Flag),
            by = "Isolate")

#Plot 1: Stacked bar — prophage family per isolate
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

# ──Plot 2: Stacked bar — BGC types per isolate
pal<-c("#771155", "#CC99BB", "#114477", "#4477AA", "#117777", "#44AAAA", "#77CCCC", "#117744", "#44AA77", "#88CCAA", "#777711", "#AAAA44", "#DDDD77", "#774411", "#AA7744", "#DDAA77", "#771122", "#AA4455", "#DD7788","#41AB5D", "#252525", "#525252", "#737373", "#969696")
bgc_type_counts %>%
  mutate(Isolate = fct_reorder(Isolate, n, .fun = sum, .desc = TRUE)) %>%
  ggplot(aes(x = Isolate, y = n, fill = Type)) +
  geom_col() +
  scale_fill_manual(values = pal)+
  facet_wrap(~ Pond, scales = "free_x") +
  labs(x = "Isolate", y = "Number of BGC Regions", fill = "BGC Type") +
  theme_bw() +
  theme(axis.text.x  = element_text(angle = 45, hjust = 1),
        legend.position = "bottom",
        legend.text  = element_text(size = 7)) +
  guides(fill = guide_legend(ncol = 3))

# Plot 3: Box/violin — BGC count by inhibition category
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
#no significant effect of weak vs. strong on the number of BSG

#Plot 4: Scatter — BGC count vs. mean inhibition
ggplot(isolate_data, aes(x = mean, y = n_bgc_total,
                         color = Pond.x, shape = Inhibition_Flag)) +
  geom_errorbar(aes(xmin = mean - se, xmax = mean + se),
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

#Plot 5: Scatter — # phages vs. mean inhibition
library(ggpubr)
ggplot(isolate_data, aes(x = mean, y = n_phages,
                         color = Pond.x, shape = Inhibition_Flag)) +
  #stat_cor(method='spearman')+
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


# Plot 6: Scatter — # phages vs. total BGCs
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

#Stat2

#Kruskal-Wallis: phage count ~ pond
kruskal.test(n_phages ~ Pond.x, data = isolate_data)
#Kruskal-Wallis chi-squared = 0.29727, df = 1, p-value = 0.5856

#Kruskal-Wallis: phage count ~ inhibition category
kruskal.test(n_phages ~ Inhibition_Flag, data = isolate_data)
#Kruskal-Wallis chi-squared = 1.7025, df = 1, p-value = 0.19

#Kruskal-Wallis: BGC count ~ inhibition category
kruskal.test(n_bgc_total ~ Inhibition_Flag, data = isolate_data)
#Kruskal-Wallis chi-squared = 1.749, df = 1, p-value = 0.186

# Kruskal-Wallis: BGC count ~ pond
kruskal.test(n_bgc_total ~ Pond.x, data = isolate_data)
#Kruskal-Wallis chi-squared = 1.6626, df = 1, p-value = 0.1973

#Dunn post-hoc tests
dunn.test(isolate_data$n_phages, isolate_data$Pond.x,
          method = "bh", kw = TRUE, label = TRUE)

dunn.test(isolate_data$n_phages, isolate_data$Inhibition_Flag,
          method = "bh", kw = TRUE, label = TRUE)

dunn.test(isolate_data$n_bgc_total, isolate_data$Inhibition_Flag,
          method = "bh", kw = TRUE, label = TRUE)

#Spearman correlations
cor.test(isolate_data$n_phages, isolate_data$mean, method = "spearman")

cor.test(isolate_data$n_bgc_total, isolate_data$mean, method = "spearman")

cor.test(isolate_data$n_phages, isolate_data$n_bgc_total, method = "spearman")
#S = 323.49, p-value = 0.002541
