##PF Pangenome viral data
#fixing up data, don't need to run after teh first time
library(tidyverse)

#read in files
files <- list.files(path = "~/Documents/GitHub/Pseudo_fluor/viral_blast_output/", 
                    pattern = "*.txt", 
                    full.names = TRUE)
final_df <- files %>%
  set_names() %>% 
  map_df(~read_tsv(.x, col_names = FALSE, col_types = cols()), .id = "source_file") %>%
  mutate(source_file = basename(source_file))

#give the file reasonable names
blast_colnames <- c(
  "Isolate",
  "query_id",        # V1: Query Seq-id
  "subject_id",      # V2: Subject Seq-id
  "perc_identity",   # V3: Percentage of identical matches
  "align_len",       # V4: Alignment length
  "mismatches",      # V5: Number of mismatches
  "gap_opens",       # V6: Number of gap openings
  "q_start",         # V7: Start of alignment in query
  "q_end",           # V8: End of alignment in query
  "s_start",         # V9: Start of alignment in subject
  "s_end",           # V10: End of alignment in subject
  "evalue",          # V11: Expect value (E-value)
  "bit_score"        # V12: Bit score
)
names(final_df)<-blast_colnames

#fix up isolate names
final_df$Isolate<-gsub('.fasta_final-viral-combined.fa_blast_out.txt', '', final_df$Isolate)

#add taxonomy of each
unique_viral_ass<-read.delim('~/Documents/GitHub/Pseudo_fluor/unique_viral_acession.txt', header=T)
finaldf2<-merge(final_df, unique_viral_ass, by.x='subject_id', by.y='Ascession', all.x=T)

#save to file
write.table(finaldf2, '~/Documents/GitHub/Pseudo_fluor/viral_blast_hits.txt', sep='\t', quote=F, row.names=F)

# ============================================================================
# Prophage Taxonomy & Chytrid Inhibition Analysis
# ============================================================================

library(tidyverse)
library(dunn.test)
library(ggpubr)   
library(ggrepel)  

# ── 0. Colour palettes ────────────────────────────────────────────────────────
pond_colors   <- c("Sixty Lake" = "#0072B2", "Conness Pond" = "#E69F00")
family_colors <- c(
  "#E69F00",  # orange
  "#56B4E9",  # sky blue
  "#009E73",  # green
  "#CC79A7",  # pink
  "#0072B2",  # blue
  "#D55E00",  # vermillion
  "#F0E442"   # yellow
)
inhibition_colors <- c("Weakly Inhibitory"    = "#56B4E9",
                       "Strongly Inhibitory"  = "#E69F00")

# ── 1. Load data ──────────────────────────────────────────────────────────────
blast <- read_tsv("~/Documents/GitHub/Pseudo_fluor/viral_blast_hits.txt")

# ── 2. Parse taxonomy ranks from FullTaxonomy column ─────────────────────────
# Format: "Viruses; Realm; Kingdom; Phylum; Class; Order; Family; Genus; Species"
# Number of levels varies; safe_nth() handles missing levels gracefully.

safe_nth <- function(lst, i) {
  map_chr(lst, ~ if (length(.x) >= i) .x[[i]] else NA_character_)
}

blast_tax <- blast %>%
  mutate(
    FullTaxonomy = str_replace_all(FullTaxonomy, "\\s+", " ") %>% str_trim(),
    tax_levels   = str_split(FullTaxonomy, ";\\s*")
  ) %>%
  mutate(
    realm        = safe_nth(tax_levels, 2),
    kingdom      = safe_nth(tax_levels, 3),
    phylum       = safe_nth(tax_levels, 4),
    class        = safe_nth(tax_levels, 5),
    order        = safe_nth(tax_levels, 6),
    family       = safe_nth(tax_levels, 7),
    genus        = safe_nth(tax_levels, 8),
    species      = safe_nth(tax_levels, 9),
    # Finest available taxonomic level (strips trailing period)
    lowest_taxon = map_chr(tax_levels, ~ str_remove(tail(.x, 1), "\\.$"))
  ) %>%
  select(-tax_levels)

# ── 3. Deduplicate: best BLAST hit per isolate × subject phage ────────────────
# Multiple contigs/HSPs can hit the same reference phage — keep highest bit_score.

blast_dedup <- blast_tax %>%
  group_by(Pond, Isolate, subject_id, Hit, family, genus, lowest_taxon,
           mean, sd, se, Inhibition_Flag) %>%
  dplyr::slice_max(bit_score, n = 1, with_ties = FALSE) %>%
  ungroup()

# ── 4. Summary tables ─────────────────────────────────────────────────────────

# 4a. Unique phage hits per isolate × family/genus
isolate_summary <- blast_dedup %>%
  group_by(Pond, Isolate, family, genus) %>%
  dplyr::summarise(
    n_hits     = dplyr::n(),
    phage_hits = paste(unique(Hit), collapse = " | "),
    .groups    = "drop"
  ) %>%
  dplyr::arrange(Pond, Isolate, family)

print(isolate_summary)

# 4b. Pond-level: how many isolates carry each taxon
pond_summary <- blast_dedup %>%
  group_by(Pond, family, genus) %>%
  dplyr::summarise(
    n_isolates   = dplyr::n_distinct(Isolate),
    isolates     = paste(sort(unique(Isolate)), collapse = ", "),
    n_phage_hits = dplyr::n_distinct(subject_id),
    .groups      = "drop"
  ) %>%
  dplyr::arrange(Pond, desc(n_isolates))

print(pond_summary)

# 4c. Wide presence/absence matrix: isolates × family
pa_matrix <- blast_dedup %>%
  dplyr::distinct(Pond, Isolate, family) %>%
  mutate(present = 1) %>%
  pivot_wider(names_from = family, values_from = present, values_fill = 0) %>%
  dplyr::arrange(Pond, Isolate)

print(pa_matrix)

# 4d. Per-isolate phage count with inhibition metadata
phage_counts <- blast_dedup %>%
  dplyr::distinct(Pond, Isolate, subject_id, mean, se, Inhibition_Flag) %>%
  dplyr::count(Pond, Isolate, mean, se, Inhibition_Flag, name = "n_phages") %>%
  mutate(Inhibition_Flag = factor(Inhibition_Flag,
                                  levels = c("Weakly Inhibitory",
                                             "Strongly Inhibitory")))

# ── 5. Plot 1: Stacked bar — phage family breakdown per isolate by pond ───────
blast_dedup %>%
  dplyr::distinct(Pond, Isolate, subject_id, family) %>%
  ggplot(aes(x = Isolate, fill = family)) +
  geom_bar() +
  facet_wrap(~ Pond, scales = "free_x") +
  scale_fill_manual(values = family_colors) +
  labs(
    title = "Unique Prophage Hits per Isolate",
    x     = "Isolate",
    y     = "Number of Unique Phage Hits",
    fill  = "Family"
  ) +
  theme_bw() +
  theme(
    axis.text.x     = element_text(angle = 45, hjust = 1),
    legend.position = "bottom"
  )

# ── 6. Plot 2: Box & whisker — # phages per pond ──────────────────────────────
ggplot(phage_counts, aes(x = Pond, y = n_phages, fill = Pond)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(aes(color = Pond), width = 0.15, size = 2.5, alpha = 0.9) +
  geom_text_repel(aes(label = Isolate), size = 2.8, max.overlaps = 15) +
  scale_fill_manual(values  = pond_colors) +
  scale_color_manual(values = pond_colors) +
  labs(
    title = "Number of Unique Prophage Hits per Isolate by Pond",
    x     = "Pond",
    y     = "Number of Unique Phage Hits"
  ) +
  theme_bw() +
  theme(legend.position = "none")

# ── 7. Plot 3: Scatter — # phages vs. mean inhibition ────────────────────────
ggplot(phage_counts, aes(x = mean, y = n_phages,
                         color = Pond, shape = Inhibition_Flag)) +
  geom_errorbarh(aes(xmin = mean - se, xmax = mean + se),
                 height = 0.2, alpha = 0.4) +
  geom_point(size = 3.5, alpha = 0.9) +
  geom_smooth(aes(group = Pond), method = "lm", se = TRUE,
              linetype = "dashed", alpha = 0.15) +
  geom_text_repel(aes(label = Isolate), size = 2.8, max.overlaps = 15) +
  scale_color_manual(values = pond_colors) +
  labs(
    title  = "Prophage Count vs. Mean Inhibition of Chytrid Fungi",
    x      = "Mean Inhibition (%)",
    y      = "Number of Unique Phage Hits",
    color  = "Pond",
    shape  = "Inhibition Category"
  ) +
  theme_bw() +
  theme(legend.position = "bottom")

# ── 8. Plot 4: Violin + box — # phages by inhibition category ────────────────
comparisons_inhib <- list(c("Weakly Inhibitory", "Strongly Inhibitory"))

ggplot(phage_counts, aes(x = Inhibition_Flag, y = n_phages, fill = Inhibition_Flag)) +
  #geom_violin(alpha = 0.4, trim = FALSE) +
  geom_boxplot(width = 0.2, outlier.shape = NA, alpha = 0.7) +
  facet_wrap(~ Pond) +
  #geom_jitter(aes(color = Pond), width = 0.1, size = 2.5, alpha = 0.9) +
  #geom_text_repel(aes(label = Isolate), size = 2.8, max.overlaps = 15) +
  #stat_compare_means(comparisons     = comparisons_inhib,
                    # method          = "wilcox.test",
                     #label           = "p.format") +
  scale_fill_manual(values  = inhibition_colors) +
  scale_color_manual(values = pond_colors) +
  
  labs(
    title = "Prophage Count by Inhibition Category",
    x     = "Inhibition Category",
    y     = "Number of Phage Hits",
    fill  = "Inhibition Category",
    color = "Pond"
  ) +
  theme_bw() +
  theme(
    legend.position = "bottom",
    axis.text.x     = element_text(angle = 25, hjust = 1)
  )

#-------------------------------
#Stats
#-------------------------------

#Kruskal-Wallis: does phage count differ between ponds?
kruskal.test(n_phages ~ Pond, data = phage_counts)
#Kruskal-Wallis chi-squared = 0.29727, df = 1, p-value = 0.5856

#Kruskal-Wallis: does phage count differ between inhibition categories?
kruskal.test(n_phages ~ Inhibition_Flag, data = phage_counts)
#Kruskal-Wallis chi-squared = 1.7025, df = 1, p-value = 0.192

#Spearman correlation: phage count vs. mean inhibition
cor.test(phage_counts$n_phages, phage_counts$mean, method = "spearman")
#S = 1289, p-value = 0.1807, rho=-0.33
