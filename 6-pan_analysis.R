# =============================================================================
# Pseudomonas Pangenome Analysis Pipeline
# =============================================================================
#load packages
library(tidyverse)
library(vegan)
library(ggplot2)
library(pheatmap)
library(ggrepel)
library(broom)
library(RColorBrewer)
library(patchwork)
library(scales)

rtab_file   <- "~/Documents/GitHub/Pseudo_fluor/panaroo_results/gene_presence_absence2.Rtab" 
meta_file   <- "~/Documents/GitHub/Pseudo_fluor/pf_genome_table.txt"         
output_dir  <- "~/Documents/GitHub/Pseudo_fluor/figures"

# Colour palette used throughout
pal_inhibitory <- c("Strongly Inhibitory" = "#D7263D", "Weakly Inhibitory" = "#4A90D9")
pal_pond       <- c("Sixty Lake" = "#2E8B57", "Conness Pond" = "#E07B39")

save_pdf <- function(plot_obj, filename, w = 8, h = 6) {
  ggsave(
    filename  = file.path(output_dir, filename),
    plot      = plot_obj,
    device    = "pdf",
    width     = w,
    height    = h,
    useDingbats = FALSE
  )
  message("Saved: ", file.path(output_dir, filename))
}

##load and cleanup data
# Load presence/absence matrix (genes x isolates)
pa_raw <- read.delim(rtab_file, header = TRUE, check.names = FALSE)

# Extract isolate names by stripping the long prefix/suffix from column names
clean_names <- colnames(pa_raw)[-1] %>%
  str_replace("^2023-08-15_256samples_", "") %>%
  str_replace("\\.bam_scaffolds\\.fasta_PROKKA_\\d+$", "") %>%
  str_replace("_scaffolds\\.fasta_PROKKA_\\d+$", "") %>%
  str_replace("\\.fasta_PROKKA_\\d+$", "")

# Build clean matrix: rows = genes, cols = isolates
pa_mat <- pa_raw[, -1] %>% as.matrix()
rownames(pa_mat) <- pa_raw$Gene
colnames(pa_mat) <- clean_names

# Remove the empty trailing row if present
pa_mat <- pa_mat[rowSums(pa_mat) > 0, ]

message(sprintf("Pangenome matrix: %d genes x %d isolates", nrow(pa_mat), ncol(pa_mat)))
#Pangenome matrix: 7778 genes x 18 isolates

# Load metadata
meta <- read.delim(meta_file, header = TRUE, stringsAsFactors = FALSE)

# Ensure isolate order matches matrix columns
meta <- meta[match(colnames(pa_mat), meta$Isolate), ]
stopifnot(all(meta$Isolate == colnames(pa_mat)))

#transpose pan matrix
pa_iso <- t(pa_mat)

# =============================================================================
# Core / Accessory / Unique gene breakdown
# =============================================================================
pan_breakdown <- function(mat) {
  mat  <- mat[rowSums(mat) > 0, , drop = FALSE]
  freq <- rowSums(mat); n <- ncol(mat)
  tibble(
    Category = factor(c("Core", "Accessory", "Unique/Singleton"),
                      levels = c("Core", "Accessory", "Unique/Singleton")),
    Count    = c(sum(freq == n), sum(freq > 1 & freq < n), sum(freq == 1))
  ) %>%
    mutate(Pct = Count / sum(Count),
           Label = sprintf("%s\n%d (%.1f%%)", Category, Count, 100 * Pct)) %>%
    structure(n_isolates = n, n_genes = nrow(mat))
}

pan_pie <- function(df, title) {
  df <- df %>% mutate(
    # small slices (<12%) get their label just outside the pie so they don't collide
    label_x   = ifelse(Pct < 0.12, 1.62, 1),
    label_col = ifelse(Category == "Core", "white", "black"),
    # slice mid-points: geom_col stacks the last factor level at the bottom
    label_y   = rev(cumsum(rev(Count)) - rev(Count) / 2)
  )
  ggplot(df, aes(x = 1, y = Count, fill = Category)) +
    geom_col(width = 1, colour = "white", linewidth = 0.5) +
    geom_text(aes(x = label_x, y = label_y, label = Label, colour = label_col),
              size = 3.2, fontface = "bold", lineheight = 0.9) +
    scale_colour_identity() +
    coord_polar(theta = "y", clip = "off") +
    xlim(0.5, 1.62) +
    scale_fill_manual(values = c("Core" = "#2166AC", "Accessory" = "#92C5DE",
                                 "Unique/Singleton" = "#F4A582")) +
    labs(title = title,
         subtitle = sprintf("%d genes | %d isolates",
                            attr(df, "n_genes"), attr(df, "n_isolates")),
         fill = NULL) +
    theme_void(base_size = 12) +
    theme(plot.title = element_text(face = "bold", hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5, colour = "grey40"),
          legend.position = "none")
}

pond_isolates <- split(meta$Isolate, meta$Pond)          # meta rows match pa_mat columns

breakdowns <- c(list("All isolates" = pan_breakdown(pa_mat)),
                lapply(pond_isolates[c("Sixty Lake", "Conness Pond")],
                       function(iso) pan_breakdown(pa_mat[, iso, drop = FALSE])))

for (nm in names(breakdowns)) {
  b <- breakdowns[[nm]]
  cat(sprintf("%-13s %4d genes | core %4d (%.1f%%) | accessory %4d (%.1f%%) | unique %4d (%.1f%%)\n",
              nm, attr(b, "n_genes"), b$Count[1], 100*b$Pct[1],
              b$Count[2], 100*b$Pct[2], b$Count[3], 100*b$Pct[3]))
}

#make figure and save
p1 <- wrap_plots(Map(pan_pie, breakdowns, names(breakdowns)), nrow = 1) +
  plot_annotation(tag_levels = "A",
                  theme = theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 14)))

save_pdf(p1, "01_pangenome_breakdown_pie.pdf", w = 15, h = 5.5)

# ── Genes that change pangenome status between ponds (Fisher's exact test) ────
# Builds on pan_breakdown(): each gene is classed within each pond as
#   Core (all 9 isolates) | Accessory (2-8) | Unique (1) | Absent (0).
# A gene "switches" if its class differs between ponds. Each gene's presence
# (Sixty Lake vs Conness Pond) is tested with Fisher's exact test.
# Genes with an identical 18-isolate presence pattern give identical tests, so
# they are collapsed into blocks and BH is applied across blocks, not genes.
# Requires: pa_mat, meta (from script 6), tidyverse, patchwork, save_pdf().

# MASS (script 2) and plyr (scripts 2, 3, 7) mask dplyr verbs when the scripts
# share one R session; use dplyr's versions explicitly inside this section.
select    <- dplyr::select
filter    <- dplyr::filter
mutate    <- dplyr::mutate
summarise <- dplyr::summarise
count     <- dplyr::count
rename    <- dplyr::rename
arrange   <- dplyr::arrange
n         <- dplyr::n

pond_a <- "Sixty Lake"; pond_b <- "Conness Pond"
iso_a  <- meta$Isolate[meta$Pond == pond_a]
iso_b  <- meta$Isolate[meta$Pond == pond_b]

classify <- function(n_present, n_iso) {
  factor(case_when(n_present == n_iso ~ "Core",
                   n_present >= 2     ~ "Accessory",
                   n_present == 1     ~ "Unique",
                   TRUE               ~ "Absent"),
         levels = c("Core", "Accessory", "Unique", "Absent"))
}

status_df <- tibble(
  Gene    = rownames(pa_mat),
  pattern = apply(pa_mat, 1, paste, collapse = ""),
  n_a     = rowSums(pa_mat[, iso_a, drop = FALSE]),
  n_b     = rowSums(pa_mat[, iso_b, drop = FALSE])
) %>%
  mutate(status_a = classify(n_a, length(iso_a)),
         status_b = classify(n_b, length(iso_b)),
         switches = status_a != status_b)

# Fisher's exact test once per presence-pattern block
block_tests <- status_df %>%
  filter(!(n_a == length(iso_a) & n_b == length(iso_b))) %>%   # core everywhere: nothing to test
  group_by(pattern, n_a, n_b) %>%
  summarise(n_genes = n(), genes = list(Gene), .groups = "drop") %>%
  mutate(p_value = map2_dbl(n_a, n_b, ~ fisher.test(matrix(
    c(.x, length(iso_a) - .x, .y, length(iso_b) - .y), 2))$p.value),
    p_adj   = p.adjust(p_value, method = "BH"),
    log2OR  = log2(((n_a + .5) * (length(iso_b) - n_b + .5)) /
                     ((length(iso_a) - n_a + .5) * (n_b + .5))))   # Haldane +0.5

gene_results <- status_df %>%
  left_join(block_tests %>% select(pattern, block_n_genes = n_genes, p_value, p_adj, log2OR),
            by = "pattern") %>%
  mutate(Direction = case_when(is.na(p_value) ~ "Core in both",
                               log2OR > 0 ~ paste("More common in", pond_a),
                               log2OR < 0 ~ paste("More common in", pond_b),
                               TRUE ~ "No difference"),
         Transition = paste0(status_a, " (", pond_a, ") -> ", status_b, " (", pond_b, ")"),
         Sig = !is.na(p_adj) & p_adj < 0.05) %>%
  arrange(p_value, desc(block_n_genes))

# ---- summaries ----
transition_tab <- gene_results %>% count(status_a, status_b, name = "n_genes", .drop = FALSE) %>%
  complete(status_a, status_b, fill = list(n_genes = 0L))   # show empty class combinations as 0
cat("\nStatus cross-tab (rows =", pond_a, ", cols =", pond_b, "):\n")
print(xtabs(n_genes ~ status_a + status_b, transition_tab))
cat(sprintf("\nGenes switching status: %d of %d\n", sum(gene_results$switches), nrow(gene_results)))
cat(sprintf("Blocks tested: %d | significant blocks (BH < 0.05): %d containing %d genes\n",
            nrow(block_tests), sum(block_tests$p_adj < 0.05),
            sum(block_tests$n_genes[block_tests$p_adj < 0.05])))
sig_genes <- gene_results %>% filter(Sig)
print(sig_genes %>% count(Transition, Direction, name = "n_genes") %>% arrange(desc(n_genes)))

write_csv(gene_results %>% select(-pattern),
          file.path(output_dir, "10_pond_status_fisher_all_genes.csv"))
write_csv(sig_genes %>% select(-pattern),
          file.path(output_dir, "10_pond_status_fisher_significant.csv"))

# ---- figure ----
# A: how genes move between pangenome classes across ponds
pA <- ggplot(transition_tab, aes(status_b, fct_rev(status_a), fill = n_genes)) +
  geom_tile(colour = "white", linewidth = 0.8) +
  geom_text(aes(label = n_genes,
                colour = case_when(n_genes > 1500 ~ "white", n_genes == 0 ~ "grey60", TRUE ~ "black")),
            size = 3.6) +
  scale_colour_identity() +
  scale_fill_gradient(low = "#F0F4FA", high = "#2166AC", trans = "log10",
                      na.value = "grey96", name = "Genes") +   # 0 -> NA on log scale -> grey
  labs(title = "Pangenome class of each gene in the two ponds",
       x = paste("Class in", pond_b), y = paste("Class in", pond_a)) +
  coord_equal() +
  theme_minimal(base_size = 11) +
  theme(panel.grid = element_blank())

# B: presence counts per pond, one point per block, significant blocks highlighted
plot_blocks <- block_tests %>%
  mutate(Sig = p_adj < 0.05,
         ex  = map_chr(genes, ~ { nm <- .x[!grepl("^group_", .x)]
         paste(head(if (length(nm)) nm else .x, 2), collapse = ", ") }))
lab_blocks <- plot_blocks %>% filter(Sig) %>% arrange(p_value) %>% slice_head(n = 6) %>%
  mutate(label = sprintf("%d gene%s (%s)", n_genes, ifelse(n_genes == 1, "", "s"), ex))

pB <- ggplot(plot_blocks, aes(n_a, n_b)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dotted", colour = "grey60") +
  geom_point(aes(size = n_genes, colour = Sig), alpha = 0.75,
             position = position_jitter(width = 0.15, height = 0.15, seed = 1)) +
  ggrepel::geom_text_repel(data = lab_blocks, aes(label = label), size = 2.7,
                           min.segment.length = 0, box.padding = 0.5, seed = 1,
                           max.overlaps = Inf) +
  scale_colour_manual(values = c(`TRUE` = "#D55E00", `FALSE` = "grey65"),
                      labels = c(`TRUE` = "BH < 0.05", `FALSE` = "n.s."), name = NULL) +
  scale_size_area(max_size = 8, breaks = c(1, 10, 100, 500), name = "Genes in block") +
  scale_x_continuous(breaks = 0:length(iso_a), limits = c(-0.5, length(iso_a) + 0.5)) +
  scale_y_continuous(breaks = 0:length(iso_b), limits = c(-0.5, length(iso_b) + 0.5)) +
  labs(title = "Gene presence per pond (Fisher's exact test per block)",
       subtitle = sprintf("%d blocks tested; %d significant (%d genes); off-diagonal = differs between ponds",
                          nrow(block_tests), sum(plot_blocks$Sig), sum(plot_blocks$n_genes[plot_blocks$Sig])),
       x = paste("Isolates carrying gene,", pond_a, "(of", length(iso_a), ")"),
       y = paste("Isolates carrying gene,", pond_b, "(of", length(iso_b), ")")) +
  coord_equal() +
  theme_classic(base_size = 11) +
  theme(plot.subtitle = element_text(size = 8.5, colour = "grey30"))

p10 <- (pA | pB) + plot_layout(widths = c(1, 1.25)) + plot_annotation(tag_levels = "A")
save_pdf(p10, "10_pond_status_fisher.pdf", w = 14, h = 6.5)

###this passes to the R script cog_enrichment_figure.R to make a more informative plot

# =============================================================================
#Gene count per isolate
# =============================================================================
gene_count_df <- tibble(
  Isolate          = colnames(pa_mat),
  GeneCount        = colSums(pa_mat),
  InhibitoryType   = meta$InhibitoryType,
  Pond             = meta$Pond
)

p2 <- ggplot(gene_count_df,
             aes(x = InhibitoryType, y = GeneCount, fill = InhibitoryType)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA, width = 0.5) +
  geom_jitter(aes(colour = Pond), width = 0.15, size = 3, stroke = 0.3) +
  geom_text_repel(aes(label = Isolate), size = 2.8, max.overlaps = 20) +
  scale_fill_manual(values = pal_inhibitory)  +
  scale_colour_manual(values = pal_pond) +
  labs(x = NULL, y = "Number of genes",
       fill = "Inhibitory Type", colour = "Pond") +
  theme_classic(base_size = 13) +
  theme(plot.title = element_text(face = "bold"),
        legend.position = "right")

save_pdf(p2, "02_gene_count_boxplot.pdf", w = 8, h = 6)

# Wilcoxon test
wtest <- wilcox.test(GeneCount ~ InhibitoryType, data = gene_count_df)
cat(sprintf("Gene count Wilcoxon test: W = %.1f, p = %.4f\n", wtest$statistic, wtest$p.value))
#Gene count Wilcoxon test: W = 26.0, p = 0.2761

wtest2 <- wilcox.test(GeneCount ~ Pond, data = gene_count_df)
cat(sprintf("Gene count Wilcoxon test: W = %.1f, p = %.4f\n", wtest2$statistic, wtest2$p.value))
#Gene count Wilcoxon test: W = 81.0, p = 0.0004

#interaction with Scheirer ray hare
library(rcompanion)
scheirerRayHare(GeneCount ~  InhibitoryType + Pond, data = gene_count_df)
#DV:  GeneCount 
#Observations:  18 
#D:  0.995872 
#MS total:  28.5 

#                    Df Sum Sq       H p.value
#InhibitoryType       1  15.16  0.5341 0.46490
#Pond                 1 343.13 12.0896 0.00051
#InhibitoryType:Pond  1  31.84  1.1219 0.28951
#Residuals           14  71.00 

#look at just Conness, plot clearly shows lower genes in inhibitory categories
gene_count_split<-split(gene_count_df, gene_count_df$Pond)

wtest3 <- wilcox.test(GeneCount ~ InhibitoryType, data = gene_count_split$`Conness Pond`)
cat(sprintf("Gene count Wilcoxon test: W = %.1f, p = %.4f\n", wtest3$statistic, wtest3$p.value))
#Gene count Wilcoxon test: W = 0.0, p = 0.0195

wtest4 <- wilcox.test(GeneCount ~ InhibitoryType, data = gene_count_split$`Sixty Lake`)
cat(sprintf("Gene count Wilcoxon test: W = %.1f, p = %.4f\n", wtest4$statistic, wtest4$p.value))
#Gene count Wilcoxon test: W = 11.0, p = 0.6949

# =============================================================================
#PCoA of pangenome
# =============================================================================
jac_dist <- vegdist(pa_iso, method = "jaccard", binary = TRUE)
pcoa_res  <- cmdscale(jac_dist, k = 4, eig = TRUE)

# Variance explained
var_exp <- pcoa_res$eig / sum(pcoa_res$eig[pcoa_res$eig > 0]) * 100

pcoa_df <- as_tibble(pcoa_res$points, .name_repair = "minimal") %>%
  setNames(paste0("PC", 1:4)) %>%
  mutate(
    Isolate        = meta$Isolate,
    InhibitoryType = meta$InhibitoryType,
    Pond           = meta$Pond,
    MeanInhib      = meta$MeanInhib
  )

# PCoA coloured by InhibitoryType, shaped by Pond
p3a <- ggplot(pcoa_df, aes(x = PC1, y = PC2,
                           colour = InhibitoryType, shape = Pond)) +
  geom_point(size = 4, stroke = 1) +
  geom_text_repel(aes(label = Isolate), size = 3, max.overlaps = 20) +
  scale_colour_manual(values = pal_inhibitory) +
  scale_shape_manual(values = c("Sixty Lake" = 16, "Conness Pond" = 17)) +
  labs(
    title    = "PCoA of Gene Presence/Absence (Jaccard Distance)",
    x        = sprintf("PC1 (%.1f%%)", var_exp[1]),
    y        = sprintf("PC2 (%.1f%%)", var_exp[2]),
    colour   = "Inhibitory Type",
    shape    = "Pond"
  ) +
  theme_classic(base_size = 13) +
  theme(plot.title = element_text(face = "bold"))

save_pdf(p3a, "03a_PCoA_inhibitory_type.pdf", w = 9, h = 7)

# PERMANOVA: does inhibitory type or pond explain variance?
set.seed(42)
perm_res <- adonis2(jac_dist ~ InhibitoryType + Pond,
                    data = meta, permutations = 999)
print(perm_res)
#         Df SumOfSqs      R2     F Pr(>F)    
#Model     2 0.103481 0.52523 8.297  0.001 ***
#Residual 15 0.093541 0.47477                 
#Total    17 0.197021 1.00000 

perm_res2 <- adonis2(jac_dist ~ InhibitoryType,
                    data = meta, permutations = 999)
print(perm_res2)
#         Df SumOfSqs      R2     F Pr(>F)   
#Model     1 0.061063 0.30993 7.186  0.008 **
#Residual 16 0.135959 0.69007                
#Total    17 0.197021 1.00000   

perm_res3 <- adonis2(jac_dist ~ Pond,
                     data = meta, permutations = 999)
print(perm_res3)
#Df SumOfSqs     R2      F Pr(>F)   
#Model     1 0.044586 0.2263 4.6798  0.005 **
#Residual 16 0.152435 0.7737                 
#Total    17 0.197021 1.0000

# =============================================================================
#Accessory genome heatmap
# =============================================================================
# Subset to accessory genes only
accessory_mat <- pa_mat[gene_freq > 1 & gene_freq < n_isolates, ]
message(sprintf("Accessory genes for heatmap: %d", nrow(accessory_mat)))
#Accessory genes for heatmap: 2146

# Annotation for columns (isolates)
col_annot <- data.frame(
  InhibitoryType = meta$InhibitoryType,
  Pond           = meta$Pond,
  row.names      = meta$Isolate
)

annot_colors <- list(
  InhibitoryType = pal_inhibitory,
  Pond           = pal_pond
)

# Use binary colour (white = absent, blue = present)
heatmap_colors <- colorRampPalette(c("white", "#2166AC"))(2)

pheatmap(
  mat               = accessory_mat,
  color             = heatmap_colors,
  annotation_col    = col_annot,
  annotation_colors = annot_colors,
  show_rownames     = nrow(accessory_mat) <= 200,  # hide row names if too many
  show_colnames     = TRUE,
  clustering_method = "ward.D2",
  fontsize_col      = 9,
  fontsize_row      = 5,
  main              = "Accessory Genome Presence/Absence",
  filename          = file.path(output_dir, "04_accessory_heatmap.pdf"),
  width             = 12,
  height            = 14
)

# =============================================================================
#Fisher's exact test — gene ~ InhibitoryType
# =============================================================================
strong_idx <- which(meta$InhibitoryType == "Strongly Inhibitory")
weak_idx   <- which(meta$InhibitoryType == "Weakly Inhibitory")
n_strong   <- length(strong_idx)
n_weak     <- length(weak_idx)

# Run Fisher's test for every gene
fisher_results <- apply(pa_mat, 1, function(row) {
  a <- sum(row[strong_idx])      # present in strong
  b <- n_strong - a              # absent in strong
  c <- sum(row[weak_idx])        # present in weak
  d <- n_weak - c                # absent in weak
  ct <- matrix(c(a, b, c, d), nrow = 2)
  ft <- fisher.test(ct, alternative = "two.sided")
  c(p_value = ft$p.value, odds_ratio = as.numeric(ft$estimate),
    n_strong_present = a, n_weak_present = c)
}) %>% t() %>% as.data.frame() %>%
  rownames_to_column("Gene") %>%
  mutate(
    p_adj       = p.adjust(p_value, method = "BH"),
    Direction   = case_when(
      n_strong_present > n_weak_present ~ "Enriched in Strong",
      n_strong_present < n_weak_present ~ "Enriched in Weak",
      TRUE ~ "No difference"
    )
  ) %>%
  arrange(p_adj)

# Significant hits at FDR < 0.05
sig_fisher <- fisher_results %>% filter(p_adj < 0.05)
cat(sprintf("Significant genes (FDR < 0.05): %d\n", nrow(sig_fisher)))
print(sig_fisher)
#at the global level, no significant effects of inhibitory type on gene abundance

write_csv(fisher_results, file.path(output_dir, "05_fisher_results_all_genes.csv"))
write_csv(sig_fisher,     file.path(output_dir, "05_fisher_results_significant.csv"))

# Volcano-style plot: -log10(p_adj) vs. log2(odds ratio)
volcano_df <- fisher_results %>%
  mutate(
    log2OR    = log2(odds_ratio + 0.01),
    neglog10p = -log10(p_adj + 1e-10),
    Sig       = p_adj < 0.05
  )

p5 <- ggplot(volcano_df, aes(x = log2OR, y = neglog10p,
                             colour = Direction, alpha = Sig)) +
  geom_point(size = 1.5) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", colour = "grey40") +
  geom_text_repel(
    data    = filter(volcano_df, Sig),
    aes(label = Gene),
    size    = 2.8,
    max.overlaps = 30
  ) +
  scale_colour_manual(values = c(
    "Enriched in Strong" = "#D7263D",
    "Enriched in Weak"   = "#4A90D9",
    "No difference"      = "grey60"
  )) +
  scale_alpha_manual(values = c(`TRUE` = 1, `FALSE` = 0.3), guide = "none") +
  labs(
    title    = "Fisher's Exact Test: Gene ~ Inhibitory Type",
    subtitle = "Dashed line = FDR 5%",
    x        = "log2(Odds Ratio)",
    y        = "-log10(adjusted p-value)",
    colour   = NULL
  ) +
  theme_classic(base_size = 13) +
  theme(plot.title = element_text(face = "bold"))

save_pdf(p5, "05_fisher_volcano.pdf", w = 9, h = 7)
#not necessary to include

# =============================================================================
#Linear regression — gene presence ~ MeanInhib
# =============================================================================
mean_inhib <- meta$MeanInhib

lm_results <- apply(pa_mat, 1, function(row) {
  df <- data.frame(presence = as.numeric(row), inhib = mean_inhib)
  fit <- lm(inhib ~ presence, data = df)
  s   <- summary(fit)
  coef_tbl <- coef(s)
  if (nrow(coef_tbl) < 2) {
    return(c(estimate = NA, std_error = NA, t_value = NA, p_value = NA, r_squared = NA))
  }
  c(
    estimate  = coef_tbl["presence", "Estimate"],
    std_error = coef_tbl["presence", "Std. Error"],
    t_value   = coef_tbl["presence", "t value"],
    p_value   = coef_tbl["presence", "Pr(>|t|)"],
    r_squared = s$r.squared
  )
}) %>% t() %>% as.data.frame() %>%
  rownames_to_column("Gene") %>%
  drop_na(p_value) %>%
  mutate(p_adj = p.adjust(p_value, method = "BH")) %>%
  arrange(p_adj)

sig_lm <- lm_results %>% filter(p_adj < 0.05)
cat(sprintf("Significant genes from regression (FDR < 0.05): %d\n", nrow(sig_lm)))
print(sig_lm)
#Significant genes from regression (FDR < 0.05): 787

write_csv(lm_results, file.path(output_dir, "06_linear_regression_all_genes.csv"))
write_csv(sig_lm,     file.path(output_dir, "06_linear_regression_significant.csv"))

# Manhattan-style plot of regression results
lm_plot_df <- lm_results %>%
  mutate(
    neglog10p = -log10(p_adj + 1e-10),
    Sig       = p_adj < 0.05,
    Direction = ifelse(estimate > 0, "Positive", "Negative")
  )

p6 <- ggplot(lm_plot_df, aes(x = estimate, y = neglog10p,
                             colour = Direction, alpha = Sig)) +
  geom_point(size = 1.5) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", colour = "grey40") +
  geom_vline(xintercept = 0, linetype = "dotted", colour = "grey60") +
  geom_text_repel(
    data    = filter(lm_plot_df, Sig),
    aes(label = Gene),
    size    = 2.8,
    max.overlaps = 30
  ) +
  scale_colour_manual(values = c("Positive" = "#D7263D", "Negative" = "#4A90D9")) +
  scale_alpha_manual(values = c(`TRUE` = 1, `FALSE` = 0.3), guide = "none") +
  labs(
    title    = "Linear Regression: Gene Presence ~ Mean Inhibition (%)",
    subtitle = "Dashed line = FDR 5% | Positive = gene presence increases inhibition",
    x        = "Regression estimate (effect on % inhibition)",
    y        = "-log10(adjusted p-value)",
    colour   = "Effect direction"
  ) +
  theme_classic(base_size = 13) +
  theme(plot.title = element_text(face = "bold"))

save_pdf(p6, "06_regression_volcano.pdf", w = 9, h = 7)

#plot inhibition values split by presence/absence for the top hit
if (nrow(sig_lm) > 0) {
  top_gene <- sig_lm$Gene[1]
  top_df   <- data.frame(
    Presence  = factor(pa_mat[top_gene, ], labels = c("Absent", "Present")),
    MeanInhib = mean_inhib,
    Isolate   = meta$Isolate,
    Pond      = meta$Pond
  )
  p6b <- ggplot(top_df, aes(x = Presence, y = MeanInhib, fill = Presence)) +
    geom_boxplot(alpha = 0.7, outlier.shape = NA, width = 0.4) +
    geom_jitter(aes(colour = Pond), width = 0.1, size = 3) +
    geom_text_repel(aes(label = Isolate), size = 2.8) +
    scale_fill_manual(values  = c("Absent" = "#AECDE8", "Present" = "#D7263D")) +
    scale_colour_manual(values = pal_pond) +
    labs(
      title  = sprintf("Top regression hit: %s", top_gene),
      x      = sprintf("%s presence", top_gene),
      y      = "Mean inhibition (%)",
      fill   = NULL, colour = "Pond"
    ) +
    theme_classic(base_size = 13) +
    theme(plot.title = element_text(face = "bold"), legend.position = "right")
  save_pdf(p6b, "06b_top_regression_gene_boxplot.pdf", w = 7, h = 6)
}

# =============================================================================
#Pond comparisons — gene content & pangenome structure per pond
# =============================================================================
ponds      <- unique(meta$Pond)
pond_idx   <- lapply(setNames(ponds, ponds), function(p) which(meta$Pond == p))

#Per-pond pangenome breakdown (core / accessory / unique within pond) -
pond_breakdown <- map_dfr(ponds, function(p) {
  idx     <- pond_idx[[p]]
  n       <- length(idx)
  sub_mat <- pa_mat[, idx]
  freq    <- rowSums(sub_mat)
  tibble(
    Pond       = p,
    N_isolates = n,
    Pan        = sum(freq > 0),
    Core       = sum(freq == n),
    Accessory  = sum(freq > 1 & freq < n),
    Unique     = sum(freq == 1)
  )
})

cat("\nPer-pond pangenome breakdown:\n")
print(pond_breakdown)
write_csv(pond_breakdown, file.path(output_dir, "07a_pond_pangenome_breakdown.csv"))

# Stacked bar comparing pond compositions
pond_bar_df <- pond_breakdown %>%
  pivot_longer(c(Core, Accessory, Unique),
               names_to = "Category", values_to = "Count") %>%
  mutate(Category = factor(Category, levels = c("Core", "Accessory", "Unique")))

p7a <- ggplot(pond_bar_df, aes(x = Pond, y = Count, fill = Category)) +
  geom_col(position = "dodge", colour = "white", linewidth = 0.4) +
  geom_text(aes(label = Count),
            position = position_dodge(width = 0.9),
            vjust = -0.4, size = 3.5) +
  scale_fill_manual(values = c(
    "Core"      = "#2166AC",
    "Accessory" = "#92C5DE",
    "Unique"    = "#F4A582"
  )) +
  labs(title = "Per-Pond Pangenome Breakdown",
       x = NULL, y = "Number of genes", fill = "Gene category") +
  theme_classic(base_size = 13) +
  theme(plot.title = element_text(face = "bold"))

save_pdf(p7a, "07a_pond_pangenome_breakdown_bar.pdf", w = 8, h = 6)

#Pond-specific genes (only present in one pond)
# For each gene present in ≥1 isolate, record which ponds carry it
gene_in_pond <- map(setNames(ponds, ponds), function(p) {
  rownames(pa_mat)[rowSums(pa_mat[, pond_idx[[p]], drop = FALSE]) > 0]
})

pond_exclusive <- map(setNames(ponds, ponds), function(p) {
  other_ponds <- setdiff(ponds, p)
  other_genes <- Reduce(union, gene_in_pond[other_ponds])
  setdiff(gene_in_pond[[p]], other_genes)
})

cat("\nPond-exclusive gene counts:\n")
walk2(names(pond_exclusive), pond_exclusive, ~cat(sprintf("  %s: %d genes\n", .x, length(.y))))
#  Sixty Lake: 949 genes
#  Conness Pond: 1075 genes

pond_excl_df <- tibble(
  Pond  = names(pond_exclusive),
  Count = lengths(pond_exclusive)
)

p7b <- ggplot(pond_excl_df, aes(x = Pond, y = Count, fill = Pond)) +
  geom_col(width = 0.5, colour = "white") +
  geom_text(aes(label = Count), vjust = -0.5, size = 4, fontface = "bold") +
  scale_fill_manual(values = pal_pond) +
  labs(title = "Pond-Exclusive Genes",
       subtitle = "Genes found only in isolates from one pond",
       x = NULL, y = "Number of genes") +
  theme_classic(base_size = 13) +
  theme(plot.title = element_text(face = "bold"), legend.position = "none")

save_pdf(p7b, "07b_pond_exclusive_genes.pdf", w = 7, h = 6)

#Fisher's exact test: gene ~ Pond
sl_idx   <- pond_idx[["Sixty Lake"]]
cp_idx   <- pond_idx[["Conness Pond"]]
n_sl     <- length(sl_idx)
n_cp     <- length(cp_idx)

pond_fisher <- apply(pa_mat, 1, function(row) {
  a <- sum(row[sl_idx]);  b <- n_sl - a
  c <- sum(row[cp_idx]);  d <- n_cp - c
  ct <- matrix(c(a, b, c, d), nrow = 2)
  ft <- fisher.test(ct, alternative = "two.sided")
  c(p_value = ft$p.value, odds_ratio = as.numeric(ft$estimate),
    n_SixtyLake = a, n_ConnessPond = c)
}) %>% t() %>% as.data.frame() %>%
  rownames_to_column("Gene") %>%
  mutate(
    p_adj     = p.adjust(p_value, method = "BH"),
    Direction = case_when(
      n_SixtyLake > n_ConnessPond ~ "Enriched in Sixty Lake",
      n_SixtyLake < n_ConnessPond ~ "Enriched in Conness Pond",
      TRUE ~ "No difference"
    )
  ) %>%
  arrange(p_adj)

sig_pond_fisher <- pond_fisher %>% filter(p_adj < 0.05)
cat(sprintf("\nPond Fisher's test — significant genes (FDR < 0.05): %d\n",
            nrow(sig_pond_fisher)))
print(sig_pond_fisher)

write_csv(pond_fisher,       file.path(output_dir, "07c_pond_fisher_all_genes.csv"))
write_csv(sig_pond_fisher,   file.path(output_dir, "07c_pond_fisher_significant.csv"))

# Volcano for pond Fisher's test
pond_volcano_df <- pond_fisher %>%
  mutate(
    log2OR    = log2(odds_ratio + 0.01),
    neglog10p = -log10(p_adj + 1e-10),
    Sig       = p_adj < 0.05
  )

p7c <- ggplot(pond_volcano_df,
              aes(x = log2OR, y = neglog10p, colour = Direction, alpha = Sig)) +
  geom_point(size = 1.5) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", colour = "grey40") +
  geom_text_repel(data = filter(pond_volcano_df, Sig),
                  aes(label = Gene), size = 2.8, max.overlaps = 30) +
  scale_colour_manual(values = c(
    "Enriched in Sixty Lake"   = "#2E8B57",
    "Enriched in Conness Pond" = "#E07B39",
    "No difference"            = "grey60"
  )) +
  scale_alpha_manual(values = c(`TRUE` = 1, `FALSE` = 0.3), guide = "none") +
  labs(title    = "Fisher's Exact Test: Gene ~ Pond",
       subtitle = "Dashed line = FDR 5%",
       x        = "log2(Odds Ratio)",
       y        = "-log10(adjusted p-value)",
       colour   = NULL) +
  theme_classic(base_size = 13) +
  theme(plot.title = element_text(face = "bold"))

save_pdf(p7c, "07c_pond_fisher_volcano.pdf", w = 9, h = 7)

#Gene count per isolate, split by pond
p7d <- ggplot(gene_count_df, aes(x = Pond, y = GeneCount, fill = Pond)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA, width = 0.45) +
  geom_jitter(aes(colour = InhibitoryType), width = 0.13, size = 3) +
  geom_text_repel(aes(label = Isolate), size = 2.8, max.overlaps = 20) +
  scale_fill_manual(values  = pal_pond) +
  scale_colour_manual(values = pal_inhibitory) +
  labs(title  = "Gene Content per Isolate by Pond",
       x = NULL, y = "Number of genes",
       fill = "Pond", colour = "Inhibitory Type") +
  theme_classic(base_size = 13) +
  theme(plot.title = element_text(face = "bold"))

save_pdf(p7d, "07d_gene_count_by_pond.pdf", w = 8, h = 6)

wtest_pond <- wilcox.test(GeneCount ~ Pond, data = gene_count_df)
cat(sprintf("Gene count by pond — Wilcoxon: W = %.1f, p = %.4f\n",wtest_pond$statistic, wtest_pond$p.value))
#Gene count by pond — Wilcoxon: W = 81.0, p = 0.0004

# =============================================================================
# Gene category breakdown — what varies across Core/Accessory/Unique
# =============================================================================
# Assign each gene to a category
gene_category <- case_when(
  gene_freq == n_isolates             ~ "Core",
  gene_freq > 1 & gene_freq < n_isolates ~ "Accessory",
  gene_freq == 1                      ~ "Unique"
)
names(gene_category) <- rownames(pa_mat)

# Build annotated gene table
gene_annot_df <- tibble(
  Gene      = rownames(pa_mat),
  Freq      = gene_freq,
  Category  = gene_category,
  # Flag named genes (not group_XXXX) as likely annotated
  IsNamed   = !str_detect(Gene, "^group_"),
  # Pond presence flags
  InSixtyLake   = rowSums(pa_mat[, sl_idx, drop = FALSE]) > 0,
  InConnessPond = rowSums(pa_mat[, cp_idx, drop = FALSE]) > 0,
  PondPattern   = case_when(
    InSixtyLake  & InConnessPond  ~ "Both ponds",
    InSixtyLake  & !InConnessPond ~ "Sixty Lake only",
    !InSixtyLake & InConnessPond  ~ "Conness Pond only",
    TRUE                          ~ "Absent"
  )
)

write_csv(gene_annot_df, file.path(output_dir, "08_gene_category_annotation.csv"))

#Named gene rate per category
named_rate <- gene_annot_df %>%
  group_by(Category) %>%
  summarise(
    Total      = n(),
    Named      = sum(IsNamed),
    Unnamed    = Total - Named,
    PctNamed   = 100 * Named / Total,
    .groups    = "drop"
  ) %>%
  mutate(Category = factor(Category, levels = c("Core", "Accessory", "Unique")))

cat("\nAnnotation rate by gene category:\n")
print(named_rate)
#  Category  Total Named Unnamed PctNamed
#1 Accessory  2146   305    1841    14.2 
#2 Core       5147  3080    2067    59.8 
#3 Unique      485    42     443     8.66

p8a <- ggplot(named_rate %>%
                pivot_longer(c(Named, Unnamed),
                             names_to = "Type", values_to = "Count"),
              aes(x = Category, y = Count, fill = Type)) +
  geom_col(position = "fill", colour = "white", linewidth = 0.3) +
  geom_text(aes(label = ifelse(Count > 50,
                               sprintf("%d\n(%.0f%%)", Count,
                                       100 * Count / (Named + Unnamed)),
                               "")),
            position = position_fill(vjust = 0.5), size = 3.2,
            data = . %>% mutate(Named = named_rate$Named[match(Category, named_rate$Category)],
                                Unnamed = named_rate$Unnamed[match(Category, named_rate$Category)])) +
  scale_fill_manual(values = c("Named" = "#2166AC", "Unnamed" = "#D9D9D9")) +
  scale_y_continuous(labels = percent_format()) +
  labs(title    = "Annotated vs. Hypothetical Genes by Category",
       subtitle = "Named = known gene symbol; Unnamed = group_XXXX (hypothetical)",
       x = NULL, y = "Proportion of genes", fill = NULL) +
  theme_classic(base_size = 13) +
  theme(plot.title = element_text(face = "bold"))

save_pdf(p8a, "08a_named_gene_rate_by_category.pdf", w = 8, h = 6)

#Pond distribution within each gene category
pond_by_cat <- gene_annot_df %>%
  filter(Category != "Unique") %>%   # unique genes are by definition single-isolate
  group_by(Category, PondPattern) %>%
  summarise(Count = n(), .groups = "drop") %>%
  mutate(Category = factor(Category, levels = c("Core", "Accessory")))

p8b <- ggplot(pond_by_cat, aes(x = Category, y = Count, fill = PondPattern)) +
  geom_col(position = "dodge", colour = "white", linewidth = 0.3) +
  geom_text(aes(label = Count),
            position = position_dodge(width = 0.9), vjust = -0.4, size = 3.2) +
  scale_fill_manual(values = c(
    "Both ponds"       = "#6A3D9A",
    "Sixty Lake only"  = "#2E8B57",
    "Conness Pond only"= "#E07B39"
  )) +
  labs(title = "Pond Distribution of Genes Within Each Category",
       x = NULL, y = "Number of genes", fill = "Pond pattern") +
  theme_classic(base_size = 13) +
  theme(plot.title = element_text(face = "bold"))

save_pdf(p8b, "08b_pond_pattern_by_category.pdf", w = 8, h = 6)

#Unique gene — which isolate carries each singleton?
singleton_df <- tibble(
  Gene     = names(gene_category[gene_category == "Unique"]),
  Isolate  = apply(pa_mat[gene_category == "Unique", , drop = FALSE], 1,
                   function(r) colnames(pa_mat)[which(r == 1)[1]])
) %>%
  left_join(meta %>% select(Isolate, Pond, InhibitoryType), by = "Isolate") %>%
  group_by(Isolate, Pond, InhibitoryType) %>%
  summarise(N_unique_genes = n(), .groups = "drop") %>%
  arrange(desc(N_unique_genes))

cat("\nUnique (singleton) gene counts per isolate:\n")
print(singleton_df)
write_csv(singleton_df, file.path(output_dir, "08c_singleton_genes_per_isolate.csv"))

p8c <- ggplot(singleton_df,
              aes(x = reorder(Isolate, N_unique_genes),
                  y = N_unique_genes,
                  fill = Pond)) +
  geom_col(colour = "white", linewidth = 0.3) +
  geom_text(aes(label = N_unique_genes), hjust = -0.2, size = 3.5) +
  scale_fill_manual(values = pal_pond) +
  coord_flip() +
  expand_limits(y = max(singleton_df$N_unique_genes) * 1.15) +
  labs(title    = "Unique (Singleton) Genes per Isolate",
       subtitle = "Genes found in only one isolate across the whole pangenome",
       x = NULL, y = "Number of unique genes", fill = "Pond") +
  theme_classic(base_size = 13) +
  theme(plot.title = element_text(face = "bold"))

save_pdf(p8c, "08c_singleton_genes_per_isolate.pdf", w = 8, h = 6)

#Heatmap of accessory genes split by pond-pattern
# Order genes by PondPattern then by frequency for interpretability
acc_annot <- gene_annot_df %>%
  filter(Category == "Accessory") %>%
  arrange(PondPattern, desc(Freq)) %>%
  column_to_rownames("Gene")

acc_mat_ordered <- accessory_mat[rownames(acc_annot), ]

row_annot <- data.frame(
  PondPattern = acc_annot$PondPattern,
  row.names   = rownames(acc_annot)
)

row_annot_colors <- list(
  PondPattern = c(
    "Both ponds"        = "#6A3D9A",
    "Sixty Lake only"   = "#2E8B57",
    "Conness Pond only" = "#E07B39"
  )
)

pheatmap(
  mat               = acc_mat_ordered,
  color             = colorRampPalette(c("white", "#2166AC"))(2),
  annotation_col    = col_annot,
  annotation_row    = row_annot,
  annotation_colors = c(annot_colors, row_annot_colors),
  cluster_rows      = FALSE,   # already sorted by pond pattern
  cluster_cols      = TRUE,
  clustering_method = "ward.D2",
  show_rownames     = nrow(acc_mat_ordered) <= 200,
  show_colnames     = TRUE,
  fontsize_col      = 9,
  fontsize_row      = 5,
  main              = "Accessory Genome — Ordered by Pond Pattern",
  filename          = file.path(output_dir, "08d_accessory_heatmap_by_pond.pdf"),
  width             = 13,
  height            = 15
)
message("Saved: ", file.path(output_dir, "08d_accessory_heatmap_by_pond.pdf"))

# =============================================================================
#Conness Pond — Fisher's exact test gene ~ InhibitoryType
# =============================================================================
# Subset matrix and metadata to Conness Pond only
cp_meta <- meta %>% filter(Pond == "Conness Pond")
cp_mat  <- pa_mat[, cp_meta$Isolate]

# Remove genes with zero variance within Conness Pond (all present or all absent)
cp_var  <- apply(cp_mat, 1, function(r) length(unique(r)) > 1)
cp_mat  <- cp_mat[cp_var, ]
message(sprintf("Conness Pond: %d isolates, %d variable genes tested",ncol(cp_mat), nrow(cp_mat)))
#Conness Pond: 9 isolates, 1526 variable genes tested

cp_strong_idx <- which(cp_meta$InhibitoryType == "Strongly Inhibitory")
cp_weak_idx   <- which(cp_meta$InhibitoryType == "Weakly Inhibitory")
cp_n_strong   <- length(cp_strong_idx)
cp_n_weak     <- length(cp_weak_idx)

cat(sprintf("  Strongly Inhibitory: %d isolates\n  Weakly Inhibitory:   %d isolates\n",
            cp_n_strong, cp_n_weak))
#Strongly Inhibitory: 5 isolates
#Weakly Inhibitory:   4 isolates

cp_fisher <- apply(cp_mat, 1, function(row) {
  a <- sum(row[cp_strong_idx]);  b <- cp_n_strong - a
  c <- sum(row[cp_weak_idx]);    d <- cp_n_weak   - c
  ct <- matrix(c(a, b, c, d), nrow = 2)
  ft <- tryCatch(
    fisher.test(ct, alternative = "two.sided"),
    error = function(e) list(p.value = NA, estimate = NA)
  )
  c(p_value          = ft$p.value,
    odds_ratio       = as.numeric(ft$estimate),
    n_strong_present = a,
    n_weak_present   = c)
}) %>% t() %>% as.data.frame() %>%
  rownames_to_column("Gene") %>%
  drop_na(p_value) %>%
  mutate(
    p_adj     = p.adjust(p_value, method = "BH"),
    Direction = case_when(
      n_strong_present > n_weak_present ~ "Enriched in Strong",
      n_strong_present < n_weak_present ~ "Enriched in Weak",
      TRUE ~ "No difference"
    ),
    Category  = gene_category[Gene]
  ) %>%
  arrange(p_adj)

sig_cp_fisher <- cp_fisher %>% filter(p_adj < 0.05)
cat(sprintf("\nConness Pond — significant genes (FDR < 0.05): %d\n",nrow(sig_cp_fisher)))
#Conness Pond — significant genes (FDR < 0.05): 1495

print(sig_cp_fisher)

write_csv(cp_fisher,     file.path(output_dir, "09_conness_fisher_all_genes.csv"))
write_csv(sig_cp_fisher, file.path(output_dir, "09_conness_fisher_significant.csv"))

#Volcano plot
cp_volcano_df <- cp_fisher %>%
  mutate(
    log2OR    = log2(odds_ratio + 0.01),
    neglog10p = -log10(p_adj + 1e-10),
    Sig       = p_adj < 0.05
  )

p9a <- ggplot(cp_volcano_df,
              aes(x = log2OR, y = neglog10p, colour = Direction, alpha = Sig)) +
  geom_point(aes(shape = Category), size = 2) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", colour = "grey40") +
  geom_text_repel(
    data    = filter(cp_volcano_df, Sig),
    aes(label = Gene),
    size    = 3, max.overlaps = 30
  ) +
  scale_colour_manual(values = c(
    "Enriched in Strong" = "#D7263D",
    "Enriched in Weak"   = "#4A90D9",
    "No difference"      = "grey60"
  )) +
  scale_shape_manual(values = c(
    "Core"      = 16,
    "Accessory" = 17,
    "Unique"    = 4
  )) +
  scale_alpha_manual(values = c(`TRUE` = 1, `FALSE` = 0.25), guide = "none") +
  labs(
    title    = "Conness Pond: Gene ~ Inhibitory Type (Fisher's Exact Test)",
    subtitle = "Dashed line = FDR 5% | Shape = gene category",
    x        = "log2(Odds Ratio)  [positive = enriched in Strongly Inhibitory]",
    y        = "-log10(adjusted p-value)",
    colour   = "Enrichment",
    shape    = "Gene category"
  ) +
  theme_classic(base_size = 13) +
  theme(plot.title = element_text(face = "bold"))

save_pdf(p9a, "09a_conness_fisher_volcano.pdf", w = 10, h = 7)

#Presence/absence heatmap for significant genes in Conness Pond
if (nrow(sig_cp_fisher) > 0) {
  sig_genes_cp <- sig_cp_fisher$Gene
  
  # Pull full matrix rows for sig genes, all isolates, ordered by pond then type
  isolate_order <- meta %>%
    arrange(Pond, InhibitoryType) %>%
    pull(Isolate)
  
  sig_full_mat <- pa_mat[sig_genes_cp, isolate_order, drop = FALSE]
  
  # Row annotation: direction + category
  row_annot_cp <- sig_cp_fisher %>%
    select(Gene, Direction, Category) %>%
    column_to_rownames("Gene")
  
  row_colors_cp <- list(
    Direction = c(
      "Enriched in Strong" = "#D7263D",
      "Enriched in Weak"   = "#4A90D9",
      "No difference"      = "grey80"
    ),
    Category = c(
      "Core"      = "#2166AC",
      "Accessory" = "#92C5DE",
      "Unique"    = "#F4A582"
    )
  )
  
  # Column annotation ordered to match isolate_order
  col_annot_ordered <- data.frame(
    InhibitoryType = meta$InhibitoryType[match(isolate_order, meta$Isolate)],
    Pond           = meta$Pond[match(isolate_order, meta$Isolate)],
    row.names      = isolate_order
  )
  
  pheatmap(
    mat               = sig_full_mat,
    color             = colorRampPalette(c("white", "#2166AC"))(2),
    annotation_col    = col_annot_ordered,
    annotation_row    = row_annot_cp,
    annotation_colors = c(annot_colors, row_colors_cp),
    cluster_rows      = TRUE,
    cluster_cols      = FALSE,   # keep pond/type grouping visible
    clustering_method = "ward.D2",
    show_rownames     = TRUE,
    show_colnames     = TRUE,
    fontsize_col      = 9,
    fontsize_row      = 8,
    main              = "Conness Pond: Significant Genes by Inhibitory Type",
    filename          = file.path(output_dir, "09b_conness_sig_genes_heatmap.pdf"),
    width             = 12,
    height            = max(6, nrow(sig_full_mat) * 0.3 + 4)
  )
  message("Saved: ", file.path(output_dir, "09b_conness_sig_genes_heatmap.pdf"))
} else {
  message("No significant genes at FDR < 0.05 in Conness Pond — heatmap skipped.")
  message("Consider inspecting 09_conness_fisher_all_genes.csv for top nominally significant genes.")
}

# Top hits — individual boxplots of inhibition by gene presence
# Plot top 6 (or fewer) significant genes showing MeanInhib split by presence
top_n  <- min(6, nrow(sig_cp_fisher))

if (top_n > 0) {
  top_genes_cp <- sig_cp_fisher$Gene[1:top_n]
  
  top_plots <- map(top_genes_cp, function(g) {
    df <- tibble(
      Presence  = factor(pa_mat[g, cp_meta$Isolate],
                         levels = c(0, 1), labels = c("Absent", "Present")),
      MeanInhib = cp_meta$MeanInhib,
      Isolate   = cp_meta$Isolate,
      Type      = cp_meta$InhibitoryType
    )
    ggplot(df, aes(x = Presence, y = MeanInhib, fill = Presence)) +
      geom_boxplot(alpha = 0.7, outlier.shape = NA, width = 0.45) +
      geom_jitter(aes(colour = Type), width = 0.12, size = 3) +
      geom_text_repel(aes(label = Isolate), size = 2.5, max.overlaps = 10) +
      scale_fill_manual(values  = c("Absent" = "#AECDE8", "Present" = "#D7263D")) +
      scale_colour_manual(values = pal_inhibitory) +
      labs(title  = g,
           x = NULL, y = "Mean inhibition (%)",
           fill = NULL, colour = NULL) +
      theme_classic(base_size = 11) +
      theme(plot.title  = element_text(face = "bold.italic", size = 10),
            legend.position = "none")
  })
  
  # Combine into one patchwork figure
  p9c <- wrap_plots(top_plots, ncol = min(3, top_n)) +
    plot_annotation(
      title    = "Conness Pond: Top Significant Genes vs. Mean Inhibition",
      subtitle = sprintf("Top %d genes by FDR — Conness Pond isolates only", top_n),
      theme    = theme(plot.title    = element_text(face = "bold", size = 13),
                       plot.subtitle = element_text(colour = "grey40"))
    )
  
  save_pdf(p9c, "09c_conness_top_genes_boxplots.pdf",
           w = min(3, top_n) * 4, h = ceiling(top_n / 3) * 4.5)
}
# =============================================================================
# Conness Pond: genes separating strongly vs weakly inhibitory isolates
# The ~1500 Fisher-significant genes fall into two presence/absence blocks, so
# the figure shows (A) the blocks across isolates and (B) what the blocks
# contain (COG categories) relative to the Conness core genome.
# =============================================================================
library(tidyverse)
library(patchwork)
select <- dplyr::select; filter <- dplyr::filter; summarise <- dplyr::summarise
mutate <- dplyr::mutate; count <- dplyr::count; arrange <- dplyr::arrange

base_dir <- "~/Documents/GitHub/Pseudo_fluor"
sig  <- read_csv( "~/Documents/GitHub/Pseudo_fluor/figures/09_conness_fisher_significant.csv", show_col_types = FALSE)
ann  <- read_tsv("~/Documents/GitHub/Pseudo_fluor/non_prokka_annotation/pangenome_annotation.tsv", show_col_types = FALSE, guess_max = 10000)
pa   <- read.delim("~/Documents/GitHub/Pseudo_fluor/panaroo_results/gene_presence_absence2.Rtab", check.names = FALSE)
meta <- read_tsv("~/Documents/GitHub/Pseudo_fluor/isolate_inhibi_category.txt", show_col_types = FALSE)

col_strong <- "#D55E00"; col_weak <- "#0072B2"; col_core <- "grey55"

cp <- meta %>% filter(Pond == "Conness Pond") %>%
  mutate(Inhibition_Flag = factor(Inhibition_Flag, c("Strongly Inhibitory", "Weakly Inhibitory"))) %>%
  arrange(Inhibition_Flag, desc(mean))

#blocks: label each significant gene by its presence pattern ------------
sig <- sig %>%
  mutate(block = case_when(n_strong_present == 5 & n_weak_present == 0 ~ "Strong-only",
                           n_strong_present == 0 & n_weak_present == 4 ~ "Weak-only",
                           TRUE ~ "Near-complete split"),
         block = factor(block, c("Strong-only", "Weak-only", "Near-complete split")))
block_n <- count(sig, block)
#  block                   n
#<fct>               <int>
#1 Strong-only           574
#2 Weak-only             916
#3 Near-complete split     5

#A: presence/absence map (isolates x genes, genes ordered by block) ------
gene_order <- sig %>% filter(block != "Near-complete split") %>% arrange(block, Gene) %>% pull(Gene)
pa_long <- pa %>% filter(Gene %in% gene_order) %>%
  select(Gene, all_of(cp$Group)) %>%
  pivot_longer(-Gene, names_to = "Isolate", values_to = "present") %>%
  left_join(sig %>% select(Gene, block), by = "Gene") %>%
  left_join(cp %>% select(Isolate = Group, Inhibition_Flag, mean), by = "Isolate") %>%
  mutate(Gene = factor(Gene, gene_order),
         Isolate = factor(Isolate, rev(cp$Group)),
         fill = ifelse(present == 1, as.character(Inhibition_Flag), "absent"))

blk_lab <- sig %>% filter(Gene %in% gene_order) %>% mutate(x = match(Gene, gene_order)) %>% group_by(block) %>%
  summarise(xmid = mean(range(x)), xmin = min(x), xmax = max(x), n = n(), .groups = "drop") %>%
  mutate(label = sprintf("%s\n%d genes", block, n))

pA <- ggplot(pa_long, aes(as.integer(Gene), Isolate, fill = fill)) +
  geom_raster() +
  geom_vline(xintercept = blk_lab$xmax[-nrow(blk_lab)] + 0.5, colour = "white", linewidth = 0.8) +
  scale_fill_manual(values = c("Strongly Inhibitory" = col_strong, "Weakly Inhibitory" = col_weak,
                               "absent" = "grey93"),
                    labels = c("absent" = "Gene absent",
                               "Strongly Inhibitory" = "Present (strongly inhibitory isolate)",
                               "Weakly Inhibitory" = "Present (weakly inhibitory isolate)"),
                    breaks = c("Strongly Inhibitory", "Weakly Inhibitory", "absent"), name = NULL) +
  scale_x_continuous(breaks = blk_lab$xmid, labels = blk_lab$label, expand = c(0, 0),
                     sec.axis = dup_axis(breaks = NULL)) +
  scale_y_discrete(labels = function(i) {
    m <- cp$mean[match(i, cp$Group)]; sprintf("%s (%.0f%%)", i, m) }) +
  labs(title = "A  Gene presence across Conness Pond isolates",
       subtitle = sprintf(paste0("%d of %d per-gene Fisher-significant genes form two complete blocks (shown); %d more follow\n",
                                 "near-complete splits (one isolate off). Isolates labelled with mean %% Bd inhibition (Jel197)."),
                          length(gene_order), nrow(sig), sum(sig$block == "Near-complete split")),
       x = NULL, y = NULL) +
  theme_minimal(base_size = 10) +
  theme(panel.grid = element_blank(), legend.position = "bottom",
        plot.title = element_text(face = "bold", size = 10.5),
        plot.subtitle = element_text(size = 8.5, colour = "grey30"),
        axis.text.x = element_text(size = 8.5, lineheight = 0.9))

#COG composition of each block vs the Conness core genome -------------
cog_long <- function(genes, group) {
  ann %>% filter(Gene %in% genes) %>%
    mutate(letters = map(COG_cat, ~ if (is.na(.x)) "No COG" else unique(strsplit(gsub(";", "", .x), "")[[1]]))) %>%
    select(Gene, letters) %>% unnest(letters) %>%
    count(letters, name = "n_genes") %>%
    mutate(group = group, pct = 100 * n_genes / length(genes))
}
core_genes <- ann$Gene[ann$class_ConnessPond == "Core"]
comp <- bind_rows(
  cog_long(sig$Gene[sig$block == "Strong-only"], "Strong-only"),
  cog_long(sig$Gene[sig$block == "Weak-only"],   "Weak-only"),
  cog_long(core_genes, "Conness core genome")) %>%
  mutate(group = factor(group, c("Strong-only", "Weak-only", "Conness core genome")))

fun <- c(X = "Mobilome: prophages, transposons", L = "Replication, recombination, repair",
         U = "Secretion, trafficking", V = "Defense mechanisms", K = "Transcription",
         T = "Signal transduction", M = "Cell wall/membrane/envelope", E = "Amino acid metabolism",
         P = "Inorganic ion transport", C = "Energy production", G = "Carbohydrate metabolism",
         H = "Coenzyme metabolism", I = "Lipid metabolism", J = "Translation, ribosomes",
         R = "General function prediction", S = "Function unknown")
keep <- c("X", "L", "U", "V", "K", "T", "M", "E", "P", "C", "G", "H", "I", "J", "R", "S")
comp_plot <- comp %>% filter(letters %in% keep) %>%
  mutate(cat = factor(paste0(letters, "  ", fun[letters]), rev(paste0(keep, "  ", fun[keep]))))
nocog <- comp %>% filter(letters == "No COG") %>% select(group, pct) %>% deframe()

pB <- ggplot(comp_plot, aes(pct, cat, fill = group)) +
  geom_col(position = position_dodge2(width = 0.8, reverse = TRUE, padding = 0.1), width = 0.75) +
  scale_fill_manual(values = c("Strong-only" = col_strong, "Weak-only" = col_weak,
                               "Conness core genome" = col_core),
                    labels = c("Strong-only" = sprintf("Strong-only block (%d genes; %.0f%% have no COG)", sum(sig$block == "Strong-only"), nocog["Strong-only"]),
                               "Weak-only"   = sprintf("Weak-only block (%d genes; %.0f%% have no COG)", sum(sig$block == "Weak-only"), nocog["Weak-only"]),
                               "Conness core genome" = sprintf("Conness core genome (%s genes; %.0f%% have no COG)", format(length(core_genes), big.mark = ","), nocog["Conness core genome"])),
                    name = NULL) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.05)), labels = function(x) paste0(x, "%")) +
  labs(title = "B  What the blocks contain: COG functional categories",
       subtitle = "Percent of all genes in each group assigned to each category (genes with several COG letters count in each)",
       x = "Genes in group", y = NULL) +
  theme_classic(base_size = 10) +
  theme(legend.position = "bottom", legend.direction = "vertical",
        plot.title = element_text(face = "bold", size = 10.5),
        plot.subtitle = element_text(size = 8.5, colour = "grey30"))

fig <- pA / pB + plot_layout(heights = c(0.9, 2)) +
  plot_annotation(caption = paste(
    "Per-gene Fisher tests are not independent here: genes sharing a presence pattern give identical results.",
    "Treat the blocks as two lineage-level differences, not ~1,500 separate associations.", sep = "\n"),
    theme = theme(plot.caption = element_text(size = 7.5, colour = "grey35", hjust = 0)))

ggsave(file.path(base_dir, "figures/09d_conness_blocks_function.pdf"), fig, width = 10, height = 10, useDingbats = FALSE)
ggsave(file.path(base_dir, "figures/09d_conness_blocks_function.png"), fig, width = 10, height = 10, dpi = 300)
write_tsv(comp, file.path(base_dir, "figures/09d_conness_blocks_cog_composition.tsv"))
