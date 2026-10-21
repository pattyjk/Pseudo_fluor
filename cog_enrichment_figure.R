# =============================================================================
# COG category composition: variable (accessory + unique) vs core genome
# Input: cog_enrichment_variable_vs_core.tsv (one row per set x COG letter:
#        n_variable, pct_variable, n_core, pct_core, log2_ratio, p_adj)
# =============================================================================
library(tidyverse)
library(patchwork)

base_dir <- "~/Documents/GitHub/Pseudo_fluor"
cog <- read_tsv(file.path(base_dir, "cog_enrichment_variable_vs_core.tsv"), show_col_types = FALSE)

sets <- c("All isolates", "Sixty Lake", "Conness Pond")
col_var  <- "#D55E00"   # enriched in variable genome
col_core <- "#0072B2"   # enriched in core genome

# order categories by the all-isolates ratio (most variable-enriched at top)
cat_order <- cog %>% filter(set == "All isolates") %>% arrange(log2_ratio) %>%
  mutate(lab = paste0(COG_letter, "  ", desc)) %>% pull(lab)

plot_df <- cog %>%
  mutate(lab   = factor(paste0(COG_letter, "  ", desc), levels = cat_order),
         set   = factor(set, levels = sets),
         Sig   = p_adj < 0.05,
         Class = case_when(Sig & log2_ratio > 0 ~ "Enriched in variable genome",
                           Sig & log2_ratio < 0 ~ "Enriched in core genome",
                           TRUE                 ~ "Not significant"))

# number of COG-annotated variable / core genes per set (for the panel headers)
totals <- plot_df %>% group_by(set) %>%
  summarise(n_var  = round(max(n_variable / (pct_variable / 100), na.rm = TRUE)),
            n_core = round(max(n_core / (pct_core / 100), na.rm = TRUE)), .groups = "drop") %>%
  mutate(strip = sprintf("%s\nCOG-annotated genes:\n%s variable | %s core", set,
                         format(n_var, big.mark = ",", trim = TRUE), format(n_core, big.mark = ",", trim = TRUE)))
plot_df <- plot_df %>% left_join(totals %>% select(set, strip), by = "set") %>%
  mutate(strip = factor(strip, levels = totals$strip[match(sets, totals$set)]))

p_main <- ggplot(plot_df, aes(log2_ratio, lab)) +
  geom_vline(xintercept = 0, colour = "grey60", linewidth = 0.4) +
  geom_segment(aes(x = 0, xend = log2_ratio, yend = lab, colour = Class), linewidth = 0.7) +
  geom_point(aes(colour = Class, fill = Class), shape = 21, size = 2.6, stroke = 0.6) +
  scale_colour_manual(values = c("Enriched in variable genome" = col_var,
                                 "Enriched in core genome"     = col_core,
                                 "Not significant"             = "grey65")) +
  scale_fill_manual(values = c("Enriched in variable genome" = col_var,
                               "Enriched in core genome"     = col_core,
                               "Not significant"             = "white")) +
  facet_wrap(~ strip, nrow = 1) +
  scale_x_continuous(limits = c(-3.6, 5.6), breaks = seq(-2, 4, 2)) +
  labs(x = expression(log[2]~"(share of variable genes / share of core genes)"),
       y = NULL, colour = NULL, fill = NULL,
       title = "Variable genome is enriched for mobile, defence and secretion functions;\ncore genome for housekeeping",
       caption = "Filled points: BH-adjusted Fisher's exact p < 0.05 within each set. Variable = accessory + unique genes.\nOnly COG-annotated genes are included; a gene with several COG letters counts once per letter.") +
  theme_classic(base_size = 10) +
  theme(strip.background = element_blank(), strip.text = element_text(size = 9, hjust = 0),
        plot.title = element_text(size = 10.5), plot.caption = element_text(size = 7.5, colour = "grey35", hjust = 0),
        panel.spacing.x = unit(0.8, "lines"), legend.position = "bottom",
        axis.text.y = element_text(size = 8.5))

# right-hand annotation: brackets spanning the variable- vs core-enriched categories
ref  <- plot_df %>% filter(set == "All isolates") %>% mutate(y = as.integer(lab))
up   <- ref %>% filter(Class == "Enriched in variable genome")
down <- ref %>% filter(Class == "Enriched in core genome")
br <- tibble(ymin = c(min(up$y), min(down$y)) - 0.35, ymax = c(max(up$y), max(down$y)) + 0.35,
             label = c("Higher share in the\nvariable genome\n(accessory + unique)",
                       "Higher share in the\ncore genome"),
             col = c(col_var, col_core))

p_brk <- ggplot(br) +
  geom_segment(aes(x = 0, xend = 0, y = ymin, yend = ymax, colour = col), linewidth = 1) +
  geom_segment(aes(x = 0, xend = -0.25, y = ymin, yend = ymin, colour = col), linewidth = 1) +
  geom_segment(aes(x = 0, xend = -0.25, y = ymax, yend = ymax, colour = col), linewidth = 1) +
  geom_text(aes(x = 0.35, y = (ymin + ymax) / 2, label = label, colour = col),
            hjust = 0, size = 3, lineheight = 0.9) +
  scale_colour_identity() +
  scale_y_continuous(limits = c(0.4, length(cat_order) + 0.6), expand = c(0, 0)) +
  scale_x_continuous(limits = c(-0.3, 3), expand = c(0, 0)) +
  coord_cartesian(clip = "off") +
  theme_void()

# match the bracket panel's y-range to the facet panels
p_main <- p_main + scale_y_discrete(expand = c(0, 0.6))
fig <- p_main + p_brk + plot_layout(widths = c(1, 0.2)) &
  theme(plot.margin = margin(5, 5, 5, 5))

ggsave(file.path(base_dir, "cog_enrichment_variable_vs_core.pdf"), fig, width = 13, height = 7.5, useDingbats = FALSE)
ggsave(file.path(base_dir, "cog_enrichment_variable_vs_core.png"), fig, width = 13, height = 7.5, dpi = 300)
