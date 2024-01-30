# Distribution of scores
for_distribution <- TilScores_BridgesPhenotype_markers %>% 
  select(BCAC_ID, percentage, BRIDGES_ID, marker) %>%
  mutate(percentage = case_when(percentage == 0.00001 ~ 0, TRUE ~ percentage),
         marker = factor(marker, levels = c("CD8", "CD20", "FOXP3", "CD163"))) %>%
  ungroup()

figure.paper.distribution.raw <- for_distribution %>%
  ggplot(mapping = aes(x = percentage, color = marker)) + 
  geom_freqpoly(binwidth = 0.01) +
  theme_classic() + 
  scale_color_manual(values=c("#01458A", "#80A2C4", "#C80000", "#E99999")) + 
  theme(legend.position = "top") + xlim(0,1) +
  labs(y = "Count", x = "Percentage")

figure.paper.distribution <- for_distribution %>%
  ggplot(mapping = aes(x = logit(percentage), color = marker)) + 
  geom_freqpoly(binwidth = 0.5) + 
  theme_classic() + 
  scale_color_manual(values=c("#01458A", "#80A2C4", "#C80000", "#E99999")) + 
  theme(legend.position = "top") + ylim(0,2500) + xlim(-12,0) +
  labs(y = "Count", x = "Logit(Percentage)")

figure.paper.distribution.comb <- figure.paper.distribution.raw + figure.paper.distribution + 
  plot_annotation(tag_levels = 'A')

# ggsave(figure.paper.distribution.comb, device = "pdf", height = 4, width = 9, dpi = 300,
#        filename = "/DATA/users/f.rojas/Projects/immune_infiltration_germline/figures/figure.paper.distribution.pdf")

# ggsave(figure.paper.distribution.comb, device = "eps", height = 4, width = 9, dpi = 300,
#        filename = "/DATA/users/f.rojas/Projects/immune_infiltration_germline/figures/figure.paper.distribution.eps")

