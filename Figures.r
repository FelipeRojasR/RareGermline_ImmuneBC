# Figures 

# Baseline Figure
figure_paper.baseline <- full_primary_results %>%
    filter(model = "baseline", cases = "overall") %>%
  mutate(ci.lb = case_when(ci.lb <= 0.3 ~ 0.3, TRUE ~ ci.lb),
         ci.ub = case_when(ci.ub >= 5 ~ 5, TRUE ~ ci.ub),
         gene = factor(gene, levels = unique(gene_and_functions %>% arrange(mechanism) %>% .$gene)), 
         variants = factor(variants, levels = c("missense", "truncating")),
         sig.pval = TRUE,
         # sig.pval = case_when(p.value <= 0.05 ~ TRUE, TRUE ~ FALSE),
         marker = factor(marker, levels = c("CD8", "CD20", "FOXP3", "CD163"))) %>% 
  ggplot(mapping = aes(x = estimate, y = fct_rev(gene), color = fct_rev(marker))) +
  theme_classic() + 
  geom_point(position = position_dodge2(width = 0.7)) + 
  facet_grid(cols = vars(variants), scales = "free") + 
  geom_linerange(aes(xmin = ci.lb, xmax = ci.ub), position = position_dodge2(width = 0.7)) + 
  geom_vline(xintercept = 1, linetype = "dashed") + 
  ####
  geom_vline(xintercept = 0.3, linetype = "solid") + 
  geom_vline(xintercept = 5, linetype = "solid") + 
  ####
  labs(x = "Odds Ratio (95% CI)", y = "") + 
  guides(color = guide_legend(reverse = TRUE)) + 
  coord_trans(x = 'log10') + 
  scale_color_manual(values=c("#E99999", "#C80000", "#80A2C4",  "#01458A")) +
  scale_x_continuous(breaks = c(0.5,1,2,5,10,20)) + 
  theme(legend.position = "top")

ggsave(figure_paper.baseline, device = "pdf", height = 10, width = 6, 
       filename = "/DATA/users/f.rojas/Projects/immune_infiltration_germline/figures/figure.paper.baseline.pdf", 
       dpi = 300)

ggsave(figure_paper.baseline, device = "eps", height = 10, width = 6, 
       filename = "/DATA/users/f.rojas/Projects/immune_infiltration_germline/figures/figure.paper.baseline.eps", 
       dpi = 600)


# Adjusted model

adjusted.joined.plot <- adjusted.joined.rbind %>%
  mutate(gene = factor(gene, levels = unique(gene_and_functions %>% arrange(mechanism) %>% .$gene)), 
         variants = factor(variants, levels = c("missense", "truncating")),
         marker = factor(marker, levels = c("CD8", "CD20", "FOXP3", "CD163")),
         model = factor(model, levels = c("unadjusted", "adjusted")),
         ci.lb = case_when(ci.lb <= 0.3 ~ 0.3, TRUE ~ ci.lb),
         ci.ub = case_when(ci.ub >= 6 ~ 6, TRUE ~ ci.ub)) %>% 
  ggplot(mapping = aes(x = estimate, y = fct_rev(gene), color = fct_rev(marker), shape = fct_rev(model))) + 
  theme_classic() + 
  geom_point(position = position_dodge2(width = 0.6)) + 
  facet_grid(rows = vars(variants), scales = "free", space = "free") +
  geom_linerange(aes(xmin = ci.lb, xmax = ci.ub), position = position_dodge2(width = 0.6)) + 
  ####
  geom_vline(xintercept = 1, linetype = "dashed") + 
  geom_vline(xintercept = 0.3, linetype = "solid") + 
  geom_vline(xintercept = 6, linetype = "solid") + 
  ####
  labs(x = "Odds Ratio (95% CI)", y = "") + 
  guides(color = guide_legend(reverse = TRUE)) + 
  coord_trans(x = 'log10') + 
  scale_color_manual(values=c("#E99999", "#C80000", "#80A2C4",  "#01458A")) +
  scale_x_continuous(breaks = c(0.5,1,2,5)) + 
  theme(legend.position = "top")

ggsave(adjusted.joined.plot, device = "pdf", height = 5.5, width = 4, 
       filename = "/DATA/users/f.rojas/Projects/immune_infiltration_germline/figures/significant.adjusted.pdf", 
       dpi = 300)

ggsave(adjusted.joined.plot, device = "eps", height = 5.5, width = 4, 
       filename = "/DATA/users/f.rojas/Projects/immune_infiltration_germline/figures/significant.adjusted.eps", 
       dpi = 600)





# Triple negative figure
figure_paper.TNBC <- full_primary_results %>%
  filter(model = "TNBC", cases = "overall")  %>%
  mutate(ci.lb = case_when(ci.lb <= 0.3 ~ 0.3, TRUE ~ ci.lb),
         ci.ub = case_when(ci.ub >= 6 ~ 6, TRUE ~ ci.ub),
         gene = factor(gene, levels = unique(gene_and_functions %>% arrange(mechanism) %>% .$gene)), 
         variants = factor(variants, levels = c("missense", "truncating")),
         marker = factor(marker, levels = c("CD8", "CD20", "FOXP3", "CD163"))) %>% 
  ggplot(mapping = aes(x = estimate, y = fct_rev(gene), color = fct_rev(marker))) +
  theme_classic() + 
  geom_point(position = position_dodge2(width = 0.7)) + 
  facet_grid(cols = vars(variants), scales = "free") + 
  geom_linerange(aes(xmin = ci.lb, xmax = ci.ub), position = position_dodge2(width = 0.7)) + 
  geom_vline(xintercept = 1, linetype = "dashed") + 
  ####
  geom_vline(xintercept = 0.3, linetype = "solid") + 
  geom_vline(xintercept = 6, linetype = "solid") + 
  ####
  labs(x = "Odds Ratio (95% CI)", y = "") + 
  guides(color = guide_legend(reverse = TRUE)) + 
  coord_trans(x = 'log10') + 
  scale_color_manual(values=c("#E99999", "#C80000", "#80A2C4",  "#01458A")) +
  scale_x_continuous(breaks = c(0.5,1,2,5,10,20)) + 
  theme(legend.position = "top")

ggsave(figure_paper.TNBC, device = "pdf", height = 10, width = 6, 
       filename = "/DATA/users/f.rojas/Projects/immune_infiltration_germline/figures/figure.paper.TNBC.pdf", 
       dpi = 300)

ggsave(figure_paper.TNBC, device = "eps", height = 10, width = 6, 
       filename = "/DATA/users/f.rojas/Projects/immune_infiltration_germline/figures/figure.paper.TNBC.eps", 
       dpi = 600)


# MECHANISMS
figure.paper.mechanisms <- full_mechanisms_results  %>%
  mutate(mechanism = factor(gene, levels = unique(gene_and_functions %>% arrange(mechanism) %>% .$mechanism)), 
         variants = factor(variants, levels = c("missense", "truncating")),
         marker = factor(marker, levels = c("CD8", "CD20", "FOXP3", "CD163"))) %>% 
  ggplot(mapping = aes(x = estimate, y = fct_rev(mechanism), color = fct_rev(marker))) +
  theme_classic() + 
  geom_point(position = position_dodge2(width = 0.7)) + 
  facet_grid(cols = vars(variants)) + 
  scale_alpha_discrete(range = c(0.55, 1)) + 
  geom_linerange(aes(xmin = ci.lb, xmax = ci.ub), position = position_dodge2(width = 0.7)) + 
  geom_vline(xintercept = 1, linetype = "dashed") + 
  labs(x = "Odds Ratio (95% CI)", y = "") + 
  guides(color = guide_legend(reverse = TRUE)) + 
  coord_trans(x = 'log10') + 
  scale_color_manual(values=c("#E99999", "#C80000", "#80A2C4",  "#01458A")) +
  scale_x_continuous(breaks = c(0.5,1,2,5,10,20)) + 
  theme(legend.position = "top") + 
  ggtitle("Mechanisms baseline")

ggsave(figure.paper.mechanisms, device = "pdf", height = 3.5, width = 7, 
       filename = "/DATA/users/f.rojas/Projects/immune_infiltration_germline/figures/figure.paper.mechanisms.pdf", 
       dpi = 300)

ggsave(figure.paper.mechanisms, device = "eps", height = 3.5, width = 7, 
       filename = "/DATA/users/f.rojas/Projects/immune_infiltration_germline/figures/figure.paper.mechanisms.eps", 
       dpi = 600)

