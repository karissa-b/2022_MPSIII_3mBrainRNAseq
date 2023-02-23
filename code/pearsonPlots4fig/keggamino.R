a <- toptabs_cqn_all %>%
      dplyr::filter(gene_id %in% KEGG$KEGG_ANTIGEN_PROCESSING_AND_PRESENTATION) %>%
      dplyr::filter(coef != "sgsh/+") %>%
      dplyr::select(gene_name, coef_mps, logFC) %>%
      spread(key = "coef_mps", value = "logFC") %>%
      ggscatter(
        x = "MPS-IIIA_A and C",
        y = "MPS-IIIC_A and C",
        alpha = 0.9,
        add = "reg.line",  # Add regressin line
        add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
        conf.int = TRUE # Add confidence interval
      ) +
      stat_cor(
        method = "pearson", label.y = 0.9
      ) +
  labs(x = "MPS-IIIA (AC)",
       y = "MPS-IIIC (AC)",
       title = "A vs C") +
  theme(aspect.ratio = 1,
       axis.title.x = element_text(color = c("#ff005e"), size = 14, face = "bold"),
       axis.title.y = element_text(color = c("#00b064"), size = 14, face = "bold")
       )

b <-toptabs_cqn_all %>%
  dplyr::filter(gene_id %in% KEGG$KEGG_ANTIGEN_PROCESSING_AND_PRESENTATION) %>%
  dplyr::filter(coef != "sgsh/+") %>%
  dplyr::select(gene_name, coef_mps, logFC) %>%
  spread(key = "coef_mps", value = "logFC") %>%
  ggscatter(
    x = "MPS-IIIA_A and C",
    y = "MPS-IIIA_A and B",
    alpha = 0.9,
    add = "reg.line",  # Add regressin line
    add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
    conf.int = TRUE # Add confidence interval
  ) +
  stat_cor(
    method = "pearson", label.y = 0.9
  ) +
  labs(x = "MPS-IIIA (AC)",
       y = "MPS-IIIA (AB)",
       title = "A vs A") +
  theme(aspect.ratio = 1,
        axis.title.x = element_text(color = c("#ff005e"), size = 14, face = "bold"),
        axis.title.y = element_text(color = c("#ff005e"), size = 14, face = "bold")
  )

c <- toptabs_cqn_all %>%
  dplyr::filter(gene_id %in% KEGG$KEGG_ANTIGEN_PROCESSING_AND_PRESENTATION) %>%
  dplyr::filter(coef != "sgsh/+") %>%
  dplyr::select(gene_name, coef_mps, logFC) %>%
  spread(key = "coef_mps", value = "logFC") %>%
  ggscatter(
    x = "MPS-IIIC_A and C",
    y = "MPS-IIIC_B and C",
    alpha = 0.9,
    add = "reg.line",  # Add regressin line
    add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
    conf.int = TRUE # Add confidence interval
  ) +
  stat_cor(
    method = "pearson", label.y = 0.9
  ) +
  labs(x = "MPS-IIIC (AC)",
       y = "MPS-IIIC (BC)",
       title = "C vs C") +
  theme(aspect.ratio = 1,
        axis.title.x = element_text(color = c("#00b064"), size = 14, face = "bold"),
        axis.title.y = element_text(color = c("#00b064"), size = 14, face = "bold")
  )

d <- toptabs_cqn_all %>%
  dplyr::filter(gene_id %in% KEGG$KEGG_GLYCOSAMINOGLYCAN_DEGRADATION) %>%
  dplyr::filter(coef != "sgsh/+") %>%
  dplyr::select(gene_name, coef_mps, logFC) %>%
  spread(key = "coef_mps", value = "logFC") %>%
  ggscatter(
    x = "MPS-IIIA_A and C",
    y = "MPS-IIIC_A and C",
    alpha = 0.9,
    add = "reg.line",  # Add regressin line
    add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
    conf.int = TRUE # Add confidence interval
  ) +
  stat_cor(
    method = "pearson", label.y = 0.9
  ) +
  labs(x = "MPS-IIIA (AC)",
       y = "MPS-IIIC (AC)",
       title = "A vs C") +
  theme(aspect.ratio = 1,
        axis.title.x = element_text(color = c("#ff005e"), size = 14, face = "bold"),
        axis.title.y = element_text(color = c("#00b064"), size = 14, face = "bold")
  )

e <-toptabs_cqn_all %>%
  dplyr::filter(gene_id %in% KEGG$KEGG_GLYCOSAMINOGLYCAN_DEGRADATION) %>%
  dplyr::filter(coef != "sgsh/+") %>%
  dplyr::select(gene_name, coef_mps, logFC) %>%
  spread(key = "coef_mps", value = "logFC") %>%
  ggscatter(
    x = "MPS-IIIA_A and C",
    y = "MPS-IIIA_A and B",
    alpha = 0.9,
    add = "reg.line",  # Add regressin line
    add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
    conf.int = TRUE # Add confidence interval
  ) +
  stat_cor(
    method = "pearson", label.y = 1.2
  ) +
  labs(x = "MPS-IIIA (AC)",
       y = "MPS-IIIA (AB)",
       title = "A vs A") +
  theme(aspect.ratio = 1,
        axis.title.x = element_text(color = c("#ff005e"), size = 14, face = "bold"),
        axis.title.y = element_text(color = c("#ff005e"), size = 14, face = "bold")
  )

f <- toptabs_cqn_all %>%
  dplyr::filter(gene_id %in% KEGG$KEGG_GLYCOSAMINOGLYCAN_DEGRADATION) %>%
  dplyr::filter(coef != "sgsh/+") %>%
  dplyr::select(gene_name, coef_mps, logFC) %>%
  spread(key = "coef_mps", value = "logFC") %>%
  ggscatter(
    x = "MPS-IIIC_A and C",
    y = "MPS-IIIC_B and C",
    alpha = 0.9,
    add = "reg.line",  # Add regressin line
    add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
    conf.int = TRUE # Add confidence interval
  ) +
  stat_cor(
    method = "pearson", label.y = 0.9
  ) +
  labs(x = "MPS-IIIC (AC)",
       y = "MPS-IIIC (BC)",
       title = "C vs C") +
  theme(aspect.ratio = 1,
        axis.title.x = element_text(color = c("#00b064"), size = 14, face = "bold"),
        axis.title.y = element_text(color = c("#00b064"), size = 14, face = "bold")
  )

g <- toptabs_cqn_all %>%
  dplyr::filter(gene_id %in% KEGG$KEGG_OTHER_GLYCAN_DEGRADATION) %>%
  dplyr::filter(coef != "sgsh/+") %>%
  dplyr::select(gene_name, coef_mps, logFC) %>%
  spread(key = "coef_mps", value = "logFC") %>%
  ggscatter(
    x = "MPS-IIIA_A and C",
    y = "MPS-IIIC_A and C",
    alpha = 0.9,
    add = "reg.line",  # Add regressin line
    add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
    conf.int = TRUE # Add confidence interval
  ) +
  stat_cor(
    method = "pearson", label.y = 0.9
  ) +
  labs(x = "MPS-IIIA (AC)",
       y = "MPS-IIIC (AC)",
       title = "A vs C") +
  theme(aspect.ratio = 1,
        axis.title.x = element_text(color = c("#ff005e"), size = 14, face = "bold"),
        axis.title.y = element_text(color = c("#00b064"), size = 14, face = "bold")
  )

h <-toptabs_cqn_all %>%
  dplyr::filter(gene_id %in% KEGG$KEGG_OTHER_GLYCAN_DEGRADATION) %>%
  dplyr::filter(coef != "sgsh/+") %>%
  dplyr::select(gene_name, coef_mps, logFC) %>%
  spread(key = "coef_mps", value = "logFC") %>%
  ggscatter(
    x = "MPS-IIIA_A and C",
    y = "MPS-IIIA_A and B",
    alpha = 0.9,
    add = "reg.line",  # Add regressin line
    add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
    conf.int = TRUE # Add confidence interval
  ) +
  stat_cor(
    method = "pearson", label.y = 0.9
  ) +
  labs(x = "MPS-IIIA (AC)",
       y = "MPS-IIIA (AB)",
       title = "A vs A") +
  theme(aspect.ratio = 1,
        axis.title.x = element_text(color = c("#ff005e"), size = 14, face = "bold"),
        axis.title.y = element_text(color = c("#ff005e"), size = 14, face = "bold")
  )

i <- toptabs_cqn_all %>%
  dplyr::filter(gene_id %in% KEGG$KEGG_OTHER_GLYCAN_DEGRADATION) %>%
  dplyr::filter(coef != "sgsh/+") %>%
  dplyr::select(gene_name, coef_mps, logFC) %>%
  spread(key = "coef_mps", value = "logFC") %>%
  ggscatter(
    x = "MPS-IIIC_A and C",
    y = "MPS-IIIC_B and C",
    alpha = 0.9,
    add = "reg.line",  # Add regressin line
    add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
    conf.int = TRUE # Add confidence interval
  ) +
  stat_cor(
    method = "pearson", label.y = 0.9
  ) +
  labs(x = "MPS-IIIC (AC)",
       y = "MPS-IIIC (BC)",
       title = "C vs C") +
  theme(aspect.ratio = 1,
        axis.title.x = element_text(color = c("#00b064"), size = 14, face = "bold"),
        axis.title.y = element_text(color = c("#00b064"), size = 14, face = "bold")
  )


ggarrange(a, b, c, d, e, f, g, h, i, labels = LETTERS[4:12] ) +

   # (a | b | c) / ( d | e | f) / (g | h | i) +
   #  plot_layout(guides= "collect", row_spacing = unit(0.5, "cm"))
  ggsave("output/plots/lysoglyootherpearsoncors.png",
         width = 5, height = 6, units = "cm", dpi = 500, scale = 5)

# KEGG_ANTIGEN_PROCESSING_AND_PRESENTATION --------------------------------
a <- toptabs_cqn_all %>%
  dplyr::filter(gene_id %in% KEGG$KEGG_ANTIGEN_PROCESSING_AND_PRESENTATION) %>%
  dplyr::filter(coef != "sgsh/+") %>%
  dplyr::select(gene_name, coef_mps, logFC) %>%
  spread(key = "coef_mps", value = "logFC") %>%
  ggscatter(
    x = "MPS-IIIA_A and C",
    y = "MPS-IIIC_A and C",
    alpha = 0.9,
    add = "reg.line",  # Add regressin line
    add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
    conf.int = TRUE # Add confidence interval
  ) +
  stat_cor(
    method = "pearson", label.y = 0.9
  ) +
  labs(x = "MPS-IIIA (AC)",
       y = "MPS-IIIC (AC)",
       title = "A vs C") +
  theme(aspect.ratio = 1,
        axis.title.x = element_text(color = c("#ff005e"), size = 14, face = "bold"),
        axis.title.y = element_text(color = c("#00b064"), size = 14, face = "bold")
  )

b <-toptabs_cqn_all %>%
  dplyr::filter(gene_id %in% KEGG$KEGG_ANTIGEN_PROCESSING_AND_PRESENTATION) %>%
  dplyr::filter(coef != "sgsh/+") %>%
  dplyr::select(gene_name, coef_mps, logFC) %>%
  spread(key = "coef_mps", value = "logFC") %>%
  ggscatter(
    x = "MPS-IIIA_A and C",
    y = "MPS-IIIA_A and B",
    alpha = 0.9,
    add = "reg.line",  # Add regressin line
    add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
    conf.int = TRUE # Add confidence interval
  ) +
  stat_cor(
    method = "pearson", label.y = 0.9
  ) +
  labs(x = "MPS-IIIA (AC)",
       y = "MPS-IIIA (AB)",
       title = "A vs A") +
  theme(aspect.ratio = 1,
        axis.title.x = element_text(color = c("#ff005e"), size = 14, face = "bold"),
        axis.title.y = element_text(color = c("#ff005e"), size = 14, face = "bold")
  )

c <- toptabs_cqn_all %>%
  dplyr::filter(gene_id %in% KEGG$KEGG_ANTIGEN_PROCESSING_AND_PRESENTATION) %>%
  dplyr::filter(coef != "sgsh/+") %>%
  dplyr::select(gene_name, coef_mps, logFC) %>%
  spread(key = "coef_mps", value = "logFC") %>%
  ggscatter(
    x = "MPS-IIIC_A and C",
    y = "MPS-IIIC_B and C",
    alpha = 0.9,
    add = "reg.line",  # Add regressin line
    add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
    conf.int = TRUE # Add confidence interval
  ) +
  stat_cor(
    method = "pearson", label.y = 0.9
  ) +
  labs(x = "MPS-IIIC (AC)",
       y = "MPS-IIIC (BC)",
       title = "C vs C") +
  theme(aspect.ratio = 1,
        axis.title.x = element_text(color = c("#00b064"), size = 14, face = "bold"),
        axis.title.y = element_text(color = c("#00b064"), size = 14, face = "bold")
  )

ggarrange(a, b, c,
          labels = LETTERS[4:7] ,
          nrow = 1) +
  ggsave("output/plots/pearsonCorKEGG_ANTIGEN_PROCESSING_AND_PRESENTATION.png",
         width = 5, height = 3, units = "cm", dpi = 500, scale = 5)


# KEGG_COMPLEMENT ---------------------------------------------------------


a <- toptabs_cqn_all %>%
  dplyr::filter(gene_id %in% KEGG$KEGG_COMPLEMENT_AND_COAGULATION_CASCADES) %>%
  dplyr::filter(coef != "sgsh/+") %>%
  dplyr::select(gene_name, coef_mps, logFC) %>%
  dplyr::distinct(gene_name, coef_mps, .keep_all = T) %>%
  spread(key = "coef_mps", value = "logFC") %>%
  ggscatter(
    x = "MPS-IIIA_A and C",
    y = "MPS-IIIC_A and C",
    alpha = 0.9,
    add = "reg.line",  # Add regressin line
    add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
    conf.int = TRUE # Add confidence interval
  ) +
  stat_cor(
    method = "pearson", label.y = 0.9
  ) +
  labs(x = "MPS-IIIA (AC)",
       y = "MPS-IIIC (AC)",
       title = "A vs C") +
  theme(aspect.ratio = 1,
        axis.title.x = element_text(color = c("#ff005e"), size = 14, face = "bold"),
        axis.title.y = element_text(color = c("#00b064"), size = 14, face = "bold")
  )


b <- toptabs_cqn_all %>%
  dplyr::filter(gene_id %in% KEGG$KEGG_COMPLEMENT_AND_COAGULATION_CASCADES) %>%
  dplyr::filter(coef != "sgsh/+") %>%
  dplyr::select(gene_name, coef_mps, logFC) %>%
  dplyr::distinct(gene_name, coef_mps, .keep_all = T) %>%
  spread(key = "coef_mps", value = "logFC") %>%
  ggscatter(
    x = "MPS-IIIA_A and C",
    y = "MPS-IIIA_A and B",
    alpha = 0.9,
    add = "reg.line",  # Add regressin line
    add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
    conf.int = TRUE # Add confidence interval
  ) +
  stat_cor(
    method = "pearson", label.y = 0.9
  ) +
  labs(x = "MPS-IIIA (AC)",
       y = "MPS-IIIA (AB)",
       title = "A vs A") +
  theme(aspect.ratio = 1,
        axis.title.x = element_text(color = c("#ff005e"), size = 14, face = "bold"),
        axis.title.y = element_text(color = c("#ff005e"), size = 14, face = "bold")
  )

c <- toptabs_cqn_all %>%
  dplyr::filter(gene_id %in% KEGG$KEGG_COMPLEMENT_AND_COAGULATION_CASCADES) %>%
  dplyr::filter(coef != "sgsh/+") %>%
  dplyr::select(gene_name, coef_mps, logFC) %>%
  dplyr::distinct(gene_name, coef_mps, .keep_all = T) %>%
  spread(key = "coef_mps", value = "logFC") %>%
  ggscatter(
    x = "MPS-IIIC_A and C",
    y = "MPS-IIIC_B and C",
    alpha = 0.9,
    add = "reg.line",  # Add regressin line
    add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
    conf.int = TRUE # Add confidence interval
  ) +
  stat_cor(
    method = "pearson", label.y = 0.9
  ) +
  labs(x = "MPS-IIIC (AC)",
       y = "MPS-IIIC (BC)",
       title = "C vs C") +
  theme(aspect.ratio = 1,
        axis.title.x = element_text(color = c("#00b064"), size = 14, face = "bold"),
        axis.title.y = element_text(color = c("#00b064"), size = 14, face = "bold")
  )

ggarrange(a, b, c,
          labels = LETTERS[4:7] ,
          nrow = 1) +
  ggsave("output/plots/pearsonCorKEGG_COMPLEMENT_AND_COAGULATION_CASCADES.png",
         width = 5, height = 3, units = "cm", dpi = 500, scale = 5)


# KEGG_AMINOSUG -----------------------------------------------------------
a <- toptabs_cqn_all %>%
  dplyr::filter(gene_id %in% KEGG$KEGG_AMINO_SUGAR_AND_NUCLEOTIDE_SUGAR_METABOLISM) %>%
  dplyr::filter(coef != "sgsh/+") %>%
  dplyr::select(gene_name, coef_mps, logFC) %>%
  dplyr::distinct(gene_name, coef_mps, .keep_all = T) %>%
  spread(key = "coef_mps", value = "logFC") %>%
  ggscatter(
    x = "MPS-IIIA_A and C",
    y = "MPS-IIIC_A and C",
    alpha = 0.9,
    add = "reg.line",  # Add regressin line
    add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
    conf.int = TRUE # Add confidence interval
  ) +
  stat_cor(
    method = "pearson", label.y = 0.9
  ) +
  labs(x = "MPS-IIIA (AC)",
       y = "MPS-IIIC (AC)",
       title = "A vs C") +
  theme(aspect.ratio = 1,
        axis.title.x = element_text(color = c("#ff005e"), size = 14, face = "bold"),
        axis.title.y = element_text(color = c("#00b064"), size = 14, face = "bold")
  )


b <- toptabs_cqn_all %>%
  dplyr::filter(gene_id %in% KEGG$KEGG_AMINO_SUGAR_AND_NUCLEOTIDE_SUGAR_METABOLISM) %>%
  dplyr::filter(coef != "sgsh/+") %>%
  dplyr::select(gene_name, coef_mps, logFC) %>%
  dplyr::distinct(gene_name, coef_mps, .keep_all = T) %>%
  spread(key = "coef_mps", value = "logFC") %>%
  ggscatter(
    x = "MPS-IIIA_A and C",
    y = "MPS-IIIA_A and B",
    alpha = 0.9,
    add = "reg.line",  # Add regressin line
    add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
    conf.int = TRUE # Add confidence interval
  ) +
  stat_cor(
    method = "pearson", label.y = 0.9
  ) +
  labs(x = "MPS-IIIA (AC)",
       y = "MPS-IIIA (AB)",
       title = "A vs A") +
  theme(aspect.ratio = 1,
        axis.title.x = element_text(color = c("#ff005e"), size = 14, face = "bold"),
        axis.title.y = element_text(color = c("#ff005e"), size = 14, face = "bold")
  )

c <- toptabs_cqn_all %>%
  dplyr::filter(gene_id %in% KEGG$KEGG_AMINO_SUGAR_AND_NUCLEOTIDE_SUGAR_METABOLISM) %>%
  dplyr::filter(coef != "sgsh/+") %>%
  dplyr::select(gene_name, coef_mps, logFC) %>%
  dplyr::distinct(gene_name, coef_mps, .keep_all = T) %>%
  spread(key = "coef_mps", value = "logFC") %>%
  ggscatter(
    x = "MPS-IIIC_A and C",
    y = "MPS-IIIC_B and C",
    alpha = 0.9,
    add = "reg.line",  # Add regressin line
    add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
    conf.int = TRUE # Add confidence interval
  ) +
  stat_cor(
    method = "pearson", label.y = 0.9
  ) +
  labs(x = "MPS-IIIC (AC)",
       y = "MPS-IIIC (BC)",
       title = "C vs C") +
  theme(aspect.ratio = 1,
        axis.title.x = element_text(color = c("#00b064"), size = 14, face = "bold"),
        axis.title.y = element_text(color = c("#00b064"), size = 14, face = "bold")
  )

ggarrange(a, b, c,
          nrow = 1) +
  ggsave("output/plots/pearsonCorKEGG_AMINO_SUGAR_AND_NUCLEOTIDE_SUGAR_METABOLISM.png",
         width = 5, height = 3, units = "cm", dpi = 500, scale = 5)

# FATTY_ACID --------------------------------------------------------------
a <- toptabs_cqn_all %>%
  dplyr::filter(gene_id %in% KEGG$KEGG_FATTY_ACID_METABOLISM) %>%
  dplyr::filter(coef != "sgsh/+") %>%
  dplyr::select(gene_name, coef_mps, logFC) %>%
  dplyr::distinct(gene_name, coef_mps, .keep_all = T) %>%
  spread(key = "coef_mps", value = "logFC") %>%
  ggscatter(
    x = "MPS-IIIA_A and C",
    y = "MPS-IIIC_A and C",
    alpha = 0.9,
    add = "reg.line",  # Add regressin line
    add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
    conf.int = TRUE # Add confidence interval
  ) +
  stat_cor(
    method = "pearson", label.y = 0.9
  ) +
  labs(x = "MPS-IIIA (AC)",
       y = "MPS-IIIC (AC)",
       title = "A vs C") +
  theme(aspect.ratio = 1,
        axis.title.x = element_text(color = c("#ff005e"), size = 14, face = "bold"),
        axis.title.y = element_text(color = c("#00b064"), size = 14, face = "bold")
  )


b <- toptabs_cqn_all %>%
  dplyr::filter(gene_id %in% KEGG$KEGG_AMINO_SUGAR_AND_NUCLEOTIDE_SUGAR_METABOLISM) %>%
  dplyr::filter(coef != "sgsh/+") %>%
  dplyr::select(gene_name, coef_mps, logFC) %>%
  dplyr::distinct(gene_name, coef_mps, .keep_all = T) %>%
  spread(key = "coef_mps", value = "logFC") %>%
  ggscatter(
    x = "MPS-IIIA_A and C",
    y = "MPS-IIIA_A and B",
    alpha = 0.9,
    add = "reg.line",  # Add regressin line
    add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
    conf.int = TRUE # Add confidence interval
  ) +
  stat_cor(
    method = "pearson", label.y = 0.9
  ) +
  labs(x = "MPS-IIIA (AC)",
       y = "MPS-IIIA (AB)",
       title = "A vs A") +
  theme(aspect.ratio = 1,
        axis.title.x = element_text(color = c("#ff005e"), size = 14, face = "bold"),
        axis.title.y = element_text(color = c("#ff005e"), size = 14, face = "bold")
  )

c <- toptabs_cqn_all %>%
  dplyr::filter(gene_id %in% KEGG$KEGG_AMINO_SUGAR_AND_NUCLEOTIDE_SUGAR_METABOLISM) %>%
  dplyr::filter(coef != "sgsh/+") %>%
  dplyr::select(gene_name, coef_mps, logFC) %>%
  dplyr::distinct(gene_name, coef_mps, .keep_all = T) %>%
  spread(key = "coef_mps", value = "logFC") %>%
  ggscatter(
    x = "MPS-IIIC_A and C",
    y = "MPS-IIIC_B and C",
    alpha = 0.9,
    add = "reg.line",  # Add regressin line
    add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
    conf.int = TRUE # Add confidence interval
  ) +
  stat_cor(
    method = "pearson", label.y = 0.9
  ) +
  labs(x = "MPS-IIIC (AC)",
       y = "MPS-IIIC (BC)",
       title = "C vs C") +
  theme(aspect.ratio = 1,
        axis.title.x = element_text(color = c("#00b064"), size = 14, face = "bold"),
        axis.title.y = element_text(color = c("#00b064"), size = 14, face = "bold")
  )

ggarrange(a, b, c,
          nrow = 1) +
  ggsave("output/plots/pearsonCorKEGG_AMINO_SUGAR_AND_NUCLEOTIDE_SUGAR_METABOLISM.png",
         width = 5, height = 3, units = "cm", dpi = 500, scale = 5)

