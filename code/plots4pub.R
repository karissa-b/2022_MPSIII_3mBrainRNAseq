# PLOTS FOR PUB
library(patchwork)
library(UpSetR)
library(pathview)

toptabs_cqn_all <- toptables_cqn.ab %>%
  bind_rows() %>%
  mutate(mps = "A and B") %>%
  bind_rows(toptables_cqn.ac %>%
              bind_rows() %>%
              mutate(mps = "A and C") ) %>%
  bind_rows(toptables_cqn.bc %>%
              bind_rows() %>%
              mutate(mps = "B and C") ) %>%
  bind_rows(toptables_cqn.b2 %>%
              bind_rows() %>%
              mutate(mps = "B2") ) %>%
  mutate(coef_mps = paste0(coef, "_", mps))

logCPm.all <- logCPM.ac %>%
  as.data.frame() %>%
  rownames_to_column("gene_id") %>%
  left_join(logCPM.ab[,(x.ab$samples %>% dplyr::filter(genotype != "sgsh/+") %>% .$sample)] %>% as.data.frame %>% rownames_to_column("gene_id") ) %>%
  left_join(logCPM.bc %>% as.data.frame %>% rownames_to_column("gene_id") ) %>%
  left_join(logCPM.b2[,(x.b2$samples %>% dplyr::filter(genotype != "EOfAD-like") %>% .$sample)] %>% as.data.frame %>% rownames_to_column("gene_id") )



# custom function to plot spearman cors
makePearsonCorrPlot <- function(geneset, plot.title) {
  plotSpearmanCor <- function(x.axis.coef, y.axis.coef, geneset) {
    toptabs_cqn_all %>%
      dplyr::filter(gene_id %in% geneset) %>%
      dplyr::filter(coef != "sgsh/+") %>%
      dplyr::select(gene_name, coef_mps, logFC) %>%
      spread(key = "coef_mps", value = "logFC") %>%
      ggscatter(
        x = paste(x.axis.coef),
        y = paste(y.axis.coef),
        alpha = 0.9,
        add = "reg.line",  # Add regressin line
        add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
        conf.int = TRUE # Add confidence interval
      ) +
      stat_cor(
        method = "pearson", label.y = 0.9
      )+
      theme(aspect.ratio = 1)
  }

  a <- plotSpearmanCor(x.axis.coef = "MPS-IIIA_A and B",
                       y.axis.coef = "MPS-IIIA_A and C",
                       geneset = geneset)

  # MPS-IIIB
  b <- plotSpearmanCor(x.axis.coef = "MPS-IIIB_A and B",
                       y.axis.coef = "MPS-IIIB_B and C",
                       geneset = geneset)

  # MPS-IIIC
  c <- plotSpearmanCor(x.axis.coef = "MPS-IIIC_A and C",
                       y.axis.coef = "MPS-IIIC_B and C",
                       geneset = geneset)

  # A v B
  d <- plotSpearmanCor(x.axis.coef = "MPS-IIIA_A and B",
                       y.axis.coef = "MPS-IIIB_A and B",
                       geneset = geneset)


  # A v C
  e <- plotSpearmanCor(x.axis.coef = "MPS-IIIA_A and C",
                       y.axis.coef = "MPS-IIIC_A and C",
                       geneset = geneset)

  # B v C
  f <- plotSpearmanCor(x.axis.coef = "MPS-IIIB_B and C",
                       y.axis.coef = "MPS-IIIC_B and C",
                       geneset = geneset)

  (a | b | c) /( d | e | f) +
    plot_layout(guides= "collect") +
    plot_annotation(title = plot.title)
}




# PCA ---------------------------------------------------------------------
#``````` A and B ```````
ggarrange(
  logCPM.ab %>%
    t() %>%
    prcomp() %>%
    autoplot(data = tibble(sample = rownames(.$x)) %>%
               left_join(x.ab$samples),
             colour = "genotype",
             shape = "Sex",
             size = 4
    ) +
    scale_color_manual(values = c("#4000ff", #WT
                                  "#ff7aab", # sgsh het
                                  "#ff005e", # A
                                  "#ffae00" #B
    )) +
    theme(aspect.ratio = 1) +
    labs(title = "A and B",
         colour = "Genotype"),

 #``````` A  and C ```````

logCPM.ac %>%
  t() %>%
  prcomp() %>%
  autoplot(data = tibble(sample = rownames(.$x)) %>%
             left_join(x.ac$samples),
           colour = "genotype",
           shape = "Sex",
           size = 4
  ) +
  scale_color_manual(values = c("#4000ff", #WT
                                "#ff005e", # A
                                "#00b064" # C
  )) +
  theme(aspect.ratio = 1) +
  labs(title = "A and C",
       colour = "Genotype"),

#``````` B and C ```````

logCPM.bc %>%
  t() %>%
  prcomp() %>%
  autoplot(data = tibble(sample = rownames(.$x)) %>%
             left_join(x.bc$samples),
           colour = "genotype",
           shape = "Sex",
           size = 4
  ) +
  scale_color_manual(values = c("#4000ff", #WT
                                "#ffae00", # B
                                "#00b064" # C
  )) +
  theme(aspect.ratio = 1) +
  labs(title = "B and C",
       colour = "Genotype"),
nrow= 1,
labels = "AUTO"
) +
  ggsave("output/plots/PCA_all.png", width = 12, height = 8, units = "cm", dpi = 500, scale = 2.5)

# NUM DE ----------------------------------------------------------------------
 ggarrange(
toptables_cqn.ab %>%
  lapply(function(x){
    x %>%
      mutate(Direction = case_when(
        sign(logFC) == "1" ~ "up",
        sign(logFC) == "-1" ~ "down"
      ))
  }) %>%
  bind_rows() %>%
  dplyr::filter(DE == TRUE) %>%
  ggplot(aes(y = `coef`)) +
  geom_bar(stat = "count",
           width = 0.9,
           aes(fill = Direction),
           position = "dodge") +
  theme(legend.position = "bottom") +
  scale_x_continuous(limits = c(0,150)) +
  labs(title = "A and B",
       x = "Number of DE genes",
       y = "Contrast",
       colour = "Direction of change") +
  theme_pubr() +
  theme(axis.title.y = element_blank()),


toptables_cqn.ac %>%
  lapply(function(x){
    x %>%
      mutate(Direction = case_when(
        sign(logFC) == "1" ~ "up",
        sign(logFC) == "-1" ~ "down"
      ))
  }) %>%
  bind_rows() %>%
  dplyr::filter(DE == TRUE) %>%
  ggplot(aes(y = `coef`)) +
  geom_bar(stat = "count",
           width = 0.9,
           aes(fill = Direction),
           position = "dodge") +

  labs(title = "A and C",
       x = "Number of DE genes",
       y = "Contrast",
       colour = "Direction of change") +
  scale_x_continuous(limits = c(0,150)) +
  theme_pubr() +
  theme(axis.title.y = element_blank()),


toptables_cqn.bc %>%
  lapply(function(x){
    x %>%
      mutate(Direction = case_when(
        sign(logFC) == "1" ~ "up",
        sign(logFC) == "-1" ~ "down"
      ))
  }) %>%
  bind_rows() %>%
  dplyr::filter(DE == TRUE) %>%
  ggplot(aes(y = `coef`)) +
  geom_bar(stat = "count",
           width = 0.9,
           aes(fill = Direction),
           position = "dodge") +
  labs(title = "B and C",
       x = "Number of DE genes",
       y = "Contrast",
       colour = "Direction of change") +
  scale_x_continuous(limits = c(0,150)) +
  theme_pubr() +
  theme(axis.title.y = element_blank()),
common.legend = T,
legend = "bottom",
nrow =1) +
  ggsave("output/plots/numDE.png", width = 12, height = 2, units = "cm", dpi = 500, scale = 2.5)




# Volcano plots -----------------------------------------------------------
ggarrange(

toptables_cqn.ab %>%
  bind_rows() %>%
  ggplot(aes(y = -log10(PValue), x = logFC, colour = DE)) +
  geom_point(
    alpha = 0.5, size = 1.25
  ) +
  facet_wrap(~coef, ncol = 1) +
  coord_cartesian(xlim = c(-2.5,2.5),
                  ylim = c(0, 10) ) +
  theme(legend.position = "bottom") +
  scale_color_manual(values = c("grey50", "red")),

toptables_cqn.ac %>%
  bind_rows() %>%
  ggplot(aes(y = -log10(PValue), x = logFC, colour = DE)) +
  geom_point(
    alpha = 0.5, size = 1.25
  ) +
  facet_wrap(~coef, ncol = 1) +
  coord_cartesian(xlim = c(-2.5,2.5),
                  ylim = c(0, 10) ) +
  theme(legend.position = "bottom") +
  scale_color_manual(values = c("grey50", "red")),

toptables_cqn.bc %>%
  bind_rows() %>%
  ggplot(aes(y = -log10(PValue), x = logFC, colour = DE)) +
  geom_point(
    alpha = 0.5, size = 1.25
  ) +
  facet_wrap(~coef, ncol = 1) +
  coord_cartesian(xlim = c(-2.5,2.5),
                  ylim = c(0, 10) ) +
  theme(legend.position = "bottom") +
  scale_color_manual(values = c("grey50", "red")),
common.legend = T,
nrow = 1,
heights = c(2,1,1)
) +
  ggsave("output/plots/volc_all.png", width = 12, height = 12, units = "cm", dpi = 500, scale = 2.4)



# GSEA heatmap summary --------------------------------------------------------------------

frykegg.all <-
  fryKEGG.ac %>%
  bind_rows() %>%
  mutate(MPS = "AC") %>%
  bind_rows(fryKEGG.ab %>%
              bind_rows() %>%
              mutate(MPS= "AB")) %>%
  bind_rows(fryKEGG.bc %>%
              bind_rows() %>%
              mutate(MPS = "BC"))

sigpaths <-
  frykegg.all %>%
  dplyr::filter(FDR.Mixed < 0.05) %>%
  .$pathway %>%
  unique()


frykegg.all %>%
  bind_rows() %>%
  dplyr::filter(pathway %in% sigpaths) %>%
  dplyr::filter(coef != "sgsh/+") %>%
  mutate(MPS = case_when(MPS == "AC" ~ "A and C",
                         MPS == "AB" ~ "A and B",
                         MPS == "BC" ~ "B and C")) %>%
  mutate(pathway = str_remove(pathway, pattern = "KEGG_")) %>%
  ggplot(aes(x = coef,
             y = pathway )) +
  geom_tile(aes(fill = -log10(PValue.Mixed),
                alpha = FDR.Mixed < 0.05)) +
  geom_label(aes(label = signif(FDR.Mixed, digits = 2)),
             fill = NA) +
  facet_wrap(~MPS, scales = "free_x") +
  scale_fill_viridis_c() +
  scale_alpha_manual(values = c(0.5,1)) +
  ggpubr::theme_pubclean() +
  theme(panel.grid.major.y = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        axis.text.y = element_text(face = "bold", size = 15),
        axis.text.x = element_text(face = "bold", size = 15, angle = 315, hjust = 0.1, vjust = 0.2),
        legend.position = "right"

  )+
  ggsave("output/plots/KEGG_pval_heatmap.png", width = 16, height = 16, units = "cm", dpi = 500, scale = 2)

# overlapping genes in the KEGGs ------------------------------------------
KEGG[sigpaths] %>%
  fromList() %>%
  upset(nsets = length(.),
        empty.intersections = "off",
        order.by = "freq")


# KEGG lysosome HEATMAP -----------------------------------------------------------

png("output/plots/keggLyso_Heatmap.png", width = 70, height = 20, units = "cm", res = 100)
toptabs_cqn_all %>%
  dplyr::filter(gene_id %in% KEGG$KEGG_LYSOSOME) %>%
  dplyr::filter(coef != "sgsh/+") %>%
  mutate(coef_mps = paste0(coef, " (", mps, ")")) %>%
  dplyr::select(gene_name, logFC, coef_mps) %>%
  dplyr::distinct(gene_name, coef_mps, .keep_all = T) %>%
  spread(key = "coef_mps", value = "logFC") %>%
  column_to_rownames("gene_name") %>%
  t() %>%
  pheatmap(color = colorRampPalette(rev(brewer.pal(n = 5,
                                                   name = "RdBu")))(100),
           main = "KEGG_LYSOSOME",
           breaks = c(seq(min(.), 0, length.out=ceiling(100/2) + 1),
                      seq(max(.)/100, max(.), length.out=50)),
           cellheight = 30, cellwidth = 16,
           angle_col = 45,
           treeheight_row = 10,
           fontsize = 15)
dev.off()

png("output/plots/keggotherGlyca_Heatmap.png", width = 70, height = 20, units = "cm", res = 100)
toptabs_cqn_all %>%
  dplyr::filter(gene_id %in% KEGG$KEGG_OTHER_GLYCAN_DEGRADATION) %>%
  dplyr::filter(coef != "sgsh/+") %>%
  mutate(coef_mps = paste0(coef, " (", mps, ")")) %>%
  dplyr::select(gene_name, logFC, coef_mps) %>%
  dplyr::distinct(gene_name, coef_mps, .keep_all = T) %>%
  spread(key = "coef_mps", value = "logFC") %>%
  column_to_rownames("gene_name") %>%
  t() %>%
  pheatmap(color = colorRampPalette(rev(brewer.pal(n = 5,
                                                   name = "RdBu")))(100),
           main = "KEGG_OTHER_GLYCAN_DEGRADATION",
           breaks = c(seq(min(.), 0, length.out=ceiling(100/2) + 1),
                      seq(max(.)/100, max(.), length.out=50)),
           cellheight = 30, cellwidth = 16,
           angle_col = 45,
           treeheight_row = 10,
           fontsize = 15)
dev.off()


png("output/plots/keggglycosdega_Heatmap.png", width = 70, height = 20, units = "cm", res = 100)
toptabs_cqn_all %>%
  dplyr::filter(gene_id %in% KEGG$KEGG_GLYCOSAMINOGLYCAN_DEGRADATION) %>%
  dplyr::filter(coef != "sgsh/+") %>%
  mutate(coef_mps = paste0(coef, " (", mps, ")")) %>%
  dplyr::select(gene_name, logFC, coef_mps) %>%
  dplyr::distinct(gene_name, coef_mps, .keep_all = T) %>%
  spread(key = "coef_mps", value = "logFC") %>%
  column_to_rownames("gene_name") %>%
  t() %>%
  pheatmap(color = colorRampPalette(rev(brewer.pal(n = 5,
                                                   name = "RdBu")))(100),
           main = "KEGG_GLYCOSAMINOGLYCAN_DEGRADATION",
           breaks = c(seq(min(.), 0, length.out=ceiling(100/2) + 1),
                      seq(max(.)/100, max(.), length.out=50)),
           cellheight = 30, cellwidth = 16,
           angle_col = 45,
           treeheight_row = 10,
           fontsize = 15)
dev.off()


# correlations of


# lyso correlations --------------------------------------------------
makePearsonCorrPlot(geneset = KEGG$KEGG_LYSOSOME, plot.title = "KEGG_LYSOSOME") +
  ggsave("output/plots/spearmancorrs_kegglyso.png", width = 5, height = 3, units = "cm", dpi = 500, scale = 5)

# glycosodeg correlations --------------------------------------------------

makePearsonCorrPlot(geneset = KEGG$KEGG_GLYCOSAMINOGLYCAN_DEGRADATION, plot.title = "KEGG_GLYCOSAMINOGLYCAN_DEGRADATION")
  ggsave("output/plots/spearmancorrs_KEGG_GLYCOSAMINOGLYCAN_DEGRADATION.png", width = 5, height = 3, units = "cm", dpi = 500, scale = 5)


# other glycan deg correlations --------------------------------------------------------
  makePearsonCorrPlot(geneset = KEGG$KEGG_OTHER_GLYCAN_DEGRADATION, plot.title = "KEGG_OTHER_GLYCAN_DEGRADATION") +
  ggsave("output/plots/spearmancorrs_KEGG_OTHER_GLYCAN_DEGRADATION.png", width = 5, height = 3, units = "cm", dpi = 500, scale = 5)

# COMPLEMENT_AND_COAGULATION_CASCADES correlations
  makePearsonCorrPlot(geneset = KEGG$KEGG_COMPLEMENT_AND_COAGULATION_CASCADES, plot.title = "KEGG_COMPLEMENT_AND_COAGULATION_CASCADES") +
    ggsave("output/plots/spearmancorrs_KEGG_COMPLEMENT_AND_COAGULATION_CASCADES.png", width = 5, height = 3, units = "cm", dpi = 500, scale = 5)

# AMINO_SUGAR_AND_NUCLEOTIDE_SUGAR_METABOLISM correlations
  makePearsonCorrPlot(geneset = KEGG$KEGG_AMINO_SUGAR_AND_NUCLEOTIDE_SUGAR_METABOLISM, plot.title = "KEGG_AMINO_SUGAR_AND_NUCLEOTIDE_SUGAR_METABOLISM") +
    ggsave("output/plots/spearmancorrs_KEGG_AMINO_SUGAR_AND_NUCLEOTIDE_SUGAR_METABOLISM.png", width = 5, height = 3, units = "cm", dpi = 500, scale = 5)


  # KEGG_VIBRIO_CHOLERAE_INFECTION correlations
  makePearsonCorrPlot(geneset = KEGG$KEGG_VIBRIO_CHOLERAE_INFECTION, plot.title = "KEGG_VIBRIO_CHOLERAE_INFECTION") +
    ggsave("output/plots/pearsoncorrsKEGG_VIBRIO_CHOLERAE_INFECTION.png", width = 5, height = 3, units = "cm", dpi = 500, scale = 5)

  # KEGG_FATTY_ACID_METABOLISM correlations
  makePearsonCorrPlot(geneset = KEGG$KEGG_FATTY_ACID_METABOLISM, plot.title = "KEGG_FATTY_ACID_METABOLISM") +
    ggsave("output/plots/pearsoncorrsKEGG_VIBRIO_CHOLERAE_INFECTION.png", width = 5, height = 3, units = "cm", dpi = 500, scale = 5)


# gsea cell type ----------------------------------------------------------
# cell type heatmap -----------------------------------------------------------------
frycell.all <-
  celltype.ac %>%
  bind_rows() %>%
  mutate(MPS = "A and C") %>%
  bind_rows(celltypes.ab %>%
              bind_rows() %>%
              mutate(MPS= "A and B")) %>%
  bind_rows(celltype.bc %>%
              bind_rows() %>%
              mutate(MPS = "B and C"))

# kegg lyso correlations --------------------------------------------------



frycell.all %>%
  dplyr::filter(coef != "sgsh/+") %>%
  ggplot(aes(x = coef,
             y = pathway)
  ) +
  geom_tile(aes(fill = -log10(PValue),
                alpha = FDR < 0.05),
            colour = "black") +
  geom_label(aes(label = signif(FDR, digits = 2)
  ),
  fill = NA) +
  facet_wrap(~MPS, scales = "free_x") +
  scale_fill_viridis_c() +
  scale_alpha_manual(values = c(0.5,1)) +
  ggpubr::theme_pubclean() +
  theme(panel.grid.major.y = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        axis.text.y = element_text(face = "bold", size = 15),
        axis.text.x = element_text(face = "bold", size = 15, angle = 315, hjust = 0.1, vjust = 0.2)
        ) +
  ggsave("output/plots/celltype_heatmap.png", width = 12, height = 16, units = "cm", dpi = 500, scale = 2)


# correlations ------------------------------------------------------------

makePearsonCorrPlot(geneset = cell_type_markers$Oligodendrocyte, plot.title = "OLIGODENDROCYTES") +
  ggsave("output/plots/Pearsoncorrs_oligoden.png", width = 5, height = 3, units = "cm", dpi = 500, scale = 5)




# Neural stem cell --------------------------------------------------------
makePearsonCorrPlot(geneset = cell_type_markers$`Neural stem cell`, plot.title = "Neural stem cell") +
    ggsave("output/plots/Pearsoncorrs_Neural stem cell.png", width = 5, height = 3, units = "cm", dpi = 500, scale = 5)


# microglia ---------------------------------------------------------------
makePearsonCorrPlot(geneset = cell_type_markers$`Microglia_apoc1 high`, plot.title = "NMicroglia_apoc1 high") +
    ggsave("output/plots/Pearsoncorrs_Microglia_apoc1 high.png", width = 5, height = 3, units = "cm", dpi = 500, scale = 5)


# heatmaps of FC ----------------------------------------------------------

# oligoden ----------------------------------------------------------------

png("output/plots/oligodenHeatmap.png", width = 50, height = 20, units = "cm", res = 100)
toptabs_cqn_all %>%
  dplyr::filter(gene_id %in% cell_type_markers$Oligodendrocyte) %>%
  dplyr::filter(coef != "sgsh/+") %>%
  mutate(coef_mps = paste0(coef, " (", mps, ")")) %>%
  dplyr::select(gene_name, logFC, coef_mps) %>%
  dplyr::distinct(gene_name, coef_mps, .keep_all = T) %>%
  spread(key = "coef_mps", value = "logFC") %>%
  column_to_rownames("gene_name") %>%
  t() %>%
  pheatmap(color = colorRampPalette(rev(brewer.pal(n = 5,
                                                   name = "RdBu")))(100),
           main = "Oligodendrocytes",
           breaks = c(seq(min(.), 0, length.out=ceiling(100/2) + 1),
                      seq(max(.)/100, max(.), length.out=50)),
           cellheight = 30, cellwidth = 12,
           angle_col = 45,
           treeheight_row = 10,
           #treeheight_col = 20
  )
dev.off()



anno_celltype <- cell_type_markers[c("Oligodendrocyte", "Neural stem cell", "Microglia_apoc1 high") ] %>%
  unlist() %>%
  as.data.frame() %>%
  rownames_to_column("pathway") %>%
  set_colnames(c("pathway", "gene_id")) %>%
  mutate(pathway = str_remove(pathway, pattern = "[0-9]+$")) %>%
  as_tibble() %>%
  pivot_wider(names_from = pathway, values_from = pathway) %>%
  left_join(x.ac$genes %>% dplyr::select(gene_id, gene_name)) %>%
  dplyr::select(-gene_id) %>%
  dplyr::distinct(gene_name, .keep_all

                  = T) %>%
  column_to_rownames("gene_name")

# neural stem cell --------------------------------------------------------

png("output/plots/oligodenHeatmap.png", width = 50, height = 20, units = "cm", res = 100)
toptabs_cqn_all %>%
  dplyr::filter(gene_id %in% cell_type_markers$Oligodendrocyte) %>%
  dplyr::filter(coef != "sgsh/+") %>%
  mutate(coef_mps = paste0(coef, " (", mps, ")")) %>%
  dplyr::select(gene_name, logFC, coef_mps) %>%
  dplyr::distinct(gene_name, coef_mps, .keep_all = T) %>%
  spread(key = "coef_mps", value = "logFC") %>%
  column_to_rownames("gene_name") %>%
  t() %>%
  pheatmap(color = colorRampPalette(rev(brewer.pal(n = 5,
                                                   name = "RdBu")))(100),
           main = "Oligodendrocytes",
           annotation_col = anno_celltype,
           breaks = c(seq(min(.), 0, length.out=ceiling(100/2) + 1),
                      seq(max(.)/100, max(.), length.out=50)),
           cellheight = 30, cellwidth = 12,
           angle_col = 45,
           treeheight_row = 10,
           #treeheight_col = 20
  )
dev.off()

# microglia heatnmap --------------------------------------------------------

png("output/plots/oligodenHeatmap.png", width = 50, height = 20, units = "cm", res = 100)
toptabs_cqn_all %>%
  dplyr::filter(gene_id %in% cell_type_markers$`Microglia_apoc1 high`) %>%
  dplyr::filter(coef != "sgsh/+") %>%
  mutate(coef_mps = paste0(coef, " (", mps, ")")) %>%
  dplyr::select(gene_name, logFC, coef_mps) %>%
  dplyr::distinct(gene_name, coef_mps, .keep_all = T) %>%
  spread(key = "coef_mps", value = "logFC") %>%
  column_to_rownames("gene_name") %>%
  t() %>%
  pheatmap(color = colorRampPalette(rev(brewer.pal(n = 5,
                                                   name = "RdBu")))(100),
           main = "Microglia_apoc1 high",
           annotation_col = anno_celltype,
           breaks = c(seq(min(.), 0, length.out=ceiling(100/2) + 1),
                      seq(max(.)/100, max(.), length.out=50)),
           cellheight = 30, cellwidth = 12,
           angle_col = 45,
           treeheight_row = 10,
           #treeheight_col = 20
  )
dev.off()

cell_type_markers[c("Oligodendrocyte", "Neural stem cell", "Microglia_apoc1 high") ] %>%
  fromList() %>%
  upset(
        order.by = "freq",
        nsets = 23,
        mb.ratio = c(0.4, 0.6)
        )

# stem cell heatmap -------------------------------------------------------

png("output/plots/stemcellHeatmap.png", width = 50, height = 20, units = "cm", res = 100)
toptabs_cqn_all %>%
  dplyr::filter(gene_id %in% cell_type_markers$`Neural stem cell`) %>%
  dplyr::filter(coef != "sgsh/+") %>%
  mutate(coef_mps = paste0(coef, " (", mps, ")")) %>%
  dplyr::select(gene_name, logFC, coef_mps) %>%
  dplyr::distinct(gene_name, coef_mps, .keep_all = T) %>%
  spread(key = "coef_mps", value = "logFC") %>%
  column_to_rownames("gene_name") %>%
  t() %>%
  pheatmap(color = colorRampPalette(rev(brewer.pal(n = 5,
                                                   name = "RdBu")))(100),
           main = "Neural stem cell",
           breaks = c(seq(min(.), 0, length.out=ceiling(100/2) + 1),
                      seq(max(.)/100, max(.), length.out=50)),
           cellheight = 30, cellwidth = 12,
           annotation_col = anno_celltype,
           angle_col = 45,
           treeheight_row = 10,
           #treeheight_col = 20
  )
dev.off()

toptables_cqn.ab %>%
  bind_rows() %>%
  dplyr::filter(gene_id %in% cell_type_markers$Oligodendrocyte) %>%
  dplyr::filter(coef != "sgsh/+") %>%
  dplyr::select(gene_name, logFC, coef) %>%
  dplyr::distinct(gene_name, coef, .keep_all = T) %>%
  spread(key = "coef", value = "logFC") %>%
  column_to_rownames("gene_name") %>%
  arrange(desc(`MPS-IIIA`)) %>%
  t() %>%
  pheatmap(color = colorRampPalette(rev(brewer.pal(n = 5,
                                                   name = "RdBu")))(100),
           main = "Oligodendrocytes in A and B",
           annotation_col = anno_celltype,
           breaks = c(seq(min(.), 0, length.out=ceiling(100/2) + 1),
                      seq(max(.)/100, max(.), length.out=50)),
           cellheight = 30, cellwidth = 12,
           angle_col = 45,
           cluster_rows = F,
           treeheight_row = 0, treeheight_col = 0
  )

toptables_cqn.ac %>%
  bind_rows() %>%
  dplyr::filter(gene_id %in% cell_type_markers$Oligodendrocyte) %>%
  dplyr::select(gene_name, logFC, coef) %>%
  dplyr::distinct(gene_name, coef, .keep_all = T) %>%
  spread(key = "coef", value = "logFC") %>%
  column_to_rownames("gene_name") %>%
  arrange(desc(`MPS-IIIC`)) %>%
  t() %>%
  pheatmap(color = colorRampPalette(rev(brewer.pal(n = 5,
                                                   name = "RdBu")))(100),
           main = "Oligodendrocytes in AC",
           breaks = c(seq(min(.), 0, length.out=ceiling(100/2) + 1),
                      seq(max(.)/100, max(.), length.out=50)),
           cellheight = 30, cellwidth = 12,
           angle_col = 45,
           cluster_rows = F,
           treeheight_row = 0, treeheight_col = 0
  )


toptables_cqn.bc %>%
  bind_rows() %>%
  dplyr::filter(gene_id %in% cell_type_markers$Oligodendrocyte) %>%
  dplyr::select(gene_name, logFC, coef) %>%
  dplyr::distinct(gene_name, coef, .keep_all = T) %>%
  spread(key = "coef", value = "logFC") %>%
  column_to_rownames("gene_name") %>%
  arrange(desc(`MPS-IIIC`)) %>%
  t() %>%
  pheatmap(color = colorRampPalette(rev(brewer.pal(n = 5,
                                                   name = "RdBu")))(100),
           main = "Oligodendrocytes in BC",
           breaks = c(seq(min(.), 0, length.out=ceiling(100/2) + 1),
                      seq(max(.)/100, max(.), length.out=50)),
           cellheight = 30, cellwidth = 12,
           angle_col = 45,
           cluster_rows = F,
           treeheight_row = 0, treeheight_col = 0
  )


# supp figs ---------------------------------------------------------------

ggarrange(
  toptable_1.ab %>%
    bind_rows() %>%
    mutate(rankstat = sign(logFC)*-log10(PValue)) %>%
    ggplot(aes(x = length, y = rankstat)) +
    geom_point(
      aes(colour = DE),
      alpha = 0.5
    ) +
    geom_smooth(se = FALSE, method = "gam") +
    facet_grid(rows = vars(coef)) +
    theme_bw() +
    theme(legend.position = "none") +
    scale_color_manual(values = c("grey50", "red")) +
    scale_x_log10()+
    labs(x = "Average transcript length per gene",
         colour = "Differentially expressed?",
         y = "sign(logFC)*-log10(PValue)") +
    coord_cartesian(ylim = c(-10, 10)),

  toptable_1.ab %>%
    bind_rows() %>%
    mutate(rankstat = sign(logFC)*-log10(PValue)) %>%
    ggplot(aes(x = gc_content, y = rankstat)) +
    geom_point(
      aes(colour = DE),
      alpha = 0.5
    ) +
    geom_smooth(se = FALSE, method = "gam") +
    facet_grid(rows = vars(coef)) +
    theme_bw() +
    theme(legend.position = "none") +
    scale_color_manual(values = c("grey50", "red")) +
    coord_cartesian(ylim = c(-10,10)) +
    labs(x = "Weighted average GC content (%) per gene",
         colour = "Differentially expressed?",
         y = "sign(logFC)*-log10(PValue)"),
  common.legend = TRUE,
  labels = "AUTO"
) +
  ggsave("output/plots/cqn/AB.png", width = 10, height = 10, units = "cm", scale = 2)

ggarrange(
  toptable_1.ac %>%
    bind_rows() %>%
    mutate(rankstat = sign(logFC)*-log10(PValue)) %>%
    ggplot(aes(x = length, y = rankstat)) +
    geom_point(
      aes(colour = DE),
      alpha = 0.5
    ) +
    geom_smooth(se = FALSE, method = "gam") +
    facet_grid(rows = vars(coef)) +
    theme_bw() +
    theme(legend.position = "none") +
    scale_color_manual(values = c("grey50", "red")) +
    scale_x_log10()+
    labs(x = "Average transcript length per gene",
         colour = "Differentially expressed?",
         y = "sign(logFC)*-log10(PValue)") +
    coord_cartesian(ylim = c(-10, 10)),

  toptable_1.ac %>%
    bind_rows() %>%
    mutate(rankstat = sign(logFC)*-log10(PValue)) %>%
    ggplot(aes(x = gc_content, y = rankstat)) +
    geom_point(
      aes(colour = DE),
      alpha = 0.5
    ) +
    geom_smooth(se = FALSE, method = "gam") +
    facet_grid(rows = vars(coef)) +
    theme_bw() +
    theme(legend.position = "none") +
    scale_color_manual(values = c("grey50", "red")) +
    coord_cartesian(ylim = c(-10,10)) +
    labs(x = "Weighted average GC content (%) per gene",
         colour = "Differentially expressed?",
         y = "sign(logFC)*-log10(PValue)"),
  common.legend = TRUE,
  labels = "AUTO"
) +
  ggsave("output/plots/cqn/AC.png", width = 10, height = 10, units = "cm", scale = 2)

ggarrange(
  toptable_1.bc %>%
    bind_rows() %>%
    mutate(rankstat = sign(logFC)*-log10(PValue)) %>%
    ggplot(aes(x = length, y = rankstat)) +
    geom_point(
      aes(colour = DE),
      alpha = 0.5
    ) +
    geom_smooth(se = FALSE, method = "gam") +
    facet_grid(rows = vars(coef)) +
    theme_bw() +
    theme(legend.position = "none") +
    scale_color_manual(values = c("grey50", "red")) +
    scale_x_log10()+
    labs(x = "Average transcript length per gene",
         colour = "Differentially expressed?",
         y = "sign(logFC)*-log10(PValue)") +
    coord_cartesian(ylim = c(-10, 10)),

  toptable_1.bc %>%
    bind_rows() %>%
    mutate(rankstat = sign(logFC)*-log10(PValue)) %>%
    ggplot(aes(x = gc_content, y = rankstat)) +
    geom_point(
      aes(colour = DE),
      alpha = 0.5
    ) +
    geom_smooth(se = FALSE, method = "gam") +
    facet_grid(rows = vars(coef)) +
    theme_bw() +
    theme(legend.position = "none") +
    scale_color_manual(values = c("grey50", "red")) +
    coord_cartesian(ylim = c(-10,10)) +
    labs(x = "Weighted average GC content (%) per gene",
         colour = "Differentially expressed?",
         y = "sign(logFC)*-log10(PValue)"),
  common.legend = TRUE,
  labels = "AUTO"
) +
  ggsave("output/plots/cqn/BC.png", width = 10, height = 10, units = "cm", scale = 2)


KEGG[sigpaths]
  fromList() %>%
  upset

# PATHVIEWS

# PATHVIEW ----------------------------------------------------------------
# Glycoaminoglycan deg ----------------------------------------------------
setwd("output/plots/pathview")
pv_input <- toptabs_cqn_all %>%
  dplyr::filter(coef != "sgsh/+") %>%
  dplyr::select(coef_mps, gene_id, logFC) %>%
  spread(key = "coef_mps", value = "logFC") %>%
  column_to_rownames("gene_id")

plotpathview <- function(coef1, coef2, keggID, outsuff) {

  pathview(species = "dre",
           pathway.id = paste(keggID),
           gene.data = pv_input[c(paste(coef1), paste(coef2))],
           gene.idtype = "ENSEMBL",
           out.suffix = outsuff,
           low = list(gene = "#2166AC"),
           high = list(gene = "#B2182B"),
           mid = list(gene= "white"),
           bins = list(gene= 20)
  )
}

plotpathview(coef1 = "MPS-IIIA_A and B",
             coef2 = "MPS-IIIB_A and B",
             outsuff = "AB Glycosaminoglycan degradation",
             keggID = "00531"
             )


plotpathview(coef1 = "MPS-IIIA_A and C",
             coef2 = "MPS-IIIC_A and C",
             outsuff = "AC Glycosaminoglycan degradation",
             keggID = "00531"
)

plotpathview(coef1 = "MPS-IIIB_B and C",
             coef2 = "MPS-IIIC_B and C",
             outsuff = "BC Glycosaminoglycan degradation",
             keggID = "00531"
)


plotpathview(coef1 = "MPS-IIIA_A and B",
             coef2 = "MPS-IIIC_A and B",
             outsuff = "AB lyso",
             keggID = "04142"
)




# IRE ---------------------------------------------------------------------

irePlot <- fryire.all %>%
  dplyr::filter(coef != "sgsh/+") %>%
  mutate(MPS = case_when(MPS == "AC" ~ "A and C",
                         MPS == "AB" ~ "A and B",
                         MPS == "BC" ~ "B and C",
                         MPS == "B2" ~ "B2") %>%
           factor(levels = c("A and C", "A and B", "B and C", "B2"))) %>%
  ggplot(aes(x = coef,
             y = pathway )) +
  geom_tile(aes(fill = -log10(PValue.Mixed),
                alpha = FDR.Mixed < 0.05),
            colour = "black") +
  geom_label(aes(label = signif(FDR.Mixed, digits = 2)),
             fill = NA) +
  facet_wrap(~MPS, scales = "free_x", nrow = 1 ) +
  scale_fill_viridis_c() +
  scale_alpha_manual(values = c(0.5,1)) +
  ggpubr::theme_pubclean() +
  theme(panel.grid.major.y = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        text = element_text(face = "bold", size = 15),
        axis.text.x = element_text(face = "bold", size = 15, angle = 315, hjust = 0.1, vjust = 0.2),
        legend.position = "right"

  )
  #ggsave("output/plots/ire_withExtraB_pvals.png", width = 5, height = 2.5, units = "cm", dpi = 500, scale = 5)

png("output/plots/ire/ire3_allHeatmap.png", width = 50, height = 20, units = "cm", res = 100)
toptabs_cqn_all %>%
  dplyr::filter(gene_id %in% ireGenes$ire3_all) %>%
  dplyr::filter(coef != "sgsh/+") %>%
  dplyr::filter(abs(logFC) < 2 ) %>%
  mutate(coef_mps = paste0(coef, " (", mps, ")")) %>%
  dplyr::select(gene_name, logFC, coef_mps) %>%
  dplyr::distinct(gene_name, coef_mps, .keep_all = T) %>%
  spread(key = "coef_mps", value = "logFC") %>%
  column_to_rownames("gene_name") %>%
  na.omit %>%
  t() %>%
  pheatmap(color = colorRampPalette(rev(brewer.pal(n = 5,
                                                   name = "RdBu")))(100),
           main = "ire3_all",
           breaks = c(seq(min(.), 0, length.out=ceiling(100/2) + 1),
                      seq(max(.)/100, max(.), length.out=50)),
           cellheight = 30,
           #cellwidth = 12,
           angle_col = 45,
           border_color = NA,
           show_colnames = F,
           treeheight_row = 10,
           treeheight_col = 20
  )
dev.off()

png("output/plots/ire/ire5_allHeatmap.png", width = 50, height = 20, units = "cm", res = 100)
toptabs_cqn_all %>%
  dplyr::filter(gene_id %in% ireGenes$ire5_all) %>%
  dplyr::filter(coef != "sgsh/+") %>%
  dplyr::filter(abs(logFC) < 2 ) %>%
  mutate(coef_mps = paste0(coef, " (", mps, ")")) %>%
  dplyr::select(gene_name, logFC, coef_mps) %>%
  dplyr::distinct(gene_name, coef_mps, .keep_all = T) %>%
  spread(key = "coef_mps", value = "logFC") %>%
  column_to_rownames("gene_name") %>%
  na.omit %>%
  t() %>%
  pheatmap(color = colorRampPalette(rev(brewer.pal(n = 5,
                                                   name = "RdBu")))(100),
           main = "ire5_all",
           breaks = c(seq(min(.), 0, length.out=ceiling(100/2) + 1),
                      seq(max(.)/100, max(.), length.out=50)),
           cellheight = 30,
           #cellwidth = 12,
           angle_col = 45,
           border_color = NA,
           show_colnames = F,
           treeheight_row = 10,
           treeheight_col = 20
  )
dev.off()

pcaMeta <- x.ab$samples %>%
  dplyr::select(sample, genotype,Experiment = MPS ) %>%
  dplyr::filter(genotype != "sgsh/+") %>%
  rbind(x.ac$samples %>%
          dplyr::select(sample, genotype, Experiment = MPS)) %>%
  rbind(x.bc$samples %>%
          dplyr::select(sample, genotype, Experiment = MPS)) %>%
  rbind(x.b2$samples %>%
          dplyr::filter(genotype %in% c("WT", "MPS-IIIB")) %>%
          dplyr::select(sample, genotype) %>%
          mutate(Experiment = "B2"))

#plot PCA fn
MPScols <- c("#4000ff", #WT
             "#ff005e", # A
             "#ffae00", #B
             "#00b064" #C
             )



plotPCA <- function(geneset, colourby, colourVec, plot.title) {
  logCPm.all %>%
    dplyr::filter(gene_id %in% geneset) %>%
    column_to_rownames("gene_id") %>%
    na.omit %>%
    t() %>%
    prcomp %>%
    autoplot(
      data = tibble(sample = rownames(.$x)) %>%
        left_join(pcaMeta),
      colour = colourby,
      size = 4,
      alpha = 0.75
    ) +
    scale_color_manual(values =  colourVec) +
    theme(aspect.ratio = 1,
          plot.background = element_rect(fill = NA),
          legend.background = element_rect(fill = NA),
          panel.background = element_rect(fill = NA),
          legend.box.background = element_rect(fill = NA),
          text = element_text(face = "bold", size = 15)

    ) +
    labs(colour = colourby,
         title = plot.title)
}


a <- plotPCA(geneset = ireGenes$ire3_all,
        colourby = "genotype",
        colourVec = MPScols,
        plot.title = "3' IRE genes")
b <- plotPCA(geneset = ireGenes$ire5_all,
          colourby = "genotype",
          colourVec = MPScols,
          plot.title = "5' IRE genes")

c <- plotPCA(geneset = ireGenes$ire3_all,
             colourby = "Experiment",
             colourVec = viridis_pal(end = 0.8)(4),
             plot.title = "3' IRE genes")
d <- plotPCA(geneset = ireGenes$ire5_all,
             colourby = "Experiment",
             colourVec = viridis_pal(end = 0.8)(4),
             plot.title = "5' IRE genes")
irePlot +
ggarrange(
  ggarrange(a, b, common.legend = T),
  ggarrange(c, d, common.legend = T)
) +
  plot_layout(ncol = 1, heights = c(0.75, 1)) +
  ggsave("output/plots/ire/irePlot.png", width = 5, height = 3, units = "cm", dpi = 500, scale = 7 )




