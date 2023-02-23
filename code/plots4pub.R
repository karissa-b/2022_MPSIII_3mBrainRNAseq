# PLOTS FOR PUB
library(patchwork)
library(UpSetR)
library(fgsea)
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
  # bind_rows(toptables_cqn.b2 %>%
  #             bind_rows() %>%
  #             mutate(mps = "B2") ) %>%
  mutate(coef_mps = paste0(coef, "_", mps))

logCPm.all <- logCPM.ac %>%
  as.data.frame() %>%
  rownames_to_column("gene_id") %>%
  left_join(logCPM.ab[,(x.ab$samples %>% dplyr::filter(genotype != "sgsh/+") %>% .$sample)] %>% as.data.frame %>% rownames_to_column("gene_id") ) %>%
  left_join(logCPM.bc %>% as.data.frame %>% rownames_to_column("gene_id") ) %>%
  left_join(logCPM.b2[,(x.b2$samples %>% dplyr::filter(genotype != "EOfAD-like") %>% .$sample)] %>% as.data.frame %>% rownames_to_column("gene_id") )



# custom function to plot pearson cors
makePearsonCorrPlot <- function(geneset, plot.title) {
  plotpearsonCor <- function(x.axis.coef, y.axis.coef, geneset) {
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

  a <- plotpearsonCor(x.axis.coef = "MPS-IIIA_A and B",
                       y.axis.coef = "MPS-IIIA_A and C",
                       geneset = geneset)

  # MPS-IIIB
  b <- plotpearsonCor(x.axis.coef = "MPS-IIIB_A and B",
                       y.axis.coef = "MPS-IIIB_B and C",
                       geneset = geneset)

  # MPS-IIIC
  c <- plotpearsonCor(x.axis.coef = "MPS-IIIC_A and C",
                       y.axis.coef = "MPS-IIIC_B and C",
                       geneset = geneset)

  # A v B
  d <- plotpearsonCor(x.axis.coef = "MPS-IIIA_A and B",
                       y.axis.coef = "MPS-IIIB_A and B",
                       geneset = geneset)


  # A v C
  e <- plotpearsonCor(x.axis.coef = "MPS-IIIA_A and C",
                       y.axis.coef = "MPS-IIIC_A and C",
                       geneset = geneset)

  # B v C
  f <- plotpearsonCor(x.axis.coef = "MPS-IIIB_B and C",
                       y.axis.coef = "MPS-IIIC_B and C",
                       geneset = geneset)

  (a | b | c) /( d | e | f) +
    plot_layout(guides= "collect") +
    plot_annotation(title = plot.title)
}

# PCA ---------------------------------------------------------------------
# PCA all samps
a <- logCPM.ab %>%
  as.data.frame() %>%
  rownames_to_column("gene_id") %>%
  left_join(logCPM.ac%>%
              as.data.frame() %>%
              rownames_to_column("gene_id")
            ) %>%
  left_join(logCPM.bc %>%
              as.data.frame() %>%
              rownames_to_column("gene_id")
  )  %>%
  na.omit() %>%
  remove_rownames() %>%
  column_to_rownames("gene_id") %>%
  t() %>%
  prcomp() %>%
  autoplot(data = tibble(sample = rownames(.$x)) %>%
             left_join(meta) %>%
             mutate(Arm = MPS),
           colour = "Arm",
           alpha = 0.75,
           size = 4
  ) +
  stat_ellipse(aes(group = MPS, colour = MPS)) +
  theme_classic()
  ggsave("output/plots/PCA_allDataOneplot.png", width = 7, height = 5, units = "cm", scale = 2)



# AB no hets

logCPM.ab %>%
  as.data.frame() %>%
  dplyr::select(c(x.ab$samples %>%
                  dplyr::filter(genotype != "sgsh/+") %>% .$sample)) %>%
  as.matrix() %>%
  t() %>%
  prcomp() %>%
  autoplot(data = tibble(sample = rownames(.$x)) %>%
             left_join(x.ab$samples),
           colour = "genotype",
           shape = "Sex",
           size = 4
  ) +
  scale_color_manual(values = c("#4000ff", #WT
                                #"#ff7aab", # sgsh het
                                "#ff005e", # A
                                "#ffae00" #B
  )) +
  theme(aspect.ratio = 1) +
  labs(title = "A and B",
       colour = "Genotype") +
  ggsave("output/plots/PCA_ab_nohet.png", width = 7, height = 5, units = "cm", scale =1.5)

# example home tank
logCPM.ac %>%
  t() %>%
  prcomp() %>%
  autoplot(data = tibble(sample = rownames(.$x)) %>%
             left_join(x.ac$samples) %>%
             mutate(Tank = str_remove(Tank, pattern = "AC_P2_lay1_tank")),
           colour = "Tank",
           shape = "genotype",
           size = 4
  ) +
  theme_classic() +
  scale_color_brewer(palette = "Dark2") +
  ggsave("output/plots/PCA_ac_hometank.png", width = 7, height = 5, units = "cm", scale = 2)

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
             left_join(x.bc$samples) %>%
             mutate(P = str_extract(Tank, pattern = "P[1|2]"),
                    lay = str_extract(Tank, pattern = "P[1|2]_lay[1-3]"),
                    DOB = as.character(DOB)),
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





  bind_rows(fryKEGG.b2 %>%
              bind_rows() %>%
              mutate(MPS= "B2"))

sigpaths <-
  frykegg.all %>%
  dplyr::filter(FDR.Mixed < 0.05) %>%
  .$pathway %>%
  unique()

sigpaths %>% str_remove("KEGG_") %>% sort

Aonly <- fryKEGG.ac$`MPS-IIIA` %>%
  dplyr::filter(FDR.Mixed < 0.05) %>%
  dplyr::filter(!pathway %in% (fryKEGG.b2$`MPS-IIIB` %>%
                  dplyr::filter(FDR.Mixed < 0.05) %>%
                  .$pathway)) %>%
  .$pathway %>%
  str_remove("KEGG_")

frykegg.all %>%
  bind_rows() %>%
  dplyr::filter(pathway %in% sigpaths) %>%
  dplyr::filter(coef != "sgsh/+") %>%
  mutate(MPS = case_when(MPS == "AC" ~ "A and C",
                         MPS == "AB" ~ "A and B",
                         MPS == "BC" ~ "B and C",
                         MPS == "B2" ~ "B_6m") %>%
           factor(levels = c("A and C", "B_6m", "A and B", "B and C"))) %>%
  mutate(pathway = str_remove(pathway, pattern = "KEGG_"),
         order = case_when(
           pathway == "LYSOSOME" ~ 1,
           pathway %in% c("OTHER_GLYCAN_DEGRADATION",
                          "GLYCOSAMINOGLYCAN_DEGRADATION")  ~ 2,
           pathway %in% c("AMINO_SUGAR_AND_NUCLEOTIDE_SUGAR_METABOLISM",
                          "ANTIGEN_PROCESSING_AND_PRESENTATION") ~ 3,
           pathway %in% c("FATTY_ACID_METABOLISM" ) ~ 4,
           pathway %in% c("SYSTEMIC_LUPUS_ERYTHEMATOSUS",
                          "SPHINGOLIPID_METABOLISM",
                          "DRUG_METABOLISM_CYTOCHROME_P450",
                          "SNARE_INTERACTIONS_IN_VESICULAR_TRANSPORT",
                          "COMPLEMENT_AND_COAGULATION_CASCADES",
                          "CELL_ADHESION_MOLECULES_CAMS") ~ 5,
           pathway %in% Aonly ~ 6,
           TRUE ~7),
         DE = FDR.Mixed < 0.05
         ) %>%
  ggplot(aes(x = coef,
             y = reorder(pathway, -order)) ) +
  geom_tile(aes(fill = -log10(PValue.Mixed),
                alpha = FDR.Mixed < 0.05)) +
  geom_label(aes(label = signif(FDR.Mixed, digits = 2)),
             fill = NA) +
  facet_wrap(~MPS, scales = "free_x", nrow = 1) +
  scale_fill_viridis_c() +
  scale_alpha_manual(values = c(0.1,1)) +
  ggpubr::theme_pubclean() +
  theme(panel.grid.major.y = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        axis.text.y = element_text(face = "bold", size = 15),
        axis.text.x = element_text(face = "bold", size = 15, angle = 315, hjust = 0.1, vjust = 0.2),
        legend.position = "right"

  ) +
    ggsave("output/plots/KEGG_pval_heatmap.png", width = 18, height = 16, units = "cm", dpi = 500, scale = 2)

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
  na.omit() %>%
  t() %>%
  pheatmap(color = colorRampPalette(rev(brewer.pal(n = 5,
                                                   name = "RdBu")))(100),
           main = "",
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
  na.omit %>%
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


png("output/plots/KEGG_AMINO_SUGAR_AND_NUCLEOTIDE_SUGAR_METABOLISM_Heatmap.png", width = 70, height = 20, units = "cm", res = 100)
toptabs_cqn_all %>%
  bind_rows() %>%
  dplyr::filter(gene_id %in% KEGG$KEGG_AMINO_SUGAR_AND_NUCLEOTIDE_SUGAR_METABOLISM) %>%
  dplyr::filter(coef != "sgsh/+") %>%
  mutate(coef_mps = paste0(coef, " (", mps, ")")) %>%
  dplyr::select(gene_name, logFC, coef_mps) %>%
  dplyr::distinct(gene_name, coef_mps, .keep_all = T) %>%
  spread(key = "coef_mps", value = "logFC") %>%
  na.omit() %>%
  column_to_rownames("gene_name") %>%
  t() %>%
  pheatmap(color = colorRampPalette(rev(brewer.pal(n = 5,
                                                   name = "RdBu")))(256),
           main = "KEGG_AMINO_SUGAR_AND_NUCLEOTIDE_SUGAR_METABOLISM",
           breaks = c(seq(min(.), 0, length.out=ceiling(256/2) + 1),
                      seq(max(.)/256, max(.), length.out=256/2)),
           cellheight = 30, cellwidth = 16,
           angle_col = 45,
           treeheight_row = 10,
           fontsize = 15)
dev.off()

png("output/plots/KEGG_ANTIGEN_PROCESSING_AND_PRESENTATION_Heatmap.png", width = 70, height = 20, units = "cm", res = 100)
toptabs_cqn_all %>%
  bind_rows() %>%
  dplyr::filter(gene_id %in% KEGG$KEGG_ANTIGEN_PROCESSING_AND_PRESENTATION) %>%
  dplyr::filter(coef != "sgsh/+") %>%
  mutate(coef_mps = paste0(coef, " (", mps, ")")) %>%
  dplyr::select(gene_name, logFC, coef_mps) %>%
  dplyr::distinct(gene_name, coef_mps, .keep_all = T) %>%
  spread(key = "coef_mps", value = "logFC") %>%
  na.omit() %>%
  column_to_rownames("gene_name") %>%
  t() %>%
  pheatmap(color = colorRampPalette(rev(brewer.pal(n = 5,
                                                   name = "RdBu")))(256),
           main = "KEGG_ANTIGEN_PROCESSING_AND_PRESENTATION",
           breaks = seq(-0.6, to = 0.6, by = 1.2/256),
           cellheight = 30, cellwidth = 16,
           angle_col = 45,
           treeheight_row = 10,
           fontsize = 15)
dev.off()

png("output/plots/KEGG_COMPLEMENT_AND_COAGULATION_CASCADES_Heatmap.png", width = 70, height = 20, units = "cm", res = 100)
toptabs_cqn_all %>%
  dplyr::filter(gene_id %in% KEGG$KEGG_COMPLEMENT_AND_COAGULATION_CASCADES) %>%
  dplyr::filter(coef != "sgsh/+") %>%
  mutate(coef_mps = paste0(coef, " (", mps, ")")) %>%
  dplyr::select(gene_name, logFC, coef_mps) %>%
  dplyr::distinct(gene_name, coef_mps, .keep_all = T) %>%
  spread(key = "coef_mps", value = "logFC") %>%
  na.omit() %>%
  column_to_rownames("gene_name") %>%
  t() %>%
  pheatmap(color = colorRampPalette(rev(brewer.pal(n = 5,
                                                   name = "RdBu")))(256),
           main = "KEGG_COMPLEMENT_AND_COAGULATION_CASCADES",
           breaks = seq(-0.6, to = 0.6, by = 1.2/256),
           cellheight = 30, cellwidth = 16,
           angle_col = 45,
           treeheight_row = 10,
           fontsize = 15)
dev.off()


png("output/plots/ac_only/KEGG_FATTY_ACID_METABOLISM_Heatmap.png", width = 70, height = 20, units = "cm", res = 100)
toptables_cqn.ac %>%
  dplyr::filter(gene_id %in% KEGG$KEGG_FATTY_ACID_METABOLISM) %>%
  dplyr::filter(coef != "sgsh/+") %>%
  mutate(coef_mps = paste0(coef, " (", mps, ")")) %>%
  dplyr::select(gene_name, logFC, coef_mps) %>%
  dplyr::distinct(gene_name, coef_mps, .keep_all = T) %>%
  spread(key = "coef_mps", value = "logFC") %>%
  na.omit() %>%
  column_to_rownames("gene_name") %>%
  t() %>%
  pheatmap(color = colorRampPalette(rev(brewer.pal(n = 5,
                                                   name = "RdBu")))(256),
           main = "KEGG_FATTY_ACID_METABOLISM",
           breaks = seq(-0.6, to = 0.6, by = 1.2/256),
           cellheight = 30, cellwidth = 16,
           angle_col = 45,
           treeheight_row = 10,
           fontsize = 15)
dev.off()

png("output/plots/ac_only/KEGG_FATTY_ACID_METABOLISM_Heatmap.png", width = 70, height = 20, units = "cm", res = 100)
toptables_cqn.ac %>%
  bind_rows() %>%
  dplyr::filter(gene_id %in% KEGG$KEGG_FATTY_ACID_METABOLISM) %>%
  dplyr::filter(coef != "sgsh/+") %>%
  dplyr::select(gene_name, logFC, coef) %>%
  dplyr::distinct(gene_name, coef, .keep_all = T) %>%
  spread(key = "coef", value = "logFC") %>%
  na.omit() %>%
  column_to_rownames("gene_name") %>% view
  t() %>%
  pheatmap(color = colorRampPalette(rev(brewer.pal(n = 5,
                                                   name = "RdBu")))(256),
           main = "KEGG_FATTY_ACID_METABOLISM",
           breaks = seq(-0.6, to = 0.6, by = 1.2/256),
           cellheight = 30, cellwidth = 16,
           angle_col = 45,
           treeheight_row = 10,
           fontsize = 15)
dev.off()

# correlations of


# correlations of


# lyso correlations --------------------------------------------------
makePearsonCorrPlot(geneset = KEGG$KEGG_LYSOSOME, plot.title = "KEGG_LYSOSOME") +
  ggsave("output/plots/pearsonCor_kegglyso.png", width = 5, height = 3, units = "cm", dpi = 500, scale = 5)

# glycosodeg correlations --------------------------------------------------

makePearsonCorrPlot(geneset = KEGG$KEGG_GLYCOSAMINOGLYCAN_DEGRADATION, plot.title = "KEGG_GLYCOSAMINOGLYCAN_DEGRADATION")
  ggsave("output/plots/pearsoncorrs_KEGG_GLYCOSAMINOGLYCAN_DEGRADATION.png", width = 5, height = 3, units = "cm", dpi = 500, scale = 5)


# other glycan deg correlations --------------------------------------------------------
  makePearsonCorrPlot(geneset = KEGG$KEGG_OTHER_GLYCAN_DEGRADATION, plot.title = "KEGG_OTHER_GLYCAN_DEGRADATION") +
  ggsave("output/plots/pearsoncorrs_KEGG_OTHER_GLYCAN_DEGRADATION.png", width = 5, height = 3, units = "cm", dpi = 500, scale = 5)

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
      %>%
    # bind_rows(celltype.b2 %>%
    #             bind_rows() %>%
    #             mutate(MPS = "B2"))



frycell.all %>%
  dplyr::filter(coef != "sgsh/+") %>%
    arrange(PValue) %>%
    mutate(order = case_when(
      pathway %in% c("Oligodendrocyte", "Oligodendrocyte_sla high",
                     "Neural stem cell", "Microglia_apoc1 high")~ 1,
        TRUE ~ 2
    )
             ) %>%
  ggplot(aes(x = coef,
             y = reorder(x = pathway, -order))
  ) +
  geom_tile(aes(fill = -log10(PValue),
                alpha = FDR < 0.05),
            colour = "black") +
  geom_label(aes(label = signif(FDR, digits = 2)
  ),
  fill = NA) +
  facet_wrap(~MPS, scales = "free_x", nrow = 1) +
  scale_fill_viridis_c() +
  scale_alpha_manual(values = c(0.1,1)) +
  ggpubr::theme_pubclean() +
  theme(panel.grid.major.y = element_blank(),
        legend.position = "right",
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        axis.text.y = element_text(face = "bold", size = 15),
        axis.text.x = element_text(face = "bold", size = 15, angle = 315, hjust = 0.1, vjust = 0.2)
        ) +
  ggsave("output/plots/celltype_heatmap.png", width = 15, height = 16, units = "cm", dpi = 500, scale = 2)


  x# correlations ------------------------------------------------------------

makePearsonCorrPlot(geneset = cell_type_markers$Oligodendrocyte, plot.title = "OLIGODENDROCYTES") +
  ggsave("output/plots/Pearsoncorrs_oligoden.png", width = 5, height = 3, units = "cm", dpi = 500, scale = 5)




# Neural stem cell --------------------------------------------------------
makePearsonCorrPlot(geneset = cell_type_markers$`Neural stem cell`, plot.title = "Neural stem cell") +
    ggsave("output/plots/Pearsoncorrs_Neural stem cell.png", width = 5, height = 3, units = "cm", dpi = 500, scale = 5)


# microglia ---------------------------------------------------------------
makePearsonCorrPlot(geneset = cell_type_markers$`Microglia_apoc1 high`, plot.title = "NMicroglia_apoc1 high") +
    ggsave("output/plots/Pearsoncorrs_Microglia_apoc1 high.png", width = 5, height = 3, units = "cm", dpi = 500, scale = 5)


# heatmaps of FC ----------------------------------------------------------#

  anno_celltype <-
    cell_type_markers[c("Oligodendrocyte", "Neural stem cell", "Microglia_apoc1 high") ] %>%
    unlist() %>%
    as.data.frame() %>%
    rownames_to_column("pathway") %>%
    set_colnames(c("pathway", "gene_id")) %>%
    mutate(pathway = str_remove(pathway, pattern = "[0-9]+$")) %>%
    as_tibble() %>%
    pivot_wider(names_from = pathway, values_from = pathway) %>%
    left_join(x.ac$genes %>% dplyr::select(gene_id, gene_name)) %>%
    dplyr::select(-gene_id) %>%
    dplyr::distinct(gene_name, .keep_all= T) %>%
    dplyr::filter(!is.na(gene_name)) %>%
    column_to_rownames("gene_name")

  # all
  toptabs_cqn_all %>%
    dplyr::filter(gene_id %in% c(cell_type_markers$Oligodendrocyte,
                                 cell_type_markers$`Microglia_apoc1 high`,
                                 cell_type_markers$`Neural stem cell`)) %>%
    dplyr::filter(coef != "sgsh/+") %>%
    mutate(coef_mps = paste0(coef, " (", mps, ")")) %>%
    dplyr::select(gene_name, logFC, coef_mps) %>%
    dplyr::distinct(gene_name, coef_mps, .keep_all = T) %>%
    spread(key = "coef_mps", value = "logFC") %>%
    column_to_rownames("gene_name") %>%
    pheatmap(color = colorRampPalette(rev(brewer.pal(n = 5,
                                                     name = "RdBu")))(100),
             main = "Oligodendrocytes",
             breaks = c(seq(-1,1, by = 2/100)),
           #  cellheight = 30, cellwidth = 12,
             angle_col = 45,
             treeheight_row = 10,
             annotation_row = anno_celltype
             #treeheight_col = 20
    )


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
           breaks = c(seq(-1,1, by = 2/100)),
           cellheight = 30, cellwidth = 12,
           angle_col = 45,
           treeheight_row = 10,
           annotation_col = anno_celltype
           #treeheight_col = 20
  )
dev.off()








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

png("output/plots/upset_celltype.png", width = 12, height = 8, units = "cm", res = 300)
cell_type_markers[c("Oligodendrocyte", "Oligodendrocyte_sla high", "Neural stem cell", "Microglia_apoc1 high" ) ] %>%
  fromList() %>%
  upset(
        order.by = "freq",
        nsets = 23,
        mb.ratio = c(0.6, 0.4),
        )
dev.off()

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
             coef2 = "MPS-IIIB_A and B",
             outsuff = "AB lyso",
             keggID = "04142"
)

plotpathview(coef1 = "MPS-IIIA_A and C",
             coef2 = "MPS-IIIC_A and C",
             outsuff = "antigenAC",
             keggID = "04612"
)


plotpathview(coef1 = "MPS-IIIA_A and C",
             coef2 = "MPS-IIIC_A and C",
             outsuff = "aminosugNucleoSug",
             keggID = "00520"
)

plotpathview(coef1 = "MPS-IIIA_A and C",
             coef2 = "MPS-IIIC_A and C",
             outsuff = "fattyacid",
             keggID = "01212"
)
01212



# IRE ---------------------------------------------------------------------

fryire.all <-
  fryIRE.ac %>%
  bind_rows() %>%
  mutate(MPS = "AC") %>%
  bind_rows(fryIRE.ab %>%
              bind_rows() %>%
              mutate(MPS= "AB")) %>%
  bind_rows(fryIRE.bc %>%
              bind_rows() %>%
              mutate(MPS = "BC"))
  #bind_rows(fryIRE.b2 %>%
             # bind_rows() %>%
             #  mutate(MPS = "B2"))

fryire.all %>%
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
  scale_alpha_manual(values = c(0.2,1)) +
  ggpubr::theme_pubclean() +
  theme(panel.grid.major.y = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        text = element_text(face = "bold", size = 16),
        axis.text.x = element_text(face = "bold", size = 15, angle = 315, hjust = 0.1, vjust = 0.2),
        legend.position = "bottom"
  ) +
  ggsave("output/plots/ire_pvals.png", width = 5, height = 2.5, units = "cm", dpi = 500, scale = 5)

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





# why is fatty acid DE in hgsnat but not sgsh ----------------------------
temp <- x.ac$genes %>%
  dplyr::filter(gene_name %in% c("acsl1b", "hadh", "adh8b", "ehhadh")) %>%
  .$gene_id

temp2 <- logCPM.ac %>%
  as.data.frame() %>%
  rownames_to_column('gene_id') %>%
  dplyr::filter(!(gene_id %in% temp)) %>%
  column_to_rownames('gene_id')
  .[-temp,]

# looking at the HM,
design.ac %>% colnames() %>% .[2:3] %>%
  sapply(function(y) {
    temp2 %>%
      fry(
        index = KEGG,
        design = design.ac,
        contrast = y,
        sort = "mixed"
      ) %>%
      rownames_to_column("pathway") %>%
      as_tibble() %>%
      mutate(coef = y)
  }, simplify = FALSE)

# leading edge genes

fgsea <- toptables_cqn.ac %>%
  sapply(function(x) {
    x %>%
      mutate(rank = sign(logFC) * log10(1/PValue)) %>%
      arrange(rank) %>%
      dplyr::select(gene_id, rank) %>% #only want the Pvalue with sign
      with(structure(rank, names = gene_id)) %>%
      rev() # reverse so the start of the list is upregulated genes
  }, simplify = FALSE
  ) %>%
  lapply(function(x){
    fgseaMultilevel(stats = x,
                    pathways = cell_type_markers) %>%
      as_tibble() %>%
      dplyr::rename(FDR = padj) %>%
      mutate(padj = p.adjust(pval, "bonferroni")) %>%
      dplyr::select(pathway, pval, FDR, padj, everything()) %>%
      arrange(pval) %>%
      mutate(sig = padj < 0.05)
  })

fgsea$`MPS-IIIA` %<>% mutate(coef = "MPS-IIIA")
fgsea$`MPS-IIIC` %<>% mutate(coef = "MPS-IIIC")


# UpSet Leading edge ------------------------------------------------------

png("output/plots/ac_only/upsetA_celltype.png", width = 12, height = 8, units = "cm", res = 300)
fgsea$`MPS-IIIA` %>%
  dplyr::filter(pathway %in% c("Oligodendrocyte",
                               "Oligodendrocyte_sla high",
                               "Neural stem cell",
                               "Microglia_apoc1 high")) %>%
  dplyr::select(pathway, leadingEdge) %>%
  unnest %>%
  split(f = .$pathway) %>%
  lapply(magrittr::extract2, "leadingEdge") %>%
  fromList() %>%
  upset(order.by = "freq",
        mb.ratio = c(0.7, 0.3))
dev.off()

png("output/plots/ac_only/upsetC_celltype.png", width = 12, height = 8, units = "cm", res = 300)
fgsea$`MPS-IIIC` %>%
  dplyr::filter(pathway %in% c("Oligodendrocyte",
                               "Oligodendrocyte_sla high",
                               "Neural stem cell",
                               "Microglia_apoc1 high")) %>%
  dplyr::select(pathway, leadingEdge) %>%
  unnest %>%
  split(f = .$pathway) %>%
  lapply(magrittr::extract2, "leadingEdge") %>%
  fromList() %>%
  upset(order.by = "freq")
dev.off()

# plot 4 bruce ------------------------------------------------------------

toptables_cqn.ac %>%
  bind_rows() %>%
  ggplot(aes(y = -log10(PValue), x = logFC, colour = DE)) +
  geom_point(
    alpha = 0.5, size = 1.25
  ) +
  facet_wrap(~coef, nrow = 1) +
  coord_cartesian(xlim = c(-2.5,2.5),
                  ylim = c(0, 10) ) +
  theme(legend.position = "bottom") +
  scale_color_manual(values = c("grey50", "red")) +
  geom_label_repel(aes(label = gene_name),
                   data = . %>%
                     dplyr::filter(FDR < 0.1) %>%
                     dplyr::filter(gene_id %in% GO$GO_OLIGODENDROCYTE_DIFFERENTIATION))

goseqac.a %>%
  dplyr::filter(FDR < 0.05) %>%
  mutate(propDE = numDEInCat/numInCat) %>%
  ggplot(aes(x = -log10(over_represented_pvalue),
             fill = propDE,
             y = reorder(category, -over_represented_pvalue))) +
  geom_col(colour = "black") +
  scale_fill_viridis_c() +
  labs(fill = "Proportion of DE\ngenes in GO term",
       x = "-log10(pvalue)" )+
  theme(panel.grid = element_blank(),
        axis.ticks = element_blank(),
        axis.title.y = element_blank(),
        text = element_text(face = "bold", size = 15),
        legend.position = "right") +
  ggsave("output/plots/GOtermssgsh_ac.png", width = 15, height = 10, units = "cm",
         dpi = 400, scale = 2.5)

goseqac.c %>%
  dplyr::filter(FDR < 0.1) %>%
  ggplot(aes(x = -log10(over_represented_pvalue),
             y = reorder(category, -over_represented_pvalue))) +
  geom_col() +
  ggtitle("Significnt GO terms in MPS-IIIC")


# Manhattenpltos ----------------------------------------------------------
plotman = function(manhatinput, axis) {
  ggplot(manhatinput, aes(x = midCum, y = -log10(PValue))) +
  geom_point(aes(color=chromosome), alpha = 0.5, size = 1) +
  geom_point(alpha = 0.5, size = 1, colour = "red",
             data = . %>%
               dplyr::filter(coef == "MPS-IIIB" & chromosome == 24)) +
  geom_point(alpha = 0.5, size = 1, colour = "red",
             data = . %>%
               dplyr::filter(coef ==  "MPS-IIIA" & chromosome == 22)) +
  geom_point(alpha = 0.5, size = 1, colour = "red",
               data = . %>%
                 dplyr::filter(coef ==  "MPS-IIIC" & chromosome == 1)) +
  geom_point(alpha = 0.5, size = 1, colour = "red",
             data = . %>%
               dplyr::filter(coef ==  "sgsh/+" & chromosome == 22)) +
  scale_color_manual(values = axis.ab$colour) +
  scale_x_continuous(label = axis.ab$chromosome, breaks = axis.ab$center) +
  scale_y_continuous(limits =c(0, 10)) + # zoom in
  facet_wrap(~coef, ncol = 1) +
  labs(x = "Chromosome", y = expression(paste(-log[10], "(p)"))) +
  theme(
    legend.position="none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  )
}

a <- plotman(manhatinput = manhat_input.ab, axis = axis.ab) +
  ggtitle("A and B")


b <- plotman(manhatinput = manhat_input.ac, axis = axis.ac) +
  ggtitle("A and C")
c <- plotman(manhatinput = manhat_input.bc, axis = axis.bc) +
  ggtitle("B and C")

a / (b + c) +
  plot_annotation(tag_levels = 'A') +
  ggsave("output/plots/stackedmanhats.png", width = 15, height = 20, units = "cm", dpi = 400, scale = 2)


