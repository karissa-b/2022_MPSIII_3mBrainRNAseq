library(tidyverse)
library(patchwork)
library(UpSetR)
library(pathview)
library(ggeasy)
library(colorspace)

x.ac <- read_rds("data/R_objects/dge_ac.rds")
toptables_cqn.ac <- read_rds("data/R_objects/toptab_ac_cqn.rds")
fryKEGG.ac <- read_rds("data/R_objects/AC_fryKEGG.rds")
fryIRE.ac <- read_rds("data/R_objects/AC_fryIRE.rds")
celltype.ac <- read_rds("data/R_objects/AC_frycellTypes.rds")

x.ab <- read_rds("data/R_objects/dge_AB.rds")
fryKEGG.ab <- read_rds("data/R_objects/AB_fryKEGG.rds")
fryIRE.ab <- read_rds("data/R_objects/AB_fryIRE.rds")
celltypes.ab <- read_rds("data/R_objects/AB_frycellTypes.rds")
toptables_cqn.ab <- read_rds("data/R_objects/AB_toptab_cqn.rds")

x.bc <- read_rds("data/R_objects/dge_bc.rds")
toptables_cqn.bc <- read_rds("data/R_objects/toptab_bc_cqn.rds")
fryKEGG.bc <- read_rds("data/R_objects/BC_fryKEGG.rds")
fryIRE.bc <- read_rds("data/R_objects/BC_fryIRE.rds")
celltype.bc <- read_rds("data/R_objects/BC_frycellTypes.rds")

x.b2 <- read_rds("data/R_objects/dge_B2.rds")
fryKEGG.b2 <- read_rds("data/R_objects/B2_fryKEGG.rds")
fryIRE.b2 <- read_rds("data/R_objects/B2_fryIRE.rds")
celltype.b2 <- read_rds("data/R_objects/B2_frycellTypes.rds")
toptables_cqn.b2 <-  read_rds("data/R_objects/B2_toptab_cqn.rds")

toptables_cqn_abc <-
  toptables_cqn.ac %>%
  bind_rows() %>%
  mutate(mps = "A and C") %>%
  bind_rows(toptables_cqn.b2 %>%
              bind_rows() %>%
              mutate(mps = "B2") ) %>%
  mutate(coef_mps = paste0(coef, "_", mps))


# PCA ---------------------------------------------------------------------
logCPM.bc %>%
  t() %>%
  prcomp() %>%
  autoplot(data = tibble(sample = rownames(.$x)) %>%
             left_join(x.bc$samples) %>%
             mutate(P = str_extract(Tank, pattern = "P[1|2]"),
                    lay = str_extract(Tank, pattern = "P[1|2]_lay[1-3]"),
                    DOB = as.character(DOB)),
           colour = "DOB",
           shape = "genotype",
           size = 4
  ) +
scale_color_brewer(palette = "Dark2", direction = -1)+
  theme(aspect.ratio = 1) +
  labs(title = "B and C",
       colour = "DOB",
       shape = "Genotype") +
  ggsave("output/plots4talk/PCA_BC_lay.png", width = 8, height =  8, units = "cm", dpi = 600, scale = 1.5)


# gsea1 -------------------------------------------------------------------

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
           factor(levels = c( "A and B", "B and C", "A and C", "B_6m"))) %>%
  mutate(pathway = str_remove(pathway, pattern = "KEGG_"),
         order = case_when(
           pathway == "LYSOSOME" ~ 1,
           pathway %in% c("OTHER_GLYCAN_DEGRADATION",
                          "GLYCOSAMINOGLYCAN_DEGRADATION")  ~ 2,
           pathway %in% c("AMINO_SUGAR_AND_NUCLEOTIDE_SUGAR_METABOLISM",
                          "ANTIGEN_PROCESSING_AND_PRESENTATION") ~ 3,
           pathway %in% c("FATTY_ACID_METABOLISM",
                          "VIBRIO_CHOLERAE_INFECTION") ~ 4,
           pathway %in% c("SYSTEMIC_LUPUS_ERYTHEMATOSUS",
                          "SPHINGOLIPID_METABOLISM",
                          "DRUG_METABOLISM_CYTOCHROME_P450",
                          "SNARE_INTERACTIONS_IN_VESICULAR_TRANSPORT",
                          "COMPLEMENT_AND_COAGULATION_CASCADES",
                          "CELL_ADHESION_MOLECULES_CAMS") ~ 5,
           pathway %in% Aonly ~ 6,
           TRUE ~7),
         DE = FDR.Mixed < 0.05,
         Direction = case_when(
           Direction == "Up" ~ 1,
           Direction == "Down" ~ -1),
         log10p = -log10(PValue.Mixed),
         dirp = log10p*Direction
         ) %>%
  ggplot(aes(y = coef,
             x = reorder(pathway, order)) ) +
  geom_tile(aes(fill = dirp,
                alpha = FDR.Mixed < 0.05)) +
  geom_label(aes(label = signif(FDR.Mixed, digits = 2)),
             fill = NA) +
  facet_wrap(~MPS, scales = "free_y",
             ncol = 1) +
  scale_fill_continuous_diverging("Blue-Red 2", mid = 0)+
  scale_alpha_manual(values = c(0.1,1)) +
  ggpubr::theme_pubclean() +
  theme(panel.grid.major.y = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        axis.text.y = element_text(face = "bold", size = 15),
        axis.text.x = element_text(face = "bold", size = 15, angle = 315, hjust = 0),
        legend.position = "top",
        strip.text = element_text(face = "bold", size = 15, colour = "black",
                                  margin = margin(t = 20, b = 10)),
        strip.background = element_rect(fill = "white"),
        plot.margin = margin(t = 0, r = 5.5, b = 0, l = 0, unit = "cm")
  ) +
  labs(alpha = "FDRp < 0.05", fill = "-log10(p-val) * sign(logFC) ")+
  ggsave("output/plots4talk/firstGSEAres.png", width = 30, height =  16, units = "cm", dpi = 600, scale = 1.5)


design.ac %>% colnames() %>% .[2:3] %>%
  sapply(function(y) {
    logCPM.ac[-4614,] %>% #without hexb
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


# gsea2 -------------------------------------------------------------------
frykegg.all <-
  fryKEGG.ac %>%
  bind_rows() %>%
  mutate(MPS = "AC") %>%
  bind_rows(fryKEGG.ab %>%
              bind_rows() %>%
              mutate(MPS= "AB")) %>%
  bind_rows(fryKEGG.bc %>%
              bind_rows() %>%
              mutate(MPS = "BC")) %>%
  bind_rows(fryKEGG.b2 %>%
              bind_rows() %>%
              mutate(MPS = "B2"))

sigpaths <-
  frykegg.acb2 %>%
  dplyr::filter(FDR.Mixed < 0.05) %>%
  .$pathway %>%
  unique()

Aonly <- fryKEGG.ac$`MPS-IIIA` %>%
  dplyr::filter(FDR.Mixed < 0.05) %>%
  dplyr::filter(!pathway %in% (fryKEGG.b2$`MPS-IIIB` %>%
                                 dplyr::filter(FDR.Mixed < 0.05) %>%
                                 .$pathway)) %>%
  .$pathway %>%
  str_remove("KEGG_")

frykegg.all %>%
  dplyr::filter(pathway %in% sigpaths) %>%
  dplyr::filter(coef != "sgsh/+") %>%

  mutate(MPS = case_when(MPS == "AC" ~ "A and C",
                         MPS == "AB" ~ "A and B",
                         MPS == "BC" ~ "B and C",
                         MPS == "B2" ~ "B_6m") %>%
           factor(levels = c( "A and B", "B and C", "A and C", "B_6m"))) %>%
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
  dplyr::filter(MPS %in% c("A and C", "B_6m")) %>%
  mutate(coef4plot = case_when(
    coef == "MPS-IIIB" ~ "MPS-IIIB (6m)",
    TRUE ~ coef
  ) %>%
    factor(levels = c("MPS-IIIA", "MPS-IIIC", "MPS-IIIB (6m)"))) %>%
  ggplot(aes(y = coef4plot,
             x = reorder(pathway, order)) ) +
  geom_tile(aes(fill = -log10(PValue.Mixed),
                alpha = FDR.Mixed < 0.05)) +
  geom_label(aes(label = signif(FDR.Mixed, digits = 2)),
             fill = NA) +
  scale_fill_viridis_c() +
  scale_alpha_manual(values = c(0.1,1)) +
  ggpubr::theme_pubclean() +
  theme(panel.grid.major.y = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        strip.text = element_text(face = "bold", size = 15),
        axis.text.y = element_text(face = "bold", size = 15),
        axis.text.x = element_text(face = "bold", size = 15, angle = 315, hjust = 0),
        legend.position = "top",
        plot.margin = margin(t = 0, r = 5.5, b = 0, l = 0, unit = "cm")


  ) +
  ggsave("output/plots4talk/3gSEAres.png", width = 40, height =  15, units = "cm", dpi = 600, scale = 1.5)


# CELLTYPE HEATMAP --------------------------------------------------------

frycell.all <-
  celltype.ac %>%
  bind_rows() %>%
  mutate(MPS = "A and C") %>%
  bind_rows(frycell.all %>%
              bind_rows() %>%
              mutate(MPS= "A and B")) %>%
  bind_rows(celltype.bc %>%
              bind_rows() %>%
              mutate(MPS = "B and C")) %>%
  bind_rows(celltype.b2 %>%
              bind_rows() %>%
              mutate(MPS = "B2"))



celltype.ac %>%
  bind_rows() %>%
  mutate(MPS = "A and C",
         Direction = case_when(Direction == "Up" ~ 1,
           Direction == "Down" ~ -1),
         log10p = -log10(PValue),
         dirp = log10p*Direction) %>%
  ggplot(aes(x = coef,
             y = reorder(pathway, PValue))
  ) +
  geom_tile(aes(fill = dirp,
                alpha = FDR < 0.05),
            colour = "black") +
  geom_label(aes(label = signif(FDR, digits = 2)),
  fill = NA) +
  scale_fill_continuous_diverging("Blue-Red 2", mid = 0) +
  scale_alpha_manual(values = c(0.1,1)) +
  ggpubr::theme_pubclean() +
  theme(panel.grid.major.y = element_blank(),
        legend.position = "right",
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        axis.text.y = element_text(face = "bold", size = 15),
        axis.text.x = element_text(face = "bold", size = 15, angle = 315, hjust = 0.1, vjust = 0.2)
  ) +
ggsave("output/plots/celltype_heatmap.png", width = 10, height = 16, units = "cm", dpi = 500, scale = 2)



# ~~~ HEATMAPS ~~~ --------------------------------------------------------
# # LYSO ------------------------------------------------------------------
png("output/plots4talk/lyso.png", width = 60, height = 20, units = "cm", res = 200)
toptables_cqn_abc %>%
  dplyr::filter(gene_id %in% KEGG$KEGG_LYSOSOME) %>%
  dplyr::filter(coef != "sgsh/+") %>%
  mutate(coef_mps = paste0(coef, " (", mps, ")")) %>%
  dplyr::select(gene_name, logFC, coef_mps) %>%
  dplyr::distinct(gene_name, coef_mps, .keep_all = T) %>%
  spread(key = "coef_mps", value = "logFC") %>%
  column_to_rownames("gene_name") %>%
  t() %>%
  pheatmap(color = colorRampPalette(rev(brewer.pal(n = 5,
                                                   name = "RdBu")))(200),
           main = "KEGG_LYSOSOME",
           breaks = c(seq(-.8,.8, by = 1.6/200)),
           cellheight = 30, cellwidth = 12,
           angle_col = 45,
           treeheight_row = 10,
           #annotation_col = anno_celltype
           #treeheight_col = 20
  )
dev.off()

plotpearson = function(geneset, # which geneset
                       x.axis,
                       y.axis,
                       xcol,
                       ycol
                       ) {
toptables_cqn_abc %>%
  dplyr::filter(gene_id %in% geneset) %>%
  dplyr::select(gene_name, coef, logFC) %>%
  spread(key = "coef", value = "logFC") %>%
  ggscatter(
    x = x.axis,
    y = y.axis,
    alpha = 0.9,
    add = "reg.line",  # Add regressin line
    add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
    conf.int = TRUE # Add confidence interval
  ) +
  stat_cor(
    method = "pearson", label.y = 0.9
  )+
  theme(aspect.ratio = 1,
        axis.title.x = element_text(colour = xcol, face = "bold"),
        axis.title.y = element_text(colour = ycol, face = "bold")
  )
}



a <- plotpearson(geneset = KEGG$KEGG_RETINOL_METABOLISM,
            x.axis = "MPS-IIIA", xcol = "#ff005e",
            y.axis = "MPS-IIIC", ycol = "#00b064")

b <- plotpearson(geneset = KEGG$KEGG_RETINOL_METABOLISM,
            x.axis = "MPS-IIIA", xcol = "#ff005e",
            y.axis = "MPS-IIIB", ycol = "#ffae00")

c <- plotpearson(geneset = KEGG$KEGG_RETINOL_METABOLISM,
            x.axis = "MPS-IIIB", xcol = "#ffae00",
            y.axis = "MPS-IIIC", ycol = "#00b064")

a + b + c +
  plot_annotation(title = 'KEGG_AMINO_SUGAR_AND_NUCLEOTIDE_SUGAR_METABOLISM' ) +
  ggsave("output/plots4talk/kegg_aminosugorrelation.png", width = 18, height =  15, units = "cm", dpi = 600, scale = 1.5)




# IRE ---------------------------------------------------------------------

fryIRE.ac %>%
  bind_rows() %>%
  ggplot(aes(x = coef,
             y = pathway )) +
  geom_tile(aes(fill = -log10(PValue.Mixed),
                alpha = FDR.Mixed < 0.05),
            colour = "black") +
  geom_label(aes(label = signif(FDR.Mixed, digits = 2)),
             fill = NA) +
  scale_fill_viridis_c() +
  scale_alpha_manual(values = c(0.2,1)) +
  ggpubr::theme_pubclean() +
  theme(panel.grid.major.y = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        text = element_text(face = "bold", size = 15),
        axis.text.x = element_text(face = "bold", size = 15, angle = 315, hjust = 0.1, vjust = 0.2),
        legend.position = "right"
  ) +
  ggsave("output/plots4talk/ire.png", width = 5, height =  5, units = "cm", dpi = 600, scale = 3)


png("output/plots4talk/ire3_all_A_Heatmap.png", width = 50, height = 20, units = "cm", res = 100)
toptables_cqn.ac$`MPS-IIIA` %>%
  bind_rows() %>%
  dplyr::filter(gene_id %in% ireGenes$ire3_all) %>%
  dplyr::select(gene_name, logFC, coef) %>%
  dplyr::distinct(gene_name, coef, .keep_all = T) %>%
  spread(key = "coef", value = "logFC") %>%
  column_to_rownames("gene_name") %>%
  na.omit %>%
  arrange(`MPS-IIIA`) %>%
  t() %>%
  pheatmap(color = colorRampPalette(rev(brewer.pal(n = 5,
                                                   name = "RdBu")))(100),
           main = "ire3_all",
          breaks = seq(-0.75,0.75, by = 1.5/100),
           cellheight = 100,
           #cellwidth = 12,
          cluster_rows = F,
          cluster_cols = F,
           angle_col = 45,
           border_color = NA,
           show_colnames = F,
           treeheight_row = 10,
           treeheight_col = 20
  )

dev.off()

design.ac %>% colnames() %>% .[2:3] %>%
  sapply(function(y) {
    logCPM.ac %>%
      roast(
        index = ireGenes,
        design = design.ac,
        contrast = y,
        sort = "mixed"
      ) %>%
      rownames_to_column("pathway") %>%
      as_tibble() %>%
      mutate(coef = y)
  }, simplify = FALSE) %>%
  bind_rows() %>%
  dplyr::select(pathway, PropDown, PropUp,coef) %>%
  left_join(fryIRE.ac %>% bind_rows(), by = c("pathway", "coef")) %>%
  mutate(numUp = PropUp*NGenes,
         numdn = PropDown*NGenes) %>%
  gather(key = "dir", value = "numberofgenes", starts_with("num")) %>%
  ggplot(aes(y = pathway)) +
  geom_col(aes(x = NGenes)) +
  geom_col(aes(x = numberofgenes, fill = dir)) +
  facet_wrap(~coef)+
  coord_cartesian(xlim = c(0, 250) )

# celltype ----------------------------------------------------------------

celltype.ac %>%
  bind_rows() %>%
  arrange(PValue) %>%
  ggplot(aes(x = coef,
             y = reorder(pathway, PValue))
  ) +
  geom_tile(aes(fill = -log10(PValue),
                alpha = FDR < 0.05),
            colour = "black") +
  geom_label(aes(label = signif(FDR, digits = 2)
  ),
  fill = NA) +
  scale_fill_viridis_c() +
  scale_alpha_manual(values = c(0.1,1)) +
  ggpubr::theme_pubclean() +
  theme(panel.grid.major.y = element_blank(),
        legend.position = "right",
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        axis.text.y = element_text(face = "bold", size = 15),
        axis.text.x = element_text(face = "bold", size = 15, angle = 315, hjust = 0.1, vjust = 0.2)
  )+
  ggsave("output/plots4talk/celltype.png", width = 5, height =  5, units = "cm", dpi = 600, scale = 5)


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



png("output/plots4talk/NSC_A_Heatmap.png", width = 50, height = 20, units = "cm", res = 100)
toptables_cqn.ac %>%
  bind_rows() %>%
  dplyr::filter(gene_id %in% cell_type_markers$`Neural stem cell`) %>%
  dplyr::select(gene_name, logFC, coef) %>%
  dplyr::distinct(gene_name, coef, .keep_all = T) %>%
  spread(key = "coef", value = "logFC") %>%
  column_to_rownames("gene_name") %>%
  na.omit %>%
  t() %>%
  pheatmap(color = colorRampPalette(rev(brewer.pal(n = 5,
                                                   name = "RdBu")))(100),
           main = "Neural stem cell",
           breaks = seq(-1,1, by = 2/100),
           cellheight = 100,
           annotation_col = anno_celltype %>% dplyr::select(-`Neural stem cell`) ,
           cluster_rows = F,
           cluster_cols = T,
           clustering_method = "ward",
           angle_col = 45,
           border_color = NA,
           show_colnames = T,
           treeheight_row = 10,
           treeheight_col = 20
  )
dev.off()


tibble(mps = c("A", "B", "C"),
       meanage = c(	15.22,	18.91,	23.43),
       stdev = c(4.22,	7.33,	9.47
       )) %>%
  ggplot(aes(x = mps, y = meanage)) +
  geom_col(aes(fill = mps)) +
  geom_errorbar(aes(ymin = meanage - stdev, ymax = meanage + stdev),
                width = 0.5,
                alpha = 0.8) +
  scale_fill_manual(values = c("#ff005e", "#ffae00", "#00b064")) +
  labs(y = "Age (years)") +
  theme(legend.position = "none") +
  ggsave("output/plots4talk/death plot.png", width = 10, height = 10, units = "cm", dpi = 600, scale = 1)
