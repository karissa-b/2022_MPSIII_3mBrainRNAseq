# plots


x.ab <-  readRDS("data/R_objects/dge_AB.rds")
fryKEGG.ab <-  readRDS("data/R_objects/AB_fryKEGG.rds")
fryIRE.ab <- readRDS("data/R_objects/AB_fryIRE.rds")
celltype.ab <- readRDS("data/R_objects/AB_frycellTypes.rds")
toptables_cqn.ab <-  readRDS("data/R_objects/AB_toptab_cqn.rds")


x.ac <- readRDS("data/R_objects/dge_ac.rds")
toptables_cqn.ac <- readRDS("data/R_objects/toptab_ac_cqn.rds")
fryKEGG.ac <- readRDS("data/R_objects/AC_fryKEGG.rds")
fryIRE.ac <- readRDS("data/R_objects/AC_fryIRE.rds")
celltype.ac <- readRDS("data/R_objects/AC_frycellTypes.rds")


x.bc <- readRDS("data/R_objects/dge_bc.rds")
toptables_cqn.bc<-readRDS("data/R_objects/toptab_bc_cqn.rds")
fryKEGG.bc <- readRDS("data/R_objects/BC_fryKEGG.rds")
fryIRE.bc<-readRDS("data/R_objects/BC_fryIRE.rds")
celltype.bc <- readRDS("data/R_objects/BC_frycellTypes.rds")

KEGG_all <- fryKEGG.ac %>%
  bind_rows() %>%
  mutate(family = "A and C") %>%
  bind_rows(fryKEGG.ab %>%
              bind_rows() %>%
              mutate(family = "A and B")
            ) %>%
  bind_rows(fryKEGG.bc %>%
              bind_rows() %>%
              mutate(family = "B and C")
            )

ire_all <- fryIRE.ac %>%
  bind_rows() %>%
  mutate(family = "A and C") %>%
  bind_rows(fryIRE.ab %>%
              bind_rows() %>%
              mutate(family = "A and B")
  ) %>%
  bind_rows(fryIRE.bc %>%
              bind_rows() %>%
              mutate(family = "B and C")
  )

celltype_all <- celltype.ac %>%
  bind_rows() %>%
  mutate(family = "A and C") %>%
  bind_rows(celltype.ab %>%
              bind_rows() %>%
              mutate(family = "A and B")
  ) %>%
  bind_rows(celltype.bc %>%
              bind_rows() %>%
              mutate(family = "B and C")
  )


sig_paths <- KEGG_all %>%
  dplyr::filter(FDR.Mixed < 0.05) %>%
  .$pathway %>%
  unique

# KEGG --------------------------------------------------------------------
KEGG_all %>%
  dplyr::filter(pathway %in% sig_paths) %>%
  ggplot(aes(x = coef, y = pathway)) +
  geom_tile(aes(size = -log10(PValue.Mixed),
                 fill = -log10(PValue.Mixed),
                 alpha = FDR.Mixed < 0.05)) +
  geom_label(aes(label = signif(FDR.Mixed, digits = 2)),
             fill= NA,
             label.size = 0) +
  facet_wrap(~family, scales = "free_x") +
  scale_alpha_manual(values = c(0.5, 1)) +
  scale_fill_viridis_c() +
  theme_pubclean() +
  theme(legend.position = "right") +
  ggsave("output/SCFplots/KEGG_pvals.png", width = 20, height = 20, units = "cm", dpi = 300, scale = 2)

#

# IRE ---------------------------------------------------------------------
ire_all %>%
  ggplot(aes(x = coef, y = pathway)) +
  geom_tile(aes(size = -log10(PValue.Mixed),
                fill = -log10(PValue.Mixed),
                alpha = FDR.Mixed < 0.05)) +
  geom_label(aes(label = signif(FDR.Mixed, digits = 2)),
             fill= NA,
             label.size = 0) +
  facet_wrap(~family, scales = "free_x") +
  scale_alpha_manual(values = c(0.5, 1)) +
  scale_fill_viridis_c() +
  theme_pubclean() +
  theme(legend.position = "right")+
  ggsave("output/SCFplots/ire_pvals.png", width = 15, height = 5, units = "cm", dpi = 300, scale = 2)

# cell type props

celltype_all %>%
  ggplot(aes(x = coef, y = pathway)) +
  geom_tile(aes(size = -log10(PValue.Mixed),
                fill = -log10(PValue.Mixed),
                alpha = FDR.Mixed < 0.05)) +
  geom_label(aes(label = signif(FDR.Mixed, digits = 2)),
             fill= NA,
             label.size = 0) +
  facet_wrap(~family, scales = "free_x") +
  scale_alpha_manual(values = c(0.5, 1)) +
  scale_fill_viridis_c() +
  theme_pubclean() +
  theme(legend.position = "right")+
  ggsave("output/SCFplots/celltype_pvals.png", width = 15, height = 5, units = "cm", dpi = 300, scale = 2)

colanno <- x.ac$samples %>%
  dplyr::select(Tank, genotype)


logCPM.ac[cell_type_markers$`Oligodendrocyte-enriched`,] %>%
  pheatmap(scale = "row",
           annotation_col = colanno,
           annotation_colors = list(
             genotype = c("WT" = "#440154FF", "MPS-IIIA" = "#21908CFF", "MPS-IIIC" = "#FDE725FF")
           ),
           color = colorRampPalette(rev(brewer.pal(n = 7, name =
                                                    "RdBu")))(100),
           show_rownames = F,

           main = "Expression of oligodendrocyte marker genes in WT, MPS-IIIA and MPS-IIIC zebrafish brains"
             )

