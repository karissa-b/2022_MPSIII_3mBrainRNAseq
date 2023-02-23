# ac only

# heatmaps ----------------------------------------------------------------

png("output/plots/ac_only/oligodenHeatmap.png", width = 50, height = 20, units = "cm", res = 100)
toptables_cqn.ac %>%
  bind_rows() %>%
  dplyr::filter(gene_id %in% cell_type_markers$Oligodendrocyte) %>%
  dplyr::select(gene_name, logFC, coef) %>%
  dplyr::distinct(gene_name, coef, .keep_all = T) %>%
  spread(key = "coef", value = "logFC") %>%
  column_to_rownames("gene_name") %>%
  t() %>%
  pheatmap(color = colorRampPalette(rev(brewer.pal(n = 5,
                                                   name = "RdBu")))(100),
           main = "Oligodendrocytes",
           breaks = c(seq(-1,1, by = 2/100)),
           cellheight = 30, cellwidth = 14,
           fontsize = 14,

           angle_col = 45,
           treeheight_row = 10,
           #annotation_col = anno_celltype
           #treeheight_col = 20
  )
dev.off()

  png("output/plots/ac_only/neuralstemHeatmap.png", width = 50, height = 20, units = "cm", res = 100)
  toptables_cqn.ac %>%
    bind_rows() %>%
    dplyr::filter(gene_id %in% cell_type_markers$`Neural stem cell`) %>%
    dplyr::select(gene_name, logFC, coef) %>%
    dplyr::distinct(gene_name, coef, .keep_all = T) %>%
    spread(key = "coef", value = "logFC") %>%
    column_to_rownames("gene_name") %>%
    t() %>%
    pheatmap(color = colorRampPalette(rev(brewer.pal(n = 5,
                                                     name = "RdBu")))(100),
             main = "Neural stem cells",
             breaks = c(seq(-1,1, by = 2/100)),
             cellheight = 30, cellwidth = 11,
             fontsize = 11,

             angle_col = 45,
             treeheight_row = 10,
             #annotation_col = anno_celltype
             #treeheight_col = 20
    )
  dev.off()

  png("output/plots/ac_only/microglia.png", width = 50, height = 20, units = "cm", res = 100)
  toptables_cqn.ac %>%
    bind_rows() %>%
    dplyr::filter(gene_id %in% cell_type_markers$`Microglia_apoc1 high`) %>%
    dplyr::select(gene_name, logFC, coef) %>%
    dplyr::distinct(gene_name, coef, .keep_all = T) %>%
    spread(key = "coef", value = "logFC") %>%
    column_to_rownames("gene_name") %>%
    t() %>%
    pheatmap(color = colorRampPalette(rev(brewer.pal(n = 5,
                                                     name = "RdBu")))(100),

             breaks = c(seq(-1,1, by = 2/100)),
             cellheight = 30, cellwidth = 11,
             fontsize = 11,

             angle_col = 45,
             treeheight_row = 10,
             #annotation_col = anno_celltype
             #treeheight_col = 20
    )
  dev.off()


  png("output/plots/ac_only/oligslahigh.png", width = 50, height = 20, units = "cm", res = 100)
  toptables_cqn.ac %>%
    bind_rows() %>%
    dplyr::filter(gene_id %in% cell_type_markers$`Oligodendrocyte_sla high`) %>%
    dplyr::select(gene_name, logFC, coef) %>%
    dplyr::distinct(gene_name, coef, .keep_all = T) %>%
    spread(key = "coef", value = "logFC") %>%
    column_to_rownames("gene_name") %>%
    t() %>%
    pheatmap(color = colorRampPalette(rev(brewer.pal(n = 5,
                                                     name = "RdBu")))(100),

             breaks = c(seq(-1,1, by = 2/100)),
             cellheight = 30, cellwidth = 11,
             fontsize = 11,

             angle_col = 45,
             treeheight_row = 10,
             #annotation_col = anno_celltype
             #treeheight_col = 20
    )
  dev.off()


# ire ---------------------------------------------------------------------
  png("output/plots/ac_only/ire_3all.png", width = 30, height = 50, units = "cm", res = 100)
  toptables_cqn.ac %>%
    bind_rows() %>%
    dplyr::filter(gene_id %in% ireGenes$ire3_all) %>%
    #dplyr::filter(abs(logFC) >= 0.3 ) %>%
    dplyr::select(gene_name, logFC, coef) %>%
    dplyr::distinct(gene_name, coef, .keep_all = T) %>%
    spread(key = "coef", value = "logFC") %>%
    column_to_rownames("gene_name") %>%
    na.omit %>%
    pheatmap(color = colorRampPalette(rev(brewer.pal(n = 5,
                                                     name = "RdBu")))(100),
             main = "ire3_all",
             breaks = seq(-1,1,by = 1*2/100),
             # breaks = c(seq(min(.), 0, length.out=ceiling(100/2) + 1),
             #            seq(max(.)/100, max(.), length.out=50)),
             #cellheight = 30,
             cellwidth = 30,
             angle_col = 45,
             border_color = NA,
             show_rownames = F,
             treeheight_row = 50,
             treeheight_col = 20
    )
dev.off()


png("output/plots/ac_only/ire_5all.png", width = 30, height = 20, units = "cm", res = 100)
toptables_cqn.ac %>%
  bind_rows() %>%
  dplyr::filter(gene_id %in% ireGenes$ire5_all) %>%
  #dplyr::filter(abs(logFC) >= 0.3 ) %>%
  dplyr::select(gene_name, logFC, coef) %>%
  dplyr::distinct(gene_name, coef, .keep_all = T) %>%
  spread(key = "coef", value = "logFC") %>%
  column_to_rownames("gene_name") %>%
  na.omit %>%
  pheatmap(color = colorRampPalette(rev(brewer.pal(n = 5,
                                                   name = "RdBu")))(100),
           main = "ire5_all",
           breaks = seq(-1,1,by = 1*2/100),
           # breaks = c(seq(min(.), 0, length.out=ceiling(100/2) + 1),
           #            seq(max(.)/100, max(.), length.out=50)),
           #cellheight = 30,
           cellwidth = 30,
           angle_col = 45,
           border_color = NA,
           show_rownames = F,
           treeheight_row = 50,
           treeheight_col = 20
  )
dev.off()


tf.ftgenes <- x.ac$genes %>%
  as_tibble() %>%
  dplyr::filter(grepl(description, pattern = "transferrin|ferritin")) %>%
  .$gene_id

logCPM.ac[tf.ftgenes,] %>%
  as.data.frame() %>%
  rownames_to_column("gene_id") %>%
  gather(key = "sample", value = "logCPM", starts_with("AC")) %>%
  left_join(x.ac$samples) %>%
  left_join(x.ac$genes) %>% view
  dplyr::filter(genotype != "MPS-IIIC") %>%
  ggplot(aes(x = genotype, y = logCPM)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(aes(shape = Sex)) +
  facet_wrap(~gene_name, scales = "free_y")




# hets --------------------------------------------------------------------

ggarrange(
  logCPM.ab %>%
    as.data.frame() %>%
    rownames_to_column("gene_id") %>%
    gather(key = "sample", value = "logCPM", starts_with("AB")) %>%
    left_join(x.ab$samples) %>%
    left_join(x.ab$genes) %>%
    dplyr::filter(genotype != "MPS-IIIB") %>%
    dplyr::select(gene_id, sample, logCPM) %>%
    spread(key = "sample", value = "logCPM") %>%
    column_to_rownames("gene_id") %>%
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
                      "#ff7aab", # sgsh het
                       "#ff005e"))


  fryKEGG.ab$`sgsh/+` %>%
    dplyr::select(pathway,pval = PValue.Mixed, FDR = FDR.Mixed) %>%
    head(10) %>%
    ggplot(aes(x = -log10(pval), y = reorder(pathway, -pval))) +
    geom_col() +
    ggsave("output/plots/ac_only/topKegghets.png", width = 10, height = 5, units = "cm",
           dpi = 400, scale = 2.5)




  ggsave("output/plots/ac_only/topKegghets.png", width = 10, height = 5, units = "cm",
           dpi = 400, scale = 2.5)


    ggplot(aes( ))
