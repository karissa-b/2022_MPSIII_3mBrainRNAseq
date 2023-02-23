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
       colour = "Genotype")
