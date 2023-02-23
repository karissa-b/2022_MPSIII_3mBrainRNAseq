makePCAplot = function(dge,
                       logCPM.file,
                       genoCols) {

  ggarrange(
logCPM.file %>%
  t() %>%
  prcomp() %>%
  autoplot(data = tibble(sample = rownames(.$x)) %>%
             left_join(dge$samples) %>%
             mutate(P = str_extract(Tank, pattern = "P[1|2]"),
                    lay = str_extract(Tank, pattern = "P[1|2]_lay[1-3]"),
                    DOB = as.character(DOB)),
           colour = "genotype",
           size = 4
  ) +
  scale_color_manual(values = genoCols) +
  labs(colour = "Genotype") +
  ggtitle("Genotype"),

logCPM.file %>%
  t() %>%
  prcomp() %>%
  autoplot(data = tibble(sample = rownames(.$x)) %>%
             left_join(dge$samples) %>%
             mutate(P = str_extract(Tank, pattern = "P[1|2]"),
                    lay = str_extract(Tank, pattern = "P[1|2]_lay[1-3]"),
                    DOB = as.character(DOB)),
           colour = "P",
           size = 4
  ) +
  scale_color_brewer(palette = "Set2") +
  labs(colour = "Parent pair") +
  ggtitle("Parent pair"),

logCPM.file %>%
  t() %>%
  prcomp() %>%
  autoplot(data = tibble(sample = rownames(.$x)) %>%
             left_join(dge$samples) %>%
             mutate(P = str_extract(Tank, pattern = "P[1|2]"),
                    lay = str_extract(Tank, pattern = "P[1|2]_lay[1-3]"),
                    DOB = as.character(DOB)),
           colour = "DOB",
           size = 4
  ) +
  scale_color_brewer(palette = "Dark2") +
  labs(colour = "Clutch") +
  ggtitle("Clutch"),

logCPM.file %>%
  t() %>%
  prcomp() %>%
  autoplot(data = tibble(sample = rownames(.$x)) %>%
             left_join(dge$samples) %>%
             mutate(P = str_extract(Tank, pattern = "P[1|2]"),
                    lay = str_extract(Tank, pattern = "P[1|2]_lay[1-3]"),
                    DOB = as.character(DOB)),
           colour = "Tank",
           size = 4
  ) +
  scale_color_brewer(palette = "Dark2") +
  labs(colour = "Home tank"),

logCPM.file %>%
  t() %>%
  prcomp() %>%
  autoplot(data = tibble(sample = rownames(.$x)) %>%
             left_join(dge$samples) %>%
             mutate(P = str_extract(Tank, pattern = "P[1|2]"),
                    lay = str_extract(Tank, pattern = "P[1|2]_lay[1-3]"),
                    DOB = as.character(DOB)),
           colour = "RIN",
           size = 4
  ) +
  scale_color_viridis_c() +
  labs(colour = "RIN"),
labels = "AUTO"
)
}

makePCAplot(dge = x.bc,
            logCPM.file = logCPM.bc,
            genoCols =  c("#4000ff", #WT
                          "#ffae00", # B
                          "#00b064")
) +
  ggsave("output/plots/PCA plots/BC.png", width = 10, height = 5, units = "cm", dpi = 400, scale = 3)

makePCAplot(dge = x.ab,
            logCPM.file = logCPM.ab,
            genoCols =  c("#4000ff", #WT
                          "#ff7aab", # sgsh het
                          "#ff005e", # A
                          "#ffae00") #B
) +
  ggsave("output/plots/PCA plots/AB.png", width = 10, height = 5, units = "cm", dpi = 400, scale = 3)

makePCAplot(dge = x.ac,
            logCPM.file = logCPM.ac,
            genoCols =  c("#4000ff", #WT
                          "#ff005e", # A
                          "#00b064") #C
) +
ggsave("output/plots/PCA plots/AC.png", width = 10, height = 5, units = "cm", dpi = 400, scale = 3)


