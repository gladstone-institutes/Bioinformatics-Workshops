cnts_casein <- cpm(y)[y$genes$Symbol == "Csn1s2b", ]
cnts_rndm <- cpm(y)[100, ]

df <- data.frame(Casein = cnts_casein,
                 Random = cnts_rndm,
                 targets[, c("CellType", "Status")])

p <- ggplot(df, aes(x=Random, y=1, 
               color=CellType, shape=Status)) + 
  geom_point(size=2.5) + 
  scale_x_log10() +
  scale_y_log10() +
  theme_classic() +
  theme(axis.title.y = element_blank(), 
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank(),
        legend.position = "top",
        legend.box = "vertical",
        plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"))
ggsave("Random.png",
       p,
       width = 3,
       height = 3)

p <- ggplot(df, aes(x=1, y=Random, 
               color=CellType, shape=Status)) + 
  geom_point(size=2.5) + 
  scale_x_log10() +
  scale_y_log10()  +
  theme(axis.title.x = element_blank(), 
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

ggsave("Random.png",
       p)

  xlab(paste("PC1:", pc1_var, "% variance")) + 
  ylab(paste("PC2:", pc2_var, "% variance"))