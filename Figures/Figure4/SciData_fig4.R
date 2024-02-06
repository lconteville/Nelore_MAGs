library("tidyverse")
library("ggtree")
library(phangorn)

tree <- read.tree("MAGs_data_input/Data/gtdbtk_archaea.tree")
tree

tree2<-phangorn::midpoint(tree)

tree2 <- reorder(tree2)

tree2$tip.label

fezes <- tree2$tip.label[grep("fecal", tree$tip.label)]
rumen <- tree2$tip.label[grep("ruminal", tree$tip.label)]

meta = data.frame(label = tree2$tip.label,site = NA)

meta$site[meta$label %in% rumen] <- "Ruminal MAGs"
meta$site[meta$label %in% fezes] <- "Fecal MAGs"

meta$site <- factor(meta$site, levels = c("Ruminal MAGs","Fecal MAGs"))


gg <- ggtree(tree2,size=.8) %<+% meta +
  geom_tippoint(aes(x=x+0.007, color = site), size=2.4)+
  scale_colour_manual(values = c("#CC0B0B","#0B7ACC"),na.translate=FALSE)+
  theme(legend.title = element_blank(), 
        legend.position = c(0.35, 0.83),
        legend.background = element_blank(),
        legend.key = element_rect(fill="#DAC8E6"),
        legend.text = element_text(size = 10.3, color = "black"),
        legend.key.size = unit(1.3, 'lines')) +
  guides(color = guide_legend(override.aes=list(size=6.1)))

gg

MRCA(tree2, "GB_GCA_017431845.1","GB_GCA_017430835.1")
MRCA(tree2, "ruminal_MAG_234","GB_GCA_900314635.1")
MRCA(tree2, "GB_GCA_017652345.1","GB_GCA_016291275.1")

gg <- gg %>% collapse(node=101)

gg <- gg %>% collapse(node=112) + 
  geom_highlight(node=97, fill="#DAC8E6", alpha = 1,to.bottom = TRUE,extend = 0.05)

gg2 <- gg +
  annotate("text", x = .37, y = 35,  size = 3.5, 
           label = "Methanobrevibacter", angle = 90)
gg2

ggsave("MAGs_data_input/Figures/Figure4/SciData_fig4.png", 
       gg2, width = 1600, height = 1600, dpi = 300, units = "px")


