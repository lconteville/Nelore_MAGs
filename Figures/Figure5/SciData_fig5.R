library(data.table)
library(dplyr)
library(ggvenn)
library(patchwork)
library(ggpubr)

#### FEZES vs RUMEN 99 ANI ####

fzlr99 <- read.csv("MAGs_data_input/Data/fzlr_Cdb_99.csv", header = TRUE)

fzlr99 <- fzlr99[,c(1,2)]

fzlr99_list <- split(fzlr99$genome,fzlr99$secondary_cluster)

fzlrshared = c()
numshared = c()
ruminal_unique = c()
fecal_unique = c()

for (i in 1:length(fzlr99_list)){
  if(any(fzlr99_list[[i]] %like% "ruminal_")){
    if(any(fzlr99_list[[i]] %like% "fecal_")){
      fzlrshared <- append(fzlrshared, fzlr99_list[[i]])
      numshared  <- append(numshared, i)
    }else{
      ruminal_unique <- append(ruminal_unique, fzlr99_list[[i]])
    }
  }else{
    fecal_unique <- append(fecal_unique, fzlr99_list[[i]])
  }
}

ruminal_shared <- fzlrshared[grepl("ruminal_",fzlrshared)]
fecal_shared <- fzlrshared[grepl("fecal_",fzlrshared)]

ruminal_list <- 1:sum(length(ruminal_unique),length(fzlrshared))
fecal_list <- (length(ruminal_list)-length(fzlrshared)+1):sum(length(ruminal_list)+length(fecal_unique))

b <- list(`Ruminal MAGs` = ruminal_list,
          `Fecal MAGs` = fecal_list)

plot_99fzlr <- ggvenn(b, c("Ruminal MAGs", "Fecal MAGs"),
                      show_percentage = FALSE,
                      stroke_size=1.5,
                      text_size = 5,
                      set_name_size = 4.5,
                      fill_color = c("#CC0B0B","#0B7ACC"),
                      fill_alpha = 0.4)  +
  geom_segment(aes(x = 0, y = -.3, xend = 0, yend =-1.2),
               arrow = arrow(length = unit(0.4, "cm")), size = .9) +
  annotate(geom="text", x=0, y=-1.5,
           label= paste(length(ruminal_shared)," Ruminal MAGs\n", length(fecal_shared)," Fecal MAGs"),
           color="black", size = 4) +
  expand_limits(y = c(-2,1)) +
  theme(plot.title = element_text(hjust = 0.5, vjust = 5, size = 13, face="bold"),
        plot.margin = margin(20, 2, 0, 2))+
  ggtitle("99% ANI")

plot_99fzlr

#### FEZES vs RUMEN 95 ANI ####

fzlr95 <- read.csv("MAGs_data_input/Data/fzlr_Cdb_95.csv", header = TRUE)

fzlr95 <- fzlr95[,c(1,2)]

fzlr95_list <- split(fzlr95$genome,fzlr95$secondary_cluster)

fzlrshared95 = c()
numshared95 = c()
ruminal_unique95 = c()
fecal_unique95 = c()

for (i in 1:length(fzlr95_list)){
  if(any(fzlr95_list[[i]] %like% "ruminal_")){
    if(any(fzlr95_list[[i]] %like% "fecal_")){
      fzlrshared95 <- append(fzlrshared95, fzlr95_list[[i]])
      numshared95  <- append(numshared95, i)
    }else{
      ruminal_unique95 <- append(ruminal_unique95, fzlr95_list[[i]])
    }
  }else{
    fecal_unique95 <- append(fecal_unique95, fzlr95_list[[i]])
  }
}

ruminal_shared95 <- fzlrshared95[grepl("ruminal_",fzlrshared95)]
fecal_shared95 <- fzlrshared95[grepl("fecal_",fzlrshared95)]

ruminal_list95 <- 1:sum(length(ruminal_unique95),length(fzlrshared95))
fecal_list95 <- (length(ruminal_list95)-length(fzlrshared95)+1):sum(length(ruminal_list95)+length(fecal_unique95))

c <- list(`Ruminal MAGs` = ruminal_list95,
          `Fecal MAGs` = fecal_list95)

plot_95fzlr <- ggvenn(c, c("Ruminal MAGs", "Fecal MAGs"),
                      show_percentage = FALSE,
                      stroke_size=1.5,
                      text_size = 5,
                      set_name_size = 4.5,
                      fill_color = c("#CC0B0B","#0B7ACC"),
                      fill_alpha = 0.4)  +
  geom_segment(aes(x = 0, y = -.3, xend = 0, yend =-1.2),
               arrow = arrow(length = unit(0.4, "cm")), size = .9) +
  annotate(geom="text", x=0, y=-1.5,
           label= paste(length(ruminal_shared95)," Ruminal MAGs\n", length(fecal_shared95)," Fecal MAGs"),
           color="black", size = 4) + 
  expand_limits(y = c(-2,1)) +
  theme(plot.title = element_text(hjust = 0.5, vjust = 5, size = 13, face="bold"),
        plot.margin = margin(20, 2, 0, 2))+
  ggtitle("95% ANI")

plot_95fzlr

#### MAGS vs SUPERSET ####

cdb <- read.csv("MAGs_data_input/Data/superset_mags_95_Cdb.csv", header = TRUE)
#View(cdb)

cdb <- cdb[,c(1,2)]

cdb$genome[cdb$genome %like% "ruminal_MAG"]<- paste0("Conteville_et_al_", cdb$genome[cdb$genome %like% "ruminal_MAG"])

cdb$genome[cdb$genome %like% "fecal_MAG"] <- paste0("Conteville_et_al_", cdb$genome[cdb$genome %like% "fecal_MAG"])

cdb$genome[!cdb$genome %like% "Conteville_et_al"] <- paste0("Others_", cdb$genome[!cdb$genome %like% "Conteville_et_al"])

cdb_list <- split(cdb$genome,cdb$secondary_cluster)
#View(cdb_list)

shared = c()
Conteville_unique = c()
Others_unique = c()

for (i in 1:length(cdb_list)){
  if(any(cdb_list[[i]] %like% "Conteville_et_al")){
    if(any(cdb_list[[i]] %like% "Others_")){
      shared <- append(shared, cdb_list[[i]])
    }else{
      Conteville_unique <- append(Conteville_unique, cdb_list[[i]])
    }
  }else{
    Others_unique <- append(Others_unique, cdb_list[[i]])
  }
}

Conteville_shared <- shared[grepl("Conteville_et_al",shared)]

Others_shared <- shared[grepl("Others_",shared)]

Others_list <- 1:sum(length(Others_unique),length(shared))

Conteville_list <- (length(Others_list)-length(shared)+1):sum(length(Others_list)+length(Conteville_unique))

a <- list(`Nelore MAGs` = Conteville_list,
          `Public Set` = Others_list)

plot_95set <- ggvenn(a, c("Nelore MAGs", "Public Set"),
                     show_percentage = FALSE,
                     stroke_size=1.5,
                     text_size = 5,
                     set_name_size = 4.5,
                     fill_color = c("#2CCFB3","#9D7ECE"),
                     fill_alpha = 0.6)  +
  geom_segment(aes(x = 0, y = -.3, xend = 0, yend =-1.2),
               arrow = arrow(length = unit(0.4, "cm")), size = .9) +
  annotate(geom="text", x=0, y=-1.5,
           label=paste(" ",length(Conteville_shared)," Nelore MAGs\n", length(Others_shared)," Public Set"),
           color="black", size = 4) +
  expand_limits(y = c(-2,1)) +
  theme(plot.title = element_text(hjust = 0.5, vjust = 5, size = 13, face="bold"),
        plot.margin = margin(20, 2, 0, 2))+
  ggtitle("95% ANI")

plot_95set

#### JUNTAR DIAGRAMAS DE VEN #####

a_b_c <- ggarrange(ggarrange(plot_99fzlr, plot_95fzlr, labels = "AUTO"),
                   plot_95set, ncol = 1, labels = c("", "C")) +
  theme(panel.background = element_rect(fill = "white",colour = "white"))



ggsave("MAGs_data_input/Figures/Figure5/SciData_fig5.png", 
       a_b_c, width = 2000, height = 2000, dpi = 300, units = "px")

