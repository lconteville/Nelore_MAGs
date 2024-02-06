library("reshape2") 
library("ggplot2")
library("tidyverse")
library("scales")
library("ggpubr")
library(tidyverse)
library(pheatmap)
library(viridis)
library(patchwork)
library(cowplot)

#### THEME SET ####
theme_set(theme_bw()+ theme(
  axis.text = element_text(size=4, color = "black"),
  axis.title = element_text(size = 4.5, color = "black"),
  axis.ticks = element_line(linewidth=0.25),
  legend.text = element_text(size = 4, color = "black"),
  legend.title = element_text(size = 4.2, color = "black"),
  strip.background = element_blank(),
  strip.text = element_text(size = 5,colour="black"),
  plot.title = element_text(size = 6,colour="black",hjust = -0.05,face="bold"),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank()))

#### MAGs QUALITY ####
quality_summary <-  read.csv("MAGs_data_input/Data/mags_quality_summary.tsv", header=TRUE, sep="\t")

quality_summary$SampleType_Quality <- paste(quality_summary$SampleType,quality_summary$Quality_assessment, sep="_")

quality_summary$SampleType_Quality[quality_summary$SampleType_Quality == "Rumen_High-quality"] <- "HQ Ruminal MAGs (n = 497)"
quality_summary$SampleType_Quality[quality_summary$SampleType_Quality == "Feces_High-quality"] <- "HQ Fecal MAGs (n = 486)"

quality_summary$SampleType_Quality[quality_summary$SampleType_Quality == "Feces_Medium-quality"] <- "Medium-quality MAGs"
quality_summary$SampleType_Quality[quality_summary$SampleType_Quality == "Rumen_Medium-quality"] <- "Medium-quality MAGs"


quality_summary$SampleType_Quality <- factor(quality_summary$SampleType_Quality, levels=c("HQ Ruminal MAGs (n = 497)","HQ Fecal MAGs (n = 486)","Medium-quality MAGs"))


points <- ggplot(quality_summary, aes(x=Completeness, y=Contamination,
                                      colour=SampleType_Quality ))+ 
  geom_point(size=.15,alpha = 0.9)+
  theme(plot.margin = unit(c(0,4,2,0), "pt"),
        legend.title = element_blank(),
        legend.margin=margin(c(0, 0, 0, 0)), 
        legend.spacing.x = unit(1, 'pt'),
        legend.key.size = unit(10, 'pt'))+
  scale_colour_manual(values=c("#CC0B0B","#0B7ACC","grey"))+
  guides(colour=guide_legend(override.aes=list(size=3, alpha =1)))+ 
  ggtitle('A')

points

leg_points <- get_legend(points)

leg_points <- as_ggplot(leg_points)

hq <- quality_summary[quality_summary$Quality_assessment == c("High-quality"),]
hq$SampleType <- factor(hq$SampleType, levels=c("Rumen","Feces"))

m5 <- hq %>%
  # create a new variable from count
  mutate(countfactor=cut(contig_bp, breaks=c( 0, 1000000, 1500000, 2000000, 2500000, 3000000, 3500000, 4000000, 4500000, 5000000, 6000000),
                         labels=c("1","1.5", "2", "2.5", "3", "3.5", "4", "4.5","5", "6"))) %>%
  # change level order
  mutate(countfactor=factor(as.character(countfactor), levels=levels(countfactor)))

histo_size <- ggplot(m5, aes(x=countfactor, fill = SampleType)) + 
  geom_histogram(stat="count")+
  facet_wrap(SampleType~.)+
  scale_fill_manual(values=c("#CC0B0B","#0B7ACC"))+
  theme(legend.position = "none",
        strip.text = element_blank(),
        plot.margin = unit(c(2,4,0,0), "pt"))+
  xlab("MAGs Size")+
  ylab("Count")+ 
  ggtitle('B')

histo_size

contigsplot <- ggplot(hq, aes(x=SampleType, y=n_contigs, fill=SampleType)) +
  geom_boxplot(lwd=.3,outlier.shape = NA, alpha = 0.9) +
  scale_fill_manual(values=c("#CC0B0B","#0B7ACC")) +
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.margin = unit(c(2,0,2,4), "pt"))+
  scale_y_continuous(breaks=pretty_breaks(n=5))+
  ylab("# contigs")+ 
  ggtitle('C')

n50 <- ggplot(hq, aes(x=SampleType, y=ctg_N50  , fill=SampleType)) +
  geom_boxplot(lwd=.3,outlier.shape = NA, alpha = 0.9) +
  scale_fill_manual(values=c("#CC0B0B","#0B7ACC")) +
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.margin = unit(c(2,0,2,4), "pt"))+
  scale_y_continuous(breaks=pretty_breaks(n=6))+
  ylab("N50")

gcplot <- ggplot(hq, aes(x=SampleType, y=gc_avg, fill=SampleType)) +
  geom_boxplot(lwd=.3,outlier.shape = NA, alpha = 0.9) +
  scale_fill_manual(values=c("#CC0B0B","#0B7ACC")) +
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.margin = unit(c(2,0,2,4), "pt"))+
  scale_y_continuous(labels = function(x) paste0(x, "%"), breaks=pretty_breaks(n=6))+
  ylab("GC Content")

points <- points + theme(legend.position = "none")

stats <- (points / histo_size / plot_layout(heights = c(1, 0.8))) |
  (contigsplot / n50 / gcplot /plot_spacer() / plot_layout(heights = c(1, 1, 1, 0.6))) |
  plot_layout(widths = c(3,1))

stats

stats_leg <- ggdraw() +
    draw_plot(stats, x = 0, y = 0, width = 1, height = 1) +
    draw_plot(leg_points, x = .6, y = .12, width = .45, height = .15)+
    theme(panel.background = element_rect(fill = "white",colour = "white"))  

stats_leg

ggsave("MAGs_data_input/Figures/Figure2/SciData_fig2.png", 
       stats_leg, width = 1000, height = 900, dpi = 300, units = "px")