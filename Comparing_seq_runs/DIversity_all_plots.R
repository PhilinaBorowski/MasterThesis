## Statistical tests Hill numbers - all sequences & fjords

library(stats)
library(ggplot2)
library(gridExtra)
library(readxl)
library(grid)
library(ggsignif)
library(dplyr)
library(stats)


setwd("/isibhv/projects/AG_John/Expeditions/HE627_Data/Comparison_all_seq_runs/Hill_numbers/")

data <- read_xlsx("Hill_numbers_all_fjords.xlsx", sheet = 3)
head(data)

names(data)[names(data)=="s.e....7"] <- "Shannon_se"
names(data)[names(data)=="s.e....9"] <- "Simpson_se"
names(data)[names(data)=="s.e....11"] <- "Richness_se"

str(data)
data$Station <- as.character(data$Station)
data$Sequencer <- as.factor(data$Sequencer)


## boxplots

## compare Sequence runs

p1 <- ggplot(data, aes(x = Sequencer, y = Richness)) + geom_boxplot(fill = c ("#33608C", "#DD708E", "#F5B355"), color = "black") + 
  labs(title = "Species Richness", subtitle = "Calculated with the iNEXT package", x = "Run", y = "Species Richness") + theme_light() + 
  scale_x_discrete(limits = c("MiSeq-18S", "NextSeq-18S", "NextSeq-16S"))
#ggsave("box_seq_rich.png", width = 5, height = 6)

p2 <- ggplot(data, aes(x = Sequencer, y = Shannon)) + geom_boxplot(fill = c ("#33608C", "#DD708E", "#F5B355"), color = "black") + 
  labs(title = "Shannon Diversity", subtitle = "Calculated with the iNEXT package", x = "Run", y = "Shannon Diversity") + theme_light() + 
  scale_x_discrete(limits = c("MiSeq-18S", "NextSeq-18S", "NextSeq-16S"))
#ggsave("box_seq_sha.png", width = 5, height = 6)

p3 <- ggplot(data, aes(x = Sequencer, y = Simpson)) + geom_boxplot(fill = c ("#33608C", "#DD708E", "#F5B355"), color = "black") + 
  labs(title = "Simpson Diversity", subtitle = "Calculated with the iNEXT package", x = "Run", y = "Simpson Diversity") + theme_light() + 
  scale_x_discrete(limits = c("MiSeq-18S", "NextSeq-18S", "NextSeq-16S"))
#ggsave("box_seq_sim.png", width = 5, height = 6)


a <- arrangeGrob(p1, p2, p3, ncol = 3, top = grid::textGrob("Diversity across all 3 sequencer runs", x = 0.02, 
                                                            hjust = 0, gp = gpar(fontsize = 15)))

ggsave("rich_sha_simp_all_seq.png", a, width = 10, height = 6)

## compare fjords


MiSeq <- subset(data, Sequencer == "MiSeq-18S")           # 96 samples
NextSeq16S <- subset(data, Sequencer == "NextSeq-16S")    # 36 samples
NextSeq18S <- subset(data, Sequencer == "NextSeq-18S")    # 133 samples


# Compare MiSeq Fjords

p4 <- ggplot(MiSeq, aes(x = Fjord, y = Richness, color = Fjord)) + geom_boxplot(fill = c("#e241ff", "#bb6fff", "#8a8bff", "#2fa1ff"), color = "black") +
  labs(title = "Species Richness", subtitle = "Calculated with the iNEXT package", x = "Fjord", y = "Species Richness") + theme_light() +
  scale_x_discrete(limits = c("van Mijenfjorden", "Kongsfjorden", "Wijdefjorden", "Rijpfjorden"))
#ggsave("mi_fjo_rich.png", width = 5, height = 6)

p5 <- ggplot(MiSeq, aes(x = Fjord, y = Shannon, color = Fjord)) + geom_boxplot(fill = c("#e241ff", "#bb6fff", "#8a8bff", "#2fa1ff"), color = "black") +
  labs(title = "Shannon Diversity", subtitle = "Calculated with the iNEXT package", x = "Fjord", y = "Shannon Diversity") + theme_light() +
  scale_x_discrete(limits = c("van Mijenfjorden", "Kongsfjorden", "Wijdefjorden", "Rijpfjorden"))
#ggsave("mi_fjo_sha.png", width = 5, height = 6)

p6 <- ggplot(MiSeq, aes(x = Fjord, y = Simpson, color = Fjord)) + geom_boxplot(fill = c("#e241ff", "#bb6fff", "#8a8bff", "#2fa1ff"), color = "black") +
  labs(title = "Simpson Diversity", subtitle = "Calculated with the iNEXT package", x = "Fjord", y = "Simpson Diversity") + theme_light() +
  scale_x_discrete(limits = c("van Mijenfjorden", "Kongsfjorden", "Wijdefjorden", "Rijpfjorden"))
#ggsave("mi_fjo_sim.png", width = 5, height = 6)

b <- arrangeGrob(p4, p5, p6, ncol = 3, top = grid::textGrob("Comparing the fjords - MiSeq 18S", x = 0.02, 
                                                            hjust = 0, gp = gpar(fontsize = 15)))

ggsave("rich_sha_simp_MiSeq.png", b, width = 12.5, height = 6)


# Compare NextSeq 16S

p7 <- ggplot(NextSeq16S, aes(x = Fjord, y = Richness, color = Fjord)) + geom_boxplot(fill = c("#e241ff", "#bb6fff", "#8a8bff", "#2fa1ff"), color = "black") +
  labs(title = "Species Richness", subtitle = "Calculated with the iNEXT package", x = "Fjord", y = "Species Richness") + theme_light() +
  scale_x_discrete(limits = c("van Mijenfjorden", "Kongsfjorden", "Wijdefjorden", "Rijpfjorden"))
#ggsave("n16_fjo_rich.png", width = 5, height = 6)

p8 <- ggplot(NextSeq16S, aes(x = Fjord, y = Shannon, color = Fjord)) + geom_boxplot(fill = c("#e241ff", "#bb6fff", "#8a8bff", "#2fa1ff"), color = "black") +
  labs(title = "Shannon Diversity", subtitle = "Calculated with the iNEXT package", x = "Fjord", y = "Shannon Diversity") + theme_light() +
  scale_x_discrete(limits = c("van Mijenfjorden", "Kongsfjorden", "Wijdefjorden", "Rijpfjorden"))
#ggsave("n16_fjo_sha.png", width = 5, height = 6)

p9 <- ggplot(NextSeq16S, aes(x = Fjord, y = Simpson, color = Fjord)) + geom_boxplot(fill = c("#e241ff", "#bb6fff", "#8a8bff", "#2fa1ff"), color = "black") +
  labs(title = "Simpson Diversity", subtitle = "Calculated with the iNEXT package", x = "Fjord", y = "Simpson Diversity") + theme_light() +
  scale_x_discrete(limits = c("van Mijenfjorden", "Kongsfjorden", "Wijdefjorden", "Rijpfjorden"))
#ggsave("n16_fjo_sim.png", width = 5, height = 6)

c <- arrangeGrob(p7, p8, p9, ncol = 3, top = grid::textGrob("Comparing the fjords - NextSeq 16S", x = 0.02, 
                                                            hjust = 0, gp = gpar(fontsize = 15)))

ggsave("rich_sha_simp_Next16.png", c, width = 12.5, height = 6)


# Compare NextSeq 18S

p10 <- ggplot(NextSeq18S, aes(x = Fjord, y = Richness, color = Fjord)) + geom_boxplot(fill = c("#e241ff", "#bb6fff", "#8a8bff", "#2fa1ff"), color = "black") +
  labs(title = "Species Richness", subtitle = "Calculated with the iNEXT package", x = "Fjord", y = "Species Richness") + theme_light() +
  scale_x_discrete(limits = c("van Mijenfjorden", "Kongsfjorden", "Wijdefjorden", "Rijpfjorden"))
#ggsave("n18_fjo_rich.png", width = 5, height = 6)

p11 <- ggplot(NextSeq18S, aes(x = Fjord, y = Shannon, color = Fjord)) + geom_boxplot(fill = c("#e241ff", "#bb6fff", "#8a8bff", "#2fa1ff"), color = "black") +
  labs(title = "Shannon Diversity", subtitle = "Calculated with the iNEXT package", x = "Fjord", y = "Shannon Diversity") + theme_light() +
  scale_x_discrete(limits = c("van Mijenfjorden", "Kongsfjorden", "Wijdefjorden", "Rijpfjorden"))
#ggsave("n18_fjo_sha.png", width = 5, height = 6)

p12 <- ggplot(NextSeq18S, aes(x = Fjord, y = Simpson, color = Fjord)) + geom_boxplot(fill = c("#e241ff", "#bb6fff", "#8a8bff", "#2fa1ff"), color = "black") +
  labs(title = "Simpson Diversity", subtitle = "Calculated with the iNEXT package", x = "Fjord", y = "Simpson Diversity") + theme_light() +
  scale_x_discrete(limits = c("van Mijenfjorden", "Kongsfjorden", "Wijdefjorden", "Rijpfjorden"))
#ggsave("n18_fjo_sim.png", width = 5, height = 6)

d <- arrangeGrob(p10, p11, p12, ncol = 3, top = grid::textGrob("Comparing the fjords - NextSeq 18S", x = 0.02, 
                                                               hjust = 0, gp = gpar(fontsize = 15)))

ggsave("rich_sha_simp_Next18.png", d, width = 12.5, height = 6)


## compare miseq vs nextseq

# richness of all seq next to each other
p4.2 <- ggplot(MiSeq, aes(x = Fjord, y = Richness, color = Fjord)) + geom_boxplot(fill = c("#e241ff", "#bb6fff", "#8a8bff", "#2fa1ff"), color = "black") +
  labs(title = "MiSeq 18S", subtitle = "Calculated with the iNEXT package", x = "Fjord", y = "Species Richness") + theme_light() +
  scale_x_discrete(limits = c("van Mijenfjorden", "Kongsfjorden", "Wijdefjorden", "Rijpfjorden")) + 
  scale_y_continuous(limits=c(0, 1300))

p10.2 <- ggplot(NextSeq18S, aes(x = Fjord, y = Richness, color = Fjord)) + geom_boxplot(fill = c("#e241ff", "#bb6fff", "#8a8bff", "#2fa1ff"), color = "black") +
  labs(title = "NextSeq 18S", subtitle = "Calculated with the iNEXT package", x = "Fjord", y = "Species Richness") + theme_light() +
  scale_x_discrete(limits = c("van Mijenfjorden", "Kongsfjorden", "Wijdefjorden", "Rijpfjorden")) + 
  scale_y_continuous(limits=c(0, 1300))


e <- arrangeGrob(p4.2, p10.2, ncol = 2, top = grid::textGrob("Comparing the fjords - Species Richness", x = 0.02, 
                                                             hjust = 0, gp = gpar(fontsize = 15)))

ggsave("rich_mi_n_across_fjords.png", e, width = 9, height = 6)

## shannon
p5.2 <- ggplot(MiSeq, aes(x = Fjord, y = Shannon, color = Fjord)) + geom_boxplot(fill = c("#e241ff", "#bb6fff", "#8a8bff", "#2fa1ff"), color = "black") +
  labs(title = "MiSeq 18S", subtitle = "Calculated with the iNEXT package", x = "Fjord", y = "Shannon Diversity") + theme_light() +
  scale_x_discrete(limits = c("van Mijenfjorden", "Kongsfjorden", "Wijdefjorden", "Rijpfjorden"))+ 
  scale_y_continuous(limits=c(0, 210))

p11.2 <- ggplot(NextSeq18S, aes(x = Fjord, y = Shannon, color = Fjord)) + geom_boxplot(fill = c("#e241ff", "#bb6fff", "#8a8bff", "#2fa1ff"), color = "black") +
  labs(title = "NextSeq 18S", subtitle = "Calculated with the iNEXT package", x = "Fjord", y = "Shannon Diversity") + theme_light() +
  scale_x_discrete(limits = c("van Mijenfjorden", "Kongsfjorden", "Wijdefjorden", "Rijpfjorden"))+ 
  scale_y_continuous(limits=c(0, 210))

f <- arrangeGrob(p5.2, p11.2, ncol = 2, top = grid::textGrob("Comparing the fjords - Shannon Diversity", x = 0.02, 
                                                             hjust = 0, gp = gpar(fontsize = 15)))

ggsave("sha_mi_n_across_fjords.png", f, width = 9, height = 6)

## simpson
p6.2 <- ggplot(MiSeq, aes(x = Fjord, y = Simpson, color = Fjord)) + geom_boxplot(fill = c("#e241ff", "#bb6fff", "#8a8bff", "#2fa1ff"), color = "black") +
  labs(title = "MiSeq 18S", subtitle = "Calculated with the iNEXT package", x = "Fjord", y = "Simpson Diversity") + theme_light() +
  scale_x_discrete(limits = c("van Mijenfjorden", "Kongsfjorden", "Wijdefjorden", "Rijpfjorden"))+ 
  scale_y_continuous(limits=c(0, 75))

p12.2 <- ggplot(NextSeq18S, aes(x = Fjord, y = Simpson, color = Fjord)) + geom_boxplot(fill = c("#e241ff", "#bb6fff", "#8a8bff", "#2fa1ff"), color = "black") +
  labs(title = "NextSeq 18S", subtitle = "Calculated with the iNEXT package", x = "Fjord", y = "Simpson Diversity") + theme_light() +
  scale_x_discrete(limits = c("van Mijenfjorden", "Kongsfjorden", "Wijdefjorden", "Rijpfjorden"))+ 
  scale_y_continuous(limits=c(0, 75))

f.2 <- arrangeGrob(p6.2, p12.2, ncol = 2, top = grid::textGrob("Comparing the fjords - Simpson Diversity", x = 0.02, 
                                                               hjust = 0, gp = gpar(fontsize = 15)))

ggsave("sim_mi_n_across_fjords.png", f.2, width = 9, height = 6)



###
## compare inner fjords to outer fjords

vMijen <- subset(data, Fjord == "van Mijenfjorden")           # 92 samples
Kongs <- subset(data, Fjord == "Kongsfjorden")                # 60 samples
Wijde <- subset(data, Fjord == "Wijdefjorden")                # 70 samples
Rijp <- subset(data, Fjord == "Rijpfjorden")                  # 43 samples


vM_Mi <- subset(vMijen, Sequencer == "MiSeq-18S")
vM_N16 <- subset(vMijen, Sequencer == "NextSeq-16S")
vM_N18 <- subset(vMijen, Sequencer == "NextSeq-18S")


## richness
pvm <- ggplot(vM_Mi, aes(x = Station, y = Richness)) + geom_boxplot(fill = c("#ED90A4", "#E59884", "#c0ab52",  "#D6a166"), color = "black") + 
  scale_x_discrete(limits = c("1", "2", "4", "3")) + labs(title = "MiSeq 18S") + theme_light() + 
  scale_y_continuous(limits=c(0, 1200)) + theme(plot.title = element_text(size = 10)) 

rvMn16 <- ggplot(vM_N16, aes(x = Station, y = Richness)) + geom_boxplot(fill = c("#ED90A4", "#E59884", "#c0ab52", "#D6a166"), color = "black") + 
  #geom_signif(comparisons = list(c("1", "2"), c("3", "4"), c("1", "3"), c("1", "4"), c("2", "3"), c("2", "4")), map_signif_level = T, step_increase = 0.2, tip_length = 0.08) + 
  scale_x_discrete(limits = c("1", "2", "4", "3")) + labs(title = "NextSeq 16S") + theme_light() + 
  scale_y_continuous(limits=c(0, 1200)) + theme(plot.title = element_text(size = 10, hjust = 0.5), axis.title.x = element_blank(), axis.title.y = element_blank()) 

rvMn18 <- ggplot(vM_N18, aes(x = Station, y = Richness)) + geom_boxplot(fill = c("#ED90A4", "#E59884", "#c0ab52", "#D6a166"), color = "black") + 
  #geom_signif(comparisons = list(c("1", "2"), c("3", "4"), c("1", "3"), c("1", "4"), c("2", "3"), c("2", "4")), map_signif_level = T, step_increase = 0.2, tip_length = 0.05) + 
  scale_x_discrete(limits = c("1", "2", "4", "3")) + labs(title = "NextSeq 18S") + theme_light() + 
  scale_y_continuous(limits=c(0, 1200)) + theme(plot.title = element_text(size = 10, hjust = 0.5), axis.title.x = element_blank()) 
ggsave("plot_test.png", rvMn18)

g <- arrangeGrob(pvm, pvn18, pvn16, ncol = 3, top = grid::textGrob("Species Richness - van Mijenfjorden", x = 0.02, 
                                                                   hjust = 0, gp = gpar(fontsize = 13)))

ggsave("vMijen_rich_stations_2.png", g, width = 8, height = 5)


## shannon
shvm <- ggplot(vM_Mi, aes(x = Station, y = Shannon)) + geom_boxplot(fill = c("#ED90A4", "#E59884", "#c0ab52", "#D6a166"), color = "black") + 
  scale_x_discrete(limits = c("1", "2", "4", "3")) + labs(title = "MiSeq 18S") + theme_light() + 
  scale_y_continuous(limits=c(0, 150)) + theme(plot.title = element_text(size = 10)) 

shvMn16 <- ggplot(vM_N16, aes(x = Station, y = Shannon)) + geom_boxplot(fill = c("#ED90A4", "#E59884", "#c0ab52", "#D6a166"), color = "black") + 
  #geom_signif(comparisons = list(c("1", "2"), c("3", "4"), c("1", "3"), c("1", "4"), c("2", "3"), c("2", "4")), map_signif_level = T, step_increase = 0.2, tip_length = 0.08) + 
  scale_x_discrete(limits = c("1", "2", "4", "3")) + #labs(title = "NextSeq 16S - Shannon Diversity") + theme_light() + 
  scale_y_continuous(limits=c(0, 150)) + #theme(plot.title = element_text(size = 10)) 
  theme_light() + theme(axis.title.x = element_blank(), axis.title.y = element_blank())

shvMn18 <- ggplot(vM_N18, aes(x = Station, y = Shannon)) + geom_boxplot(fill = c("#ED90A4", "#E59884", "#c0ab52", "#D6a166"), color = "black") + 
  #geom_signif(comparisons = list(c("1", "2"), c("3", "4"), c("1", "3"), c("1", "4"), c("2", "3"), c("2", "4")), map_signif_level = T, step_increase = 0.2, tip_length = 0.08) + 
  scale_x_discrete(limits = c("1", "2", "4", "3")) + #labs(title = "NextSeq 18S - Shannon Diversity") + theme_light() + 
  scale_y_continuous(limits=c(0, 150)) + #theme(plot.title = element_text(size = 10)) 
  theme_light() + theme(axis.title.x = element_blank())

h <- arrangeGrob(shvm, shvn18, shvn16, ncol = 3, top = grid::textGrob("Shannon Diversity - van Mijenfjorden", x = 0.02, 
                                                                      hjust = 0, gp = gpar(fontsize = 13)))
#top = grid::textGrob("Title", x = 0, hjust = 0))
#?textGrob
#?gpar
ggsave("vMijen_shannon_stations.png", h, width = 8, height = 5)


## simpson 
sivm <- ggplot(vM_Mi, aes(x = Station, y = Simpson)) + geom_boxplot(fill = c("#ED90A4", "#E59884", "#c0ab52", "#D6a166"), color = "black") + 
  scale_x_discrete(limits = c("1", "2", "4", "3")) + labs(title = "MiSeq 18S") + theme_light() + 
  scale_y_continuous(limits=c(0, 60)) + theme(plot.title = element_text(size = 10)) 

sivn16 <- ggplot(vM_N16, aes(x = Station, y = Simpson)) + geom_boxplot(fill = c("#ED90A4", "#E59884", "#c0ab52", "#D6a166"), color = "black") + 
  #geom_signif(comparisons = list(c("1", "2"), c("3", "4"), c("1", "3"), c("1", "4"), c("2", "3"), c("2", "4")), map_signif_level = T, step_increase = 0.2, tip_length = 0.08) + 
  scale_x_discrete(limits = c("1", "2", "4", "3")) + #labs(title = "NextSeq 16S - Simpson Diversity") + theme_light() + 
  scale_y_continuous(limits=c(0, 60)) + #theme(plot.title = element_text(size = 10)) 
  theme_light() + theme(axis.title.y = element_blank())

sivn18 <- ggplot(vM_N18, aes(x = Station, y = Simpson)) + geom_boxplot(fill = c("#ED90A4", "#E59884", "#c0ab52", "#D6a166"), color = "black") + 
  #geom_signif(comparisons = list(c("1", "2"), c("3", "4"), c("1", "3"), c("1", "4"), c("2", "3"), c("2", "4")), map_signif_level = T, step_increase = 0.2, tip_length = 0.08) + 
  scale_x_discrete(limits = c("1", "2", "4", "3")) + #labs(title = "NextSeq 18S - Simpson Diversity") + theme_light() + 
  scale_y_continuous(limits=c(0, 60)) + #theme(plot.title = element_text(size = 10)) 
  theme_light() 

i <- arrangeGrob(sivm, sivn18, sivn16, ncol = 3, top = grid::textGrob("Simpson Diversity - van Mijenfjorden", x = 0.02, 
                                                                      hjust = 0, gp = gpar(fontsize = 13)))
#top = grid::textGrob("Title", x = 0, hjust = 0))
#?textGrob
#?gpar
ggsave("vMijen_simp_stations.png", i, width = 8, height = 5)


## Kongsfjorden 

Ko_Mi <- subset(Kongs, Sequencer == "MiSeq-18S")
Ko_N16 <- subset(Kongs, Sequencer == "NextSeq-16S")
Ko_N18 <- subset(Kongs, Sequencer == "NextSeq-18S")


## richness
kmi <- ggplot(Ko_Mi, aes(x = Station, y = Richness)) + geom_boxplot(fill = c("#a1a5ec", "#c699e7", "#6fb1e7"), color = "black") + 
  scale_x_discrete(limits = c("12", "13", "11")) + labs(title = "MiSeq 18S") + theme_light() + 
  scale_y_continuous(limits=c(0, 1650)) + theme(plot.title = element_text(size = 10)) 

rkn18 <- ggplot(Ko_N18, aes(x = Station, y = Richness)) + geom_boxplot(fill = c("#a1a5ec", "#c699e7", "#6fb1e7"), color = "black") + 
  #geom_signif(comparisons = list(c("11", "12"), c("11", "13"), c("12", "13")), map_signif_level = T, tip_length = 0.08) + 
  scale_x_discrete(limits = c("12", "13", "11")) + labs(title = "NextSeq 18S") + theme_light() + 
  scale_y_continuous(limits=c(0, 1650)) + theme(plot.title = element_text(size = 10, hjust = 0.5), axis.title.x = element_blank()) 

rkn16 <- ggplot(Ko_N16, aes(x = Station, y = Richness)) + geom_boxplot(fill = c("#a1a5ec", "#c699e7", "#6fb1e7"), color = "black") + 
  #geom_signif(comparisons = list(c("11", "12"), c("11", "13"), c("12", "13")), map_signif_level = T, tip_length = 0.08) + 
  scale_x_discrete(limits = c("12", "13", "11")) + labs(title = "NextSeq 16S") + theme_light() + 
  scale_y_continuous(limits=c(0, 1650)) + theme(plot.title = element_text(size = 10, hjust = 0.5), axis.title.x = element_blank(), axis.title.y = element_blank()) 

j <- arrangeGrob(kmi, kn18, kn16, ncol = 3, top = grid::textGrob("Species richness - Kongsfjorden", x = 0.02, 
                                                                 hjust = 0, gp = gpar(fontsize = 13)))

ggsave("kongs_rich_stat.png", j, width = 8, height = 5)


## shannon 

skmi <- ggplot(Ko_Mi, aes(x = Station, y = Shannon)) + geom_boxplot(fill = c("#a1a5ec", "#c699e7", "#6fb1e7"), color = "black") + 
  scale_x_discrete(limits = c("12", "13", "11")) + labs(title = "MiSeq 18S") + theme_light() + 
  scale_y_continuous(limits=c(0, 220)) + theme(plot.title = element_text(size = 10)) 

skn18 <- ggplot(Ko_N18, aes(x = Station, y = Shannon)) + geom_boxplot(fill = c("#a1a5ec", "#c699e7", "#6fb1e7"), color = "black") + 
  #geom_signif(comparisons = list(c("11", "12"), c("11", "13"), c("12", "13")), map_signif_level = T, tip_length = 0.08) + 
  scale_x_discrete(limits = c("12", "13", "11")) + #labs(title = "NextSeq 18S - Shannon Diversity") + theme_light() + 
  scale_y_continuous(limits=c(0, 220)) + #theme(plot.title = element_text(size = 10)) 
  theme_light() + theme(axis.title.x = element_blank())

skn16 <- ggplot(Ko_N16, aes(x = Station, y = Shannon)) + geom_boxplot(fill = c("#a1a5ec", "#c699e7", "#6fb1e7"), color = "black") + 
  #geom_signif(comparisons = list(c("11", "12"), c("11", "13"), c("12", "13")), map_signif_level = T, tip_length = 0.08) + 
  scale_x_discrete(limits = c("12", "13", "11")) + #labs(title = "NextSeq 16S - Shannon Diversity") + theme_light() + 
  scale_y_continuous(limits=c(0, 220)) + #theme(plot.title = element_text(size = 10)) 
  theme_light() + theme(axis.title.x = element_blank(), axis.title.y = element_blank())

k <- arrangeGrob(skmi, skn18, skn16, ncol = 3, top = grid::textGrob("Shannon Diversity - Kongsfjorden", x = 0.02, 
                                                                    hjust = 0, gp = gpar(fontsize = 13)))

ggsave("kongs_shan_stat.png", k, width = 8, height = 5)


## simpson

sikmi <- ggplot(Ko_Mi, aes(x = Station, y = Simpson)) + geom_boxplot(fill = c("#a1a5ec", "#c699e7", "#6fb1e7"), color = "black") + 
  scale_x_discrete(limits = c("12", "13", "11")) + labs(title = "MiSeq 18S") + theme_light() + 
  scale_y_continuous(limits=c(0, 75)) + theme(plot.title = element_text(size = 10)) 

sikn18 <- ggplot(Ko_N18, aes(x = Station, y = Simpson)) + geom_boxplot(fill = c("#a1a5ec", "#c699e7", "#6fb1e7"), color = "black") + 
  #geom_signif(comparisons = list(c("11", "12"), c("11", "13"), c("12", "13")), map_signif_level = T, tip_length = 0.08) + 
  scale_x_discrete(limits = c("12", "13", "11")) + #labs(title = "NextSeq 18S - Simpson diversity") + theme_light() + 
  scale_y_continuous(limits=c(0, 75)) + #theme(plot.title = element_text(size = 10)) 
  theme_light()

sikn16 <- ggplot(Ko_N16, aes(x = Station, y = Simpson)) + geom_boxplot(fill = c("#a1a5ec", "#c699e7", "#6fb1e7"), color = "black") + 
  #geom_signif(comparisons = list(c("11", "12"), c("11", "13"), c("12", "13")), map_signif_level = T, tip_length = 0.08) + 
  scale_x_discrete(limits = c("12", "13", "11")) + #labs(title = "NextSeq 16S - Simpson diversity") + theme_light() + 
  scale_y_continuous(limits=c(0, 75)) + #theme(plot.title = element_text(size = 10)) 
  theme_light()

l <- arrangeGrob(sikmi, sikn18, sikn16, ncol = 3, top = grid::textGrob("Simpson Diversity - Kongsfjorden", x = 0.02, 
                                                                       hjust = 0, gp = gpar(fontsize = 13)))

ggsave("kongs_simp_stat.png", l, width = 8, height = 5)


## Wijdefjorden

Wi_Mi <- subset(Wijde, Sequencer == "MiSeq-18S")
Wi_N16 <- subset(Wijde, Sequencer == "NextSeq-16S")
Wi_N18 <- subset(Wijde, Sequencer == "NextSeq-18S")

## richness

wmi <- ggplot(Wi_Mi, aes(x = Station, y = Richness)) + geom_boxplot(fill = c("#a3b353", "#4fbf85", "#7eba68"), color = "black") + 
  scale_x_discrete(limits = c("5", "7", "6")) + labs(title = "MiSeq 18S") + theme_light() + 
  scale_y_continuous(limits=c(0, 1100)) + 
  theme(plot.title = element_text(size = 10)) 

wn18 <- ggplot(Wi_N18, aes(x = Station, y = Richness)) + geom_boxplot(fill = c("#a3b353", "#4fbf85", "#7eba68"), color = "black") + 
  #geom_signif(comparisons = list(c("5", "6"), c("5", "7"), c("6", "7"), map_signif_level = T, tip_length = 0.08) + 
  scale_x_discrete(limits = c("5", "7", "6")) + labs(title = "NextSeq 18S") + theme_light() + 
  scale_y_continuous(limits=c(0, 1100)) + 
  theme(plot.title = element_text(size = 10, hjust = 0.5), axis.title.x = element_blank())

wn16 <- ggplot(Wi_N16, aes(x = Station, y = Richness)) + geom_boxplot(fill = c("#a3b353", "#4fbf85", "#7eba68"), color = "black") + 
  #geom_signif(comparisons = list(c("5", "6"), c("5", "7"), c("6", "7"), map_signif_level = T, tip_length = 0.08) + 
  scale_x_discrete(limits = c("5", "7", "6")) + labs(title = "NextSeq 16S") + theme_light() + 
  scale_y_continuous(limits=c(0, 1100)) + 
  theme(plot.title = element_text(size = 10, hjust = 0.5), axis.title.x = element_blank(), axis.title.y = element_blank()) 

m <- arrangeGrob(wmi, wn18, wn16, ncol = 3, top = grid::textGrob("Species richness - Wijdefjorden", x = 0.02, 
                                                                 hjust = 0, gp = gpar(fontsize = 13)))

ggsave("wijde_rich_stat.png", m, width = 8, height = 5)


## shannon 

swmi <- ggplot(Wi_Mi, aes(x = Station, y = Shannon)) + geom_boxplot(fill = c("#a3b353", "#4fbf85", "#7eba68"), color = "black") + 
  scale_x_discrete(limits = c("5", "7", "6")) + labs(title = "MiSeq 18S") + theme_light() + 
  scale_y_continuous(limits=c(0, 160)) + 
  theme(plot.title = element_text(size = 10)) 

swn18 <- ggplot(Wi_N18, aes(x = Station, y = Shannon)) + geom_boxplot(fill = c("#a3b353", "#4fbf85", "#7eba68"), color = "black") + 
  #geom_signif(comparisons = list(c("5", "6"), c("5", "7"), c("6", "7")), map_signif_level = T, tip_length = 0) + 
  scale_x_discrete(limits = c("5", "7", "6")) + #labs(title = "NextSeq 18S") + theme_light() + 
  scale_y_continuous(limits=c(0, 170)) + theme_light() +
  #theme(plot.title = element_text(size = 10)) 
  theme(axis.title.x = element_blank())

swn16 <- ggplot(Wi_N16, aes(x = Station, y = Shannon)) + geom_boxplot(fill = c("#a3b353", "#4fbf85", "#7eba68"), color = "black") + 
  #geom_signif(comparisons = list(c("5", "6"), c("5", "7"), c("6", "7")), map_signif_level = T, tip_length = 0) + 
  scale_x_discrete(limits = c("5", "7", "6")) + #labs(title = "NextSeq 16S") + theme_light() + 
  scale_y_continuous(limits=c(0, 160)) + theme_light() + 
  #theme(plot.title = element_text(size = 10)) 
  theme(axis.title.x = element_blank(), axis.title.y = element_blank())

n <- arrangeGrob(swmi, swn18, swn16, ncol = 3, top = grid::textGrob("Shannon Diversity - Wijdefjorden", x = 0.02, 
                                                                    hjust = 0, gp = gpar(fontsize = 13)))

ggsave("wijde_sha_stat.png", n, width = 8, height = 5)


## simpson 

siwmi <- ggplot(Wi_Mi, aes(x = Station, y = Simpson)) + geom_boxplot(fill = c("#a3b353", "#4fbf85", "#7eba68"), color = "black") + 
  scale_x_discrete(limits = c("5", "7", "6")) + labs(title = "MiSeq 18S") + theme_light() + 
  scale_y_continuous(limits=c(0, 70)) + 
  theme(plot.title = element_text(size = 10)) 

siwn18 <- ggplot(Wi_N18, aes(x = Station, y = Simpson)) + geom_boxplot(fill = c("#a3b353", "#4fbf85", "#7eba68"), color = "black") + 
  #geom_signif(comparisons = list(c("5", "6"), c("5", "7"), c("6", "7")), map_signif_level = T, tip_length = 0) + 
  geom_signif(comparisons = list(c("5", "7")), map_signif_level = F, tip_length = 0, textsize = 8, test = t.test) + 
  scale_x_discrete(limits = c("5", "7", "6")) + #labs(title = "NextSeq 18S") + theme_light() + 
  scale_y_continuous(limits=c(0, 70)) + theme_light()
  #theme(plot.title = element_text(size = 10)) 
  

# ------

group_by(Wi_N18, Station == "5", Station == "7") %>%
  summarise(
    count = n(),
    mean = mean(Simpson, na.rm = TRUE),
    sd = sd(Simpson, na.rm = TRUE)
  )

# Compute the analysis of variance
Wres.aov <- aov(Simpson ~ Station, data = Wi_N18)
# Summary of the analysis
summary(Wres.aov)

TukeyHSD(Wres.aov)


# subset, so only the 2 statons
st57 <- subset(Wi_N18, Station == "5" | Station == "7")

# Compute the analysis of variance
res.aov <- aov(Simpson ~ Station, data = st57)
# Summary of the analysis
summary(res.aov)

TukeyHSD(res.aov)

# or 

st57 %>%
  group_by(Station) %>%
  summarise(mean = mean(Simpson),
            sd = sd(Simpson))

#fit the one-way ANOVA model
model <- aov(Simpson ~ Station, data = st57)

#view the model output
summary(model)

res <- t.test(Simpson ~ Station, data = st57)
res

## ------

siwn16 <- ggplot(Wi_N16, aes(x = Station, y = Simpson)) + geom_boxplot(fill = c("#a3b353", "#4fbf85", "#7eba68"), color = "black") + 
  #geom_signif(comparisons = list(c("5", "6"), c("5", "7"), c("6", "7")), map_signif_level = T, tip_length = 0) + 
  scale_x_discrete(limits = c("5", "7", "6")) + #labs(title = "NextSeq 16S") + theme_light() + 
  scale_y_continuous(limits=c(0, 70)) + theme_light() 
  #theme(plot.title = element_text(size = 10)) 

o <- arrangeGrob(siwmi, siwn18, siwn16, ncol = 3, top = grid::textGrob("Simpson Diversity - Wijdefjorden", x = 0.02, 
                                                                       hjust = 0, gp = gpar(fontsize = 13)))

ggsave("wijde_sim_stat.png", o, width = 8, height = 5)


## rijpfjorden 

Ri_Mi <- subset(Rijp, Sequencer == "MiSeq-18S")
Ri_N16 <- subset(Rijp, Sequencer == "NextSeq-16S")
Ri_N18 <- subset(Rijp, Sequencer == "NextSeq-18S")

## richness

rmi <- ggplot(Ri_Mi, aes(x = Station, y = Richness)) + geom_boxplot(fill = c("#28bbd7", "#00c1a3"), color = "black") + 
  scale_x_discrete(limits = c("10", "8")) + labs(title = "MiSeq 18S") + theme_light() + 
  scale_y_continuous(limits=c(0, 840)) + 
  theme(plot.title = element_text(size = 10)) 

rrn18 <- ggplot(Ri_N18, aes(x = Station, y = Richness)) + geom_boxplot(fill = c("#28bbd7", "#00c1a3"), color = "black") + 
  #geom_signif(comparisons = list(c("8", "10")), map_signif_level = T, tip_length = 0) + 
  scale_x_discrete(limits = c("10", "8")) + labs(title = "NextSeq 18S") + theme_light() + 
  scale_y_continuous(limits=c(0, 900)) + 
  theme(plot.title = element_text(size = 10, hjust = 0.5), axis.title.x = element_blank()) 

rrn16 <- ggplot(Ri_N16, aes(x = Station, y = Richness)) + geom_boxplot(fill = c("#28bbd7", "#00c1a3"), color = "black") + 
  #geom_signif(comparisons = list(c("8", "10")), map_signif_level = T, tip_length = 0) + 
  scale_x_discrete(limits = c("10", "8")) + labs(title = "NextSeq 16S") + theme_light() + 
  scale_y_continuous(limits=c(0, 840)) +
  theme(plot.title = element_text(size = 10), axis.title.x = element_blank(), axis.title.y = element_blank())

p <- arrangeGrob(rmi, rn18, rn16, ncol = 3, top = grid::textGrob("Species richness - Rijpfjorden", x = 0.02, 
                                                                 hjust = 0, gp = gpar(fontsize = 13)))

ggsave("rijp_rich_stat.png", p, width = 8, height = 5)
#Removed 1 row containing missing values or values outside the scale range (`stat_boxplot()`). 
# no data for station 9


## shannon
srmi <- ggplot(Ri_Mi, aes(x = Station, y = Shannon)) + geom_boxplot(fill = c("#28bbd7", "#00c1a3"), color = "black") + 
  scale_x_discrete(limits = c("10", "8")) + labs(title = "MiSeq 18S") + theme_light() + 
  scale_y_continuous(limits=c(0, 170)) + 
  theme(plot.title = element_text(size = 10)) 

srn18 <- ggplot(Ri_N18, aes(x = Station, y = Shannon)) + geom_boxplot(fill = c("#28bbd7", "#00c1a3"), color = "black") + 
  #geom_signif(comparisons = list(c("8", "10")), map_signif_level = T, tip_length = 0) + 
  scale_x_discrete(limits = c("10", "8")) + #labs(title = "NextSeq 18S") + theme_light() + 
  scale_y_continuous(limits=c(0, 170)) + theme_light() + 
  #theme(plot.title = element_text(size = 10)) 
  theme(axis.title.x = element_blank())

srn16 <- ggplot(Ri_N16, aes(x = Station, y = Shannon)) + geom_boxplot(fill = c("#28bbd7", "#00c1a3"), color = "black") + 
  #geom_signif(comparisons = list(c("8", "10")), map_signif_level = T, tip_length = 0) + 
  scale_x_discrete(limits = c("10", "8")) + #labs(title = "NextSeq 16S") + theme_light() + 
  scale_y_continuous(limits=c(0, 170)) + theme_light() +
  #theme(plot.title = element_text(size = 10) 
  theme(axis.title.x = element_blank(), axis.title.y = element_blank()) 
  
q <- arrangeGrob(srmi, srn18, srn16, ncol = 3, top = grid::textGrob("Shannon Diversity - Rijpfjorden", x = 0.02, 
                                                                    hjust = 0, gp = gpar(fontsize = 13)))

ggsave("rijp_sha_stat.png", q, width = 8, height = 5)
#Removed 1 row containing missing values or values outside the scale range (`stat_boxplot()`). 
# no data for station 9


## simpson
sirmi <- ggplot(Ri_Mi, aes(x = Station, y = Simpson)) + geom_boxplot(fill = c("#28bbd7", "#00c1a3"), color = "black") + 
  scale_x_discrete(limits = c("10", "8")) + labs(title = "MiSeq 18S") + theme_light() + 
  scale_y_continuous(limits=c(0, 75)) + 
  theme(plot.title = element_text(size = 10)) 

sirn18 <- ggplot(Ri_N18, aes(x = Station, y = Simpson)) + geom_boxplot(fill = c("#28bbd7", "#00c1a3"), color = "black") + 
  #geom_signif(comparisons = list(c("8", "10")), map_signif_level = T, tip_length = 0) + 
  scale_x_discrete(limits = c("10", "8")) + #labs(title = "NextSeq 18S") + theme_light() + 
  scale_y_continuous(limits=c(0, 75)) + theme_light() 
  #theme(plot.title = element_text(size = 10)) 

sirn16 <- ggplot(Ri_N16, aes(x = Station, y = Simpson)) + geom_boxplot(fill = c("#28bbd7", "#00c1a3"), color = "black") + 
  #geom_signif(comparisons = list(c("8", "10")), map_signif_level = T, tip_length = 0) + 
  scale_x_discrete(limits = c("10", "8")) + #labs(title = "NextSeq 16S") + theme_light() + 
  scale_y_continuous(limits=c(0, 75)) + theme_light() + 
  #theme(plot.title = element_text(size = 10)) 
  theme(axis.title.y = element_blank())

q <- arrangeGrob(sirmi, sirn18, sirn16, ncol = 3, top = grid::textGrob("Simpson Diversity - Rijpfjorden", x = 0.02, 
                                                                       hjust = 0, gp = gpar(fontsize = 13)))

ggsave("rijp_sim_stat.png", q, width = 8, height = 5)
#Removed 1 row containing missing values or values outside the scale range (`stat_boxplot()`). 
# no data for station 9





## compare NextSeq 18S FO2 filters with NextSeq 16S

# create subset of NextSeq 18S with only F02 stations

# compare richness

n18_f02 <- subset(NextSeq18S, Sample_method == "F02")

nf02r <- ggplot(n18_f02, aes(x = Fjord, y = Richness)) + geom_boxplot(fill = c("#e241ff", "#bb6fff", "#8a8bff", "#2fa1ff"), color = "black") +
  labs(title = "NextSeq 18S (F02 filters)", subtitle = "Calculated with the iNEXT package", x = "Fjord", y = "Species Richness") + theme_light() +
  scale_x_discrete(limits = c("van Mijenfjorden", "Kongsfjorden", "Wijdefjorden", "Rijpfjorden")) + 
  scale_y_continuous(limits=c(0, 1260))

n16r <- ggplot(NextSeq16S, aes(x = Fjord, y = Richness)) + geom_boxplot(fill = c("#e241ff", "#bb6fff", "#8a8bff", "#2fa1ff"), color = "black") +
  labs(title = "NextSeq 16S", subtitle = "Calculated with the iNEXT package", x = "Fjord", y = "Species Richness") + theme_light() +
  scale_x_discrete(limits = c("van Mijenfjorden", "Kongsfjorden", "Wijdefjorden", "Rijpfjorden")) + 
  scale_y_continuous(limits=c(0, 1260))


r <- arrangeGrob(nf02r, n16r, ncol = 2, top = grid::textGrob("Comparing picoplankton - Species Richness", x = 0.02, 
                                                             hjust = 0, gp = gpar(fontsize = 15)))

ggsave("pico_rich.png", r, width = 9, height = 6)  
  

# shannon 

nf02s <- ggplot(n18_f02, aes(x = Fjord, y = Shannon)) + geom_boxplot(fill = c("#e241ff", "#bb6fff", "#8a8bff", "#2fa1ff"), color = "black") +
  labs(title = "NextSeq 18S (F02 filters)", subtitle = "Calculated with the iNEXT package", x = "Fjord", y = "Species Richness") + theme_light() +
  scale_x_discrete(limits = c("van Mijenfjorden", "Kongsfjorden", "Wijdefjorden", "Rijpfjorden")) + 
  scale_y_continuous(limits=c(0, 200))

n16s <- ggplot(NextSeq16S, aes(x = Fjord, y = Shannon)) + geom_boxplot(fill = c("#e241ff", "#bb6fff", "#8a8bff", "#2fa1ff"), color = "black") +
  labs(title = "NextSeq 16S", subtitle = "Calculated with the iNEXT package", x = "Fjord", y = "Species Richness") + theme_light() +
  scale_x_discrete(limits = c("van Mijenfjorden", "Kongsfjorden", "Wijdefjorden", "Rijpfjorden")) + 
  scale_y_continuous(limits=c(0, 200))


s <- arrangeGrob(nf02s, n16s, ncol = 2, top = grid::textGrob("Comparing picoplankton - Shannon diversity", x = 0.02, 
                                                             hjust = 0, gp = gpar(fontsize = 15)))

ggsave("pico_sha.png", s, width = 9, height = 6)  



# simpson 

nf02si <- ggplot(n18_f02, aes(x = Fjord, y = Simpson)) + geom_boxplot(fill = c("#e241ff", "#bb6fff", "#8a8bff", "#2fa1ff"), color = "black") +
  labs(title = "NextSeq 18S (F02 filters)", subtitle = "Calculated with the iNEXT package", x = "Fjord", y = "Species Richness") + theme_light() +
  scale_x_discrete(limits = c("van Mijenfjorden", "Kongsfjorden", "Wijdefjorden", "Rijpfjorden")) + 
  scale_y_continuous(limits=c(0, 75))

n16si <- ggplot(NextSeq16S, aes(x = Fjord, y = Simpson)) + geom_boxplot(fill = c("#e241ff", "#bb6fff", "#8a8bff", "#2fa1ff"), color = "black") +
  labs(title = "NextSeq 16S", subtitle = "Calculated with the iNEXT package", x = "Fjord", y = "Species Richness") + theme_light() +
  scale_x_discrete(limits = c("van Mijenfjorden", "Kongsfjorden", "Wijdefjorden", "Rijpfjorden")) + 
  scale_y_continuous(limits=c(0, 75))


t <- arrangeGrob(nf02si, n16si, ncol = 2, top = grid::textGrob("Comparing picoplankton - Simpson diversity", x = 0.02, 
                                                             hjust = 0, gp = gpar(fontsize = 15)))

ggsave("pico_simp.png", t, width = 9, height = 6)  






## ordering plots for thesis

#vm_ko_page <- 

vm <- grid.arrange(rvMn18, rvMn16,
          shvMn18, shvMn16,
          sivn18, sivn16,
          nrow = 3, ncol = 2,
          top = textGrob("van Mijenforden", gp = gpar (fontsize = 14)))
ggsave("van_Mi_plots.png", vm, width = 17, height = 13, units = "cm")          
          
          
ko <- grid.arrange(rkn18, rkn16,
                  skn18, skn16,
                  sikn18, sikn16, 
                  nrow = 3, ncol = 2,
                  top = textGrob("Kongsfjorden", gp = gpar (fontsize = 14)))
ggsave("Ko_plots.png", ko, width = 17, height = 13, units = "cm")          
         

wij <- grid.arrange(wn18, wn16,
                    swn18, swn16,
                    siwn18, siwn16,
                    nrow = 3, ncol = 2,
                    top = textGrob("Wijdefjorden", gp = gpar (fontsize = 14)))
ggsave("Wi_plots.png", wij , width = 17, height = 13, units = "cm")


rijp <- grid.arrange(rrn18, rrn16,
                     srn18, srn16,
                     sirn18, sirn16,
                     nrow = 3, ncol = 2,
                     top = textGrob("Rijpfjorden", gp = gpar (fontsize = 14)))
ggsave("Ri_plots.png", rijp, width = 17, height = 13, units = "cm")



### comparing mouths + inner fjords 

#mouths: 1, 12, 5, 10
#inner: 3, 11, 6, 9
#centre: 4, 13, 7, 8


sub_mouth_16S <- subset(NextSeq16S, Station == "1" | Station == "12" | Station == "5" | Station == "10")   # 12
sub_mouth_18S <- subset(NextSeq18S, Station == "1" | Station == "12" | Station == "5" | Station == "10")   # 44

inner_16S <- subset(NextSeq16S, Station == "3" | Station == "11" | Station == "6" | Station == "9")        # 9
inner_18S <- subset(NextSeq18S, Station == "3" | Station == "11" | Station == "6" | Station == "9")        # 33

centre_16S <- subset(NextSeq16S, Station == "4" | Station == "13" | Station == "7" | Station == "8")       # 12
centre_18S <- subset(NextSeq18S, Station == "4" | Station == "13" | Station == "7" | Station == "8")       # 41
  
  
## looking at richness
  
mouth16 <- ggplot(sub_mouth_16S, aes(x = Station, y = Richness)) + geom_boxplot(fill = c("#ED90A4", "#a1a5ec", "#a3b353", "#28bbd7"), color = "black") + 
  geom_signif(comparisons = list(c("1", "12"),  c("12", "5"), c("5", "10"), c("12", "10")), map_signif_level = T, step_increase = 0.1, tip_length = 0.02, test = t.test, textsize = 8) + 
  scale_x_discrete(limits = c("1", "12", "5", "10")) + labs(title = "Species Richness - Fjord Mouths - NextSeq 16S") + theme_light() + 
  scale_y_continuous(limits=c(300, 710)) + theme(plot.title = element_text(size = 16, hjust = 0.5), axis.title = element_text(size = 14))

ggsave("16S_fjord_mouths.png", mouth16, width = 6, height = 6)  


mouth18 <- ggplot(sub_mouth_18S, aes(x = Station, y = Richness)) + geom_boxplot(fill = c("#ED90A4", "#a1a5ec", "#a3b353", "#28bbd7"), color = "black") + 
  #geom_signif(comparisons = list(c("1", "12"),  c("12", "5"), c("5", "10"), c("1", "10"), c("1", "5"), c("12", "10")), map_signif_level = T, step_increase = 0.1, tip_length = 0.02, test = t.test, textsize = 8) + 
  scale_x_discrete(limits = c("1", "12", "5", "10")) + labs(title = "Species Richness - Fjord Mouths - NextSeq 18S") + theme_light() + 
  scale_y_continuous(limits=c(0, 1250)) + theme(plot.title = element_text(size = 16, hjust = 0.5), axis.title = element_text(size = 14))

ggsave("18S_fjord_mouths.png", mouth18, width = 6, height = 6)  


inner16 <- ggplot(inner_16S, aes(x = Station, y = Richness)) + geom_boxplot(fill = c("#D6a166", "#6fb1e7", "#7eba68"), color = "black") + 
  geom_signif(comparisons = list(c("11", "6")), map_signif_level = T, step_increase = 0.1, tip_length = 0.02, test = t.test, textsize = 8) + 
  scale_x_discrete(limits = c("3", "11", "6")) + labs(title = "Species Richness - Inner fjord station - NextSeq 16S") + theme_light() + 
  scale_y_continuous(limits=c(300, 610)) + theme(plot.title = element_text(size = 16, hjust = 0.5), axis.title = element_text(size = 14))

ggsave("16S_inner_stats.png", inner16, width = 6, height = 6)


inner18 <- ggplot(inner_18S, aes(x = Station, y = Richness)) + geom_boxplot(fill = c("#D6a166", "#6fb1e7", "#7eba68"), color = "black") + 
  #geom_signif(comparisons = list(c("11", "6"), c("3", "11"), c("3", "6")), map_signif_level = T, step_increase = 0.1, tip_length = 0.02, test = t.test, textsize = 8) + 
  scale_x_discrete(limits = c("3", "11", "6")) + labs(title = "Species Richness - Inner fjord station - NextSeq 18S") + theme_light() + 
  scale_y_continuous(limits=c(0, 1200)) + theme(plot.title = element_text(size = 16, hjust = 0.5), axis.title = element_text(size = 14))

ggsave("18S_inner_stats.png", inner18, width = 6, height = 6)


centre16 <- ggplot(centre_16S, aes(x = Station, y = Richness)) + geom_boxplot(fill = c("#c0ab52", "#c699e7", "#4fbf85", "#00c1a3"), color = "black") + 
  geom_signif(comparisons = list(c("13", "7"), c("13", "8")), map_signif_level = T, step_increase = 0.1, tip_length = 0.02, test = t.test, textsize = 8) + 
  scale_x_discrete(limits = c("4", "13", "7", "8")) + labs(title = "Species Richness - Central fjord station - NextSeq 16S") + theme_light() + 
  scale_y_continuous(limits=c(300, 560)) + theme(plot.title = element_text(size = 16, hjust = 0.5), axis.title = element_text(size = 14))

ggsave("16S_centre.png", centre16, width = 6, height = 6)


centre18 <- ggplot(centre_18S, aes(x = Station, y = Richness)) + geom_boxplot(fill = c("#c0ab52", "#c699e7", "#4fbf85", "#00c1a3"), color = "black") + 
  #geom_signif(comparisons = list(c("4", "13"), c("4", "7"), c("4", "8"), c("13", "7"), c("13", "8"), c("7", "8")), map_signif_level = T, step_increase = 0.1, tip_length = 0.02, test = t.test, textsize = 8) + 
  scale_x_discrete(limits = c("4", "13", "7", "8")) + labs(title = "Species Richness - Central fjord station - NextSeq 18S") + theme_light() + 
  scale_y_continuous(limits=c(0, 1300)) + theme(plot.title = element_text(size = 16, hjust = 0.5), axis.title = element_text(size = 14))

ggsave("18S_centre.png", centre18, width = 6, height = 6)


# shannon 

sha_mouth16 <- ggplot(sub_mouth_16S, aes(x = Station, y = Shannon)) + geom_boxplot(fill = c("#ED90A4", "#a1a5ec", "#a3b353", "#28bbd7"), color = "black") + 
  #geom_signif(comparisons = list(c("1", "12"),  c("1", "5"), c("1", "10"), c("12", "5"), c("5", "10"), c("12", "10")), map_signif_level = T, step_increase = 0.1, tip_length = 0.02, test = t.test, textsize = 8) + 
  scale_x_discrete(limits = c("1", "12", "5", "10")) + labs(title = "Shannon Diversity - Fjord Mouths - NextSeq 16S") + theme_light() + 
  scale_y_continuous(limits=c(40, 100)) + theme(plot.title = element_text(size = 16, hjust = 0.5), axis.title = element_text(size = 14))

ggsave("16S_fjord_sha_mouths.png", sha_mouth16, width = 6, height = 6)  


sha_mouth18 <- ggplot(sub_mouth_18S, aes(x = Station, y = Shannon)) + geom_boxplot(fill = c("#ED90A4", "#a1a5ec", "#a3b353", "#28bbd7"), color = "black") + 
  #geom_signif(comparisons = list(c("1", "12"),  c("12", "5"), c("5", "10"), c("1", "10"), c("1", "5"), c("12", "10")), map_signif_level = T, step_increase = 0.1, tip_length = 0.02, test = t.test, textsize = 8) + 
  scale_x_discrete(limits = c("1", "12", "5", "10")) + labs(title = "Shannon Diversity - Fjord Mouths - NextSeq 18S") + theme_light() + 
  scale_y_continuous(limits=c(0, 200)) + theme(plot.title = element_text(size = 16, hjust = 0.5), axis.title = element_text(size = 14))

ggsave("18S_fjord_sha_mouths.png", sha_mouth18, width = 6, height = 6)  


sha_inner16 <- ggplot(inner_16S, aes(x = Station, y = Shannon)) + geom_boxplot(fill = c("#D6a166", "#6fb1e7", "#7eba68"), color = "black") + 
  #geom_signif(comparisons = list(c("3", "11"), c("3", "6"), c("11", "6")), map_signif_level = T, step_increase = 0.1, tip_length = 0.02, test = t.test, textsize = 8) + 
  scale_x_discrete(limits = c("3", "11", "6")) + labs(title = "Shannon Diversity - Inner fjord station - NextSeq 16S") + theme_light() + 
  scale_y_continuous(limits=c(60, 110)) + theme(plot.title = element_text(size = 16, hjust = 0.5), axis.title = element_text(size = 14))

ggsave("16S_sha_inner_stats.png", sha_inner16, width = 6, height = 6)


sha_inner18 <- ggplot(inner_18S, aes(x = Station, y = Shannon)) + geom_boxplot(fill = c("#D6a166", "#6fb1e7", "#7eba68"), color = "black") + 
  #geom_signif(comparisons = list(c("11", "6"), c("3", "11"), c("3", "6")), map_signif_level = T, step_increase = 0.1, tip_length = 0.02, test = t.test, textsize = 8) + 
  scale_x_discrete(limits = c("3", "11", "6")) + labs(title = "Shannon Diversity - Inner fjord station - NextSeq 18S") + theme_light() + 
  scale_y_continuous(limits=c(0, 165)) + theme(plot.title = element_text(size = 16, hjust = 0.5), axis.title = element_text(size = 14))

ggsave("18S_sha_inner_stats.png", sha_inner18, width = 6, height = 6)


sha_centre16 <- ggplot(centre_16S, aes(x = Station, y = Shannon)) + geom_boxplot(fill = c("#c0ab52", "#c699e7", "#4fbf85", "#00c1a3"), color = "black") + 
  geom_signif(comparisons = list(c("13", "8"), c("7", "8")), map_signif_level = T, step_increase = 0.2, tip_length = 0.02, test = t.test, textsize = 8) + 
  scale_x_discrete(limits = c("4", "13", "7", "8")) + labs(title = "Shannon Diversity - Central fjord station - NextSeq 16S") + theme_light() + 
  scale_y_continuous(limits=c(50, 120)) + theme(plot.title = element_text(size = 16, hjust = 0.5), axis.title = element_text(size = 14))

ggsave("16S_sha_centre.png", sha_centre16, width = 6, height = 6)


sha_cnetre18 <- ggplot(centre_18S, aes(x = Station, y = Shannon)) + geom_boxplot(fill = c("#c0ab52", "#c699e7", "#4fbf85", "#00c1a3"), color = "black") + 
  #geom_signif(comparisons = list(c("4", "13"), c("4", "7"), c("4", "8"), c("13", "7"), c("13", "8"), c("7", "8")), map_signif_level = T, step_increase = 0.1, tip_length = 0.02, test = t.test, textsize = 8) + 
  scale_x_discrete(limits = c("4", "13", "7", "8")) + labs(title = "Shannon Diversity - Central fjord station - NextSeq 18S") + theme_light() + 
  scale_y_continuous(limits=c(0, 200)) + theme(plot.title = element_text(size = 16, hjust = 0.5), axis.title = element_text(size = 14))

ggsave("18S_sha_centre.png", sha_cnetre18, width = 6, height = 6)


## simpson

si_mouth16 <- ggplot(sub_mouth_16S, aes(x = Station, y = Simpson)) + geom_boxplot(fill = c("#ED90A4", "#a1a5ec", "#a3b353", "#28bbd7"), color = "black") + 
  #geom_signif(comparisons = list(c("1", "12"),  c("1", "5"), c("1", "10"), c("12", "5"), c("5", "10"), c("12", "10")), map_signif_level = T, step_increase = 0.1, tip_length = 0.02, test = t.test, textsize = 8) + 
  scale_x_discrete(limits = c("1", "12", "5", "10")) + labs(title = "Simpson Diversity - Fjord Mouths - NextSeq 16S") + theme_light() + 
  scale_y_continuous(limits=c(10, 45)) + theme(plot.title = element_text(size = 16, hjust = 0.5), axis.title = element_text(size = 14))

ggsave("16S_fjord_si_mouths.png", si_mouth16, width = 6, height = 6)  


si_mouth18 <- ggplot(sub_mouth_18S, aes(x = Station, y = Simpson)) + geom_boxplot(fill = c("#ED90A4", "#a1a5ec", "#a3b353", "#28bbd7"), color = "black") + 
  geom_signif(comparisons = list(c("1", "12"), c("10", "1")), map_signif_level = T, step_increase = 0.2, tip_length = 0.02, test = t.test, textsize = 8) + 
  scale_x_discrete(limits = c("1", "12", "5", "10")) + labs(title = "Simpson Diversity - Fjord Mouths - NextSeq 18S") + theme_light() + 
  scale_y_continuous(limits=c(0, 100)) + theme(plot.title = element_text(size = 16, hjust = 0.5), axis.title = element_text(size = 14))

ggsave("18S_fjord_si_mouths.png", si_mouth18, width = 6, height = 6)  


si_inner16 <- ggplot(inner_16S, aes(x = Station, y = Simpson)) + geom_boxplot(fill = c("#D6a166", "#6fb1e7", "#7eba68"), color = "black") + 
  #geom_signif(comparisons = list(c("3", "11"), c("3", "6"), c("11", "6")), map_signif_level = T, step_increase = 0.1, tip_length = 0.02, test = t.test, textsize = 8) + 
  scale_x_discrete(limits = c("3", "11", "6")) + labs(title = "Simpson Diversity - Inner fjord station - NextSeq 16S") + theme_light() + 
  scale_y_continuous(limits=c(0, 60)) + theme(plot.title = element_text(size = 16, hjust = 0.5), axis.title = element_text(size = 14))

ggsave("16S_si_inner_stats.png", si_inner16, width = 6, height = 6)


si_inner18 <- ggplot(inner_18S, aes(x = Station, y = Simpson)) + geom_boxplot(fill = c("#D6a166", "#6fb1e7", "#7eba68"), color = "black") + 
  #geom_signif(comparisons = list(c("11", "6"), c("3", "11"), c("3", "6")), map_signif_level = T, step_increase = 0.1, tip_length = 0.02, test = t.test, textsize = 8) + 
  scale_x_discrete(limits = c("3", "11", "6")) + labs(title = "Simpson Diversity - Inner fjord station - NextSeq 18S") + theme_light() + 
  scale_y_continuous(limits=c(0, 75)) + theme(plot.title = element_text(size = 16, hjust = 0.5), axis.title = element_text(size = 14))

ggsave("18S_si_inner_stats.png", si_inner18, width = 6, height = 6)


si_centre16 <- ggplot(centre_16S, aes(x = Station, y = Simpson)) + geom_boxplot(fill = c("#c0ab52", "#c699e7", "#4fbf85", "#00c1a3"), color = "black") + 
  geom_signif(comparisons = list(c("4", "7"), c("7", "8")), map_signif_level = T, step_increase = 0.2, tip_length = 0.02, test = t.test, textsize = 8) +  
  scale_x_discrete(limits = c("4", "13", "7", "8")) + labs(title = "Simpson Diversity - Central fjord station - NextSeq 16S") + theme_light() + 
  scale_y_continuous(limits=c(20, 50)) + theme(plot.title = element_text(size = 16, hjust = 0.5), axis.title = element_text(size = 14))

ggsave("16S_si_centre.png", si_centre16, width = 6, height = 6)


si_cnetre18 <- ggplot(centre_18S, aes(x = Station, y = Simpson)) + geom_boxplot(fill = c("#c0ab52", "#c699e7", "#4fbf85", "#00c1a3"), color = "black") + 
  #geom_signif(comparisons = list(c("4", "13"), c("4", "7"), c("4", "8"), c("13", "7"), c("13", "8"), c("7", "8")), map_signif_level = T, step_increase = 0.1, tip_length = 0.02, test = t.test, textsize = 8) + 
  scale_x_discrete(limits = c("4", "13", "7", "8")) + labs(title = "Simpson Diversity - Central fjord station - NextSeq 18S") + theme_light() + 
  scale_y_continuous(limits=c(0, 75)) + theme(plot.title = element_text(size = 16, hjust = 0.5), axis.title = element_text(size = 14))

ggsave("18S_si_centre.png", si_cnetre18, width = 6, height = 6)





# t tests

st57 <- subset(Wi_N18, Station == "5" | Station == "7")

res <- t.test(Simpson ~ Station, data = st57)
res



st112 <- subset(sub_mouth_16S, Station == "1" | Station == "12")
t112 <- t.test(Richness ~ Station, data = st112)
t112

st125 <- subset(sub_mouth_16S, Station == "12" | Station == "5")
t125 <- t.test(Richness ~ Station, data = st125)
t125

st510 <- subset(sub_mouth_16S, Station =="5" | Station == "10")
t510 <- t.test(Richness ~ Station, data = st510)
t510

st1210 <- subset(sub_mouth_16S, Station == "12" | Station == "10")
t1210 <- t.test(Richness ~ Station, data = st1210)
t1210

st137 <- subset(centre_16S, Station == "13" | Station == "7")
t137 <- t.test(Richness ~ Station, data = st137)
t137

st138 <- subset(centre_16S, Station == "13" | Station == "8")
t138 <- t.test(Richness ~ Station, data = st138)
t138

st116 <- subset(inner_16S, Station == "11" | Station == "6")
t116 <- t.test(Richness ~ Station, data = st116)
t116

st138.2 <- subset(centre_16S, Station == "13" | Station == "8")
t138.2 <- t.test(Shannon ~ Station, data = st138.2)
t138.2

st78 <- subset(centre_16S, Station == "7" | Station == "8")
t78 <- t.test(Shannon ~ Station, data = st78)
t78

st47 <- subset(centre_16S, Station == "4" | Station == "7")
t47 <- t.test(Simpson ~ Station, data = st47)
t47

t78.2 <- t.test(Simpson ~ Station, data = st78)
t78.2

st112.2 <- subset(sub_mouth_18S, Station == "1" | Station == "12")
t112.2 <- t.test(Simpson ~ Station, data = st112.2)
t112.2

st110 <- subset(sub_mouth_18S, Station == "1" | Station == "10")
t110 <- t.test(Simpson ~ Station, data = st110)
t110


### look into the means of each fjord

# NextSeq 18S

# vMijen
mean(vM_N18$Richness)
sd(vM_N18$Richness)
mean(vM_N18$Shannon)
sd(vM_N18$Shannon)
mean(vM_N18$Simpson)
sd(vM_N18$Simpson)

aggregate(vM_N18$Richness, list(vM_N18$Station), FUN=mean)
aggregate(vM_N18$Richness, list(vM_N18$Station), FUN=sd)

aggregate(vM_N18$Shannon, list(vM_N18$Station), FUN=mean)
aggregate(vM_N18$Shannon, list(vM_N18$Station), FUN=sd)

aggregate(vM_N18$Simpson, list(vM_N18$Station), FUN=mean)
aggregate(vM_N18$Simpson, list(vM_N18$Station), FUN=sd)

# Kongs
mean(Ko_N18$Richness)
sd(Ko_N18$Richness)
mean(Ko_N18$Shannon)
sd(Ko_N18$Shannon)
mean(Ko_N18$Simpson)
sd(Ko_N18$Simpson)

aggregate(Ko_N18$Richness, list(Ko_N18$Station), FUN=mean)
aggregate(Ko_N18$Richness, list(Ko_N18$Station), FUN=sd)

aggregate(Ko_N18$Shannon, list(Ko_N18$Station), FUN=mean)
aggregate(Ko_N18$Shannon, list(Ko_N18$Station), FUN=sd)

aggregate(Ko_N18$Simpson, list(Ko_N18$Station), FUN=mean)
aggregate(Ko_N18$Simpson, list(Ko_N18$Station), FUN=sd)

# Wijde
mean(Wi_N18$Richness)
sd(Wi_N18$Richness)
mean(Wi_N18$Shannon)
sd(Wi_N18$Shannon)
mean(Wi_N18$Simpson)
sd(Wi_N18$Simpson)

aggregate(Wi_N18$Richness, list(Wi_N18$Station), FUN=mean)
aggregate(Wi_N18$Richness, list(Wi_N18$Station), FUN=sd)

aggregate(Wi_N18$Shannon, list(Wi_N18$Station), FUN=mean)
aggregate(Wi_N18$Shannon, list(Wi_N18$Station), FUN=sd)

aggregate(Wi_N18$Simpson, list(Wi_N18$Station), FUN=mean)
aggregate(Wi_N18$Simpson, list(Wi_N18$Station), FUN=sd)

# Rijp 
mean(Ri_N18$Richness)
sd(Ri_N18$Richness)
mean(Ri_N18$Shannon)
sd(Ri_N18$Shannon)
mean(Ri_N18$Simpson)
sd(Ri_N18$Simpson)

aggregate(Ri_N18$Richness, list(Ri_N18$Station), FUN=mean)
aggregate(Ri_N18$Richness, list(Ri_N18$Station), FUN=sd)

aggregate(Ri_N18$Shannon, list(Ri_N18$Station), FUN=mean)
aggregate(Ri_N18$Shannon, list(Ri_N18$Station), FUN=sd)

aggregate(Ri_N18$Simpson, list(Ri_N18$Station), FUN=mean)
aggregate(Ri_N18$Simpson, list(Ri_N18$Station), FUN=sd)



# NextSeq 16S

# vMijen
mean(vM_N16$Richness)
sd(vM_N16$Richness)
mean(vM_N16$Shannon)
sd(vM_N16$Shannon)
mean(vM_N16$Simpson)
sd(vM_N16$Simpson)

aggregate(vM_N16$Richness, list(vM_N16$Station), FUN=mean)
aggregate(vM_N16$Richness, list(vM_N16$Station), FUN=sd)

aggregate(vM_N16$Shannon, list(vM_N16$Station), FUN=mean)
aggregate(vM_N16$Shannon, list(vM_N16$Station), FUN=sd)

aggregate(vM_N16$Simpson, list(vM_N16$Station), FUN=mean)
aggregate(vM_N16$Simpson, list(vM_N16$Station), FUN=sd)

# Kongs
mean(Ko_N16$Richness)
sd(Ko_N16$Richness)
mean(Ko_N16$Shannon)
sd(Ko_N16$Shannon)
mean(Ko_N16$Simpson)
sd(Ko_N16$Simpson)

aggregate(Ko_N16$Richness, list(Ko_N16$Station), FUN=mean)
aggregate(Ko_N16$Richness, list(Ko_N16$Station), FUN=sd)

aggregate(Ko_N16$Shannon, list(Ko_N16$Station), FUN=mean)
aggregate(Ko_N16$Shannon, list(Ko_N16$Station), FUN=sd)

aggregate(Ko_N16$Simpson, list(Ko_N16$Station), FUN=mean)
aggregate(Ko_N16$Simpson, list(Ko_N16$Station), FUN=sd)

# Wijde
mean(Wi_N16$Richness)
sd(Wi_N16$Richness)
mean(Wi_N16$Shannon)
sd(Wi_N16$Shannon)
mean(Wi_N16$Simpson)
sd(Wi_N16$Simpson)

aggregate(Wi_N16$Richness, list(Wi_N16$Station), FUN=mean)
aggregate(Wi_N16$Richness, list(Wi_N16$Station), FUN=sd)

aggregate(Wi_N16$Shannon, list(Wi_N16$Station), FUN=mean)
aggregate(Wi_N16$Shannon, list(Wi_N16$Station), FUN=sd)

aggregate(Wi_N16$Simpson, list(Wi_N16$Station), FUN=mean)
aggregate(Wi_N16$Simpson, list(Wi_N16$Station), FUN=sd)

# Rijp 
mean(Ri_N16$Richness)
sd(Ri_N16$Richness)
mean(Ri_N16$Shannon)
sd(Ri_N16$Shannon)
mean(Ri_N16$Simpson)
sd(Ri_N16$Simpson)

aggregate(Ri_N16$Richness, list(Ri_N16$Station), FUN=mean)
aggregate(Ri_N16$Richness, list(Ri_N16$Station), FUN=sd)

aggregate(Ri_N16$Shannon, list(Ri_N16$Station), FUN=mean)
aggregate(Ri_N16$Shannon, list(Ri_N16$Station), FUN=sd)

aggregate(Ri_N16$Simpson, list(Ri_N16$Station), FUN=mean)
aggregate(Ri_N16$Simpson, list(Ri_N16$Station), FUN=sd)





























