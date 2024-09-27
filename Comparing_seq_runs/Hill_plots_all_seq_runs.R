## Statistical tests Hill numbers - all sequences & fjords

library(stats)
library(ggplot2)
library(gridExtra)
library(readxl)
library(grid)

setwd("Z:/AG_John/Expeditions/HE627_Data/Comparison_all_seq_runs/Hill_numbers/")

data <- read_xlsx("Hill_numbers_all_fjords.xlsx", sheet = 3)
head(data)


names(data)[names(data)=="s.e....7"] <- "Shannon_se"
names(data)[names(data)=="s.e....9"] <- "Simpson_se"
names(data)[names(data)=="s.e....11"] <- "Richness_se"

str(data)
data$Station <- as.character(data$Station)
data$Sequencer <- as.factor(data$Sequencer)

#manova(cbind(rv1, rv2, …) ~ iv, data)
#manova(cbinf(response1, response2) ~factor, data)
try <- manova(cbind(Shannon, Simpson) ~ Sample, data)
summary(try)
summary.aov(try)

try2 <- manova(cbind(Shannon, Simpson) ~Fjord, data)
summary(try2)

#try3 <- manova(cbind(Fjord$Wijdefjorden, Fjord$Kongsfjorden), ~ Shannon, data)


try4 <- aov(Shannon ~ Fjord, data = data)
summary(try4)


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
## compare stations? 

ggplot(data, aes(x = Station, y = Richness, color = Station)) + geom_boxplot(fill = c("#ED90A4", "#E59884", "#D6a166", "#c0ab52", "#a3b353", "#7eba68", "#4fbf85",
                                                                                               "#00c1a3", "#00c0c0", "#28bbd7", "#6fb1e7", "#a1a5ec", "#c699e7"), color = "black") + 
                                                                                                 scale_x_discrete(limits = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13"))
ggsave("stat_rich_pretty.png")


ggplot(data, aes(x = Station, y = Shannon)) + geom_boxplot()
ggsave("stat_sha.png")

ggplot(data, aes(x = Station, y = Simpson)) + geom_boxplot()
ggsave("stat_simp.png")



###
## compare methods?

# compare methods via PCA later

#ggplot(data, aes(x = Sample_method, y = Richness)) + geom_boxplot()
#ggsave("meth_rich.png")

#ggplot(data, aes(x = Sample_method, y = Shannon)) + geom_boxplot()
#ggsave("meth_sha.png")

#ggplot(data, aes(x = Sample_method, y = Simpson)) + geom_boxplot()
#ggsave("meth_simp.png")




###
## compare inner fjords to outer fjords

vMijen <- subset(data, Fjord == "van Mijenfjorden")           # 92 samples
Kongs <- subset(data, Fjord == "Kongsfjorden")                # 60 samples
Wijde <- subset(data, Fjord == "Wijdefjorden")                # 70 samples
Rijp <- subset(data, Fjord == "Rijpfjorden")                  # 43 samples



## van Mijenfjorden

#ggplot(vMijen, aes(x = Station, y = Richness)) + geom_boxplot(fill = c("#ED90A4", "#E59884", "#D6a166", "#c0ab52"), color = "black") + 
#  scale_x_discrete(limits = c("1", "2", "3", "4")) + labs(title = "Species Richness van Mijenfjorden") + 
#  facet_grid(cols = vars(Seqeuencer))

#vMijen$Sequencer_ID <- as.factor(vMijen$Sequencer)
#a <- ggplot(vMijen, aes(x = Station, y = Richness)) + geom_boxplot(fill = c("#ED90A4", "#E59884", "#D6a166", "#c0ab52"), color = "black") + 
#  scale_x_discrete(limits = c("1", "2", "3", "4")) + 
#  labs(title = "Species Richness van Mijenfjorden") #+ 

#a + ggplot2::facet_grid(Seqeuencer_ID~.)
#ggsave("richness_vmijen.png")

vM_Mi <- subset(vMijen, Sequencer == "MiSeq-18S")
vM_N16 <- subset(vMijen, Sequencer == "NextSeq-16S")
vM_N18 <- subset(vMijen, Sequencer == "NextSeq-18S")


## richness
pvm <- ggplot(vM_Mi, aes(x = Station, y = Richness)) + geom_boxplot(fill = c("#ED90A4", "#E59884", "#c0ab52",  "#D6a166"), color = "black") + 
  scale_x_discrete(limits = c("1", "2", "4", "3")) + labs(title = "MiSeq 18S") + theme_light() + 
  scale_y_continuous(limits=c(0, 1200)) + theme(plot.title = element_text(size = 10)) 

pvn16 <- ggplot(vM_N16, aes(x = Station, y = Richness)) + geom_boxplot(fill = c("#ED90A4", "#E59884", "#c0ab52", "#D6a166"), color = "black") + 
  scale_x_discrete(limits = c("1", "2", "4", "3")) + labs(title = "NextSeq 16S") + theme_light() + 
  scale_y_continuous(limits=c(0, 1200)) + theme(plot.title = element_text(size = 10)) 

pvn18 <- ggplot(vM_N18, aes(x = Station, y = Richness)) + geom_boxplot(fill = c("#ED90A4", "#E59884", "#c0ab52", "#D6a166"), color = "black") + 
  scale_x_discrete(limits = c("1", "2", "4", "3")) + labs(title = "NextSeq 18S") + theme_light() + 
  scale_y_continuous(limits=c(0, 1200)) + theme(plot.title = element_text(size = 10)) 

#tg <- textGrob('Species Richness - van Mijenfjorden', gp = gpar(fontsize = 14, fontface = 'bold'))
#sg <- textGrob('From outer fjord (St.1) to inner fjord (St.4)', gp = gpar(fontsize = 12))

#lt <- list(tg, sg)
#heights <- do.call(unit.c, lapply(lt, function(.g) 1.5*grobHeight(.g)))
#titles <- gtable::gtable_matrix('title', 
#                                grobs = matrix(lt, ncol=1), 
#                                widths = unit(1,'npc'),
#                                heights = heights)

g <- arrangeGrob(pvm, pvn18, pvn16, ncol = 3, top = grid::textGrob("Species Richness - van Mijenfjorden", x = 0.02, 
                                                                   hjust = 0, gp = gpar(fontsize = 13)))
#g2 <- arrangeGrob(pvm, pvn18, pvn16, ncol = 3, top = textGrob(lt), x = 0.02, hjust = 0)
#g3 <- arrangeGrob(pvm, pvn18, pvn16, ncol = 3, top = c(grid::textGrob("Species Richness - van Mijenfjorden", x = 0.02, 
#                                                                    hjust = 0, gp = gpar(fontface = "bold", fontsize = 13)), grid::textGrob("From outer fjord (St.1) to inner fjord (St.4)",
#                                                                                                                               x = 0.02, hjust = 0, gp = gpar(fontsize = 10))))
#g4 <- arrangeGrob(pvm, pvn18, pvn16, ncol = 3, top = grid::textGrob(lt), x = 0.02, hjust = 0)
#g5 <- arrangeGrob(pvm, pvn18, pvn16, ncol = 3, top = titles)

#top = grid::textGrob("Title", x = 0, hjust = 0))
#?textGrob
#?gpar
ggsave("vMijen_rich_stations_2.png", g, width = 8, height = 5)


## shannon
shvm <- ggplot(vM_Mi, aes(x = Station, y = Shannon)) + geom_boxplot(fill = c("#ED90A4", "#E59884", "#c0ab52", "#D6a166"), color = "black") + 
  scale_x_discrete(limits = c("1", "2", "4", "3")) + labs(title = "MiSeq 18S") + theme_light() + 
  scale_y_continuous(limits=c(0, 150)) + theme(plot.title = element_text(size = 10)) 

shvn16 <- ggplot(vM_N16, aes(x = Station, y = Shannon)) + geom_boxplot(fill = c("#ED90A4", "#E59884", "#c0ab52", "#D6a166"), color = "black") + 
  scale_x_discrete(limits = c("1", "2", "4", "3")) + labs(title = "NextSeq 16S") + theme_light() + 
  scale_y_continuous(limits=c(0, 150)) + theme(plot.title = element_text(size = 10)) 

shvn18 <- ggplot(vM_N18, aes(x = Station, y = Shannon)) + geom_boxplot(fill = c("#ED90A4", "#E59884", "#c0ab52", "#D6a166"), color = "black") + 
  scale_x_discrete(limits = c("1", "2", "4", "3")) + labs(title = "NextSeq 18S") + theme_light() + 
  scale_y_continuous(limits=c(0, 150)) + theme(plot.title = element_text(size = 10)) 

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
  scale_x_discrete(limits = c("1", "2", "4", "3")) + labs(title = "NextSeq 16S") + theme_light() + 
  scale_y_continuous(limits=c(0, 60)) + theme(plot.title = element_text(size = 10)) 

sivn18 <- ggplot(vM_N18, aes(x = Station, y = Simpson)) + geom_boxplot(fill = c("#ED90A4", "#E59884", "#c0ab52", "#D6a166"), color = "black") + 
  scale_x_discrete(limits = c("1", "2", "4", "3")) + labs(title = "NextSeq 18S") + theme_light() + 
  scale_y_continuous(limits=c(0, 60)) + theme(plot.title = element_text(size = 10)) 

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

kn18 <- ggplot(Ko_N18, aes(x = Station, y = Richness)) + geom_boxplot(fill = c("#a1a5ec", "#c699e7", "#6fb1e7"), color = "black") + 
  scale_x_discrete(limits = c("12", "13", "11")) + labs(title = "NextSeq 18S") + theme_light() + 
  scale_y_continuous(limits=c(0, 1650)) + theme(plot.title = element_text(size = 10)) 

kn16 <- ggplot(Ko_N16, aes(x = Station, y = Richness)) + geom_boxplot(fill = c("#a1a5ec", "#c699e7", "#6fb1e7"), color = "black") + 
  scale_x_discrete(limits = c("12", "13", "11")) + labs(title = "NextSeq 16S") + theme_light() + 
  scale_y_continuous(limits=c(0, 1650)) + theme(plot.title = element_text(size = 10)) 

j <- arrangeGrob(kmi, kn18, kn16, ncol = 3, top = grid::textGrob("Species richness - Kongsfjorden", x = 0.02, 
                                                                 hjust = 0, gp = gpar(fontsize = 13)))

ggsave("kongs_rich_stat.png", j, width = 8, height = 5)


## shannon 

skmi <- ggplot(Ko_Mi, aes(x = Station, y = Shannon)) + geom_boxplot(fill = c("#a1a5ec", "#c699e7", "#6fb1e7"), color = "black") + 
  scale_x_discrete(limits = c("12", "13", "11")) + labs(title = "MiSeq 18S") + theme_light() + 
  scale_y_continuous(limits=c(0, 220)) + theme(plot.title = element_text(size = 10)) 

skn18 <- ggplot(Ko_N18, aes(x = Station, y = Shannon)) + geom_boxplot(fill = c("#a1a5ec", "#c699e7", "#6fb1e7"), color = "black") + 
  scale_x_discrete(limits = c("12", "13", "11")) + labs(title = "NextSeq 18S") + theme_light() + 
  scale_y_continuous(limits=c(0, 220)) + theme(plot.title = element_text(size = 10)) 

skn16 <- ggplot(Ko_N16, aes(x = Station, y = Shannon)) + geom_boxplot(fill = c("#a1a5ec", "#c699e7", "#6fb1e7"), color = "black") + 
  scale_x_discrete(limits = c("12", "13", "11")) + labs(title = "NextSeq 16S") + theme_light() + 
  scale_y_continuous(limits=c(0, 220)) + theme(plot.title = element_text(size = 10)) 

k <- arrangeGrob(skmi, skn18, skn16, ncol = 3, top = grid::textGrob("Shannon Diversity - Kongsfjorden", x = 0.02, 
                                                                    hjust = 0, gp = gpar(fontsize = 13)))

ggsave("kongs_shan_stat.png", k, width = 8, height = 5)


## simpson

sikmi <- ggplot(Ko_Mi, aes(x = Station, y = Simpson)) + geom_boxplot(fill = c("#a1a5ec", "#c699e7", "#6fb1e7"), color = "black") + 
  scale_x_discrete(limits = c("12", "13", "11")) + labs(title = "MiSeq 18S") + theme_light() + 
  scale_y_continuous(limits=c(0, 75)) + theme(plot.title = element_text(size = 10)) 

sikn18 <- ggplot(Ko_N18, aes(x = Station, y = Simpson)) + geom_boxplot(fill = c("#a1a5ec", "#c699e7", "#6fb1e7"), color = "black") + 
  scale_x_discrete(limits = c("12", "13", "11")) + labs(title = "NextSeq 18S") + theme_light() + 
  scale_y_continuous(limits=c(0, 75)) + theme(plot.title = element_text(size = 10)) 

sikn16 <- ggplot(Ko_N16, aes(x = Station, y = Simpson)) + geom_boxplot(fill = c("#a1a5ec", "#c699e7", "#6fb1e7"), color = "black") + 
  scale_x_discrete(limits = c("12", "13", "11")) + labs(title = "NextSeq 16S") + theme_light() + 
  scale_y_continuous(limits=c(0, 75)) + theme(plot.title = element_text(size = 10)) 

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
  scale_x_discrete(limits = c("5", "7", "6")) + labs(title = "NextSeq 18S") + theme_light() + 
  scale_y_continuous(limits=c(0, 1100)) + 
  theme(plot.title = element_text(size = 10)) 

wn16 <- ggplot(Wi_N16, aes(x = Station, y = Richness)) + geom_boxplot(fill = c("#a3b353", "#4fbf85", "#7eba68"), color = "black") + 
  scale_x_discrete(limits = c("5", "7", "6")) + labs(title = "NextSeq 16S") + theme_light() + 
  scale_y_continuous(limits=c(0, 1100)) + 
  theme(plot.title = element_text(size = 10)) 

m <- arrangeGrob(wmi, wn18, wn16, ncol = 3, top = grid::textGrob("Species richness - Wijdefjorden", x = 0.02, 
                                                                 hjust = 0, gp = gpar(fontsize = 13)))

ggsave("wijde_rich_stat.png", m, width = 8, height = 5)


## shannon 

swmi <- ggplot(Wi_Mi, aes(x = Station, y = Shannon)) + geom_boxplot(fill = c("#a3b353", "#4fbf85", "#7eba68"), color = "black") + 
  scale_x_discrete(limits = c("5", "7", "6")) + labs(title = "MiSeq 18S") + theme_light() + 
  scale_y_continuous(limits=c(0, 160)) + 
  theme(plot.title = element_text(size = 10)) 

swn18 <- ggplot(Wi_N18, aes(x = Station, y = Shannon)) + geom_boxplot(fill = c("#a3b353", "#4fbf85", "#7eba68"), color = "black") + 
  scale_x_discrete(limits = c("5", "7", "6")) + labs(title = "NextSeq 18S") + theme_light() + 
  scale_y_continuous(limits=c(0, 160)) + 
  theme(plot.title = element_text(size = 10)) 

swn16 <- ggplot(Wi_N16, aes(x = Station, y = Shannon)) + geom_boxplot(fill = c("#a3b353", "#4fbf85", "#7eba68"), color = "black") + 
  scale_x_discrete(limits = c("5", "7", "6")) + labs(title = "NextSeq 16S") + theme_light() + 
  scale_y_continuous(limits=c(0, 160)) + 
  theme(plot.title = element_text(size = 10)) 

n <- arrangeGrob(swmi, swn18, swn16, ncol = 3, top = grid::textGrob("Shannon Diversity - Wijdefjorden", x = 0.02, 
                                                                    hjust = 0, gp = gpar(fontsize = 13)))

ggsave("wijde_sha_stat.png", n, width = 8, height = 5)


## simpson 

siwmi <- ggplot(Wi_Mi, aes(x = Station, y = Simpson)) + geom_boxplot(fill = c("#a3b353", "#4fbf85", "#7eba68"), color = "black") + 
  scale_x_discrete(limits = c("5", "7", "6")) + labs(title = "MiSeq 18S") + theme_light() + 
  scale_y_continuous(limits=c(0, 70)) + 
  theme(plot.title = element_text(size = 10)) 

siwn18 <- ggplot(Wi_N18, aes(x = Station, y = Simpson)) + geom_boxplot(fill = c("#a3b353", "#4fbf85", "#7eba68"), color = "black") + 
  scale_x_discrete(limits = c("5", "7", "6")) + labs(title = "NextSeq 18S") + theme_light() + 
  scale_y_continuous(limits=c(0, 70)) + 
  theme(plot.title = element_text(size = 10)) 

siwn16 <- ggplot(Wi_N16, aes(x = Station, y = Simpson)) + geom_boxplot(fill = c("#a3b353", "#4fbf85", "#7eba68"), color = "black") + 
  scale_x_discrete(limits = c("5", "7", "6")) + labs(title = "NextSeq 16S") + theme_light() + 
  scale_y_continuous(limits=c(0, 70)) + 
  theme(plot.title = element_text(size = 10)) 

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

rn18 <- ggplot(Ri_N18, aes(x = Station, y = Richness)) + geom_boxplot(fill = c("#28bbd7", "#00c1a3"), color = "black") + 
  scale_x_discrete(limits = c("10", "8")) + labs(title = "NextSeq 18S") + theme_light() + 
  scale_y_continuous(limits=c(0, 840)) + 
  theme(plot.title = element_text(size = 10)) 

rn16 <- ggplot(Ri_N16, aes(x = Station, y = Richness)) + geom_boxplot(fill = c("#28bbd7", "#00c1a3"), color = "black") + 
  scale_x_discrete(limits = c("10", "8")) + labs(title = "NextSeq 16S") + theme_light() + 
  scale_y_continuous(limits=c(0, 840)) + 
  theme(plot.title = element_text(size = 10)) 

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
  scale_x_discrete(limits = c("10", "8")) + labs(title = "NextSeq 18S") + theme_light() + 
  scale_y_continuous(limits=c(0, 170)) + 
  theme(plot.title = element_text(size = 10)) 

srn16 <- ggplot(Ri_N16, aes(x = Station, y = Shannon)) + geom_boxplot(fill = c("#28bbd7", "#00c1a3"), color = "black") + 
  scale_x_discrete(limits = c("10", "8")) + labs(title = "NextSeq 16S") + theme_light() + 
  scale_y_continuous(limits=c(0, 170)) + 
  theme(plot.title = element_text(size = 10)) 

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
  scale_x_discrete(limits = c("10", "8")) + labs(title = "NextSeq 18S") + theme_light() + 
  scale_y_continuous(limits=c(0, 75)) + 
  theme(plot.title = element_text(size = 10)) 

sirn16 <- ggplot(Ri_N16, aes(x = Station, y = Simpson)) + geom_boxplot(fill = c("#28bbd7", "#00c1a3"), color = "black") + 
  scale_x_discrete(limits = c("10", "8")) + labs(title = "NextSeq 16S") + theme_light() + 
  scale_y_continuous(limits=c(0, 75)) + 
  theme(plot.title = element_text(size = 10)) 

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
