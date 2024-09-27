## Statistical tests Hill numbers - all sequences & fjords

library(stats)
library(ggplot2)

setwd("/isibhv/projects/AG_John/Expeditions/HE627_Data/")

data <- read_xlsx("/isibhv/projects/AG_John/Expeditions/HE627_Data/Hill_numbers_all_fjords.xlsx", sheet = 2)
head(data)


names(data)[names(data)=="s.e....7"] <- "Shannon_se"
names(data)[names(data)=="s.e....9"] <- "Simpson_se"
names(data)[names(data)=="s.e....11"] <- "Richness_se"


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



?boxplot()
boxplot(data$Fjord ~ data$Shannon)
