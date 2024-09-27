## upsetr 

#install.packages("UpSetR")
library(UpSetR)
library(vegan)
library(ComplexHeatmap)

setwd("/isibhv/projects/AG_John/Expeditions/HE627_Data/Sequences_MiSeq/UpSetR/")
wd <- setwd("/isibhv/projects/AG_John/Expeditions/HE627_Data/Sequences_MiSeq/UpSetR/")


up_miseq <- read.csv("/isibhv/projects/AG_John/Expeditions/HE627_Data/Sequences_MiSeq/UpSetR/MiSeq_18S_ASVs_no_only_numbers.csv")
up_miseq <- as.data.frame(up_miseq)

#upset(fromExpression(up_miseq))

## create a binary table of presence + absence
up_pa <- decostand(up_miseq, method = "pa")

# into dataframe
up_pa <- as.data.frame(up_pa)

# merge columns into fjords

# Operators & and | perform element-wise operation producing result having length of the longer operand
# but different outocme if i use &
up_pa$van_Mijenfjorden <- as.integer(up_pa$St1.F02.1_S1|up_pa$St1.F02.2_S2|up_pa$St1.F02.3_S3|
                                       up_pa$St1.F3.1_S4|up_pa$St1.F3.2_S5|up_pa$St1.F3.3_S6|
                                       up_pa$St1.N_S10|up_pa$St1.P1_S7|                                 # deleted St1 P2
                                       up_pa$St1.P3_S9|up_pa$St2.F02.1_S11|up_pa$St2.F02.2_S12|
                                       up_pa$St2.F02.3_S13|up_pa$St2.F3.1_S44|up_pa$St2.F3.2_S45|
                                       up_pa$St2.F3.3_S46|up_pa$St2.P1_S77|up_pa$St2.P2_S78|
                                       up_pa$St2.P3_S79|up_pa$St3.F02.1_S14|up_pa$St3.F02.2_S15|
                                       up_pa$St3.F02.3_S16|up_pa$St3.F3.1_S47|up_pa$St3.F3.2_S48|
                                       up_pa$St3.F3.3_S49|up_pa$St3.P1_S80|up_pa$St3.P2_S81|
                                       up_pa$St3.P3_S82|up_pa$St4.F02.1_S17|up_pa$St4.F02.2_S18|
                                       up_pa$St4.F02.3_S19|up_pa$St4.F3.1_S50|up_pa$St4.F3.2_S51|
                                       up_pa$St4.F3.3_S52|up_pa$St4.P1_S83|up_pa$St4.P2_S84|
                                       up_pa$St4.P3_S85)

up_pa$Wijdefjorden <- as.integer(up_pa$St5.F02.1_S20|up_pa$St5.F02.2_S21|up_pa$St5.F02.3_S22|
                                   up_pa$St5.F3.1_S53|up_pa$St5.F3.2_S54|up_pa$St5.F3.3_S55|
                                   up_pa$St5.P1_S86|up_pa$St5.P2_S87|up_pa$St5.P3_S88|
                                   up_pa$St6.F02.1_S23|up_pa$St6.F02.2_S24|up_pa$St6.F02.3_S25|
                                   up_pa$St6.F3.1_S56|up_pa$St6.F3.2_S57|up_pa$St6.F3.3_S58|
                                   up_pa$St6.P1_S89|up_pa$St6.P2_S90|up_pa$St6.P3_S91|
                                   up_pa$St7.F02.1_S26|up_pa$St7.F02.2_S27|up_pa$St7.F02.3_S28|
                                   up_pa$St7.F3.1_S59|up_pa$St7.F3.2_S60|up_pa$St7.F3.3_S61|
                                   up_pa$St7.P1_S92|up_pa$St7.P2_S93|up_pa$St7.P3_S94)

up_pa$Rijpfjorden <- as.integer(up_pa$St8.F02.1_S29|up_pa$St8.F02.2_S30|up_pa$St8.F02.3_S31|
                                  up_pa$St8.F3.1_S62|up_pa$St1.F3.2_S5|up_pa$St8.F3.3_S64|
                                  up_pa$St8.P1_S95|up_pa$St8.P2_S96|up_pa$St10.F02.1_S32|
                                  up_pa$St10.F02.2_S33|up_pa$St10.F02.3_S34|up_pa$St10.F3.1_S65|
                                  up_pa$St10.F3.2_S66|up_pa$St10.F3.3_S67)

up_pa$Kongsfjorden <- as.integer(up_pa$St11.F02.1_S35|up_pa$St11.F02.2_S36|up_pa$St11.F02.3_S37|
                                   up_pa$St11.F3.1_S68|up_pa$St11.F3.2_S69|up_pa$St11.F3.3_S70|
                                   up_pa$St12.F02.1_S38|up_pa$St12.F02.2_S39|up_pa$St12.F02.3_S40|
                                   up_pa$St12.F3.1_S71|up_pa$St12.F3.2_S72|up_pa$St12.F3.3_S73|
                                   up_pa$St13.F02.1_S41|up_pa$St13.F02.2_S42|up_pa$St13.F02.3_S43|
                                   up_pa$St13.F3.1_S74|up_pa$St13.F3.2_S75|up_pa$St13.F3.3_S76)

head(up_pa)


## set seed ? 
##


# define the fjords and filter out the 0s
fjords <- up_pa[,c(96:99)]
fjords <- fjords %>% filter_all(any_vars(. !=0))
filter_all(fjords, any_vars(. != 0))
is.na(fjords)


## create a combination matrix for the plot
co_dis_fjords <- make_comb_mat(fjords)
co_dis_fjords
# 4 sets and 15 combos
# mode: distinct


UpSet(co_dis_fjords, set_order = c("van_Mijenfjorden", "Kongsfjorden", "Wijdefjorden", "Rijpfjorden"), 
                                         comb_order = order(comb_size(co_dis_fjords)))
ggsave("p1_1.png")
# doesnt work # blank


p2 <- print(UpSet(co_dis_fjords, nsets = 4, set_order = c("van_Mijenfjorden", "Kongsfjorden", "Wijdefjorden", "Rijpfjorden"), 
                  comb_order = order(comb_size(co_dis_fjords))))
ggsave("p1.png")
# doesnt work # blank


UpSet(fjords)
UpSet(up_pa)
# both dont work

UpSet(up_miseq)
#doesnt work

UpSet(co_dis_fjords, set_order = fjords, comb_order = order(comb_size(co_dis_fjords)))
# doesnt work

UpSet(co_dis_fjords, sets = c("van_Mijenfjorden", "Kongsfjorden", "Wijdefjorden", "Rijpfjorden"))



upset1 <- UpSet(co_dis_fjords)
ggsave("upset1.png") 
## nothing if i run it again 

upset2 <- UpSet(co_dis_fjords, comb_order = order(comb_size(co_dis_fjords)))
ggsave("upset2.png") ## both look the same

upset3 <- UpSet(co_dis_fjords, order_by = "freq")
# doesnt work 


co_inter_fjords <- make_comb_mat(fjords, mode = "intersect")
co_inter_fjords

upset4 <- UpSet(co_inter_fjords)
ggsave("upset4.png")    # looks just like the other two






co_uni_fjords <- make_comb_mat(fjords, mode = "union")
co_uni_fjords

upset5 <- UpSet(co_uni_fjords)
ggsave("upset5.png")

## all look the same


upset_try <- UpSet(fjords)
upset_try2 <- UpSet(up_pa)

# example of list input (list of named vectors)
listInput <- list(one = c(1, 2, 3, 5, 7, 8, 11, 12, 13), two = c(1, 2, 4, 5, 
                                                                 10), three = c(1, 5, 6, 7, 8, 9, 10, 12, 13))

# example of expression input
expressionInput <- c(one = 2, two = 1, three = 2, `one&two` = 1, `one&three` = 4, 
                     `two&three` = 1, `one&two&three` = 2)

upset(fromList(listInput), order.by = "freq")
ggsave("try.png")


