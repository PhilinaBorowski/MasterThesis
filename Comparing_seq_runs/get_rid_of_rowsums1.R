## get rid of ASVs w only 1 count
## merged NextSeq + MiSeq dataset + 16S

library(dplyr)

setwd("/isibhv/projects/AG_John/Expeditions/HE627_Data/Comparison_all_seq_runs/Test_w_normalisation_of_mega_values/")

## miseq 18s

data_mi <- read.csv("MiSeq_18S_normalisation.csv")

str(data_mi)
data_mi[2:95] <- lapply(data_mi[, c(2:95)], as.numeric)

# figure out, which rows only have 1 ASV
?rowsum()
?rowSums()

rowSums(data_mi[2:95])
rowSums(data_mi[, c(2:95)])   # same data output


data_mi$row_sum <- rowSums(data_mi[ , c(2:95)])
str(data_mi)

# get rid of all rows w rowsum = 1
sub_test <- subset(data_mi, row_sum > 1)

rowSums(sub_test[2:95])

sub_row_sum <- subset(sub_test, row_sum == 1)

write.csv(sub_test, "MiSeq_18S_reduced_no_1rowsum.csv")


4111-2523


### nextseq 18s

data_ne <- read.csv("NextSeq_18S_normalisation.csv")

str(data_ne)
data_ne[2:134] <- lapply(data_ne[, c(2:134)], as.numeric)

rowSums(data_ne[2:134])


data_ne$row_sum <- rowSums(data_ne[ , c(2:134)])
str(data_ne)

# get rid of all rows w rowsum = 1
sub_test_ne <- subset(data_ne, row_sum > 1)

rowSums(sub_test_ne[2:134])

sub_row_sum_ne <- subset(sub_test_ne, row_sum == 1)

write.csv(sub_test_ne, "NextSeq_18S_reduced_no_1rowsum.csv")

15738-10425


## nextseq 16s

setwd("/isibhv/projects/AG_John/Expeditions/HE627_Data/Sequences_NextSeq_Neu/Run_53_UJOHN_007_Prokaryotes_HE627_Dedmar_Jakob/HE627_16S/")

data_ne16 <- read.csv("NextSeq_16S_ampli_cleaned.csv")

str(data_ne16)
data_ne16[2:37] <- lapply(data_ne16[, c(2:37)], as.numeric)

rowSums(data_ne16[2:37])


data_ne16$row_sum <- rowSums(data_ne16[ , c(2:37)])
str(data_ne16)

# get rid of all rows w rowsum = 1
sub_test_ne16 <- subset(data_ne16, row_sum > 1)

rowSums(sub_test_ne16[2:37])

sub_row_sum_ne16 <- subset(sub_test_ne16, row_sum == 1)

write.csv(sub_test_ne16, "NextSeq_16S_reduced_no_1rowsum.csv")

2044-2044

