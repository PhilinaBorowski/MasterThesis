## work w phyloseq - MiSeq data
##

install.packages("viridis")  # Install
library(viridis)           # Load

library(dplyr)
library(ggplot2)
library(phyloseq)
library(tibble)
#library(microbiomeMarker)
library(ape)
library(vegan)
library(microViz)
library(tidyr)
library(zCompositions)
library(compositions)


# is there a function to see which packages were actually used? 


## set wd
setwd("/isibhv/projects/AG_John/Expeditions/HE627_Data/Sequences_MiSeq/Phyloseq_analysis/")
wd <- ("/isibhv/projects/AG_John/Expeditions/HE627_Data/Sequences_MiSeq/Phyloseq_analysis/")

##### get rid of Kingdom = NA #####
cleaning <- read.csv("MiSeq_18S_reduced_no_1rowsum_94.csv", row.names = 1)
str(cleaning)

totallyclean <- cleaning %>% drop_na(Kingdom) # got rid of 21 ASVs
2523-2502

write.csv(totallyclean, "MiSeq_18S_totally_cleaned_ready_for_phyloseq_94.csv")

## manually create asvs + taxa csvs

## read them in
asvs <- read.csv("MiSeq_18S_Abundances_wo_singl_metaz_for_phyloseq_analysis.csv")           ## need the ones wo metazoa + singletons !!!
print(asvs)
taxa <- read.csv("MiSeq_18S_ASV_taxonomy_wo_singl_metaz_for_phyloseq_analysis.csv")                 ###ASV0001 - ASV2502
meta_samples <- read.csv("meta_phyloseq_18S_miseq.csv")  

## add row names

asvs <- asvs %>% tibble::column_to_rownames("otu")
taxa <- taxa %>% tibble::column_to_rownames("otu")
meta_samples <- meta_samples %>% tibble::column_to_rownames("sample")


#### --------------
##preparation for phyloseq / normalisation / pca later



## prep for CLR + CLR 
#d.czm <- cmultRepl(t(d.1),  label=0, method="CZM", frac= 0.99, adjust=FALSE, z.warning = 0.99)
#d.clr <- t(apply(d.czm, 1, function(x){log(x) - mean(log(x))})) #centered log ratio transformation

d.czm <- cmultRepl(t(asvs),  label=0, method="CZM", frac= 0.999, adjust=FALSE, z.warning = 0.999)
d.clr <- (apply(d.czm, 1, function(x){log(x) - mean(log(x))})) #centered log ratio transformation
#lapply(d.clr, class)%>%unlist() # but still some nas
print(d.clr)

clr2 <- clr(d.czm)
print(clr2)
clr2 <- t(as.data.frame(clr2))
print(clr2)


#move all into positive values
asvs_norm <- apply(d.clr, 2, function(x) x - min(x[x < 0]))

asvs_norm2 <- apply(clr2, 2, function(x) x - min(x[x < 0]))

Warning message:
In min(x[x < 0]) : no non-missing arguments to min; returning Inf

#make it relative
asvs_norm_P <- apply(asvs_norm, 2, function(x) {x/sum(x)}) 

asvs_norm_P2 <- apply(asvs_norm2, 2, function(x) {x/sum(x)})






# create a physloseq object for "normal" values and for % values

## creating phyloseq objects

asvs_norm <- as.matrix(asvs_norm)
print(asvs_norm)
# St 6 P1 + P3 & St 8 P1 & P2 now NAs
taxa <- as.matrix(taxa)
#taxa <- as.data.frame(taxa)
str(asvs_norm)

OTU = otu_table(asvs_norm, taxa_are_rows = T)
TAX = tax_table(taxa)
samples = sample_data(meta_samples)

## do the names match? 
sample_names(OTU)           # all - are now a . 
sample_names(samples)        

MiSeq_18S <- phyloseq(OTU, TAX, samples)
MiSeq_18S
# 2502 taxa + 94 samples
# 94 samples w 19 variables
# 2502 taxa w 9 tax ranks 

#summarize_phyloseq(MiSeq_18S)   # which package? 
#Total number of reads = 3826300"
#Min. number of reads = 682"
#Max. number of reads = 997727"

sample_names(MiSeq_18S)
rank_names(MiSeq_18S)
sample_variables(MiSeq_18S)

str(MiSeq_18S)


### #w / clr 2

asvs_norm2 <- as.matrix(asvs_norm2)
print(asvs_norm2)
# St8 P1 now -inf

OTU.clr = otu_table(asvs_norm2, taxa_are_rows = T)

## do the names match? 
sample_names(OTU)           # all - are now a . 
sample_names(samples)        

MiSeq_18S2 <- phyloseq(OTU.clr, TAX, samples)
MiSeq_18S2


## create a phyloseq pbject for the normalised values

asvs_norm_P <- as.matrix(asvs_norm_P)
str(asvs_norm_P)

OTU2 = otu_table(asvs_norm_P, taxa_are_rows = T)

MiSeq_18SP <- phyloseq(OTU2, TAX, samples)
MiSeq_18SP




## rarefying? 
# basic curves

#rarecurve(t(otu_table(MiSeq_18S)), step=100)
# doenst work, cause logical subscript too long





# graphic summary

#png("Mi_plot.png")
#plot_bar(MiSeq_18S)
#dev.off()

#png("sampl_tech.png")
#plot_bar(MiSeq_18S, "sampling_method", "Abundance", title= "Abundances per sampling method")
#dev.off()

#png("stat.png")
#plot_bar(MiSeq_18S, "stat", "Abundance", title= "Abundances per station")
#dev.off()


## try to plot 

#png("MiSeq_familiy.png", width = 3800, height = 800)
#plot_bar(MiSeq_18S, fill = "Family")   # takes long, cant see shit
#dev.off()

#png("MiSeq_class.png", width = 1800, height = 800)
#plot_bar(MiSeq_18S, fill = "Class")    # takes long, cant see shit 
#dev.off()

#png("MiSeq_phylum.png", width = 800, height = 800)
#plot_bar(MiSeq_18S, fill = "Phylum")
#dev.off()


#png("rarecurves.png")
#rarecurve(t(otu_table(MiSeq_18S)), step=100, cex=0.5)


## remove sample St7.F3.3, St13.F02.1, St2.F3.1
#MiSeq_18S_cleaned <- subset_samples(MiSeq_18S, sample_names(MiSeq_18S) !="St7.F3.3" & 
#                                      sample_names(MiSeq_18S) != "St13.F02.1" & 
#                                      sample_names(MiSeq_18S) != "St2.F3.1")



## plot 

png("Mi_plot_norm.png")
plot_bar(MiSeq_18S)
dev.off()
# 2 samples w weird values below 0

# "the more negative a value is, the more likely that it was zero, or very small, 
# in the original “raw” count matrix. For most distances and hypotheses, 
# these values are probably not very important, or even negligible"




png("MiSeq_phylum_try2.png", width = 1500)
plot_bar(MiSeq_18S, fill = "Phylum") + geom_bar(aes(fill=Phylum), stat="identity", position="stack") + labs(title = "MiSeq 18S - Phylum per sample") +
                                        theme_light() + scale_fill_viridis(discrete = T, option = "C")
#Removed 10008 rows containing missing values or values outside the scale range (`geom_bar()`). 
dev.off()

png("MiSeq_phylum_try3.png", width = 1000, height = 1500)
plot_bar(MiSeq_18S, fill = "Phylum") + geom_bar(aes(fill=Phylum), stat="identity", position="stack") + coord_flip() + theme_light()
#Removed 10008 rows containing missing values or values outside the scale range (`geom_bar()`). 
dev.off()

png("MiSeq_phylum_try2_P.png", width = 1500)
plot_bar(MiSeq_18SP, fill = "Phylum") + geom_bar(aes(fill=Phylum), stat="identity", position="stack") 
#Removed 10008 rows containing missing values or values outside the scale range (`geom_bar()`). 
dev.off()

print(MiSeq_18S@otu_table)
# St6 P1 + P3 St8 P1 + P2 NAs!!

png("MiSeq_phylum_try_w_clr2.png", width = 1500)
plot_bar(MiSeq_18S2, fill = "Phylum") + geom_bar(aes(fill=Phylum), stat="identity", position="stack")
#Removed 10008 rows containing missing values or values outside the scale range (`geom_bar()`). 
dev.off()

print(MiSeq_18S2@otu_table)


## try to merge and only look at phylum
mi.phylum = tax_glom(MiSeq_18S_norm, taxrank="Phylum", NArm=FALSE)
#Found more than one class "phylo" in cache; using the first, from namespace 'phyloseq'
#Also defined by ‘tidytree’
mi.phylum

png("mi_phylum_sub.png", width = 1500)
plot_bar(mi.phylum, fill="Phylum") + theme_light()
dev.off()
# like how it looks lol 
# but ofc loss of info 

png("mi_phylum_sub2.png", width = 1500)
plot_bar(mi.phylum, fill="Phylum") 
dev.off()
# sample names much nicer!!




## try to merge triplicates together? and look at them together? 
## try mean of sample values

Mi_18_merged <- merge_samples(MiSeq_18S_norm, "sample_merg", fun = mean)
#sample_names(MiSeq_18S_norm)
#sample_variables(MiSeq_18S_norm)

#In asMethod(object) : NAs introduced by coercion
Mi_18_merged

png("Mi_bar_merged.png")
plot_bar(Mi_18_merged, fill="Phylum", title = "Merged samples (mean of replicates)") + geom_bar(aes(color=Phylum, fill=Phylum), stat="identity", position="stack")
dev.off()



## or try to sum up samples

Mi_18_merged_add <- merge_samples(MiSeq_18S_norm, "sample_merg", fun = sum)
#In asMethod(object) : NAs introduced by coercion
Mi_18_merged_add

png("Mi_bar_merged_sum.png")
plot_bar(Mi_18_merged_add, fill="Phylum", title = "Merged samples (sum of replicates)") + geom_bar(aes(color=Phylum, fill=Phylum), stat="identity", position="stack")
dev.off()

## look the same? 


png("Mi_bar.png")
plot_bar(MiSeq_18S_norm, fill="Phylum", title = "CLR transformed") + geom_bar(aes(color=Phylum, fill=Phylum), stat="identity", position="stack")
dev.off()

png("Mi_bar2.png", width = 1500)
plot_bar(MiSeq_18S_norm, fill="Phylum", title = "CLR transformed") + geom_bar(aes(color=Phylum, fill=Phylum), stat = "identity")
dev.off()





## ------------------

## pca for beta diversity 
?distance
# no aitchison
# but https://search.r-project.org/CRAN/refmans/vegan/html/vegdist.html here
# vegdist w aitchison distance 

# > wunifrac_dist = phyloseq::distance(ps.rarefied, method="unifrac", weighted=F)
# > ordination = ordinate(ps.rarefied, method="PCoA", distance=wunifrac_dist)
# > plot_ordination(ps.rarefied, ordination, color="Season") + theme(aspect.ratio=1)


MiSeq_18S_trans_clr <- transform(MiSeq_18S, transform = "clr")
# Note that small pseudocount is added if data contains zeroes

#Found more than one class "phylo" in cache; using the first, from namespace 'phyloseq'
#Also defined by ‘tidytree’

png("Mi_plot_transform_clr.png")
plot_bar(MiSeq_18S_trans_clr)
dev.off()

MiSeq_18S_clr <- decostand(MiSeq_18S, method = "clr")
#  no method for coercing this S4 class to a vector

mi_norm_ait_dist <- vegdist(MiSeq_18S_norm, method = "aitchison")
mi_norm_ait_dist <- vegdist(MiSeq_18S_trans_clr, method = "aitchison", pseudocount > 0)


## try pca + noramlisation w/ microViz package


## whats pca? 
# see any clustering or other patterns of microbiota (dis)similarity in (many) samples
ord_explore()
# too use shiny app


Mi_clr_mViz <- MiSeq_18S %>% tax_transform(trans = "clr", rank = "Family")


MiSeq_18S %>% tax_fix()
# some ASVs without any taxa; not even Eukaryotes

tax_fix_interactive(MiSeq_18S)

MiSeq_18S_update <- MiSeq_18S %>%
tax_fix(
  min_length = 4,
  unknowns = c(""),
  sep = " ", anon_unique = FALSE,
  suffix_rank = "current"
)









## tree ? 
## maybe do trees at the end to place the rares somewhere

#phy_tree(MiSeq_18S)

#myTaxa = names(sort(taxa_sums(MiSeq_18S), decreasing = TRUE)[1:10])
#ex1 = prune_taxa(myTaxa, MiSeq_18S)
#plot(phy_tree(ex1), show.node.label = TRUE)

#plot_tree(ex1, color = "SampleType", label.tips = "Phylum", ladderize = "left", justify = "left" , size = "Abundance")

Mi_tree = rtree(ntaxa(MiSeq_18S), rooted=TRUE, tip.label=taxa_names(MiSeq_18S))

png("mi_tree2.png", height = 20000)
plot(Mi_tree)
dev.off()
# can nearly read all the ASV names 


MiSeq_18S <- merge_phyloseq(MiSeq_18S, Mi_tree)
#Found more than one class "phylo" in cache; using the first, from namespace 'phyloseq'
#Also defined by ‘tidytree’
MiSeq_18S
# now w tree


head(phy_tree(MiSeq_18S)$node.label, 10)

png("try_tree.png", height = 2000)
plot_tree(MiSeq_18S)
dev.off()

png("tree_colour9_w_species.png", height = 24000, width = 1500)
plot_tree(MiSeq_18S, color = "Phylum", label.tips = "Species", ladderize = "left" , size = "Abundance", text.size = 3) #shape = "stat")
dev.off()


png("tree_colour10_min_abund.png", height = 24000, width = 1500)
plot_tree(MiSeq_18S, color = "Phylum", label.tips = "Species", ladderize = "left" , size = "Abundance", text.size = 3, base.spacing = 0)
dev.off()
# sizebase doesnt change anything ? 
# min.abundance = show numbr of total ASVs on node
# base.spacing = 0, only 1 dot




