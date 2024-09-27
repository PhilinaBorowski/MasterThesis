### phyloseq - NextSeq prep
## create a metadata file with the same sample names as in the asv table

setwd("/isibhv/projects/AG_John/Expeditions/HE627_Data/Comparison_all_seq_runs/Phyloseq")
# copy the meta data file + asv table in that folder

## used the "reduced_no_1rowsum" files + NextSeq 16S for ASVs tables

meta <- read.csv("HE627_CTD_only_first_casts.csv", header = T, dec = ".")


## 
library(dplyr)
##

## rename
meta <- meta %>% rename(
  event = Event,
  cruise = Cruise,
  stat = Station,
  cast = Cast,
  stat.cast = Station.Cast,
  lat = Latitude,
  long = Longitude,
  depth = Depth.water..m.,
  pres = Press..dbar.,
  temp = Temp...C.,
  cond = Conductivity..mS.cm.,
  sal = Sal,
  tpot = Tpot...C.,
  dens = Density..kg.m..3.,
  oxyg = O2..µmol.l.,
  oxyg_sat = O2.sat....,
  atten = Attenuation....,
  chla = Chl.a..µg.l.,
  nobs = NOBS....
)
str(meta)
meta$stat <- as.factor(meta$stat)
meta <- na.omit(meta)

## create means of the meta data
## 3-40m depth

mean(meta$temp)

meta %>% filter(depth >= 3) %>% filter(depth <= 40) %>% summarise(temp=mean(temp),
                                                                  sal=mean(sal),
                                                                  pres=mean(pres),
                                                                  cond=mean(cond),
                                                                  tpot=mean(tpot),
                                                                  dens=mean(dens),
                                                                  oxyg=mean(oxyg),
                                                                  oxyg_sat=mean(oxyg_sat),
                                                                  atten=mean(atten),
                                                                  chla=mean(chla),
                                                                  nobs=mean(nobs))

xy <- meta %>% group_by(stat) %>% filter(between(depth, 3, 40)) %>% summarise(temp=mean(temp),
                                                                             sal=mean(sal),
                                                                             pres=mean(pres),
                                                                             cond=mean(cond),
                                                                             tpot=mean(tpot),
                                                                             dens=mean(dens),
                                                                             oxyg=mean(oxyg),
                                                                             oxyg_sat=mean(oxyg_sat),
                                                                             atten=mean(atten),
                                                                             chla=mean(chla),
                                                                             nobs=mean(nobs),
                                                                             lat=mean(lat),
                                                                             long=mean(long))

## copied that in a new excel 
## incl sequencer, DNA, depth with sample names, method of sample taking
## meta_phyloseq_18S_miseq.csv

write.csv(xy, "metadata_upper_layers.csv")

# for deep DNA samples: copied metadata from deepest depth 








