### phyloseq
## create a metadata file with the same sample names as in the asv table

setwd("C:/Users/phili/OneDrive/Desktop/Master_Thesis/Sequencing/Phyloseq analysis")
# copy the meta data file + asv table in that folder

meta <- read.csv("HE627_CTD_only_first_casts.csv", header = T, dec = ".")
asvs <- read.csv("assigned_all_18S_MiSeq_analysis.csv")

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
meta$stat <- as.numeric(meta$stat)
meta <- na.omit(meta)

## create means of the meta data
## 3-40m depth

mean(meta$temp)
mean(meta$temp, trim = meta$depth(3-40))
# doesnt work

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

xy <-meta %>% group_by(stat) %>% filter(between(depth, 3, 40)) %>% summarise(temp=mean(temp),
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










