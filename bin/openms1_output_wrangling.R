### openms to tidy csv

library(tidyverse)
library(data.table)
library(here)

##1. 60 min drug data
#setwd("~/iclouddrive/Documents/Manuscripts/figeys/clusterMethods/manuscript/MetaProClust-MS1/")

features60mins <-  fread(here::here("data", "60mins_ms1", "60minsdrug-ms1features_noheader.txt"),
                         quote="",sep='\t')
ms1_openms_samples <- read.csv(here::here("data", "60mins_ms1", "sampleinfo.csv"), header=F)
ms1_openms_metadata <- read.csv(here::here("data", "meta_data.csv")) %>%
    mutate(plate_cell = paste0("P2", Cell))

consensus_features <- features60mins %>% dplyr::select(.,contains("cf"))
ms1 <- features60mins %>% dplyr::select(.,-contains("cf"))
rt <- ms1 %>% dplyr::select(.,contains("rt"))
mz <- ms1 %>% dplyr::select(.,contains("mz"))
charge <-  ms1 %>% dplyr::select(.,contains("charge"))
intensity <- ms1 %>% dplyr::select(., contains("intensity")) %>% as.data.frame()
colnames(intensity) <- ms1_openms_samples$V3

#write.csv(intensity, file=here::here("data", "60mins_ms1","ms1_60min_intensity.csv"))


##2. IBD MS1-only data
#header <- read.table(here::here("ibd", "data", "ibd-nc-ms1featuresNOHEADER.txt"), header = TRUE, nrow = 1)
#featuresibd <-  read.csv(here::here("ibd", "data", "ibd-nc_ms1features_20mins-forinput.csv"),
#                      sep=',')
featuresibd <-  read.csv(here::here("ibd", "data", "ibd-nc_ms1features_20mins-120ppm-forinput.csv"),
                         sep=',')
#featuresibd <-  read.csv(here::here("ibd", "data", "ibd-nc_ms1features_forinput.csv"),
 #                        sep=',')
ibd_samples <- read.csv(here::here("ibd", "data", "ibd-nc-samples.csv"), header=F)
consensus_features <- featuresibd %>% dplyr::select(.,contains("cf")) %>% 
    mutate(feat_num = 1:nrow(featuresibd))


## trypsin?
consensus_features %>% filter(between(mz_cf, 576.0, 577.99)) %>% arrange(desc(intensity_cf))
consensus_features %>% filter(between(mz_cf, 1081.0, 1082.99)) %>% arrange(desc(intensity_cf))


ms1 <- featuresibd %>% dplyr::select(.,-contains("cf"))
rt <- ms1 %>% dplyr::select(.,contains("rt"))
mz <- ms1 %>% dplyr::select(.,contains("mz"))
charge <-  ms1 %>% dplyr::select(.,contains("charge"))
intensity <- ms1 %>% dplyr::select(., contains("intensity")) %>% as.data.frame()
colnames(intensity) <- trimws(ibd_samples$V6, which="both")

intensity_withcf <- data.frame(feature_num = consensus_features$feat_num, RT = consensus_features$rt_cf, 
                               MZ = consensus_features$mz_cf, intensity)

#write.csv(intensity_withcf, file=here::here("ibd", "data", "ibd-nc_ms1intensity-20mins.csv"))
#write.csv(intensity_withcf, file=here::here("ibd", "data", "ibd-nc_ms1intensity-20mins-100ppm.csv"))
write.csv(intensity_withcf, file=here::here("ibd", "data", "ibd-nc_ms1intensity-20mins-120ppm.csv"))


# checking if IBD peptides.txt gives similar distribution
