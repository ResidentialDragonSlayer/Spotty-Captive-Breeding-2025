##PCA for 28 clouded leopards and 2 tigers
##no UK11, UK13, US24, or US25
##so n=25
##reseq analysis for SMSC

library(tidyverse)
setwd("C:/Users/hrwil/OneDrive - George Mason University - O365 Production/Studies/CloudedLeopardSMSC/pca/")
run <- "pca_take3"
dropped <- c("UK11", "UK13", "US24", "US25")

##Reading in eigen files

pca<-read_table(paste0(run, ".eigenvec"), col_names = FALSE)

eigenval<-scan(paste0(run, ".eigenval"))

##Cleaning up table

pca <- pca[,-1]
colnames(pca)[1] <- "ID in WGS VCF"

sample_ids <- read_csv("../metadata/sample_ids_30_samples.csv")

pca <- sample_ids %>%
  select(`This Study ID`, `ID in WGS VCF`) %>%
  right_join(pca) %>%
  filter(!(`This Study ID` %in% dropped))
colnames(pca)[1] <- "ID"

# get population info

pop <- str_sub(pca$ID, 1, 2)
pca$Population <- pop
pca <- pca %>% filter(!(Population == "SU"))

##PLOT PCA
p=ggplot(pca, aes(x=X3,y=X4, color=Population, label=ID))  +
  xlab("PC 1") + ylab("PC 2") +
  geom_point(size=4) +
  scale_color_manual(breaks = c("US", "UK","CH"),
                     values=c("#C77CFF","#00BFC4","#F8766D")) +
  theme_bw()+ theme(axis.text=element_text(size=10), legend.text = element_text(size=10),
                    legend.title=element_text(size=10))
library(ggrepel)
##FIX OVERLAPPING LABELS
p+ geom_label_repel(max.overlaps = 50, show.legend = FALSE, size=4 ) + 
  theme(axis.text=element_text(size=10),
        axis.title=element_text(size=10))


