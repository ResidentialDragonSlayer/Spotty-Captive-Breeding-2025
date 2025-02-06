# visualize admixture and evalAdmix output from clouded leopards 
# Should be good to use with any sample number
# Uncomment the ggsave line in for loop the first time

# "Sakaar is the collection point for all lost and unloved things.
# Like you."

setwd("C:/Users/hrwil/OneDrive - George Mason University - O365 Production/Studies/CloudedLeopardSMSC/admixture/sakaar")
library(tidyverse)
source("visFuns.R")

sample_ids = "../../metadata/sample_ids_28_samples.csv"
high_k = 8
run = "no_sunda_take3" # change this to whatever prefix was used
evalAdmixName = "evalAdmix_no_sunda"

# load k values into a list of data frames
files <- list.files()
runs <- str_subset(files, run)
qs <- str_subset(runs, ".Q")
evals <- str_subset(files, evalAdmixName)
mods <- lapply(qs, read.table)

# change sample names to study IDs from vcf IDs
labs <- read.csv(sample_ids)
labs2 <- filter(labs[, c(1, 3)], str_sub(Sample.IDs, 1, 2) != "TH")
colnames(labs2) <- c("newID", "oldID")
inds <- read.table(paste0(run, ".list"))
colnames(inds) <- "oldID"
new_labs <- inner_join(inds, labs2)
n <- nrow(new_labs)

# visualize evalAdmix outputs
for (f in evals) {
  r <- as.matrix(read.table(f))
  plotCorRes(cor_mat = r, pop = pop, title = f, max_z=0.25, min_z=-0.25)
}

# generate all ADMIXTURE plots
plotz <- as.list(2:high_k)
cols <- rainbow(8)
for (tbl in mods) {
  i = ncol(tbl)
  indc <- rep(new_labs$newID, i)
  vals <- c(as.matrix(tbl))
  k <- as.character(rep(1:i, each=n))
  plot <- data.frame(indc, vals, k)
  plot1 <- ggplot(plot, aes(fill=k, y=vals, x=indc)) +
    geom_bar(position="fill", stat="identity", col="black", width = 0.8) +
    theme_bw() +
    labs(y = "Ancestry", x = paste0("Individual (total = ", nrow(tbl), ")")) +
    theme(axis.text.x = element_text(angle=90)) +
    scale_fill_manual(values=cols) + 
    scale_y_continuous(expand = c(0, 0), limits = c(0, 1.03))
  #ggsave(paste0(run, "_k", i, ".jpg"), plot = plot1, path = "../admixture_plots")
  plotz[[i - 1]] <- plot1
}

plotz

# this line for creating manuscript plot to match RADseq only
tbl <- mods[[2]][, c(3, 2, 1)]
