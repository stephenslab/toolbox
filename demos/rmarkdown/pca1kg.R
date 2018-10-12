# I estimate that about 12.5 GB of memory will be needed to run this
# analysis.
#
# Download the 1kg.rds file from Dropbox from this URL:
#
#   https://www.dropbox.com/s/e6nddjkd73gpg93/1kg.rds?dl=0
#
library(rsvd)
library(ggplot2)
library(cowplot)
source("../code/geno.utils.R")

# Load the 1000 Genomes genotype data.
cat("Loading data.\n")
geno <- readRDS("../data/1kg.rds")

# load the 1000 Genomes population labels.
labels.1kg <- read.table("../data/omni_samples.20141118.panel",
                         sep = " ",header = TRUE,as.is = "id")

# Use the rpca function to compute the first 10 PCs---that is, the 10
# components that explain the most variation in the genotypes.
cat("Computing PCs.\n")

start = proc.time()
out.pca       <- rpca(geno,k = 10,center = TRUE,scale = FALSE,retx = TRUE)
run_time = proc.time() - start
print(paste0("PCA finished after: ", run_time))

pcs           <- out.pca$x
colnames(pcs) <- paste0("PC",1:10)

# Let's take a quick look at the PCA results.
summary(out.pca)

# This adds a "label" column to the `pcs` table.
pcs <- add.poplabels(pcs,labels.1kg)

# Create a PC plot with population labels.
p <- labeled.pc.plot(pcs,x = "PC1",y = "PC2",
                     label = "label",size = 2)
print(p)
