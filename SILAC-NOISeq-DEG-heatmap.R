#### Generating a heatmap of SILAC protein DEG's

## Install and load packages
#install.packages("pheatmap")
#install.packages("RColorBrewer")
library(pheatmap)
library(RColorBrewer)

## Set working directory
setwd("~/Desktop/Analysis") 

## Setting up colourpalette for pheatmap
my_palette <- colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100) # RdYlBu , PiYG, Spectral, PRGn
col_breaks = c(seq(-4,-0.3,length=35),
                seq(-0.29,0.29,length=30),
                seq(0.3,4,length=35))

## Load dataset and hitlist, and merge to subset silac data for genes in the hitlist. Also ajust rownames.  
dataset <- read.table("~/Desktop/Analysis/silacheatmap24h.txt", header=TRUE, sep="\t", comment="", quote="", stringsAsFactors=FALSE)
hits24h <- read.table("~/Desktop/Analysis/hits-24h.txt", header=FALSE, sep="\t", comment="", quote="", stringsAsFactors=FALSE)
colnames(hits24h) <- c("Symbol")
combined <- merge(hits24h, dataset, by="Symbol")
rownames(combined) <- combined[,1]
combined <- combined[,c(-1)]

## Plot pheatmap
pheatmap(combined,main ="Silac analysis hits combined - 24h", fontsize_row = 1, cellwidth = 50, color = my_palette, breaks=col_breaks, filename="pm-Silac-hits-clustered-24h.pdf")