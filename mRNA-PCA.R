#### PCA analysis of mRNASeq expression data jof 5 cell lines

## Set working directory
setwd("~/Desktop/PCA")

## Install (if required) and load packages
# install.packages("rgl")
# install.packages("ggplot2")
library(rgl)
library(ggplot2)

## Load mRNASeq dataset
dataset <- read.table("~/Desktop/PCA/log4v1.txt", header=TRUE, sep="\t", comment="", quote="", stringsAsFactors=FALSE)

## Arrange data and add rownames/colnames and remove NA's
x <- data.matrix(dataset[,2:21])
rownames(x) <- dataset[,1]
colnames(x) <- c("COLO320DM-C1","COLO320DM-C2","COLO320DM-T1","COLO320DM-T2",
                 "ABC-1-C1","ABC-1-C2","ABC-1-T1","ABC-1-T2",
                 "UO-31-C1","UO-31-C2","UO-31-T1","UO-31-T2",
                 "OVCAR-4-C1","OVCAR-4-C2","OVCAR-4-T1","OVCAR-4-T2",
                 "RKO-C1","RKO-C2","RKO-T1","RKO-T2")
x <- na.omit(log2(x))

## PCA-calulation
xx <- prcomp(t(x))
pca_data=prcomp(t(x))
pca_data_perc=round(100*pca_data$sdev^2/sum(pca_data$sdev^2),1)

## Set factor for plotting PCA, and set up df
condition = rep(c("Control","Control","G007-LK","G007-LK"), times=5)
samples = rep(c("COLO320DM","ABC-1","UO-31","OVCAR-4","RKO"), each=4)
df_pca_data=data.frame(PC1 = pca_data$x[,1], PC2 = pca_data$x[,2], sample = samples, condition=condition)

## Print out eps file
setEPS()
postscript(file="PCA-2D.eps", width=6, height=4)
ggplot(df_pca_data, aes(PC1,PC2, color = sample))+
  geom_point(size=3, aes(shape=condition))+
  labs(x=paste0("PC1 (",pca_data_perc[1],")"), y=paste0("PC2 (",pca_data_perc[2],")"))
dev.off()

