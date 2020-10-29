### Vulcanoplot from NOISeq DEG his across 5 cell lines
## 0) Install/load required packages
## Repeat step 1:3 for each of the cell lines
## 1) Load dataset into dataframe and prepare data for plotting
## 2) Edit datalabels to GENEID
## 3) Generate vulcanplot
## 4) Combine all vulcanoplots in a esi file

# Set working directory:
setwd("~/Desktop/Analysis")

## 0)
# source("https://bioconductor.org/biocLite.R")
# biocLite("org.Hs.eg.db")
library(org.Hs.eg.db)


## 1) COLO320DM vulcanoplot
# Load txt file:
res1 <- read.table("~/Desktop/Analysis/1.txt", header=TRUE, sep="\t", comment="", quote="", stringsAsFactors=FALSE)
head(res1)
colnames(res1) <- c("GeneID","log2FoldChange","Probability")
## 2)
# Make a vector containing the entrez genecodes, and use it to get out a vector of geneIDs and replace it in the datafile res
a1 <- res1[,1]
a.symbol1 <- as.vector(unlist(mget(a1, envir=org.Hs.egSYMBOL, ifnotfound=NA)))
res1[,1] <- a.symbol1
## 3)
## Generate pdf file of resulting plot
pdf(file="COLO320DM-vulcanplota.pdf", width=6, height=10)
with(res1, plot(log2FoldChange, -log10(1-Probability), pch=20, main="Volcano plot COLO320DM TvC", xlim=c(-5,5)))
with(subset(res1, Probability>.8 ), points(log2FoldChange,-log10(1-Probability), pch=20, col="red"))
dev.off()


## 1) ABC-1 vulcano plot
# Load txt file:
res2 <- read.table("~/Desktop/Analysis/2.txt", header=TRUE, sep="\t", comment="", quote="", stringsAsFactors=FALSE)
head(res2)
colnames(res2) <- c("GeneID","log2FoldChange","Probability")
## 2)
# Make a vector containing the entrez genecodes, and use it to get out a vector of geneIDs and replace it in the datafile res
a2 <- res2[,1]
a.symbol2 <- as.vector(unlist(mget(a2, envir=org.Hs.egSYMBOL, ifnotfound=NA)))
res2[,1] <- a.symbol2
## 3)
## Generate pdf file of resulting plot
pdf(file="ABC-1-vulcanplot.pdf", width=6, height=10)
with(res2, plot(log2FoldChange, -log10(1-Probability), pch=20, main="Volcano plot ABC-1 TvC", xlim=c(-5,5)))
with(subset(res2, Probability>.8 ), points(log2FoldChange, -log10(1-Probability), pch=20, col="red"))
dev.off()


## 1) UO-31 vulcano plot 
# Load txt file:
res3 <- read.table("~/Desktop/Analysis/3.txt", header=TRUE, sep="\t", comment="", quote="", stringsAsFactors=FALSE)
head(res3)
colnames(res3) <- c("GeneID","log2FoldChange","Probability")
## 2)
# Make a vecotor containing the entrez genecodes, and use it to get out a vector of geneIDs and replace it in the datafile res
a3 <- res3[,1]
a.symbol3 <- as.vector(unlist(mget(a3, envir=org.Hs.egSYMBOL, ifnotfound=NA)))
res3[,1] <- a.symbol3
## 3)
## Generate pdf file of resulting plot
pdf(file="UO-31-vulcanplot.pdf", width=6, height=10)
with(res3, plot(log2FoldChange, -log10(1-Probability), pch=20, main="Volcano plot UO-31 TvC", xlim=c(-5,5)))
with(subset(res3, Probability>.8 ), points(log2FoldChange, -log10(1-Probability), pch=20, col="red"))
dev.off()


## 1) OVCAR4 vulcano plot
# Load txt file:
res4 <- read.table("~/Desktop/Analysis/4.txt", header=TRUE, sep="\t", comment="", quote="", stringsAsFactors=FALSE)
head(res4)
colnames(res4) <- c("GeneID","log2FoldChange","Probability")
## 2)
# Make a vecotor containing the entrez genecodes, and use it to get out a vector of geneIDs and replace it in the datafile res
a4 <- res4[,1]
a.symbol4 <- as.vector(unlist(mget(a4, envir=org.Hs.egSYMBOL, ifnotfound=NA)))
res4[,1] <- a.symbol4
## 3)
## Generate pdf file of resulting plot
pdf(file="OVCAR4-vulcanplot.pdf", width=6, height=10)
with(res4, plot(log2FoldChange, -log10(1-Probability), pch=20, main="Volcano plot OVCAR4 TvC", xlim=c(-5,5)))
with(subset(res4, Probability>.8 ), points(log2FoldChange, -log10(1-Probability), pch=20, col="red"))
dev.off()



## 1) RKO vulcano plot
# Load txt file:
res8 <- read.table("~/Desktop/Analysis/8.txt", header=TRUE, sep="\t", comment="", quote="", stringsAsFactors=FALSE)
head(res8)
colnames(res8) <- c("GeneID","log2FoldChange","Probability")
## 2)
# Make a vecotor containing the entrez genecodes, and use it to get out a vector of geneIDs and replace it in the datafile res
a8 <- res8[,1]
a.symbol8 <- as.vector(unlist(mget(a8, envir=org.Hs.egSYMBOL, ifnotfound=NA)))
res8[,1] <- a.symbol8
## 3)
## Generate pdf file of resulting plot
pdf(file="RKO-vulcanplot.pdf", width=6, height=10)
with(res8, plot(log2FoldChange, -log10(1-Probability), pch=20, main="Volcano plot RKO TvC", xlim=c(-5,5)))
with(subset(res8, Probability>.8 ), points(log2FoldChange, -log10(1-Probability), pch=20, col="red"))
dev.off()

##4) Generating a combined plot of all the plots in eps format, using postscript
setEPS()
postscript(file="CombinedVulcanoplots.eps", width=30, height=10)
par(mfrow= c(1,5))
with(res1, plot(log2FoldChange, -log10(1-Probability), pch=20, main="COLO320DM", ylim=c(0.1,1.5), xlim=c(-5,5), cex.lab=2, cex.axis=2, cex.main=3))
with(subset(res1, Probability>.8 ), points(log2FoldChange,-log10(1-Probability), pch=20, col="red"))

with(res3, plot(log2FoldChange, -log10(1-Probability), pch=20, main="UO-31", ylim=c(0.1,1.5),xlim=c(-5,5), cex.lab=2, cex.axis=2, cex.main=3))
with(subset(res3, Probability>.8 ), points(log2FoldChange, -log10(1-Probability), pch=20, col="red"))

with(res4, plot(log2FoldChange, -log10(1-Probability), pch=20, main="OVCAR-4", ylim=c(0.1,1.5),xlim=c(-5,5), cex.lab=2, cex.axis=2, cex.main=3))
with(subset(res4, Probability>.8 ), points(log2FoldChange, -log10(1-Probability), pch=20, col="red"))

with(res2, plot(log2FoldChange, -log10(1-Probability), pch=20, main="ABC-1", ylim=c(0.1,1.5),xlim=c(-5,5), cex.lab=2, cex.axis=2, cex.main=3))
with(subset(res2, Probability>.8 ), points(log2FoldChange, -log10(1-Probability), pch=20, col="red"))

with(res8, plot(log2FoldChange, -log10(1-Probability), pch=20, main="RKO", ylim=c(0.1,1.5),xlim=c(-5,5), cex.lab=2, cex.axis=2, cex.main=3))
with(subset(res8, Probability>.8 ), points(log2FoldChange, -log10(1-Probability), pch=20, col="red"))

dev.off()



