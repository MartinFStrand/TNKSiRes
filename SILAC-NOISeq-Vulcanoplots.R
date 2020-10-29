### NOISEQ DEG's analysis on SILAC dataset and generation of vulcanoplots of DEG's

### 1. Installing needed packages (if not installed)
# install.packages("BiocManager")
# BiocManager::install("NOISeq") # Tutorial pdf: http://www.bioconductor.org/packages/release/bioc/vignettes/NOISeq/inst/doc/NOISeq.pdf
# install.packages("devtools")
# devtools::install_github("stephenturner/annotables") # https://www.r-bloggers.com/annotables-r-data-package-for-annotatingconverting-gene-ids/
# install.packages("dplyr")

### 2. Loading needed packages
library(NOISeq)
library(annotables)
library(dplyr)


### 3. Generating datafiles needed for NOISeq analysis of DEGs
# 3.1 Import the human gene library grch38 from annotables and store as a table
a <- (grch38 %>% select(ensgene, entrez, symbol, chr, start, end, description, biotype)) 


# Loading dataset and merging with geneIDs from the human gene library
setwd("~/Desktop/Analysis") 
dataset <- read.table("~/Desktop/Analysis/Silacfull.txt", header=TRUE, sep="\t", comment="", quote="", stringsAsFactors=FALSE)
dataset.a <- merge(dataset, a, by="symbol")
dataset.u <- dataset.a[!duplicated(dataset.a[,1]),] # Removing symbol duplicates
rownames(dataset.u)<-dataset.u[,1]
dataset.u[is.na(dataset.u)] <- 0.03 #setting NA's to background "expression" for analysis

# 3.4 Make a chromosomal list mychroms in list format
mychroms <- dataset.u[,c(52:54)]


### COLO320DM analysis (24 hours)

# 3.5 Make proteomics counts list... 
mycounts <- dataset.u[,c(2:5)]

# 3.6 Make Factor list myfactors in list format
treatment <- c("Control","Control","Treatment","Treatment")
replicates <- c("1","2","1","2")
myfactors <- data.frame(treatment, replicates)
rownames(myfactors) <- c("COLO320DMC1","COLO320DMC2","COLO320DMT1","COLO320DMT2")

# 3.7 Make a double: mybiotypes
biotypes <- dataset.u[56]
mybiotypes <- unlist(biotypes)
names(mybiotypes) <- rownames(dataset.u)
head(mybiotypes)

# 3.9 Make a double: mylength 
lengde <- dataset.u[54]-dataset.u[53]
mylength <- unlist(lengde)
names(mylength) <- rownames(dataset.u)
head(mylength)


## 4 Generating a NOISeq object for analysis
mydata <- readData(data = mycounts, length = mylength, biotype = mybiotypes, chromosome = mychroms, factors = myfactors) # addData can be used to add more data to the dataset. 
str(mydata)

# 4.1 Investigating the dataset
head(assayData(mydata)$exprs)
head(pData(mydata))
head(featureData(mydata)@data)

# 4.2 Exploring biodetection in samples
myexplodata <- dat(mydata, type = "biodetection")
explo.plot(myexplodata, plottype = "persample")
mynicedata <- dat2save(myexplodata)
head(mynicedata)

# 4.3 Checking data before DEG analysis
mysaturation = dat(mydata, k = 0, ndepth = 7, type = "saturation")
explo.plot(mysaturation, toplot = 1, samples = 1:4, yleftlim = NULL, yrightlim = NULL)
explo.plot(mysaturation, toplot = "protein_coding", samples = 1:4)
mycountsbio = dat(mydata, factor = NULL, type = "countsbio")
explo.plot(mycountsbio, toplot = "protein_coding", samples = NULL, plottype = "boxplot")
explo.plot(mycountsbio, toplot = 1, samples = NULL, plottype = "barplot")

mylengthbias = dat(mydata, factor = "treatment", type = "lengthbias")
explo.plot(mylengthbias, samples = NULL, toplot = "global")

# 4.4 Filtering
myfilt = filtered.data(mycounts, factor = myfactors$treatment, norm = FALSE,depth = NULL, method = 1, cv.cutoff = 100, cpm = 1, p.adj = "fdr")

# 4.5 NOISeq analysis
mynoiseq1 = noiseq(mydata, k = 0.5, norm = "n", factor = "treatment", pnr = 0.2, nss = 5, v = 0.02, lc = 1, replicates = "biological") # set norm as n if normalization is prev done. (allready have rpkm values)
head(mynoiseq@results[[1]])

mynoiseq.deg = degenes(mynoiseq1, q = 0.8, M = NULL)
mynoiseq.deg1 = degenes(mynoiseq1, q = 0.8, M = "up")
mynoiseq.deg2 = degenes(mynoiseq1, q = 0.8, M = "down")

DE.plot(mynoiseq1, q = 0.8, graphic = "expr", log.scale = TRUE)
DE.plot(mynoiseq1, q = 0.8, graphic = "MD")

# 4.6 Write results to file: 
write.table(mynoiseq1@results[[1]], "~/Desktop/Analysis/COLO320DMNOISeqResults-Silac24.txt", sep="\t")


res1 <- mynoiseq1@results[[1]]
resu1 <- res1[,c(3,5)]
resu1$Symbol <- rownames(res1)
resu1 <- resu1[,c(3,1,2)]
colnames(resu1) <- c("GeneID","log2FoldChange","Probability")

pdf(file="COLO320DM-24-vulcanplota-silac.pdf", width=6, height=10)
with(resu1, plot(log2FoldChange, -log10(1-Probability), pch=20, main="Volcano plot COLO320DM TvC 24h", xlim=c(-3,3)))
with(subset(resu1, Probability>.8 ), points(log2FoldChange,-log10(1-Probability), pch=20, col="red"))
dev.off()




#### ABC-1 analysis (24 hours)

# 3.5 Make proteomics counts list... 
mycounts <- dataset.u[,c(10:13)]

# 3.6 Make Factor list myfactors in list format
treatment <- c("Control","Control","Treatment72","Treatment72")
replicates <- c("1","2","1","2")
myfactors <- data.frame(treatment, replicates)
rownames(myfactors) <- c("ABC-1C1","ABC-1C2","ABC-1T1","ABC-1T2")

### 4 Generating a NOISeq object for analysis
mydata <- readData(data = mycounts, length = mylength, biotype = mybiotypes, chromosome = mychroms, factors = myfactors) # addData can be used to add more data to the dataset. 
str(mydata)

# 4.5 NOISeq analysis
mynoiseq2 = noiseq(mydata, k = 0.5, norm = "n", factor = "treatment", pnr = 0.2, nss = 5, v = 0.02, lc = 1, replicates = "biological") # set norm as n if normalization is prev done. (allready have rpkm values)
head(mynoiseq@results[[1]])

mynoiseq.deg = degenes(mynoiseq2, q = 0.8, M = NULL)
mynoiseq.deg1 = degenes(mynoiseq2, q = 0.8, M = "up")
mynoiseq.deg2 = degenes(mynoiseq2, q = 0.8, M = "down")

DE.plot(mynoiseq2, q = 0.8, graphic = "expr", log.scale = TRUE)
DE.plot(mynoiseq2, q = 0.8, graphic = "MD")

# 4.6 Write results to file: 
write.table(mynoiseq2@results[[1]], "~/Desktop/Analysis/ABC1-Results-Silac 24 hours.txt", sep="\t")
res2 <- mynoiseq2@results[[1]]
resu2 <- res2[,c(3,5)]
resu2$Symbol <- rownames(res2)
resu2 <- resu2[,c(3,1,2)]
colnames(resu2) <- c("GeneID","log2FoldChange","Probability")

pdf(file="ABC-1-24-vulcanplota-silac.pdf", width=6, height=10)
with(resu2, plot(log2FoldChange, -log10(1-Probability), pch=20, main="Volcano plot ABC-1TvC 24h", xlim=c(-3,3)))
with(subset(resu2, Probability>.8 ), points(log2FoldChange,-log10(1-Probability), pch=20, col="red"))
dev.off()



##### UO-31 analysis (24 hours)

# 3.5 Make proteomics counts list... 
mycounts <- dataset.u[,c(18:21)]

# 3.6 Make Factor list myfactors in list format
treatment <- c("Control","Control","Treatment72","Treatment72")
replicates <- c("1","2","1","2")
myfactors <- data.frame(treatment, replicates)
rownames(myfactors) <- c("UO-31C1","UO-31C2","UO-31T1","UO-31T2")

### 4 Generating a NOISeq object for analysis
mydata <- readData(data = mycounts, length = mylength, biotype = mybiotypes, chromosome = mychroms, factors = myfactors) # addData can be used to add more data to the dataset. 
str(mydata)

# 4.5 NOISeq analysis
mynoiseq3 = noiseq(mydata, k = 0.5, norm = "n", factor = "treatment", pnr = 0.2, nss = 5, v = 0.02, lc = 1, replicates = "biological") # set norm as n if normalization is prev done. (allready have rpkm values)
head(mynoiseq@results[[1]])

mynoiseq.deg = degenes(mynoiseq3, q = 0.8, M = NULL)
mynoiseq.deg1 = degenes(mynoiseq3, q = 0.8, M = "up")
mynoiseq.deg2 = degenes(mynoiseq3, q = 0.8, M = "down")

DE.plot(mynoiseq3, q = 0.8, graphic = "expr", log.scale = TRUE)
DE.plot(mynoiseq3, q = 0.8, graphic = "MD")

# 4.6 Write results to file: 
write.table(mynoiseq3@results[[1]], "~/Desktop/Analysis/UO-31-Results-Silac 24 hours.txt", sep="\t")
res3 <- mynoiseq3@results[[1]]
resu3 <- res3[,c(3,5)]
resu3$Symbol <- rownames(res3)
resu3 <- resu3[,c(3,1,2)]
colnames(resu3) <- c("GeneID","log2FoldChange","Probability")

pdf(file="UO-31-24-vulcanplota-silac.pdf", width=6, height=10)
with(resu3, plot(log2FoldChange, -log10(1-Probability), pch=20, main="Volcano plot UO-31 TvC 24h", xlim=c(-3,3)))
with(subset(resu3, Probability>.8 ), points(log2FoldChange,-log10(1-Probability), pch=20, col="red"))
dev.off()


#### OVCAR-4 analysis (24 hours)

# 3.5 Make proteomics counts list... 
mycounts <- dataset.u[,c(26:29)]

# 3.6 Make Factor list myfactors in list format
treatment <- c("Control","Control","Treatment72","Treatment72")
replicates <- c("1","2","1","2")
myfactors <- data.frame(treatment, replicates)
rownames(myfactors) <- c("OVCAR-4C1","OVCAR-4C2","OVCAR-4T1","OVCAR-4T2")

### 4 Generating a NOISeq object for analysis
mydata <- readData(data = mycounts, length = mylength, biotype = mybiotypes, chromosome = mychroms, factors = myfactors) # addData can be used to add more data to the dataset. 
str(mydata)

# 4.5 NOISeq analysis
mynoiseq4 = noiseq(mydata, k = 0.5, norm = "n", factor = "treatment", pnr = 0.2, nss = 5, v = 0.02, lc = 1, replicates = "biological") # set norm as n if normalization is prev done. (allready have rpkm values)
head(mynoiseq@results[[1]])

mynoiseq.deg = degenes(mynoiseq4, q = 0.8, M = NULL)
mynoiseq.deg1 = degenes(mynoiseq4, q = 0.8, M = "up")
mynoiseq.deg2 = degenes(mynoiseq4, q = 0.8, M = "down")

DE.plot(mynoiseq4, q = 0.8, graphic = "expr", log.scale = TRUE)
DE.plot(mynoiseq4, q = 0.8, graphic = "MD")

# 4.6 Write results to file: 
write.table(mynoiseq4@results[[1]], "~/Desktop/Analysis/OVCAR-4-Results-Silac 24 hours.txt", sep="\t")
res4 <- mynoiseq4@results[[1]]
resu4 <- res4[,c(3,5)]
resu4$Symbol <- rownames(res4)
resu4 <- resu4[,c(3,1,2)]
colnames(resu4) <- c("GeneID","log2FoldChange","Probability")

pdf(file="OVCAR-4-24-vulcanplota-silac.pdf", width=6, height=10)
with(resu4, plot(log2FoldChange, -log10(1-Probability), pch=20, main="Volcano plot OVCAR-4 TvC 24h", xlim=c(-3,3)))
with(subset(resu4, Probability>.8 ), points(log2FoldChange,-log10(1-Probability), pch=20, col="red"))
dev.off()


#### RKO analysis (24 hours)

# 3.5 Make proteomics counts list... 
mycounts <- dataset.u[,c(34:37)]

# 3.6 Make Factor list myfactors in list format
treatment <- c("Control","Control","Treatment72","Treatment72")
replicates <- c("1","2","1","2")
myfactors <- data.frame(treatment, replicates)
rownames(myfactors) <- c("RKOC1","RKOC2","RKOT1","RKOT2")

### 4 Generating a NOISeq object for analysis
mydata <- readData(data = mycounts, length = mylength, biotype = mybiotypes, chromosome = mychroms, factors = myfactors) # addData can be used to add more data to the dataset. 
str(mydata)

# 4.5 NOISeq analysis
mynoiseq8 = noiseq(mydata, k = 0.5, norm = "n", factor = "treatment", pnr = 0.2, nss = 5, v = 0.02, lc = 1, replicates = "biological") # set norm as n if normalization is prev done. (allready have rpkm values)
head(mynoiseq@results[[1]])

mynoiseq.deg = degenes(mynoiseq8, q = 0.8, M = NULL)
mynoiseq.deg1 = degenes(mynoiseq8, q = 0.8, M = "up")
mynoiseq.deg2 = degenes(mynoiseq8, q = 0.8, M = "down")

DE.plot(mynoiseq8, q = 0.8, graphic = "expr", log.scale = TRUE)
DE.plot(mynoiseq8, q = 0.8, graphic = "MD")

# 4.6 Write results to file: 
write.table(mynoiseq8@results[[1]], "~/Desktop/Analysis/RKO-Results-Silac 24 hours.txt", sep="\t")
res8 <- mynoiseq8@results[[1]]
resu8 <- res8[,c(3,5)]
resu8$Symbol <- rownames(res8)
resu8 <- resu8[,c(3,1,2)]
colnames(resu8) <- c("GeneID","log2FoldChange","Probability")

pdf(file="RKO-24h vulcanplota-silac.pdf", width=6, height=10)
with(resu8, plot(log2FoldChange, -log10(1-Probability), pch=20, main="Volcano plot RKO TvC 24h", xlim=c(-3,2)))
with(subset(resu8, Probability>.8 ), points(log2FoldChange,-log10(1-Probability), pch=20, col="red"))
dev.off()




## Making a combined vulcanoplot for all cell lines
setEPS()
postscript(file="CombinedVulcanoplots.eps", width=35, height=15)
par(mfrow= c(1,5))

with(resu1, plot(log2FoldChange, -log10(1-Probability), pch=20, main="COLO320DM TvC 24h", ylim=c(0.0, 3.5), xlim=c(-2.5,2)))
with(subset(resu1, Probability>.8 ), points(log2FoldChange,-log10(1-Probability), pch=20, col="red"))

with(resu3, plot(log2FoldChange, -log10(1-Probability), pch=20, main="UO-31 TvC 24h", ylim=c(0.0, 3.5), xlim=c(-2.5,2)))
with(subset(resu3, Probability>.8 ), points(log2FoldChange,-log10(1-Probability), pch=20, col="red"))

with(resu4, plot(log2FoldChange, -log10(1-Probability), pch=20, main="OVCAR-4 TvC 24h", ylim=c(0.0, 3.5), xlim=c(-2.5,2)))
with(subset(resu4, Probability>.8 ), points(log2FoldChange,-log10(1-Probability), pch=20, col="red"))

with(resu2, plot(log2FoldChange, -log10(1-Probability), pch=20, main="ABC-1TvC 24h", ylim=c(0.0, 3.5), xlim=c(-2.5,2)))
with(subset(resu2, Probability>.8 ), points(log2FoldChange,-log10(1-Probability), pch=20, col="red"))

with(resu8, plot(log2FoldChange, -log10(1-Probability), pch=20, main="RKO TvC 24h", ylim=c(0.0, 3.5), xlim=c(-2.5,2)))
with(subset(resu8, Probability>.8 ), points(log2FoldChange,-log10(1-Probability), pch=20, col="red"))


dev.off()



