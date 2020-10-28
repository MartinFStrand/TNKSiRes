#### Waterfall plot of mutation data from CCLE, COSMIC and CANSAR:
# Method references:  
# https://genviz.org/module-03-genvisr/0003/02/01/waterfall_GenVisR/
# https://bioconductor.org/packages/release/bioc/vignettes/GenVisR/inst/doc/waterfall_introduction.html

### Install and load packages
# library("devtools")
# devtools::install_github("griffithlab/GenVisR")
# install.packages("dplyr")
library(GenVisR)
library(dplyr)

### Set working directory
setwd("~/Desktop/mutAnalysis")

### Load driver mutation dataset (from this paper: https://pubmed.ncbi.nlm.nih.gov/29625053/)
DRIVER <- read.table("~/Desktop/mutAnalysis/DRIVER2.txt", header=TRUE, sep="\t", comment="", quote="", stringsAsFactors=FALSE)

### Load in mutation datasets for each cell line and prepare data for waterfall-plot. 
## COLO320DM:
# Load mRNASeq-mutation call data and database mutation data: 
coloCOSMIC <- read.delim("~/Desktop/mutAnalysis/1-COSMIC.csv", sep = ",")
coloCCLE <- read.table("~/Desktop/mutAnalysis/1-CCLE.txt", header=TRUE, sep="\t", comment="", quote="", stringsAsFactors=FALSE)
coloCANSAR <- read.table("~/Desktop/mutAnalysis/1-CANSAR.txt", header=FALSE, sep="\t", comment="", quote="", stringsAsFactors=FALSE)

# Change categories for CCLEdataset to match COSMIC and CANSAR
coloCCLE <- coloCCLE %>% mutate(Variant.Classification = if_else(Variant.Classification == 'Missense_Mutation', 'Substitution - Missense', Variant.Classification))
coloCCLE <- coloCCLE %>% mutate(Variant.Classification = if_else(Variant.Classification == 'Silent', 'Substitution - coding silent', Variant.Classification))
coloCCLE <- coloCCLE %>% mutate(Variant.Classification = if_else(Variant.Classification == 'Nonsense_Mutation', 'Substitution - Nonsense', Variant.Classification))
coloCCLE <- coloCCLE %>% mutate(Variant.Classification = if_else(Variant.Classification == 'Frame_Shift_Ins', 'Insertion - Frameshift', Variant.Classification))
coloCCLE <- coloCCLE %>% mutate(Variant.Classification = if_else(Variant.Classification == 'In_Frame_Del', 'Deletion - In frame', Variant.Classification))
coloCCLE <- coloCCLE %>% mutate(Variant.Classification = if_else(Variant.Classification == 'Frame_Shift_Del', 'Deletion - Frameshift', Variant.Classification))

# Filter out mutations in genes that are identified as drivers for the 4 mutation datasets: 
coloCCLE_sel <- merge(DRIVER, coloCCLE, by.x="Gene", by.y="Hugo.Symbol")
coloCCLE_sel <- select(coloCCLE_sel, Gene,Variant.Classification,Protein.Change)
sel <- which(coloCCLE_sel$Variant.Classification== "Substitution - coding silent")
coloCCLE_sel <- coloCCLE_sel[-sel,]

coloCOSMIC_sel <- merge(DRIVER, coloCOSMIC, by="Gene")
coloCOSMIC_sel <- select(coloCOSMIC_sel, Gene, AA.Mutation, Type)
sel <- which(coloCOSMIC_sel$Type== "Substitution - coding silent")
coloCOSMIC_sel <- coloCOSMIC_sel[-sel,]
sel <- which(coloCOSMIC_sel$Type== "Unknown")
coloCOSMIC_sel <- coloCOSMIC_sel[-sel,]

coloCANSAR_sel <- merge(DRIVER, coloCANSAR, by.x="Gene", by.y="V2")
coloCANSAR_sel <- select(coloCANSAR_sel, Gene, V5, V7)
sel <- which(coloCANSAR_sel$V7== "Substitution - coding silent")

# Add cell-line and source to filtered datasets:
coloCCLE_sel$cell.line <- c("COLO320DM")
coloCCLE_sel <- coloCCLE_sel[,c(4,1,2,3)]
colnames(coloCCLE_sel) <- c("sample", "gene", "variant_class", "amino.acid.change")

coloCOSMIC_sel$cell.line <- c("COLO320DM")
coloCOSMIC_sel <- coloCOSMIC_sel[,c(4,1,3,2)]
colnames(coloCOSMIC_sel) <- c("sample", "gene", "variant_class", "amino.acid.change")

coloCANSAR_sel$cell.line <- c("COLO320DM")
coloCANSAR_sel <- coloCANSAR_sel[,c(4,1,3,2)]
colnames(coloCANSAR_sel) <- c("sample", "gene", "variant_class", "amino.acid.change")


## RKO:
rkoCOSMIC <- read.delim("~/Desktop/mutAnalysis/8-COSMIC.csv", sep = ",")
rkoCCLE <- read.table("~/Desktop/mutAnalysis/8-CCLE.txt", header=TRUE, sep="\t", comment="", quote="", stringsAsFactors=FALSE)
rkoCANSAR <- read.table("~/Desktop/mutAnalysis/8-CANSAR.txt", header=FALSE, sep="\t", comment="", quote="", stringsAsFactors=FALSE)

# Change categories for CCLEdataset to match COSMIC and CANSAR
rkoCCLE <- rkoCCLE %>% mutate(Variant.Classification = if_else(Variant.Classification == 'Missense_Mutation', 'Substitution - Missense', Variant.Classification))
rkoCCLE <- rkoCCLE %>% mutate(Variant.Classification = if_else(Variant.Classification == 'Silent', 'Substitution - coding silent', Variant.Classification))
rkoCCLE <- rkoCCLE %>% mutate(Variant.Classification = if_else(Variant.Classification == 'Nonsense_Mutation', 'Substitution - Nonsense', Variant.Classification))
rkoCCLE <- rkoCCLE %>% mutate(Variant.Classification = if_else(Variant.Classification == 'Frame_Shift_Ins', 'Insertion - Frameshift', Variant.Classification))
rkoCCLE <- rkoCCLE %>% mutate(Variant.Classification = if_else(Variant.Classification == 'In_Frame_Del', 'Deletion - In frame', Variant.Classification))
rkoCCLE <- rkoCCLE %>% mutate(Variant.Classification = if_else(Variant.Classification == 'Frame_Shift_Del', 'Deletion - Frameshift', Variant.Classification))

# Filter out mutations in genes that are identified as drivers for the 4 mutation datasets: 
rkoCCLE_sel <- merge(DRIVER, rkoCCLE, by.x="Gene", by.y="Hugo.Symbol")
rkoCCLE_sel <- select(rkoCCLE_sel, Gene,Variant.Classification,Protein.Change)
sel <- which(rkoCCLE_sel$Variant.Classification== "Substitution - coding silent")
rkoCCLE_sel <- rkoCCLE_sel[-sel,]

rkoCOSMIC_sel <- merge(DRIVER, rkoCOSMIC, by="Gene")
rkoCOSMIC_sel <- select(rkoCOSMIC_sel, Gene, AA.Mutation, Type)
sel <- which(rkoCOSMIC_sel$Type== "Substitution - coding silent")
rkoCOSMIC_sel <- rkoCOSMIC_sel[-sel,]
sel <- which(rkoCOSMIC_sel$Type== "Unknown")
rkoCOSMIC_sel <- rkoCOSMIC_sel[-sel,]

rkoCANSAR_sel <- merge(DRIVER, rkoCANSAR, by.x="Gene", by.y="V1")
rkoCANSAR_sel <- select(rkoCANSAR_sel, Gene, V4, V6)
sel <- which(rkoCANSAR_sel$V6== "Substitution - coding silent")
rkoCANSAR_sel <- rkoCANSAR_sel[-sel,]
sel <- which(rkoCANSAR_sel$V6== "Unknown")
rkoCANSAR_sel <- rkoCANSAR_sel[-sel,]

# Add cell-line and source to filtered datasets:
rkoCCLE_sel$cell.line <- c("RKO")
rkoCCLE_sel <- rkoCCLE_sel[,c(4,1,2,3)]
colnames(rkoCCLE_sel) <- c("sample", "gene", "variant_class", "amino.acid.change")

rkoCOSMIC_sel$cell.line <- c("RKO")
rkoCOSMIC_sel <- rkoCOSMIC_sel[,c(4,1,3,2)]
colnames(rkoCOSMIC_sel) <- c("sample", "gene", "variant_class", "amino.acid.change")

rkoCANSAR_sel$cell.line <- c("RKO")
rkoCANSAR_sel <- rkoCANSAR_sel[,c(4,1,3,2)]
colnames(rkoCANSAR_sel) <- c("sample", "gene", "variant_class", "amino.acid.change")


## ABC1:
abcCOSMIC <- read.delim("~/Desktop/mutAnalysis/2-COSMIC.csv", sep = ",")
abcCCLE <- read.table("~/Desktop/mutAnalysis/2-CCLE.txt", header=TRUE, sep="\t", comment="", quote="", stringsAsFactors=FALSE)
abcCANSAR <- read.table("~/Desktop/mutAnalysis/2-CANSAR.txt", header=FALSE, sep="\t", comment="", quote="", stringsAsFactors=FALSE)

# Change categories for CCLEdataset to match COSMIC and CANSAR
abcCCLE <- abcCCLE %>% mutate(Variant.Classification = if_else(Variant.Classification == 'Missense_Mutation', 'Substitution - Missense', Variant.Classification))
abcCCLE <- abcCCLE %>% mutate(Variant.Classification = if_else(Variant.Classification == 'Silent', 'Substitution - coding silent', Variant.Classification))
abcCCLE <- abcCCLE %>% mutate(Variant.Classification = if_else(Variant.Classification == 'Nonsense_Mutation', 'Substitution - Nonsense', Variant.Classification))
abcCCLE <- abcCCLE %>% mutate(Variant.Classification = if_else(Variant.Classification == 'Frame_Shift_Ins', 'Insertion - Frameshift', Variant.Classification))
abcCCLE <- abcCCLE %>% mutate(Variant.Classification = if_else(Variant.Classification == 'In_Frame_Del', 'Deletion - In frame', Variant.Classification))
abcCCLE <- abcCCLE %>% mutate(Variant.Classification = if_else(Variant.Classification == 'Frame_Shift_Del', 'Deletion - Frameshift', Variant.Classification))

# Filter out mutations in genes that are identified as drivers for the 4 mutation datasets: 
abcCCLE_sel <- merge(DRIVER, abcCCLE, by.x="Gene", by.y="Hugo.Symbol")
abcCCLE_sel <- select(abcCCLE_sel, Gene,Variant.Classification,Protein.Change)
sel <- which(abcCCLE_sel$Variant.Classification== "Substitution - coding silent")
abcCCLE_sel <- abcCCLE_sel[-sel,]

abcCOSMIC_sel <- merge(DRIVER, abcCOSMIC, by="Gene")
abcCOSMIC_sel <- select(abcCOSMIC_sel, Gene, AA.Mutation, Type)
sel <- which(abcCOSMIC_sel$Type== "Substitution - coding silent")
abcCOSMIC_sel <- abcCOSMIC_sel[-sel,]
sel <- which(abcCOSMIC_sel$Type== "Unknown")
abcCOSMIC_sel <- abcCOSMIC_sel[-sel,]

abcCANSAR_sel <- merge(DRIVER, abcCANSAR, by.x="Gene", by.y="V1")
abcCANSAR_sel <- select(abcCANSAR_sel, Gene, V4, V6)
sel <- which(abcCANSAR_sel$V6== "Substitution - coding silent")
abcCANSAR_sel <- abcCANSAR_sel[-sel,]
sel <- which(abcCANSAR_sel$V6== "Unknown")
abcCANSAR_sel <- abcCANSAR_sel[-sel,]


# Add cell-line and source to filtered datasets:
abcCCLE_sel$cell.line <- c("ABC-1")
abcCCLE_sel <- abcCCLE_sel[,c(4,1,2,3)]
colnames(abcCCLE_sel) <- c("sample", "gene", "variant_class", "amino.acid.change")

abcCOSMIC_sel$cell.line <- c("ABC-1")
abcCOSMIC_sel <- abcCOSMIC_sel[,c(4,1,3,2)]
colnames(abcCOSMIC_sel) <- c("sample", "gene", "variant_class", "amino.acid.change")

abcCANSAR_sel$cell.line <- c("ABC-1")
abcCANSAR_sel <- abcCANSAR_sel[,c(4,1,3,2)]
colnames(abcCANSAR_sel) <- c("sample", "gene", "variant_class", "amino.acid.change")


## UO31:
UOCOSMIC <- read.delim("~/Desktop/mutAnalysis/3-COSMIC.csv", sep = ",")
UOCCLE <- read.table("~/Desktop/mutAnalysis/3-CCLE.txt", header=TRUE, sep="\t", comment="", quote="", stringsAsFactors=FALSE)
UOCANSAR <- read.table("~/Desktop/mutAnalysis/3-CANSAR.txt", header=FALSE, sep="\t", comment="", quote="", stringsAsFactors=FALSE)

# Change categories for CCLEdataset to match COSMIC and CANSAR
UOCCLE <- UOCCLE %>% mutate(Variant.Classification = if_else(Variant.Classification == 'Missense_Mutation', 'Substitution - Missense', Variant.Classification))
UOCCLE <- UOCCLE %>% mutate(Variant.Classification = if_else(Variant.Classification == 'Silent', 'Substitution - coding silent', Variant.Classification))
UOCCLE <- UOCCLE %>% mutate(Variant.Classification = if_else(Variant.Classification == 'Nonsense_Mutation', 'Substitution - Nonsense', Variant.Classification))
UOCCLE <- UOCCLE %>% mutate(Variant.Classification = if_else(Variant.Classification == 'Frame_Shift_Ins', 'Insertion - Frameshift', Variant.Classification))
UOCCLE <- UOCCLE %>% mutate(Variant.Classification = if_else(Variant.Classification == 'In_Frame_Del', 'Deletion - In frame', Variant.Classification))
UOCCLE <- UOCCLE %>% mutate(Variant.Classification = if_else(Variant.Classification == 'Frame_Shift_Del', 'Deletion - Frameshift', Variant.Classification))

# Filter out mutations in genes that are identified as drivers for the 4 mutation datasets: 
UOCCLE_sel <- merge(DRIVER, UOCCLE, by.x="Gene", by.y="Hugo.Symbol")
UOCCLE_sel <- select(UOCCLE_sel, Gene,Variant.Classification,Protein.Change)
sel <- which(UOCCLE_sel$Variant.Classification== "Substitution - coding silent")
UOCCLE_sel <- UOCCLE_sel[-sel,]

UOCOSMIC_sel <- merge(DRIVER, UOCOSMIC, by="Gene")
UOCOSMIC_sel <- select(UOCOSMIC_sel, Gene, AA.Mutation, Type)
sel <- which(UOCOSMIC_sel$Type== "Substitution - coding silent")
UOCOSMIC_sel <- UOCOSMIC_sel[-sel,]
sel <- which(UOCOSMIC_sel$Type== "Unknown")
UOCOSMIC_sel <- UOCOSMIC_sel[-sel,]

UOCANSAR_sel <- merge(DRIVER, UOCANSAR, by.x="Gene", by.y="V1")
UOCANSAR_sel <- select(UOCANSAR_sel, Gene, V4, V6)
sel <- which(UOCANSAR_sel$V6== "Substitution - coding silent")
UOCANSAR_sel <- UOCANSAR_sel[-sel,]


# Add cell-line and source to filtered datasets:
UOCCLE_sel$cell.line <- c("UO-31")
UOCCLE_sel <- UOCCLE_sel[,c(4,1,2,3)]
colnames(UOCCLE_sel) <- c("sample", "gene", "variant_class", "amino.acid.change")

UOCOSMIC_sel$cell.line <- c("UO-31")
UOCOSMIC_sel <- UOCOSMIC_sel[,c(4,1,3,2)]
colnames(UOCOSMIC_sel) <- c("sample", "gene", "variant_class", "amino.acid.change")

UOCANSAR_sel$cell.line <- c("UO-31")
UOCANSAR_sel <- UOCANSAR_sel[,c(4,1,3,2)]
colnames(UOCANSAR_sel) <- c("sample", "gene", "variant_class", "amino.acid.change")


## OVCAR-4:
OVCARCOSMIC <- read.delim("~/Desktop/mutAnalysis/4-COSMIC.csv", sep = ",")
OVCARCCLE <- read.table("~/Desktop/mutAnalysis/4-CCLE.txt", header=TRUE, sep="\t", comment="", quote="", stringsAsFactors=FALSE)
OVCARCANSAR <- read.table("~/Desktop/mutAnalysis/4-CANSAR.txt", header=FALSE, sep="\t", comment="", quote="", stringsAsFactors=FALSE)

# Change categories for CCLEdataset to match COSMIC and CANSAR
OVCARCCLE <- OVCARCCLE %>% mutate(Variant.Classification = if_else(Variant.Classification == 'Missense_Mutation', 'Substitution - Missense', Variant.Classification))
OVCARCCLE <- OVCARCCLE %>% mutate(Variant.Classification = if_else(Variant.Classification == 'Silent', 'Substitution - coding silent', Variant.Classification))
OVCARCCLE <- OVCARCCLE %>% mutate(Variant.Classification = if_else(Variant.Classification == 'Nonsense_Mutation', 'Substitution - Nonsense', Variant.Classification))
OVCARCCLE <- OVCARCCLE %>% mutate(Variant.Classification = if_else(Variant.Classification == 'Frame_Shift_Ins', 'Insertion - Frameshift', Variant.Classification))
OVCARCCLE <- OVCARCCLE %>% mutate(Variant.Classification = if_else(Variant.Classification == 'In_Frame_Del', 'Deletion - In frame', Variant.Classification))
OVCARCCLE <- OVCARCCLE %>% mutate(Variant.Classification = if_else(Variant.Classification == 'Frame_Shift_Del', 'Deletion - Frameshift', Variant.Classification))

# Filter out mutations in genes that are identified as drivers for the 4 mutation datasets: 
OVCARCCLE_sel <- merge(DRIVER, OVCARCCLE, by.x="Gene", by.y="Hugo.Symbol")
OVCARCCLE_sel <- select(OVCARCCLE_sel, Gene,Variant.Classification,Protein.Change)
sel <- which(OVCARCCLE_sel$Variant.Classification== "Substitution - coding silent")
OVCARCCLE_sel <- OVCARCCLE_sel[-sel,]

OVCARCOSMIC_sel <- merge(DRIVER, OVCARCOSMIC, by="Gene")
OVCARCOSMIC_sel <- select(OVCARCOSMIC_sel, Gene, AA.Mutation, Type)
sel <- which(OVCARCOSMIC_sel$Type== "Substitution - coding silent")
OVCARCOSMIC_sel <- OVCARCOSMIC_sel[-sel,]
sel <- which(OVCARCOSMIC_sel$Type== "Unknown")
OVCARCOSMIC_sel <- OVCARCOSMIC_sel[-sel,]

OVCARCANSAR_sel <- merge(DRIVER, OVCARCANSAR, by.x="Gene", by.y="V1")
OVCARCANSAR_sel <- select(OVCARCANSAR_sel, Gene, V4, V6)
sel <- which(OVCARCANSAR_sel$V6== "Substitution - coding silent")
OVCARCANSAR_sel <- OVCARCANSAR_sel[-sel,]
sel <- which(OVCARCANSAR_sel$V6== "Unknown")
OVCARCANSAR_sel <- OVCARCANSAR_sel[-sel,]


# Add cell-line and source to filtered datasets:
OVCARCCLE_sel$cell.line <- c("OVCAR-4")
OVCARCCLE_sel <- OVCARCCLE_sel[,c(4,1,2,3)]
colnames(OVCARCCLE_sel) <- c("sample", "gene", "variant_class", "amino.acid.change")

OVCARCOSMIC_sel$cell.line <- c("OVCAR-4")
OVCARCOSMIC_sel <- OVCARCOSMIC_sel[,c(4,1,3,2)]
colnames(OVCARCOSMIC_sel) <- c("sample", "gene", "variant_class", "amino.acid.change")

OVCARCANSAR_sel$cell.line <- c("OVCAR-4")
OVCARCANSAR_sel <- OVCARCANSAR_sel[,c(4,1,3,2)]
colnames(OVCARCANSAR_sel) <- c("sample", "gene", "variant_class", "amino.acid.change")



### Combine all datasets for waterfall plots
allmut = rbind(coloCCLE_sel, coloCOSMIC_sel, coloCANSAR_sel, 
               UOCCLE_sel, UOCOSMIC_sel, UOCANSAR_sel, 
              OVCARCCLE_sel, OVCARCOSMIC_sel, OVCARCANSAR_sel, 
               abcCCLE_sel, abcCOSMIC_sel, abcCANSAR_sel, 
               rkoCCLE_sel, rkoCOSMIC_sel, rkoCANSAR_sel)

### Generating waterfall plot
mutation_priority <- as.character(unique(allmut$variant_class))
pdf("waterfall-allcells-mutationdatabases-combined.pdf")
waterfall(allmut, fileType = "Custom", 
          mainXlabel = TRUE, main_geneLabSize = 2.5, 
          mainDropMut = TRUE, variant_class_order=mutation_priority, 
          sampOrder = c("COLO320DM",
                        "UO-31",
                        "OVCAR-4",
                        "ABC-1",
                        "RKO"))
dev.off()


