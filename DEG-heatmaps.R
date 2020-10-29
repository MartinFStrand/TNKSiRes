#### Generating heatmaps and clustering based on top DEGs across the cell line panel mRNASeq BGI

# Set working directory:
setwd("~/Desktop/Analysis")

# Load txt file containing top hits (above 0.8 prob from all cell lines)
DEGtop <- read.table("~/Desktop/Analysis/DEGtop5.txt", header=TRUE, sep="\t", comment="", quote="", stringsAsFactors=FALSE)
head(DEGtop)

# Load txt file containing expression data across all cell lines:
dataset <- read.table("~/Desktop/Analysis/log.txt", header=TRUE, sep="\t", comment="", quote="", stringsAsFactors=FALSE)
x <- data.matrix(dataset[,2:37])
rownames(x) <- dataset[,1]
colnames(x) <- names(dataset)[-1]
rownames(x) <- paste(1:nrow(x),rownames(x),sep="_")

# Process raw data into log2 values, including a sort of data to get as many log2 ratios as possible from the dataset. 
y <- log2(x)
y1 <- rowSums(is.na(y)) == 0 ## what rows dont have NA's 


ga <- subset(y, y1, drop = FALSE) ## subset all rows not containing NA's to a new data.matrix
g1 <- cbind(ga[,3]-ga[,1],ga[,4]-ga[,2], ## combine the data to generate log2(Treated/control) for all pairs.
            ga[,7]-ga[,5],ga[,8]-ga[,6],
            ga[,11]-ga[,9],ga[,12]-ga[,10],
            ga[,15]-ga[,13],ga[,16]-ga[,14],
            ga[,19]-ga[,17],ga[,20]-ga[,18],
            ga[,23]-ga[,21],ga[,24]-ga[,22],
            ga[,27]-ga[,25],ga[,28]-ga[,26],
            ga[,31]-ga[,29],ga[,32]-ga[,30],
            ga[,35]-ga[,33],ga[,36]-ga[,34])


y2 <- y[!y1,] # comtains at least one NA.

g2 <- matrix(NA,nrow(y2),ncol(g1))


# Rearranging data to get all possible log2-ratios. 
for(i in 1:nrow(g2)){
  for(j in 1:(ncol(g2)/2)){
    pick <- ((j-1)*4 +1):(j*4)
    if(sum(is.na(y2[i,pick]))<2 ){ # contains at least 1 NA 
      g2[i,c(2*j -1,2*j)] <- y2[i,pick[c(3,4)]] - y2[i,pick[c(1,2)]] ## One T-C log2 ratio calculated, and one ends up as NA
    }else if(sum(is.na(y2[i,pick])) == 2 ){ ## contains 2 NA
      tmp <- y2[i,pick]
      if(sum(is.na(tmp[1:2]))==2 | (sum(is.na(tmp[3:4]))==2)){
        g2[i,c(2*j -1,2*j)] <- NA # NA for both T or both C.
      }
      g2[i,2*j-1] <- y2[i,pick[which(!is.na(tmp))[2]]] - y2[i,pick[which(!is.na(tmp))[1]]]
      g2[i,2*j] <- NA
    }else{ # 3 or more NA's
      g2[i,c(2*j -1,2*j)] <- NA # NA as result
    }
    print(g2[i,c(2*j -1,2*j)])
  }
}

# Combine data withouth NA with rearraged data:
rownames(g2) <- rownames(x)[!y1]
rr <- rownames(x)
gene <- rbind(g1,g2)[rr,]
rownames(gene) <- dataset[,1]
head(gene)

# Remove cell lines not to be used in the study, and rearrange the remaining cell lines.
gene <- gene[,c(1:2,5:6,7:8,3:4,17:18)]

# Preparing list of DEG's
DEGtop <- na.omit(DEGtop)
DEGlist <- DEGtop[,1]
rownames(DEGtop) <- DEGlist


# Sorting by GeneID
gene <- as.data.frame(gene) # convderting to df to enable filtering
matches <- na.omit(match(rownames(DEGtop),rownames(gene)))
sel <- gene[matches,]
head(sel)
sel[is.na(sel)] <- 0 # set NA's as 0 for heatmap-plot

# Prepping data for heatmap: 
sel_data <- as.matrix(cbind(rowMeans(sel[,c(1:2)], na.rm=TRUE),
            rowMeans(sel[,c(3:4)], na.rm=TRUE),
            rowMeans(sel[,c(5:6)], na.rm=TRUE),
            rowMeans(sel[,c(7:8)], na.rm=TRUE),
            rowMeans(sel[,c(9:10)], na.rm=TRUE)))

colnames(sel_data) <- c("COLO320DM",
                        "UO-31",
                        "OVCAR-4",
                        "ABC-1",
                        "RKO")

# Generating pheatmap:
library(pheatmap)
library(RColorBrewer)
my_palette <- colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100) 
col_breaks = c(seq(-4,-0.3,length=35),
               seq(-0.29,0.29,length=30),
               seq(0.3,4,length=35))

pheatmap(sel_data, cluster_cols = TRUE, cellwidth=15, color = my_palette, breaks=col_breaks, show_rownames = FALSE, file= "DEGtop_newFS.pdf")
