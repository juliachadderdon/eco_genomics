library(DESeq2)
library(ggplot2)

#install.packages("DESeq2")
setwd("~/projects/eco_genomics/transcriptomics/")

#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")

#BiocManager::install("DESeq2")

#BiocManager::install("Biobase", dependencies = T, force = T)

#Import counts matrix

countsTable <- read.table("/gpfs1/cl/pbio3990/Transcriptomics/tonsa_counts.txt",
                          header = T, row.names = 1)

countsTableRound <- round(countsTable) #Because DESeq2 does not like decimals
tail(countsTable)

#Next is the conditions file
conds <- read.delim("/gpfs1/cl/pbio3990/Transcriptomics/experimental_details.txt",
                    header = T, stringsAsFactors = T, row.names = 1)

head(conds)


##explore counts matrix

#see how many reads from each sample with colSums

colSums(countsTableRound)
mean(colSums(countsTableRound))
#good #of reads --> start w 20m=great, ~18m we hav after cleaning=great

barplot(colSums(countsTableRound), names.arg = colnames(countsTableRound), 
        cex.names = 0.5, las = 2, ylim = c(0,30000000))
abline(h=mean(colSums(countsTableRound)), col = "blue4", lwd=2)
  #abline makes line across the barplot

#the average number of counts per gene 
rowSums(countsTableRound)
mean(rowSums(countsTableRound)) #3244.739 
median(rowSums(countsTableRound)) #64

apply(countsTableRound,2,mean) # 2 tells it to go across rows, 
#gives a sense of variation in sequencing effort across samples

apply(countsTableRound,1,mean) # 1 tells it to go with columns


## DESeq2 time

dds <- DESeqDataSetFromMatrix(countData = countsTableRound, colData = conds,
                              design = ~DevTemp + FinalTemp)
dim(dds)

#Filtering

dds <- dds[rowSums(counts(dds) >= 10) >=15, ]
nrow(dds) #35,527 = number of transcripts with >10 reads and >or= to 15 samples

#run deseq model to test for global differential gene expression
dds <- DESeq(dds)

#list results generated
resultsNames(dds)

# visualize global gene expression patterns using PCA
# first- transform data for plotting using variance stabilization

vsd <- vst(dds, blind=FALSE)

pcaData <- plotPCA(vsd, intgroup=c("DevTemp", "FinalTemp"), returnData=TRUE)
percentVar <- round(100*attr(pcaData, "percentVar"))

final_temp_colors <- c("BASE"="lightpink","a28"="hotpink","a33"="red")
shapes_choose <- c("D18"=16, "D22"=18)

p <- ggplot(pcaData, aes(PC1, PC2, color=FinalTemp, shape=DevTemp)) +
  geom_point(size=5) + 
  scale_shape_manual(values = shapes_choose) + 
  scale_color_manual(values = final_temp_colors) +
  labs(x = paste0('PC1: ', percentVar[1], ' %'),
       y = paste0("PC2: ", percentVar[2], ' %')) +
  theme_bw(base_size = 16)
  

p 





