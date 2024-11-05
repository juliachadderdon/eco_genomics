# Script for analyzing & visualizing gene correlation networks

library(DESeq2)
library(ggplot2)
library(WGCNA); options(stringsAsFactors = FALSE);
library(tidyverse)
library(CorLevelPlot)
library(gridExtra)
library(Rmisc)

setwd("~/projects/eco_genomics/transcriptomics/")


# step 1: import counts data

countsTable <- read.table("/gpfs1/cl/pbio3990/Transcriptomics/tonsa_counts.txt",
                          header = T, row.names = 1)

countsTableRound <- round(countsTable) #Because DESeq2 does not like decimals
tail(countsTable)

#Next is the conditions file
conds <- read.delim("/gpfs1/cl/pbio3990/Transcriptomics/experimental_details.txt",
                    header = T, stringsAsFactors = T, row.names = 1)
head(conds)


traitData = read.table("/gpfs1/cl/pbio3990/Transcriptomics/Trait_Data.txt",
                       header=T, row.names=1)

# filter matrix to just BASE data 

filtered_count_matrix_BASEonly <- countsTable[, conds$FinalTemp == "BASE"]
filtered_sample_metadata_BASEonly <- conds[conds$FinalTemp == "BASE", ]
rounded_filtered_count_matrix <- round(filtered_count_matrix_BASEonly)


# step 2: detecting outliers

# detect outlier genes
gsg <- goodSamplesGenes(t(rounded_filtered_count_matrix))
summary(gsg)

table(gsg$goodGenes) # 37235 bad genes, 82203 good genes

table(gsg$goodSamples) # all good

# filter out bad genes
data_WGCNA <- rounded_filtered_count_matrix[gsg$goodGenes==TRUE, ]
dim(data_WGCNA) # only the good 82203 left

#use clustering with tree dendrogram to identify outlier samples
htree <- hclust(dist(t(data_WGCNA)), method = 'average')
plot(htree)

# PCA - outlier detection method
pca <- prcomp(data_WGCNA)
pca_data <- pca$x
#make a data frame
pca_data <- as.data.frame(pca_data)

pca.var <- pca$sdev^2
pca.var.percent <- round(pca.var/sum(pca.var)*100, digits = 2)

ggplot(pca_data, aes(PC1, PC2))+
  geom_point()+
  geom_text(label = rownames(pca_data))+
  labs(x = paste0('PC1: ', pca.var.percent[1], '%'),
       y = paste0('PC2: ', pca.var.percent[2]))


# step 3: normalization

colData <- row.names(filtered_sample_metadata_BASEonly)

dds_WGCNA <- DESeqDataSetFromMatrix(countData = data_WGCNA, 
                                    colData = filtered_sample_metadata_BASEonly,
                                    design = ~1) # there r no specified groups

# filter
dds_WGCNA_75 <- dds_WGCNA[rowSums(counts(dds_WGCNA) >= 15) >= 6, ]
nrow(dds_WGCNA_75) # filtered down to 29,559 transcripts

# variance normalization/ stabilization
dds_norm <- vst(dds_WGCNA_75) 

# get and save normalized counts to use below
norm.counts <- assay(dds_norm) %>% t()


# step 4: network construction

# choose set of soft-thresholding powers
power <- c(c(1:10), seq(from=12, to=50, by=2))

# call the network topology analysis function
sft <- pickSoftThreshold(norm.counts,
                         powerVector = power,
                         networkType = "signed", # to focus on transcripts that are positively correlated
                         verbose = 5)

sft.data <- sft$fitIndices

# plot to pick power

a1 <- ggplot(sft.data, aes(Power, SFT.R.sq, label = Power))+
  geom_point()+
  geom_text(nudge_y = 0.1)+
  geom_hline(yintercept = 0.8, color = "red")+
  labs(x = "Power", y = "Scale free topology model fit, signed R^2")+
  theme_classic()

a2 <- ggplot(sft.data, aes(Power, mean.k., label = Power))+
  geom_point()+
  geom_text(nudge_y = 0.1)+
  geom_hline(yintercept = 0.8, color = "red")+
  labs(x = "Power", y = "Mean Connectivity")+
  theme_classic()

grid.arrange(a1, a2, nrow = 2)





