library(vcfR)

setwd("/gpfs1/cl/pbio3990/PopulationGenomics/")
list.files()

vcf <- read.vcfR("variants/Centaurea_filtered.vcf.gz")
vcf
head(vcf)

dna <- ape::read.dna("reference/GCA_030169165.1_ASM3016916v1_genomic.fa.gz", format="fasta")

gff <- read.table("reference/GCA_030169165.1_ASM3016916v1_genomic_genes.gff.gz", sep="\t", quote="")

chr1 <- create.chromR(name="Chromosome 1", vcf=vcf, seq=dna, ann=gff)

plot(chr1)

pdf(file="~/projects/eco_genomics/population_genomics/figures/ChromoPlot.pdf")
chromoqc(chr1, xlim=c(1e1, 1.1e8))
dev.off()

DP <- extract.gt(vcf, element="DP", as.numeric=T)

DP[1:5,1:10]

quantile(DP)

DP[DP==0] <- NA

quantile(DP, na.rm=T)

heatmap.bp(DP[1:1000,], rlabels = F, clabels = F)

library(SNPfiltR)

vcf.filt <- hard_filter(vcf, depth=3)
# filtering out low depths (lower than 30)

max_depth(vcf.filt, maxdepth = 60)
#filter out depths 2x mean read (>60 reads/SNPs)

meta <- read.csv("metadata/meta4vcf.csv", header=T)

meta2 <- meta[,c(1,4)]
# all rows, columns 1-4 : [rows,columns]

names(meta2) <- c("id","pop")

meta2$id = as.factor(meta2$id)
meta2$pop = as.factor(meta2$pop)

vcf.filt.indMiss_0.5 <- missing_by_sample(vcf.filt,
                                      popmap=meta2,
                                      cutoff=0.50)

vcf.filt.indMiss_0.5 <- filter_biallelic(vcf.filt.indMiss_0.5)
vcf.filt.indMiss_0.5 <- min_mac(vcf.filt.indMiss_0.5, min.mac = 1)


vcf.filt.indSNPMiss_0.5 <- missing_by_snp(vcf.filt.indMiss_0.5, cutoff=0.5)

DP2 <- extract.gt(vcf.filt.indSNPMiss_0.5,
                  element = "DP",
                  as.numeric = T)
#heatmap.bp(as.matrix(DP2[1:5000,]),
           rlabels = F, clabels = F)

write.vcf(vcf.filt.indSNPMiss_0.5,
          "~/projects/eco_genomics/population_genomics/outputs/vcf_final.filtered.0.5.vcf.gz")


vcf.filt.indMiss_0.9 <- missing_by_sample(vcf.filt,
                                          popmap=meta2,
                                          cutoff=.9)

vcf.filt.indMiss_0.9 <- filter_biallelic(vcf.filt.indMiss_0.9)
vcf.filt.indMiss_0.9 <- min_mac(vcf.filt.indMiss, min.mac = 1)

vcf.filt.indSNPMiss_0.9 <- missing_by_snp(vcf.filt.indMiss_0.9, cutoff=.5)

DP2 <- extract.gt(vcf.filt.indSNPMiss_0.9,
                  element = "DP",
                  as.numeric = T)
heatmap.bp(as.matrix(DP2[1:5000,]),
           rlabels = F, clabels = F)

write.vcf(vcf.filt.indSNPMiss_0.9,
          "~/projects/eco_genomics/population_genomics/outputs/vcf_final.filtered.0.9.vcf.gz")




