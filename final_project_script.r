library(tidyverse)
library(vcfR)
library(SNPfiltR)
library(LEA)

capituladata <- read.csv ("/gpfs1/cl/pbio3990/PopulationGenomics/traits/capitula/PNW_EU_NE_capitulummeasurements.csv")


setwd("/gpfs1/cl/pbio3990/PopulationGenomics/traits/capitula/")


dat <- read.csv("PNW_EU_NE_capitulummeasurements.csv", header=T)

str(dat)

dat = dat %>%
  mutate(id = paste(Pop, IndID, sep="_"), .before=5) %>%
  group_by(id) %>%
  slice_head(n=1)

head(dat)

summary(dat)


pca1 <- prcomp(dat[,c(7:20)], center=T, scale=T)

pca2 <- prcomp(dat[,c(16,18,19)], center=T, scale=T)


autoplot(pca1, data = dat, colour = 'Region',
         loadings = TRUE, loadings.colour = 'blue',
         loadings.label = TRUE, loadings.label.size = 3)

autoplot(pca2, data = dat, colour = 'Region',
         loadings = TRUE, loadings.colour = 'blue',
         loadings.label = TRUE, loadings.label.size = 3)


# pcadapt !!!

library(tidyverse)
library(ggplot2)
library(vcfR)
library(qqman)
library(pcadapt)


fullvcf <- read.vcfR("/gpfs1/cl/pbio3990/PopulationGenomics/variants/vcf_final.filtered.vcf")

vcfR <- read.vcfR("/gpfs1/cl/pbio3990/PopulationGenomics/variants/vcf_final.filtered.vcf.gz")

meta <- read.csv("/gpfs1/cl/pbio3990/PopulationGenomics/metadata/meta4vcf.csv")
meta2 <- meta[meta$id %in% colnames(vcfR@gt[,-1]),]


#just europe data in pcadapt

europe <- c("WEU", "CEU", "NEU","SEU")

europemeta <- meta2 %>% filter(region %in% europe)
europeids <- europemeta$id

europe_indices <- which(colnames(fullvcf@gt) %in% europeids)

europe.vcf <- fullvcf[, c(1,europe_indices)]


write.vcf(europe.vcf, "final_project/europe.vcf") 

europe.vcf <- read.pcadapt("final_project/europe.vcf", type="vcf")

head(europe.vcf)

tail(europe.vcf)

pcadapt.europe <- pcadapt(europe.vcf,
                       K=2,
                       method="componentwise",
                       min.maf=0.01,
                       LD.clumping = list(size=500, thr=0.2))


summary(pcadapt.europe) 

# Manhattan plot,  without chromosomal info
plot(pcadapt.europe, K=2) 



vcfR.fix <- as.data.frame(vcfR@fix[,1:2])

chr.main <- unique(vcfR.fix$CHROM)[1:8]

chrnum <- as.data.frame(cbind(chr.main, seq(1,8,1)))

Pval <- pcadapt.europe$pvalues
pcadapt.europe.MHplot <- cbind(vcfR.fix, Pval)

pcadapt.europe.MHplot <- left_join(chrnum, pcadapt.europe.MHplot, join_by(chr.main==CHROM))

pcadapt.europe.MHplot <- pcadapt.europe.MHplot %>%
  mutate(SNP=paste0(chr.main,"_",POS))

# Set variables as numeric that contain numbers
pcadapt.europe.MHplot$V2 = as.numeric(pcadapt.europe.MHplot$V2)
pcadapt.europe.MHplot$POS = as.numeric(pcadapt.europe.MHplot$POS)
pcadapt.europe.MHplot$pPC1 = as.numeric(pcadapt.europe.MHplot[,4])
pcadapt.europe.MHplot$pPC2 = as.numeric(pcadapt.europe.MHplot[,5])

# Drop loci with "NA" values -- these got filtered out by the min.maf option in the pcadapt() call
pcadapt.europe.MHplot <- pcadapt.europe.MHplot %>% drop_na(pPC1)

# Let's make separate ones for selection outleirs on PC1 and then do one for PC2 next 
manhattan(pcadapt.europe.MHplot,
          chr="V2",
          bp="POS",
          p="pPC1",
          col=c("blue4","orange3"),
          logP=T,
          ylab="-log10 p-value",
          genomewideline=F,
          main="PCAdapt genome scan for selection in Europe (PC1)")

manhattan(pcadapt.MHplot,
          chr="V2",
          bp="POS",
          p="pPC2",
          col=c("blue4","orange3"),
          logP=T,
          ylab="-log10 p-value",
          genomewideline=F,
          main="PCAdapt genome scan for selection (PC2)")

# There are more sophisticated ways to pull out the identities of the outliers
# but a quick way we did in class was to filter based on the quantile of the p-values,
# here, choosing a quantile of 0.001, corresponding to the smallest 0.1% of p-values as outliers

View(pcadapt.europe.MHplot %>%
       filter(pPC1<quantile(pcadapt.europe.MHplot$pPC1,0.001)))

eu_pcadapt <- pcadapt.europe.MHplot %>% filter(pPC1<quantile(pcadapt.europe.MHplot$pPC1,0.001))

View(eu_pcadapt)



# PNW

pnw <-"PNW"

pnwmeta <- meta2 %>% filter(region %in% pnw)
pnwids <- pnwmeta$id

pnw_indices <- which(colnames(fullvcf@gt) %in% pnwids)

pnw.vcf <- fullvcf[, c(1,pnw_indices)]

library(SNPfiltR)

pnw.vcf <- min_mac(pnw.vcf, min.mac=1)

vcfR.pnw.fix <- as.data.frame(pnw.vcf@fix[,1:2])

write.vcf(pnw.vcf, "final_project/pnw.vcf") 

read.vcfR("final_project/pnw.vcf") 


pnw.vcf <- read.pcadapt("final_project/pnw.vcf", type="vcf")

pcadapt.pnw <- pcadapt(pnw.vcf,
                          K=2,
                          method="componentwise",
                          min.maf=0.01,
                          LD.clumping = list(size=500, thr=0.2))


summary(pcadapt.pnw) 

# Manhattan plot,  without chromosomal info
plot(pcadapt.pnw, K=2) 




vcfR.pnw.fix <- as.data.frame(pnw.vcf@fix[,1:2])

chr.main <- unique(vcfR.pnw.fix$CHROM)[1:8]

chrnum <- as.data.frame(cbind(chr.main, seq(1,8,1)))

Pval.pnw <- pcadapt.pnw$pvalues
pcadapt.pnw.MHplot <- cbind(vcfR.pnw.fix, Pval.pnw)

View(pcadapt.pnw)
dim(vcfR.fix)
dim(Pval.pnw)

pcadapt.pnw.MHplot <- left_join(chrnum, pcadapt.pnw.MHplot, join_by(chr.main==CHROM))

pcadapt.pnw.MHplot <- pcadapt.pnw.MHplot %>%
  mutate(SNP=paste0(chr.main,"_",POS))

# Set variables as numeric that contain numbers
pcadapt.pnw.MHplot$V2 = as.numeric(pcadapt.pnw.MHplot$V2)
pcadapt.pnw.MHplot$POS = as.numeric(pcadapt.pnw.MHplot$POS)
pcadapt.pnw.MHplot$pPC1 = as.numeric(pcadapt.pnw.MHplot[,4])
pcadapt.pnw.MHplot$pPC2 = as.numeric(pcadapt.pnw.MHplot[,5])

# Drop loci with "NA" values -- these got filtered out by the min.maf option in the pcadapt() call
pcadapt.pnw.MHplot <- pcadapt.pnw.MHplot %>% drop_na(pPC1)

# Let's make separate opnws for selection outleirs on PC1 and then do one for PC2 next 
manhattan(pcadapt.pnw.MHplot,
          chr="V2",
          bp="POS",
          p="pPC1",
          col=c("blue4","orange3"),
          logP=T,
          ylab="-log10 p-value",
          genomewideline =F,
          main="PCAdapt genome scan for selection in PNW (PC1)")


View(pcadapt.pnw.MHplot %>%
       filter(pPC1<quantile(pcadapt.pnw.MHplot$pPC1,0.001)))

pnw_pcadapt <- pcadapt.pnw.MHplot %>% filter(pPC1<quantile(pcadapt.pnw.MHplot$pPC1,0.001))

View(pnw_pcadapt)


# NE US

ne <-"NE"

nemeta <- meta2 %>% filter(region %in% ne)
neids <- nemeta$id

ne_indices <- which(colnames(fullvcf@gt) %in% neids)

ne.vcf <- fullvcf[, c(1,ne_indices)]


write.vcf(ne.vcf, "final_project/ne.vcf") 


ne.vcf <- read.pcadapt("final_project/ne.vcf", type="vcf")

pcadapt.ne <- pcadapt(ne.vcf,
                       K=2,
                       method="componentwise",
                       min.maf=0.01,
                       LD.clumping = list(size=500, thr=0.2))


summary(pcadapt.ne) 

# Manhattan plot,  without chromosomal info
plot(pcadapt.ne, K=2)


vcfR.fix <- as.data.frame(vcfR@fix[,1:2])

chr.main <- unique(vcfR.fix$CHROM)[1:8]

chrnum <- as.data.frame(cbind(chr.main, seq(1,8,1)))

Pval.ne <- pcadapt.ne$pvalues
pcadapt.ne.MHplot <- cbind(vcfR.fix, Pval.ne)

dim(Pval.ne)

pcadapt.ne.MHplot <- left_join(chrnum, pcadapt.ne.MHplot, join_by(chr.main==CHROM))

pcadapt.ne.MHplot <- pcadapt.ne.MHplot %>%
  mutate(SNP=paste0(chr.main,"_",POS))

# Set variables as numeric that contain numbers
pcadapt.ne.MHplot$V2 = as.numeric(pcadapt.ne.MHplot$V2)
pcadapt.ne.MHplot$POS = as.numeric(pcadapt.ne.MHplot$POS)
pcadapt.ne.MHplot$pPC1 = as.numeric(pcadapt.ne.MHplot[,4])
pcadapt.ne.MHplot$pPC2 = as.numeric(pcadapt.ne.MHplot[,5])

# Drop loci with "NA" values -- these got filtered out by the min.maf option in the pcadapt() call
pcadapt.ne.MHplot <- pcadapt.ne.MHplot %>% drop_na(pPC1)

# Let's make separate ones for selection outleirs on PC1 and then do one for PC2 next 
manhattan(pcadapt.ne.MHplot,
          chr="V2",
          bp="POS",
          p="pPC1",
          col=c("blue4","orange3"),
          logP=T,
          ylab="-log10 p-value",
          genomewideline=F,
          main="PCAdapt genome scan for selection in NE (PC1)")

View(pcadapt.ne.MHplot %>%
       filter(pPC1<quantile(pcadapt.ne.MHplot$pPC1,0.001)))

ne_pcadapt <- pcadapt.ne.MHplot %>% filter(pPC1<quantile(pcadapt.ne.MHplot$pPC1,0.001))

View(ne_pcadapt)

# analyze overlapping sig SNPS! 

# ne & eu
ne_eu_shared_snps <- eu_pcadapt %>%
  select(SNP) %>% intersect(., ne_pcadapt %>% select(SNP)) # 1 in common
view(ne_eu_shared_snps)

# ne & pnw
ne_pnw_shared_snps <- pnw_pcadapt %>%
  select(SNP) %>% intersect(., ne_pcadapt %>% select(SNP)) # 1 in common
view(ne_pnw_shared_snps)

# pnw & eu
pnw_eu_shared_snps <- eu_pcadapt %>%
  select(SNP) %>% intersect(., pnw_pcadapt %>% select(SNP)) # 0 in common
view(pnw_eu_shared_snps)


# lowering the significance for selected SNPs since only 0-1 in common btwn groups!

ne_pcadapt.01 <- pcadapt.ne.MHplot %>% filter(pPC1<quantile(pcadapt.ne.MHplot$pPC1,0.01))
View(ne_pcadapt.01)

eu_pcadapt.01 <- pcadapt.europe.MHplot %>% filter(pPC1<quantile(pcadapt.europe.MHplot$pPC1,0.01))
View(eu_pcadapt.01)

pnw_pcadapt.01 <- pcadapt.pnw.MHplot %>% filter(pPC1<quantile(pcadapt.pnw.MHplot$pPC1,0.01))
View(pnw_pcadapt.01)



# ne & eu
ne_eu_shared_snps.01 <- eu_pcadapt.01 %>%
  select(SNP) %>% intersect(., ne_pcadapt.01 %>% select(SNP)) # 21 in common
view(ne_eu_shared_snps.01)

# ne & pnw
ne_pnw_shared_snps.01 <- pnw_pcadapt.01 %>%
  select(SNP) %>% intersect(., ne_pcadapt.01 %>% select(SNP)) # 10 in common
view(ne_pnw_shared_snps.01)

# pnw & eu
pnw_eu_shared_snps.01 <- eu_pcadapt.01 %>%
  select(SNP) %>% intersect(., pnw_pcadapt.01 %>% select(SNP)) # 8 in common
view(pnw_eu_shared_snps.01)


## Making euler plots!

# MATH

# ne & eu
92-21 # eu = 71
100-21 # ne = 79
# combined = 21

# ne & pnw
100-10 # ne = 90
86-10 # pnw = 76
# combined = 10

# pnw & eu
86-8 # pnw = 78
91-8 # eu = 83
# combined = 8

# okay actually making euler plots now

library(eulerr)

ne_eu_Euler <- euler(c("Northeast"=79, "Europe"=71, "Northeast&Europe"=21))
plot(ne_eu_Euler, quantities=TRUE, fill=c("orange3", "cadetblue", "pink3"))

ne_pnw_Euler <- euler(c("Northeast"=90, "PNW"=76, "Northeast&PNW"=10))
plot(ne_pnw_Euler, quantities=TRUE, fill=c("orange3", "darkseagreen", "pink3"))

pnw_eu_Euler <- euler(c("PNW"=78, "Europe"=83, "PNW&Europe"=8))
plot(pnw_eu_Euler, quantities=TRUE, fill=c("darkseagreen", "cadetblue", "pink3"))









