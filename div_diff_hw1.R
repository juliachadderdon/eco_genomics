
###0.5 missingness

library(vcfR)
library(tidyverse)
library(qqman)

x11.options(type="cairo")

# read in vcf file from repo
vcf.5 <- read.vcfR("~/projects/eco_genomics/population_genomics/outputs/vcf_final.filtered.0.5.vcf.gz")

#read in meta data - info on population of origin, what region
meta <- read.csv("/gpfs1/cl/pbio3990/PopulationGenomics/metadata/meta4vcf.csv")

head(meta)
#vcf file has 595 samples (type vcf into console)
dim(meta)
#meta has 629 indivs --> need to make vcf and meta compatable

meta2.5 <- meta[meta$id %in% colnames(vcf.5@gt[,-1]),]
# %in% keeps only data that both files have a counterpart to

dim(meta2.5)

#calculate diversity stats using the genetic_diff fxn in vcfR

vcf.5.dif <- genetic_diff(vcf,
                        pops = as.factor(meta2.5$region),
                        method = "nei")
# pop = factor we're interested in

str(vcf.5.dif)

chr.main <- unique(vcf.5.dif$CHROM)[1:8]

chrnum <- as.data.frame(cbind(chr.main, seq(1, 8, 1)))

#merge this and vcf
vcf.dif.5.MHplot <- left_join(chrnum, vcf.5.dif, join_by(chr.main==CHROM))

head(vcf.dif.5.MHplot)
dim(vcf.dif.5.MHplot)

vcf.dif.5.MHplot <- vcf.dif.5.MHplot %>%
  filter(Gst>0) %>%
  mutate(SNP=paste0(chr.main, "_", POS))  

vcf.dif.5.MHplot$V2=as.numeric(vcf.dif.5.MHplot$V2)

vcf.dif.5.MHplot$POS=as.numeric(vcf.dif.5.MHplot$POS)

manhattan(vcf.dif.5.MHplot, 
          chr="V2",
          bp="POS",
          p="Gst",
          col=c("blue4","orange3"),
          logp = F,
          ylab="Fst among regions",
          suggestiveline = quantile(vcf.dif.MHplot$Gst, 0.999))

write.csv(vcf.dif.5.MHplot, "~/projects/eco_genomics/population_genomics/outputs/Genetic_Diff_byRegion.5.csv",
          quote=F,
          row.names=F)

names(vcf.dif.5.MHplot)

vcf.dif.5.MHplot %>%
  as_tibble() %>%
  pivot_longer(c(4:9)) %>%
  ggplot(aes(x=value, fill=name)) +
  geom_histogram(position = "identity", alpha=1, bins=50) +
  labs(title = "Genome-wide expected Heterozygosity (Hs)", 
       fill="Regions",
       x="Gene Diveristy within Regions",
       y="Counts of SNPs")
#tidy operations use %>% --> pipe results
#use + to add more features to ggplot (title, axes lables, etc)

ggsave("Histogram_GenomDiversity_byRegion.pdf",
       path = "~/projects/eco_genomics/population_genomics/figures/")
#ggsave = saves last plot

vcf.dif.5.MHplot %>%
  as_tibble() %>%
  pivot_longer(c(4:9)) %>%
  group_by(name) %>%
  filter(value!=0) %>%
  summarise(avg_Hs=mean(value), StdDev_Hs=sd(value), N_Hs=n())




###0.9 missingness 


library(vcfR)
library(tidyverse)
library(qqman)

x11.options(type="cairo")

# read in vcf file from repo
vcf.9 <- read.vcfR("~/projects/eco_genomics/population_genomics/outputs/vcf_final.filtered.0.9.vcf.gz")

#read in meta data - info on population of origin, what region
meta <- read.csv("/gpfs1/cl/pbio3990/PopulationGenomics/metadata/meta4vcf.csv")

head(meta)
#vcf file has 595 samples (type vcf into console)
dim(meta)
#meta has 629 indivs --> need to make vcf and meta compatable

meta2.9 <- meta[meta$id %in% colnames(vcf.9@gt[,-1]),]
# %in% keeps only data that both files have a counterpart to

dim(meta2.9)

#calculate diversity stats using the genetic_diff fxn in vcfR

vcf.9.dif <- genetic_diff(vcf.9,
                          pops = as.factor(meta2.9$region),
                          method = "nei")
# pop = factor we're interested in

str(vcf.9.dif)

chr.main <- unique(vcf.9.dif$CHROM)[1:8]

chrnum <- as.data.frame(cbind(chr.main, seq(1, 8, 1)))

#merge this and vcf
vcf.dif.9.MHplot <- left_join(chrnum, vcf.9.dif, join_by(chr.main==CHROM))

head(vcf.dif.9.MHplot)
dim(vcf.dif.9.MHplot)

vcf.dif.9.MHplot <- vcf.dif.9.MHplot %>%
  filter(Gst>0) %>%
  mutate(SNP=paste0(chr.main, "_", POS))  

vcf.dif.9.MHplot$V2=as.numeric(vcf.dif.9.MHplot$V2)

vcf.dif.9.MHplot$POS=as.numeric(vcf.dif.9.MHplot$POS)

manhattan(vcf.dif.9.MHplot, 
          chr="V2",
          bp="POS",
          p="Gst",
          col=c("blue4","orange3"),
          logp = F,
          ylab="Fst among regions",
          suggestiveline = quantile(vcf.dif.MHplot$Gst, 0.999))

write.csv(vcf.dif.9.MHplot, "~/projects/eco_genomics/population_genomics/outputs/Genetic_Diff_byRegion.9.csv",
          quote=F,
          row.names=F)

names(vcf.dif.9.MHplot)

vcf.dif.9.MHplot %>%
  as_tibble() %>%
  pivot_longer(c(4:9)) %>%
  ggplot(aes(x=value, fill=name)) +
  geom_histogram(position = "identity", alpha=1, bins=50) +
  labs(title = "Genome-wide expected Heterozygosity (Hs)", 
       fill="Regions",
       x="Gene Diveristy within Regions",
       y="Counts of SNPs")
#tidy operations use %>% --> pipe results
#use + to add more features to ggplot (title, axes lables, etc)

ggsave("Histogram_GenomDiversity_byRegion.pdf",
       path = "~/projects/eco_genomics/population_genomics/figures/")
#ggsave = saves last plot

vcf.dif.9.MHplot %>%
  as_tibble() %>%
  pivot_longer(c(4:9)) %>%
  group_by(name) %>%
  filter(value!=0) %>%
  summarise(avg_Hs=mean(value), StdDev_Hs=sd(value), N_Hs=n())
