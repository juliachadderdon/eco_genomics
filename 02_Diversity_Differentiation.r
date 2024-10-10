# estimating diversity and genetic differentation in the filtered Centaurea data

library(vcfR)
library(tidyverse)
library(qqman)

x11.options(type="cairo")

# read in vcf file from repo
vcf <- read.vcfR("~/projects/eco_genomics/population_genomics/outputs/vcf_final.filtered.vcf.gz")

#read in meta data - info on population of origin, what region
meta <- read.csv("/gpfs1/cl/pbio3990/PopulationGenomics/metadata/meta4vcf.csv")

head(meta)
#vcf file has 595 samples (type vcf into console)
dim(meta)
#meta has 629 indivs --> need to make vcf and meta compatable

meta2 <- meta[meta$id %in% colnames(vcf@gt[,-1]),]
  # %in% keeps only data that both files have a counterpart to

dim(meta2)

#calculate diversity stats using the genetic_diff fxn in vcfR

vcf.dif <- genetic_diff(vcf,
                        pops = as.factor(meta2$region),
                        method = "nei")
# pop = factor we're interested in

str(vcf.dif)

chr.main <- unique(vcf.dif$CHROM)[1:8]

chrnum <- as.data.frame(cbind(chr.main, seq(1, 8, 1)))

#merge this and vcf
vcf.dif.MHplot <- left_join(chrnum, vcf.dif, join_by(chr.main==CHROM))

head(vcf.dif.MHplot)
dim(vcf.dif.MHplot)

vcf.dif.MHplot <- vcf.dif.MHplot %>%
                  filter(Gst>0) %>%
                  mutate(SNP=paste0(chr.main, "_", POS))  

vcf.dif.MHplot$V2=as.numeric(vcf.dif.MHplot$V2)

vcf.dif.MHplot$POS=as.numeric(vcf.dif.MHplot$POS)

manhattan(vcf.dif.MHplot, 
          chr="V2",
          bp="POS",
          p="Gst",
          col=c("blue4","orange3"),
          logp = F,
          ylab="Fst among regions",
          suggestiveline = quantile(vcf.dif.MHplot$Gst, 0.999))

write.csv(vcf.dif.MHplot, "~/projects/eco_genomics/population_genomics/outputs/Genetic_Diff_byRegion.csv",
          quote=F,
          row.names=F)

names(vcf.dif.MHplot)

vcf.dif.MHplot %>%
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

vcf.dif.MHplot %>%
  as_tibble() %>%
  pivot_longer(c(4:9)) %>%
  group_by(name) %>%
  filter(value!=0) %>%
  summarise(avg_Hs=mean(value), StdDev_Hs=sd(value), N_Hs=n())

# can ask for different statistics with summerise
# != means does not equal




