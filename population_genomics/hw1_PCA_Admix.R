library(tidyverse)
library(vcfR)
library(SNPfiltR)
library(LEA)

options(bitmapType="cairo")

setwd("~/projects/eco_genomics/population_genomics/")

### .5


vcf.5 <- read.vcfR("~/projects/eco_genomics/population_genomics/outputs/vcf_final.filtered.0.5.vcf.gz")

# need to thin SNPs for LD b4 we run PCA and Admixture 
# to satisfy assumptions of independence among loci

vcf.thin.5 <- distance_thin(vcf.5, min.distance = 500)

meta <- read.csv("/gpfs1/cl/pbio3990/PopulationGenomics/metadata/meta4vcf.csv")
dim(meta)

meta2.5 <- meta[meta$id %in% colnames(vcf.5@gt[, -1]) , ]
dim(meta2.5)

write.vcf(vcf.thin.5, "~/projects/eco_genomics/population_genomics/outputs/vcf_final.filtered.thinned.5.vcf.gz")
dim(vcf.thin.5)
# hide uncompressed vcf file, bc too big for github, outside of repo

system("gunzip -c ~/projects/eco_genomics/population_genomics/outputs/vcf_final.filtered.thinned.5.vcf.gz > ~/vcf_final.filtered.thinned.5.vcf")

geno <- vcf2geno(input.file="/gpfs1/home/j/c/jchadder/vcf_final.filtered.thinned.5.vcf",
                 output.file = "outputs/vcf_final.filtered.thinned.5.vcf.geno")

mygeno <- read.geno("outputs/vcf_final.filtered.thinned.5.vcf.geno")

CentPCA.5 <- LEA::pca(mygeno, scale=TRUE)

CentPCA.5 <- LEA::pca("outputs/vcf_final.filtered.thinned.5.vcf.geno", scale=TRUE)

CentPCA.5 <- load.pcaProject("vcf_final.filtered.thinned.5.vcf.pcaProject")
#load in previous pca

show(CentPCA.5)

plot(CentPCA.5)

plot(CentPCA.5$projections,
     col=as.factor(meta2.5$region))
legend("bottomright", legend=as.factor(unique(meta2.5$region)), 
                                       fill=as.factor(unique(meta2.5$region)))

ggplot(as.data.frame(CentPCA.5$projections),
       aes(x=V1, y=V2, color=meta2.5$region, shape=meta2.5$continent))+
  geom_point(alpha=1)+
  labs(title="Centaurea genetic PCA", x="PC1", y="PC2", color="Region", shape="Continent")
# close parentheses at end of aes line for ggplot, add + to add more
# could set xlim and ylim to zoom in  on specific stuff
# if want to see PC2 & PC3, etc just change x= and y=

ggsave("figures/CentPCA.5_PC1vPC2.pdf", width=6, height=6, units="in")
# ggsave saves last plot



### .9 

vcf.9 <- read.vcfR("~/projects/eco_genomics/population_genomics/outputs/vcf_final.filtered.0.9.vcf.gz")

# need to thin SNPs for LD b4 we run PCA and Admixture 
# to satisfy assumptions of independence among loci

vcf.thin.9 <- distance_thin(vcf.9, min.distance = 500)
vcf.thin.9 <- min_mac(vcf.9, min.mac=1)

meta <- read.csv("/gpfs1/cl/pbio3990/PopulationGenomics/metadata/meta4vcf.csv")
dim(meta)

meta2.9 <- meta[meta$id %in% colnames(vcf.9@gt[, -1]) , ]
dim(meta2.9)

write.vcf(vcf.thin.9, "~/projects/eco_genomics/population_genomics/outputs/vcf_final.filtered.thinned.9.vcf.gz")
dim(vcf.thin.9)
# hide uncompressed vcf file, bc too big for github, outside of repo


system("gunzip -c ~/projects/eco_genomics/population_genomics/outputs/vcf_final.filtered.thinned.9.vcf.gz > ~/vcf_final.filtered.thinned.9.vcf")

geno <- vcf2geno(input.file="/gpfs1/home/j/c/jchadder/vcf_final.filtered.thinned.9.vcf",
                 output.file = "~/projects/eco_genomics/population_genomics/outputs/vcf_final.filtered.thinned.9.vcf.geno")

CentPCA.9 <- LEA::pca("~/projects/eco_genomics/population_genomics/outputs/vcf_final.filtered.thinned.9.vcf.geno", scale=TRUE)

#CentPCA.9 <- load.pcaProject("vcf_final.filtered.thinned.9.vcf.pcaProject")
#load in previous pca

show(CentPCA.9)

plot(CentPCA.9)

plot(CentPCA.9$projections,
     col=as.factor(meta2.9$region))
legend("bottomright", legend=as.factor(unique(meta2.9$region)), 
       fill=as.factor(unique(meta2.9$region)))

ggplot(as.data.frame(CentPCA.9$projections),
       aes(x=V1, y=V2, color=meta2.9$region, shape=meta2.9$continent))+
  geom_point(alpha=1)+
  labs(title="Centaurea genetic PCA", x="PC1", y="PC2", color="Region", shape="Continent")

