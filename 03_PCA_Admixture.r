library(tidyverse)
library(vcfR)
library(SNPfiltR)
library(LEA)

options(bitmapType="cairo")

setwd("~/projects/eco_genomics/population_genomics/")

vcf <- read.vcfR("outputs/vcf_final.filtered.vcf.gz")

# need to thin SNPs for LD b4 we run PCA and Admixture 
# to satisfy assumptions of independence among loci

vcf.thin <- distance_thin(vcf, min.distance = 500)

meta <- read.csv("/gpfs1/cl/pbio3990/PopulationGenomics/metadata/meta4vcf.csv")
dim(meta)

meta2 <- meta[meta$id %in% colnames(vcf@gt[, -1]) , ]
dim(meta2)

write.vcf(vcf.thin, "outputs/vcf_final.filtered.thinned.vcf.gz")

# hide uncompressed vcf file, bc too big for github, outside of repo

system("gunzip -c ~/projects/eco_genomics/population_genomics/outputs/vcf_final.filtered.thinned.vcf.gz > ~/vcf_final.filtered.thinned.vcf")

geno <- vcf2geno(input.file="/gpfs1/home/j/c/jchadder/vcf_final.filtered.thinned.vcf",
                 output.file = "outputs/vcf_final.filtered.thinned.vcf.geno")

CentPCA <- LEA::pca("outputs/vcf_final.filtered.thinned.vcf.geno", scale=TRUE)


CentPCA <- load.pcaProject("vcf_final.filtered.thinned.vcf.pcaProject")
  #load in previous pca

show(CentPCA)

plot(CentPCA)

plot(CentPCA$projections,
     col=as.factor(meta2$region))
legend("bottomright", legend=as.factor(unique(meta2$region)), 
                                       fill=as.factor(unique(meta2$region)))

ggplot(as.data.frame(CentPCA$projections),
      aes(x=V1, y=V2, color=meta2$region, shape=meta2$continent))+
      geom_point(alpha=1)+
      labs(title="Centaurea genetic PCA", x="PC1", y="PC2", color="Region", shape="Continent")
# close parentheses at end of aes line for ggplot, add + to add more
# could set xlim and ylim to zoom in  on specific stuff
# if want to see PC2 & PC3, etc just change x= and y=

ggsave("figures/CentPCA_PC1vPC2.pdf", width=6, height=6, units="in")
# ggsave saves last plot


# running admix analysis and create plots
# use LEA R package with function "snmf"

CentAdmix <- snmf("outputs/vcf_final.filtered.thinned.vcf.geno",
                  K=1:10,
                  entropy = T,
                  repetitions = 3,
                  project = "new")  
  #if adding to analysis later, could choose project= "continue"

plot(CentAdmix)
  #can play around with K to fit data
  # cross-entropy indicates how well data fits K value (lower values better)

par(mfrow=c(2,1))
plot(CentAdmix, col="blue4", main="SNMF")
plot(CentPCA$eigenvalues[1:10], 
     ylab="Eigenvalues", xlab="Number of PCs",
     col="blue4", main="PCA")
# PCA and SNMF match up, yay!

dev.off()
#resets plotting

myK=5
 # so you don't change whole script's K while we play around
CE = cross.entropy(CentAdmix, K=myK)
best = which.min(CE)
  # see which CE is best

myKQ = Q(CentAdmix, K=myK, run=best)

myKQmeta = cbind(myKQ, meta2)

my.colors = c("blue4", "pink2", "lightblue", "olivedrab", "orange2")

myKQmeta = as_tibble(myKQmeta) %>%
  group_by(continent) %>%
  arrange(region, pop, .by_group = TRUE)

barplot(as.matrix(t(myKQmeta[ , 1:myK])),
        border = NA,
        space = 0,
        col = my.colors[1:myK],
        xlab = "Geographic Regions", ylab = "Ancestry Proportions",
        main = paste0("Ancestry matrix K=",myK))
a


