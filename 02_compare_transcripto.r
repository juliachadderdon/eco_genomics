# load in 01_DeSeq2.r code before this !

library(pheatmap)
resultsNames(dds)

options(bitmapType="cairo")

# pull out results for developmental temperature 22 vs 18
res_D22vsD18 <- results(dds, name = "DevTemp_D22_vs_D18", alpha=.05)
  
# order by significance
res_D22vsD18 <-res_D22vsD18[order(res_D22vsD18$padj),]
head(res_D22vsD18)

summary(res_D22vsD18)

#make counts plot
d <- plotCounts(dds, gene="TRINITY_DN140854_c0_g5_i2", int=(c("DevTemp", "FinalTemp")), returnData=TRUE)
d
#plot
p <- ggplot(d, aes(x=DevTemp, y=count, color=DevTemp, shape=FinalTemp))+
  theme_minimal() + theme(text=element_text(size=20), 
                          panel.grid.major = element_line(color = "grey"))

p <- p + geom_point(position=position_jitter(w=0.2, h=0), size=3)
p

#plot MA 
plotMA(res_D22vsD18, ylim=c(-4,4))

# volcano plot
# --> convert deseq results object into a data frame to plot
res_df <- as.data.frame(res_D22vsD18)

# add column to data frame to denote if gene is sig differentially exp
res_df$Significant <- ifelse(res_df$padj < 0.05 & abs(res_df$log2FoldChange) >1,
                             "Significant", "Not Significant")
#plot
ggplot(res_df, aes(x = log2FoldChange, y = -log(padj), color=Significant))+
  geom_point(alpha = 0.8)+
  scale_color_manual(values = c("slateblue", "darkgreen"))+
  labs(x = "Log2 Fold Change", y = "log10 Adjusted P-value", title = "Volcano Plot")+
  theme_minimal()+
  theme(legend.position = "top")+
  geom_hline(yintercept=-log10(0.05), linetype = "dashed", color = "orange")+
  geom_vline(xintercept = c(-1,1), linetype="dashed", color = "orange")

## heat map time!

# variance stabilization on data
vsd <- vst(dds, blind = FALSE)

#look at 20 top (most significant) genes in file- all is too many visually
topgenes <- head(rownames(res_D22vsD18), 20)
matrix <- assay(vsd)[topgenes, ] #uses counts that had variance stabilization
df <- as.data.frame(colData(dds)[,c("DevTemp", "FinalTemp")])

pheatmap(matrix, annotation_col=df, show_rownames=FALSE, cluster_col=T, cluster_rows=T)
# columns are samples, rows are genes, color denotes magnitide






