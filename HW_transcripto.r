# transcriptomics homework script- question 3

library(DESeq2)
library(ggplot2)
library(eulerr)

# set up groups within DESeq object
dds$group <- factor(paste0(dds$DevTemp, dds$FinalTemp))
design(dds) <- ~ group  
dds <- DESeq(dds)  
dim(dds)  

resultsNames((dds))
# 1] "Intercept"               "group_D18A33_vs_D18A28" 
# [3] "group_D18BASE_vs_D18A28" "group_D22A28_vs_D18A28" 
# [5] "group_D22A33_vs_D18A28"  "group_D22BASE_vs_D18A28"




# 1. compare gene expression in dev18 at BASE vs 28

res_D18_BASE_D18_A28 <- results(dds, contrast = c("group", "D18BASE", "D18A28"),
                                 alpha = 0.05)
res_D18_BASE_D18_A28 <- res_D18_BASE_D18_A28[!is.na(res_D18_BASE_D18_A28$padj),]
res_D18_BASE_D18_A28 <- res_D18_BASE_D18_A28[order(res_D18_BASE_D18_A28$padj),]

head(res_D18_BASE_D18_A28)
summary(res_D18_BASE_D18_A28)

# make a list of DEGs

degs_D18_BASE_D18_A28 <- row.names(res_D18_BASE_D18_A28[res_D18_BASE_D18_A28$padj < 0.05,])
plotMA(res_D18_BASE_D18_A28, ylim = c(-4,4))


# 2. compare gene expression in dev18 at BASE vs 33

res_D18_BASE_D18_A33 <- results(dds, contrast = c("group", "D18BASE", "D22A33"),
                                alpha = 0.05)
res_D18_BASE_D18_A33 <- res_D18_BASE_D18_A33[!is.na(res_D18_BASE_D18_A33$padj),]
res_D18_BASE_D18_A33 <- res_D18_BASE_D18_A33[order(res_D18_BASE_D18_A33$padj),]

head(res_D18_BASE_D18_A33)
summary(res_D18_BASE_D18_A33)

# make a list of DEGs

degs_D18_BASE_D18_A33 <- row.names(res_D18_BASE_D18_A33[res_D18_BASE_D18_A33$padj < 0.05,])
plotMA(res_D18_BASE_D18_A33, ylim = c(-4,4))




# 3. compare gene expression in dev22 at BASE vs 28

res_D22_BASE_D22_A28 <- results(dds, contrast = c("group", "D22BASE", "D22A28"),
                                alpha = 0.05)
res_D22_BASE_D22_A28 <- res_D22_BASE_D22_A28[!is.na(res_D22_BASE_D22_A28$padj),]
res_D22_BASE_D22_A28 <- res_D22_BASE_D22_A28[order(res_D22_BASE_D22_A28$padj),]

head(res_D22_BASE_D22_A28)
summary(res_D22_BASE_D22_A28)

# make a list of DEGs

degs_D22_BASE_D22_A28 <- row.names(res_D22_BASE_D22_A28[res_D22_BASE_D22_A28$padj < 0.05,])
plotMA(res_D22_BASE_D22_A28, ylim = c(-4,4))


# 4. compare gene expression in dev22 at BASE vs 33

res_D22_BASE_D22_A33 <- results(dds, contrast = c("group", "D22BASE", "D22A33"),
                                alpha = 0.05)
res_D22_BASE_D22_A33 <- res_D22_BASE_D22_A33[!is.na(res_D22_BASE_D22_A33$padj),]
res_D22_BASE_D22_A33 <- res_D22_BASE_D22_A33[order(res_D22_BASE_D22_A33$padj),]

head(res_D22_BASE_D22_A33)
summary(res_D22_BASE_D22_A33)

# make a list of DEGs

degs_D22_BASE_D22_A33 <- row.names(res_D22_BASE_D22_A33[res_D22_BASE_D22_A33$padj < 0.05,])
plotMA(res_D22_BASE_D22_A33, ylim = c(-4,4))



# Make Euler Plot! 

# btwn 18 s
length(intersect(degs_D18_BASE_D18_A28, degs_D18_BASE_D18_A33)) # 37

length(degs_D18_BASE_D18_A28) # 41
length(degs_D18_BASE_D18_A33) # 559

41-37 # 4
559-37 # 522

D18_Euler <-euler(c("A28"=4, "A33"=522, "A28&A33"=37))

plot(D18_Euler, lty=1:3, quantities=TRUE, fill=c("lavender", "darkseagreen", "pink"))


# btwn 22s
length(intersect(degs_D22_BASE_D22_A28, degs_D22_BASE_D22_A33)) # 144

length(degs_D22_BASE_D22_A28) # 289
length(degs_D22_BASE_D22_A33) # 1564

289-144 # 145
1564-144 # 1420


D22_Euler <-euler(c("A28"=145, "A33"=1420, "A28&A33"=144))

plot(D22_Euler, lty=1:3, quantities=TRUE, fill=c("lavender", "darkseagreen", "pink"))




# Make Scatter Plot - 18

res_D18_BASEvsA28 <- as.data.frame(results(dds, 
                                           contrast = c("group", "D18BASE", "D18A28"), alpha = .05))

res_D18_BASEvsA33 <- as.data.frame(results(dds, 
                                           contrast = c("group", "D18BASE", "D18A33"), alpha = .05))

# merge data frames
res_df18 <- merge(res_D18_BASEvsA28, res_D18_BASEvsA33, by = "row.names", 
                  suffixes = c(".BASEvs28", ".BASEvs33"))
rownames(res_df18) <- res_df18$Row.names

res_df18 <- res_df18[,-1]

library(dplyr)
library(tidyr)

# define color mapping logic with mutate function

res_df18 <- res_df18 %>% 
  mutate(fill = case_when(padj.BASEvs28 < 0.05 & stat.BASEvs28 < 0 ~ "darkseagreen4", 
                          padj.BASEvs28 < 0.05 & stat.BASEvs28 > 0 ~ "hotpink2",
                          padj.BASEvs33 < 0.05 & stat.BASEvs33 < 0 ~ "cornflowerblue",
                          padj.BASEvs33 < 0.05 & stat.BASEvs33 > 0 ~ "orange2"))

# count number of points per fill color
color_counts <- res_df18 %>%
  group_by(fill) %>%
  summarise(count = n())

label_positions <- data.frame(fill = c("darkseagreen4", "hotpink2","cornflowerblue","orange2"),
                              x_pos = c(1, 5, 0, -7.5), y_pos = c(-5, 0, 9, 3))


label_data <- merge(color_counts, label_positions, by="fill")

plot18 <- ggplot(res_df18, aes(x = log2FoldChange.BASEvs28, y = log2FoldChange.BASEvs33, color = fill)) +
  geom_point(alpha = 1) +
  scale_color_identity() +
  geom_text(data = label_data, aes(x=x_pos, y=y_pos, label=count, color=fill, size=5)) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black")+
  geom_abline(intercept = 0, slope = -1, linetype = "dashed", color = "grey")+
  xlim(-10,10) + ylim(-10,10)+
  labs(x = "Log2FoldChange BASE vs 28", 
       y = "Log2FoldChange BASE vs 33",
       title = "How does responce in Dev 18 vary by temp treatment?") +
  theme_minimal()

plot18





# Make Scatter Plot - 22

res_D22_BASEvsA28 <- as.data.frame(results(dds, 
                                           contrast = c("group", "D22BASE", "D22A28"), alpha = .05))

res_D22_BASEvsA33 <- as.data.frame(results(dds, 
                                          contrast = c("group", "D22BASE", "D22A33"), alpha = .05))

# merge data frames
res_df22 <- merge(res_D22_BASEvsA28, res_D22_BASEvsA33, by = "row.names", 
                  suffixes = c(".BASEvs28", ".BASEvs33"))
rownames(res_df22) <- res_df22$Row.names

res_df22 <- res_df22[,-1]

# define color mapping logic with mutate function

res_df22 <- res_df22 %>% 
  mutate(fill = case_when(padj.BASEvs28 < 0.05 & stat.BASEvs28 < 0 ~ "darkseagreen4", 
                          padj.BASEvs28 < 0.05 & stat.BASEvs28 > 0 ~ "hotpink2",
                          padj.BASEvs33 < 0.05 & stat.BASEvs33 < 0 ~ "cornflowerblue",
                          padj.BASEvs33 < 0.05 & stat.BASEvs33 > 0 ~ "orange2"))

# count number of points per fill color
color_counts <- res_df22 %>%
  group_by(fill) %>%
  summarise(count = n())

label_positions <- data.frame(fill = c("darkseagreen4", "hotpink2","cornflowerblue","orange2"),
                              x_pos = c(1, 5, 0, -7.5), y_pos = c(-5, 0, 9, 3))


label_data <- merge(color_counts, label_positions, by="fill")

plot22 <- ggplot(res_df22, aes(x = log2FoldChange.BASEvs28, y = log2FoldChange.BASEvs33, color = fill)) +
  geom_point(alpha = 1) +
  scale_color_identity() +
  geom_text(data = label_data, aes(x=x_pos, y=y_pos, label=count, color=fill, size=5)) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black")+
  geom_abline(intercept = 0, slope = -1, linetype = "dashed", color = "grey")+
  xlim(-10,10) + ylim(-10,10)+
  labs(x = "Log2FoldChange BASE vs 28", 
       y = "Log2FoldChange BASE vs 33",
       title = "How does responce in Dev 22 vary by temp treatment?") +
  theme_minimal()

plot22
 



