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


# 1. compare baseline gene expression between developmental treatment groups
res_D18_BASE_D22_BASE <- results(dds, contrast = c("group", "D18BASE", "D22BASE"),
                            alpha = 0.05)
res_D18_BASE_D22_BASE <- res_D18_BASE_D22_BASE[!is.na(res_D18_BASE_D22_BASE$padj),]
res_D18_BASE_D22_BASE <- res_D18_BASE_D22_BASE[order(res_D18_BASE_D22_BASE$padj),]

head(res_D18_BASE_D22_BASE)
summary(res_D18_BASE_D22_BASE)

# make a list of which genes in comparisons of interest are differentailly expressed
  # (list of DEGs)

degs_D18_BASE_D22_BASE <- row.names(res_D18_BASE_D22_BASE[res_D18_BASE_D22_BASE$padj < 0.05,])
plotMA(res_D18_BASE_D22_BASE, ylim = c(-4,4)) #plot to see!

# 2. compare gene expression btwn developmental temp treatment groups

res_D18_A28_D22_A28 <- results(dds, contrast = c("group", "D18A28", "D22A28"),
                                 alpha = 0.05)
res_D18_A28_D22_A28 <- res_D18_A28_D22_A28[!is.na(res_D18_A28_D22_A28$padj),]
res_D18_A28_D22_A28 <- res_D18_A28_D22_A28[order(res_D18_A28_D22_A28$padj),]

head(res_D18_A28_D22_A28)
summary(res_D18_A28_D22_A28)

degs_D18_A28_D22_A28 <- row.names(res_D18_A28_D22_A28[res_D18_A28_D22_A28$padj < 0.05,])
plotMA(res_D18_A28_D22_A28, ylim = c(-4,4))


# 3. compare gene exp btwn dev temp and treat groups at A28

res_D18_A33_D22_A33 <- results(dds, contrast = c("group", "D18A33", "D22A33"),
                               alpha = 0.05)
res_D18_A33_D22_A33 <- res_D18_A33_D22_A33[!is.na(res_D18_A33_D22_A33$padj),]
res_D18_A33_D22_A33 <- res_D18_A33_D22_A33[order(res_D18_A33_D22_A33$padj),]

head(res_D18_A33_D22_A33)
summary(res_D18_A33_D22_A33)

degs_D18_A33_D22_A33 <- row.names(res_D18_A33_D22_A33[res_D18_A33_D22_A33$padj < 0.05,])
plotMA(res_D18_A33_D22_A33, ylim = c(-4,4))


length(degs_D18_BASE_D22_BASE) # 1935
length(degs_D18_A28_D22_A28) # 296
length(degs_D18_A33_D22_A33) # 78

# look at overlaps in which genes are differentially expressed im multiple contrasts
  # lookat  dif treatments  expressing the same gene exp: how many

length(intersect(degs_D18_BASE_D22_BASE, degs_D18_A28_D22_A28)) # 107
length(intersect(degs_D18_BASE_D22_BASE, degs_D18_A33_D22_A33)) # 44
length(intersect(degs_D18_A33_D22_A33, degs_D18_A28_D22_A28))   # 29

# overlap btwn all three
nested_intersection <- intersect(degs_D18_BASE_D22_BASE, degs_D18_A28_D22_A28)
length(intersect(degs_D18_A33_D22_A33, nested_intersection)) #23

### 10/22 - making euler plots

# calculate number of unique genes in each portion of the euler plot

1935-107-44+23 # = 1807 genes diff expressed uniquely at baseline btwn 18 vs 22
296-107-29+23 # = 183 genes uniquely expressed when exposed to 28
78-44-29+23 # = 28 uniquely expressed when exposed to 33

107-23 # = 84 genes uniq to BASE ans A28
44-23 # = 21 genes uniq to BASE and A33
29-23 # = 6 genes uniq to A28 and A33


myEuler <-euler(c("BASE"=1807, "A28"=183, "A33"=28, "BASE&A28"=84, 
                  "BASE&A33"=21, "A28&A33"=6, "BASE&A28&A33"=23))

plot(myEuler, lty=1:3, quantities=TRUE)


## make scatter plot of responses to A28/A33 when copepods at 18 vs 22


# contrast D18_A28vsBASE
res_D18_BASEvsA28 <- as.data.frame(results(dds, 
                    contrast = c("group", "D18BASE", "D18A28"), alpha = .05))

# contrast D22_A28vsBASE

res_D22_BASEvsA28 <- as.data.frame(results(dds, 
                    contrast = c("group", "D22BASE", "D22A28"), alpha = .05))

# merge data frames
res_df28 <- merge(res_D18_BASEvsA28, res_D22_BASEvsA28, by = "row.names", 
                  suffixes = c(".18", ".22"))
rownames(res_df28) <- res_df28$Row.names

res_df28 <- res_df28[,-1]

library(dplyr)
library(tidyr)

# define color mapping logic with mutate function

res_df28 <- res_df28 %>% 
  mutate(fill = case_when(padj.18 < 0.05 & stat.18 < 0 ~ "turquoise", 
                          padj.18 < 0.05 & stat.18 > 0 ~ "magenta",
                          padj.22 < 0.05 & stat.22 < 0 ~ "blue2",
                          padj.22 < 0.05 & stat.22 > 0 ~ "orange"))

# count number of points per fill color
color_counts <- res_df28 %>%
  group_by(fill) %>%
  summarise(count = n())

label_positions <- data.frame(fill = c("turquoise", "magenta","blue2","orange"),
  x_pos = c(1, 5, 0, -7.5), y_pos = c(-5, 0, 9, 3))


label_data <- merge(color_counts, label_positions, by="fill")

plot28 <- ggplot(res_df28, aes(x = log2FoldChange.18, y = log2FoldChange.22, color = fill)) +
  geom_point(alpha = 1) +
  scale_color_identity() +
  geom_text(data = label_data, aes(x=x_pos, y=y_pos, label=count, color=fill, size=5)) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black")+
  geom_abline(intercept = 0, slope = -1, linetype = "dashed", color = "grey")+
  xlim(-10,10) + ylim(-10,10)+
  labs(x = "Log2FoldChange 28 vs base at 18", 
       y = "Log2FoldChange 28 vs base at 22",
       title = "How does responce to 28 C vary by DevTemp?") +
  theme_minimal()

plot28


## Repeat for A33

# contrast D18_A33vsBASE
res_D18_BASEvsA33 <- as.data.frame(results(dds, 
                                           contrast = c("group", "D18BASE", "D18A33"), alpha = .05))

# contrast D22_A33vsBASE

res_D22_BASEvsA33 <- as.data.frame(results(dds, 
                                           contrast = c("group", "D22BASE", "D22A33"), alpha = .05))

# merge data frames
res_df33 <- merge(res_D18_BASEvsA33, res_D22_BASEvsA33, by = "row.names", 
                  suffixes = c(".18", ".22"))
rownames(res_df33) <- res_df33$Row.names

res_df33 <- res_df33[,-1]

library(dplyr)
library(tidyr)

# define color mapping logic with mutate function

res_df33 <- res_df33 %>% 
  mutate(fill = case_when(padj.18 < 0.05 & stat.18 < 0 ~ "turquoise", 
                          padj.18 < 0.05 & stat.18 > 0 ~ "magenta",
                          padj.22 < 0.05 & stat.22 < 0 ~ "blue2",
                          padj.22 < 0.05 & stat.22 > 0 ~ "orange"))

# count number of points per fill color
color_counts <- res_df33 %>%
  group_by(fill) %>%
  summarise(count = n())

label_positions <- data.frame(fill = c("turquoise", "magenta","blue2","orange"),
                              x_pos = c(1, 5, 0, -7.5), y_pos = c(-5, 0, 9, 3))


label_data <- merge(color_counts, label_positions, by="fill")

plot33 <- ggplot(res_df33, aes(x = log2FoldChange.18, y = log2FoldChange.22, color = fill)) +
  geom_point(alpha = 1) +
  scale_color_identity() +
  geom_text(data = label_data, aes(x=x_pos, y=y_pos, label=count, color=fill, size=5)) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black")+
  geom_abline(intercept = 0, slope = -1, linetype = "dashed", color = "grey")+
  labs(x = "Log2FoldChange 33 vs base at 18", 
       y = "Log2FoldChange 33 vs base at 22",
       title = "How does responce to 33 C vary by DevTemp?") +
  theme_minimal()

plot33

ggsave("~/projects/eco_genomics2024/transcriptomics/figures/combined_scatter_plot.png", 
       combined_plot, width = 12, hight = 6)

## put the 2 plots together in a 2 panel plot

library(gridExtra)

combined_plot <- grid.arrange(plot28, plot33, ncol=2)

# notes about plot: many more diff. exp. genes in response to 33c



