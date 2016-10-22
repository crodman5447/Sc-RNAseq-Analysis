# Set the working directory
setwd("/Users/mcuoco/Google Drive/Projects/OvCa Patient Ascites")

# Load required packages
library(ggplot2)
library(reshape2)
library(stats)
#library(plyr)
#library(NMF)
#library(tsne)
#library(scde)
#library(corrplot)
#library(devtools)
#library(Rtsne)

# Read in the TPM (transcript-per-million) matrices (TPM matrics are calculated by RSEM)
ov_rsem1 = read.table("ov_rsem1", header=T, sep="\t", row.names=1)
ov_rsem2 = read.table("ov_rsem2", header=T, sep="\t", row.names=1)
ov_rsem3 = read.table("ov_rsem3", header=T, sep="\t", row.names=1)

# Keep samples that only belong to this ovarian project (sample names begin with "GL"); some samples are from different projects
ov1 = ov_rsem1[ ,grep("GL", colnames(ov_rsem1))] # note: can also do grep("^GL", vector) to find GL's at the beginning of sample name
ov2 = ov_rsem2[ ,grep("GL", colnames(ov_rsem2))]
ov3 = ov_rsem3[ ,grep("GL", colnames(ov_rsem3))]

# Combine ovarian samples into one data frame
ov_tpm = data.frame(ov1, ov2, ov3)

# Convert TPM into something more appropriate (from Itay Tirosh)
# We expect there to be about 100,000 transcripts per cell, so we divide the TPM
# by 10 to make it out of 100,000 instead. This number is more representative
# of the actual mount of transcripts expressed per cell. We also add 1 to each value
# so that we end up not taking the log(0), which is undefined.
E = log2(ov_tpm/10 + 1) # I am following Itay's nomenclature as described in Melanoma Science 2016

# Calculate the pooled expression for each gene from original expression matrix (this calculation lacks the divison by 10)
G_avg = apply(ov_tpm, 1, function(x) mean(x))
Ep = log2(G_avg + 1)
# We now exclude genes from our expression matrix with an aggregate expression below 4 (from Itay)
# This helps reduce a lot of the noise from the original 20,000 genes, especially because lowly expressed genes are not
# as "accurate" as more highly expressed genes
Ep_logical = (Ep >= 4)
Ef = E[Ep_logical, ] # omitted genes with very small frequency of detection

# For each cell, determine the number of genes detected (at least 1 read was mapped)
# This is the best estimate we have for cell complexity
G = apply(Ef, 2, function(x) sum(x > 0)) # where G is the number of unique genes detected for a set of cells
# I want to create a vector with patient ID information so that I use this to color-code
# a plot containing each cell and the corresponding genes detected 
patient_id = sapply(names(G), function(x) strsplit(x, "_")[[1]][1]) # create a vector of patient IDs
G_ID = data.frame(G, patient_id) # create a data frame with genes detected and patient ID as variables

# How many cells do we have from each patient?
T = table(G_ID$patient_id)
total_samples_gl15 = T[1]
total_samples_gl17 = T[2]
total_samples_gl20 = T[3]

# Create our first qc-plot showing the genes detected across all cells
plot_G_ID = ggplot(data=G_ID, aes(x=patient_id, y=G, color=patient_id)) + geom_jitter() + xlab("") + ylab("Genes detected") + ggtitle("Genes detected across all cells before filtering") + theme_bw()
legend_names = c(paste(c("GL15 (", total_samples_gl15, " cells)"), collapse=""), paste(c("GL17 (", total_samples_gl17, " cells)"), collapse=""), paste(c("GL20 (", total_samples_gl20, " cells)"), collapse=""))
plot_G_ID = plot_G_ID + scale_color_discrete(name="Patient ID", labels=legend_names) + geom_hline(yintercept=1700) # rename legend title and draw horizontal line at 1700 (cutoff)
plot_G_ID
ggsave("plot_G_ID.png", width=6, height=5)

# We now import Itay/Rajul's curated list of housekeeping genes (Melanoma Science 2016)
HK = read.table("aad0501_Table_S16.txt", header=TRUE, sep="\t") # where HK represents a set of housekeeping genes
HK = as.character(HK[,1]) # convert to character vector so that we can call the length() function

# Get the average expression of each HK across all cells
HK_set = row.names(Ef) %in% HK # Mark rows in E as TRUE if the row name is also present in the vector 'HK'
HK_tpm = Ef[HK_set, ] # Subset our original expression matrix to contain only genes present in the vector 'HK'
HK_avg = apply(HK_tpm, 2, function(x) mean(x)) # get the average housekeeping gene expression, E, for each cell
# We are restructuring the data frame for downstream plotting utility. The wide to long format conversion is from the
# reshape2 package (developed by Hadley Wickham, author of ggplot2).
HK_wide = t(HK_avg) # transform data frame to switch rows with columns
HK_long = melt(HK_wide)
HK_long = HK_long[ ,2:3] # omit the first column that contains "1" (not sure why this shows up)
names(HK_long) = c("Cell", "Ef") # rename the variables of this data frame

# We now only select the cells that have either more than 1700 detected genes OR an average housekeeping expression, E, above 3
Ef_logical = (G >= 1700) & (HK_long$Ef >= 3) # where Ef represents the filtered expression matrix (this returns a logical vector)
Ef = Ef[ ,Ef_logical]

# For comparison, I want to generate a plot of genes detected after filtering (code is similar to that described above)
G_filter = apply(Ef, 2, function(x) sum(x > 0))
id_filter = sapply(names(G_filter), function(x) strsplit(x, "_")[[1]][1])
G_id_filter = data.frame(G_filter, id_filter)
# How many cells do we have from each patient after filtering?
T_filter = table(G_id_filter$id_filter)
gl15_filter = T_filter[1]
gl17_filter = T_filter[2]
gl20_filter = T_filter[3]
# Create our second qc-plot showing the genes detected across all cells after filtering
plot_G_ID_filter = ggplot(data=G_id_filter, aes(x=id_filter, y=G_filter, color=id_filter)) + geom_jitter() + xlab("") + ylab("Genes detected") + ggtitle("Genes detected across all cells after filtering") + theme_bw()
legend_names_filter = c(paste(c("GL15 (", gl15_filter, " cells)"), collapse=""), paste(c("GL17 (", gl17_filter, " cells)"), collapse=""), paste(c("GL20 (", gl20_filter, " cells)"), collapse=""))
plot_G_ID_filter = plot_G_ID_filter + scale_color_discrete(name="Patient ID", labels=legend_names_filter) + geom_hline(yintercept=1700)
plot_G_ID_filter
ggsave("plot_G_ID_filter.png", width=6, height=5)

# For the remaining cells and genes, we define relative expression by centering the data
# We center the data around 0 in log2 space because it is more intuitive to see fold-change relationships
Ef_avg = apply(Ef, 1, function(x) mean(x))
Er = Ef - Ef_avg # where Er represents the relative (i.e. data is centered around 0) expression matrix
# At this point:
# 0 means average expression
# 1 means 2-fold higher (2^1)
# 2 means 4-fold higher (2^2)
# 3 means 8-fold higher (2^3)

# Using principal component analysis for dimensionality reduction
base_pca = prcomp(t(Er)) # the prcomp() is in the stats package
# We want to visualize the first two principal components in 2D and color the points by patient ID
id_factor=factor(id_filter)
# source: http://hms-dbmi.github.io/scw/analysis-of-heterogeneity-and-subpopulations.html
#id_color = rainbow(length(levels(id_factor)))[id_factor] # the rainbow(n) creates a vector of n contiguous colors
#plot(base_pca$x[,1], base_pca$x[,2], pch=16, col=id_color, main='PCA for GL15, GL17, and GL20', xlab="PC1", ylab="PC2")
P = as.data.frame(base_pca$x)
pca_plot = ggplot(data=P, aes(x=P$PC1, y=P$PC2, color=id_factor)) + geom_point(size=1) + theme_bw() + xlab("PC1") + ylab("PC2") + ggtitle("PCA for GL15, GL17, GL20 Human Ovarian Ascites")
pca_plot = pca_plot + scale_color_discrete(name="Patient ID", labels=legend_names_filter)
pca_plot = pca_plot + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) # remove gridlines from plot: http://felixfan.github.io/ggplot2-remove-grid-background-margin/
pca_plot
ggsave("plot_pc1pc2.png", width=7, height=5)

# Graph PCA plot by cell complexity
pca_complexity_plot = ggplot(data=P, aes(x=P$PC1, y=P$PC2, color=G_filter)) + geom_point(size=1) + theme_bw() + xlab("PC1") + ylab("PC2") + ggtitle("PCA for GL15, GL17, GL20 Human Ovarian Ascites")
pca_complexity_plot = pca_complexity_plot + scale_color_continuous(name="Genes detected")
pca_complexity_plot = pca_complexity_plot + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
pca_complexity_plot
ggsave("plot_pc1pc2_complexity.png", width=7, height=5)

