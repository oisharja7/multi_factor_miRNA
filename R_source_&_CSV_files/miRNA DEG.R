library(limma)

# working directory
setwd("D:/DEG analysis_oisharja/multi_factor_miRNA")

# list agilengt text files 
exp_files <- list.files(pattern = "*.txt")

# read raw agilent files
data_raw <- read.maimages(exp_files, source = "agilent", green.only = TRUE)

# background correction
data_bgc <- backgroundCorrect(data_raw, method = "normexp", offset = 16)

# normalizing between array 
data_norm <- normalizeBetweenArrays(data_bgc, method = "quantile")

# grouping
group <- factor(c(
  rep("Normal", 4),
  rep("Primary", 4),
  rep("Recurrent", 4),
  rep("Primary", 4),
  rep("Recurrent", 4)
))

length(group)

# design matrix
design <- model.matrix(~0 + group)
colnames(design) <- levels(group)

# fitting linear model 
fit <- lmFit(data_norm, design)

# contrast
contrast.matrix <- makeContrasts(
  Prim_vs_Norm = Primary - Normal,
  Rec_vs_Norm = Recurrent - Normal,
  Rec_vs_Prim = Recurrent - Primary,
  levels = group
)

fit2 <- contrasts.fit(fit, contrast.matrix)

# eBayes moderation
fit2 <- eBayes(fit2)

# top tables for each contrast 
results_prim_norm <- topTable(fit2, coef = "Prim_vs_Norm", adjust.method = "BH", number = Inf)
results_rec_norm <- topTable(fit2, coef = "Rec_vs_Norm", adjust.method = "BH", number = Inf)
results_rec_prim <- topTable(fit2, coef = "Rec_vs_Prim", adjust.method = "BH", number = Inf)

# signifiant miRNAs 
sig_prim_norm <- results_prim_norm[results_prim_norm$adj.P.Val < 0.05 & abs(results_prim_norm$logFC) > 1, ]

sig_rec_norm <- results_rec_norm[results_rec_norm$adj.P.Val < 0.05 & abs(results_rec_norm$logFC) > 1, ]

sig_rec_prim <- results_rec_prim[results_rec_prim$adj.P.Val < 0.05 & abs(results_rec_prim$logFC) > 1, ]

# export 
write.csv(sig_prim_norm, "sig_prim_norm_miRNA.csv")
write.csv(sig_rec_norm, "sig_rec_norm_miRNA.csv")
write.csv(sig_rec_prim, "sig_rec_prim_miRNA.csv")

# visualisation 

library(ggplot2)

# Primary vs Normal 
# conversion to data frame
df_prim_norm <- as.data.frame(results_prim_norm)
# label significant 
df_prim_norm$Significant <- ifelse(df_prim_norm$adj.P.Val < 0.05 & abs(df_prim_norm$logFC) > 1, "Significant", "Not Significant")

# volcano plot 
ggplot(df_prim_norm, aes(x = logFC, y = -log10(adj.P.Val))) +
         geom_point(aes(color = Significant), size = 1) +
         scale_color_manual(values = c("Not Significant" = "orange", "Significant" = "blue")) +
         theme_minimal() +
         labs(
           title = "Volcano Plot of Differentially Expressed miRNAs: Primary vs Normal",
           x = "Log2 Fold Change",
           y = "-Log10(Adjusted P. Value)"
         )


# Recurrent vs Normal 
# conversion to data frame
df_rec_norm <- as.data.frame(results_rec_norm)
# label significant 
df_rec_norm$Significant <- ifelse(df_rec_norm$adj.P.Val < 0.05 & abs(df_rec_norm$logFC) > 1, "Significant", "Not Significant")

# volcano plot 
ggplot(df_rec_norm, aes(x = logFC, y = -log10(adj.P.Val))) +
  geom_point(aes(color = Significant), size = 1) +
  scale_color_manual(values = c("Not Significant" = "green", "Significant" = "brown")) +
  theme_minimal() +
  labs(
    title = "Volcano Plot of Differentially Expressed miRNAs: Recurrent vs Normal",
    x = "Log2 Fold Change",
    y = "-Log10(Adjusted P. Value)"
  )


# Recurrent vs Primary 
# conversion to data frame
df_rec_prim <- as.data.frame(results_rec_prim)
# label significant 
df_rec_prim$Significant <- ifelse(df_rec_prim$adj.P.Val < 0.05 & abs(df_rec_prim$logFC) > 1, "Significant", "Not Significant")

# volcano plot 
ggplot(df_rec_prim, aes(x = logFC, y = -log10(adj.P.Val))) +
  geom_point(aes(color = Significant), size = 1) +
  scale_color_manual(values = c("Not Significant" = "black", "Significant" = "yellow")) +
  theme_minimal() +
  labs(
    title = "Volcano Plot of Differentially Expressed miRNAs: Recurrent vs Primary",
    x = "Log2 Fold Change",
    y = "-Log10(Adjusted P. Value)"
  )


head(data_norm$genes)

head(rownames(results_prim_norm))

head(rownames(data_norm$E))


# Heatmap: Primary vs Normal

# Ensure Row Names are Assigned to the Expression Matrix
# assigned the SystematicName from data_norm$genes to data_norm$E.
# (This step assumes your data_norm$genes contains the SystematicName column.)
rownames(data_norm$E) <- data_norm$genes$SystematicName

# Top 50 Significant miRNAs Using Numeric Indices
# Order the results by adjusted p-value and extract the row names.
# In your results, the row names are numeric indices (as characters).
ordered_indices_1 <- as.numeric(rownames(results_prim_norm[order(results_prim_norm$adj.P.Val), ]))
top_indices_primvsnorm <- head(ordered_indices, 50)

# Subset the Expression Matrix Using These Numeric Indices
# assuming that the numeric indices in results_prim_norm refer to row positions in data_norm$E.
sig_expr_primvsnorm <- data_norm$E[top_indices_primvsnorm, ]

# Make sure the row names of the subset are unique.
rownames(sig_expr_primvsnorm) <- make.unique(rownames(sig_expr_primvsnorm))

# Plot the Heatmap
# Make sure sample annotation data frame, sample_annotation, has row names matching the column names of data_norm$E)
library(pheatmap)
pheatmap(sig_expr_primvsnorm,
         scale = "row",                       # Standardize each row (miRNA)
         annotation_col = sample_annotation,  # Sample annotation (e.g., group labels)
         main = "Heatmap of Top 50 Significant miRNAs (Primary vs Normal)",
         angle_col = 45,
         fontsize_col = 8,
         cellwidth = 30)




# Heatmap: Recurrent vs Normal 

# Top 50 Significant miRNAs Using Numeric Indices
ordered_indices_2 <- as.numeric(rownames(results_rec_norm[order(results_rec_norm$adj.P.Val), ]))
top_indices_recvsnorm <- head(ordered_indices_2, 50)

# Subset the Expression Matrix Using These Numeric Indices
sig_expr_recvsnorm <- data_norm$E[top_indices_recvsnorm, ]

# Make sure the row names of the subset are unique.
rownames(sig_expr_recvsnorm) <- make.unique(rownames(sig_expr_recvsnorm))

# Plot the Heatmap
library(pheatmap)
pheatmap(sig_expr_recvsnorm,
         scale = "row",                       # Standardize each row (miRNA)
         annotation_col = sample_annotation,  # Sample annotation (e.g., group labels)
         main = "Heatmap of Top 50 Significant miRNAs (Recurrent vs Normal)",
         angle_col = 45,
         fontsize_col = 8,
         cellwidth = 30)



# Top 50 Significant miRNAs Using Numeric Indices
ordered_indices_3 <- as.numeric(rownames(results_rec_prim[order(results_rec_prim$adj.P.Val), ]))
top_indices_recvsprim <- head(ordered_indices_3, 50)

# Subset the Expression Matrix Using These Numeric Indices
sig_expr_recvsprim <- data_norm$E[top_indices_recvsprim, ]

# Make sure the row names of the subset are unique.
rownames(sig_expr_recvsprim) <- make.unique(rownames(sig_expr_recvsprim))

# Plot the Heatmap
library(pheatmap)
pheatmap(sig_expr_recvsprim,
         scale = "row",                       # Standardize each row (miRNA)
         annotation_col = sample_annotation,  # Sample annotation (e.g., group labels)
         main = "Heatmap of Top 50 Significant miRNAs (Recurrent vs Primary)",
         angle_col = 45,
         fontsize_col = 8,
         cellwidth = 30)


