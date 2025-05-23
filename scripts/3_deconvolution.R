# Import libraries
library(DeconRNASeq)
library(ggpubr)
library(reshape2)
library(vegan)

### Pancreas dataset
# Upload pseudobulk samples
dat_pancreas=read.csv("data/pancreas_wilcoxon_datasets_df_top_100_genes.csv", row.names = 1, check.names = F)
# Upload signatures for tissue types and remove cancer signature
sig_pancreas=read.csv("data/pancreas_wilcoxon_signatures_df_top_100_genes.csv",row.names = 1, check.names = F)
sig_pancreas = sig_pancreas[,-1]

# Run Deconvolution
decon=DeconRNASeq(dat_pancreas, sig_pancreas,proportions, checksig=FALSE, known.prop = F, use.scale = TRUE, fig = TRUE)
# Extract the results ("fractions" of each cell type per sample)
fractions=as.data.frame(decon$out.all)
rownames(fractions)=colnames(dat_pancreas)
# Check that fractions of each cell type are equal to 1 for each sample
rowSums(fractions)
# Visualisation
fr=melt(t(fractions))
cancer_samples <- rownames(fractions)[grepl("Cancer", rownames(fractions))]
fractions_cancer <- fractions[cancer_samples, ]
fr_cancer <- melt(t(fractions_cancer))
ggbarplot(fr_cancer, x="Var2", y="value", fill="Var1", palette = "Set2", size = 0.5) + 
  rotate_x_text(90) + 
  theme(axis.text.x = element_text(size = 5))

### BRCA dataset
# Upload pseudobulk samples
dat_brca=read.csv("data/breast_wilcoxon_datasets_df_top_100_genes.csv", row.names = 1, check.names = F)
# Upload signatures for tissue types and remove cancer signatures
sig_brca=read.csv("D:/Обучение в ИБ/Диплом/sc-RNA_data/deconvolution/breast_wilcoxon_signatures_df_top_100_genes.csv",row.names = 1, check.names = F)
sig_brca=sig_brca[,c(3,4,5,7)]

# Run Deconvolution
decon=DeconRNASeq(dat_brca, sig_brca, proportions, checksig=FALSE, known.prop = F, use.scale = TRUE, fig = FALSE)
# Extract the results ("fractions" of each cell type per sample)
fractions=as.data.frame(decon$out.all)
rownames(fractions)=colnames(dat_brca)
# Check that fractions of each cell type are equal to 1 for each sample
rowSums(fractions)
# Visualisation
fr=melt(t(fractions))
cancer_samples <- rownames(fractions)[grepl("Cancer|HER2\\+|ER\\+|PR\\+|TNBC", rownames(fractions))]
fractions_cancer <- fractions[cancer_samples, ]
fr_cancer <- melt(t(fractions_cancer))
ggbarplot(fr_cancer, x="Var2", y="value", fill="Var1", palette = "Set2", size = 0.5) + 
  rotate_x_text(90) + 
  theme(axis.text.x = element_text(size = 5))

### Ovarian dataset
dat_ovary=read.csv("data/ovary_wilcoxon_datasets_df_top_100_genes.csv", row.names = 1, check.names = F)
sig_ovary=read.csv("data/ovary_wilcoxon_signatures_df_top_100_genes.csv",row.names = 1, check.names = F)
sig_ovary = sig_ovary[,-1]

# Run Deconvolution
decon=DeconRNASeq(dat_ovary, sig_ovary, proportions, checksig=FALSE, known.prop = F, use.scale = TRUE, fig = TRUE)
# Extract the results ("fractions" of each cell type per sample)
fractions=as.data.frame(decon$out.all)
rownames(fractions)=colnames(dat_ovary)
# Check that fractions of each cell type are equal to 1 for each sample
rowSums(fractions)
# Visualisation
fr=melt(t(fractions))
cancer_samples <- rownames(fractions)[grepl("Cancer", rownames(fractions))]
fractions_cancer <- fractions[cancer_samples, ]
fr_cancer <- melt(t(fractions_cancer))
ggbarplot(fr_cancer, x="Var2", y="value", fill="Var1", palette = "Set2", size = 0.5) + 
  rotate_x_text(90) + 
  theme(axis.text.x = element_text(size = 5))
