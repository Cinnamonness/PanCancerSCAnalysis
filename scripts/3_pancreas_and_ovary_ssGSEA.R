# Import libraries
library(GSVA)
library(vegan)
library(ggplot2)
library(ggrepel)
library(dplyr)
library(tidyr)
library(patchwork)
library(FactoMineR)  
library(factoextra)  

# Load gene expression dataset
datasets <- read.csv("data/pancreas_wilcoxon_datasets_df_top_100_genes.csv", 
                     row.names = 1, 
                     check.names = FALSE,
                     stringsAsFactors = FALSE)
# Define gene signatures for different cell types
genelist <- list(
  "Pancreas Normal_Acinar cell" = c('GPHA2',
                                    'SLC38A5',
                                    'NIBAN1',
                                    'DEF6',
                                    'CELA3B',
                                    'CELA3A',
                                    'SPINK1',
                                    'CLPSL1',
                                    'DPEP1',
                                    'PLA2G1B',
                                    'REG3G',
                                    'EGF',
                                    'CELA2B',
                                    'AQP8',
                                    'REG1A',
                                    'PRSS3',
                                    'MT1G',
                                    'SLC39A5',
                                    'MT1H',
                                    'MYRIP',
                                    'CEL',
                                    'PNLIP',
                                    'TRAF5',
                                    'PNLIPRP1',
                                    'CUZD1',
                                    'PM20D1',
                                    'AZGP1',
                                    'FKBP11',
                                    'TEX11',
                                    'LDHB',
                                    'SLC16A10',
                                    'SLC39A8',
                                    'PTF1A',
                                    'TMEM97',
                                    'PAIP2B',
                                    'SLC22A31',
                                    'AQP12B',
                                    'LMO3',
                                    'BCAT1',
                                    'BCL2L14',
                                    'MYH7',
                                    'SPX',
                                    'NRG4',
                                    'CTRL',
                                    'PDIA2',
                                    'C2orf88',
                                    'SLC43A1',
                                    'IL22RA1',
                                    'CLPS',
                                    'CPB1',
                                    'RNASE1',
                                    'ANPEP',
                                    'SYCN',
                                    'TPST2',
                                    'CPA1',
                                    'CPA4',
                                    'FBP2',
                                    'RARRES2',
                                    'CPA2',
                                    'CTRC',
                                    'GATA4',
                                    'KLK1',
                                    'CELA2A',
                                    'SERPINI2',
                                    'MXRA5',
                                    'QPRT',
                                    'CTRB1',
                                    'CTRB2',
                                    'ZG16',
                                    'PGGHG',
                                    'LPAR3',
                                    'PRSS1',
                                    'HEYL',
                                    'MAP3K5',
                                    'TMEM131L',
                                    'HHIP-AS1',
                                    'LFNG',
                                    'GSTA2',
                                    'AADAC',
                                    'TRHDE',
                                    'MECOM',
                                    'PLTP',
                                    'ANKRD62',
                                    'FITM1',
                                    'ARHGDIG',
                                    'KCTD15',
                                    'RAB3D',
                                    'CYP2E1',
                                    'CXCL12',
                                    'NUPR1',
                                    'GATA2',
                                    'PODXL',
                                    'PIWIL2',
                                    'PSAT1',
                                    'CYTIP',
                                    'PABPC4',
                                    'AKR7A3',
                                    'GP2',
                                    'P2RX1',
                                    'TMEM52'),
  "Pancreas Normal_Trunk" = c('SNRPD3',
                              'NCOA7',
                              'TINAGL1',
                              'IL17RD',
                              'FLNB',
                              'PTPRG',
                              'COL16A1',
                              'ECHDC1',
                              'SERINC2',
                              'AKAP7',
                              'EPB41L2',
                              'PTPRK',
                              'KIAA1671',
                              'TBC1D10A',
                              'LIF',
                              'CCN2',
                              'VNN1',
                              'STX7',
                              'PRKCD',
                              'PTP4A2',
                              'TWF2',
                              'PERP',
                              'TNFAIP3',
                              'CGAS',
                              'KDM6B',
                              'SAT2',
                              'KIFC3',
                              'PLLP',
                              'RIPOR1',
                              'RO60',
                              'TRADD',
                              'SLC34A2',
                              'SOD3',
                              'DHX15',
                              'NFATC3',
                              'EIF5A',
                              'IL34',
                              'CDH3',
                              'PPARGC1A',
                              'KLF3',
                              'SEL1L3',
                              'ZNF267',
                              'GORAB',
                              'BCL7C',
                              'ARL8A',
                              'HSD3B7',
                              'ELF3',
                              'RNPEP',
                              'STX4',
                              'DNAJA2',
                              'TENT4B',
                              'IRX3',
                              'PYCARD',
                              'VPS35',
                              'LAD1',
                              'TGFB1I1',
                              'ZNF281',
                              'CNGA1',
                              'NIPAL1',
                              'CHIC2',
                              'BRCC3',
                              'OCIAD2',
                              'FLNA',
                              'MELTF',
                              'DOK4',
                              'CES1',
                              'LPCAT2',
                              'COQ9',
                              'INAVA',
                              'NEK7',
                              'SCNN1G',
                              'FLII',
                              'ALDH3A2',
                              'GPRC5B',
                              'IQCK',
                              'SHMT1',
                              'ARL6IP1',
                              'IL4R',
                              'NUAK2',
                              'BCKDK',
                              'ADORA2B',
                              'ZNF286A',
                              'KDM5B',
                              'ADORA1',
                              'PRSS8',
                              'SCO1',
                              'GINM1',
                              'ADGRG6',
                              'SFT2D1',
                              'AFDN',
                              'EZR',
                              'SLC22A3',
                              'WTAP',
                              'MAPK1',
                              'MCM3',
                              'APOL1',
                              'ADGRF1',
                              'C6orf141',
                              'MARCKSL1',
                              'GNAI2'),
  "Pancreas Normal_α" = c('DZIP3',
                          'NXPE3',
                          'ANXA6',
                          'GRM7-AS3',
                          'GRM7',
                          'KCNIP1',
                          'FBLL1',
                          'GABRG2',
                          'ABCG1',
                          'DSCAM',
                          'CHL1',
                          'NSG2',
                          'CNTN4',
                          'LDLRAP1',
                          'RCAN3',
                          'PCDHB9',
                          'PCDHB16',
                          'CCDC181',
                          'GRAMD2A',
                          'LOXL1',
                          'DRAIC',
                          'EMP2',
                          'LINC01290',
                          'ROGDI',
                          'CDIP1',
                          'NMNAT2',
                          'GAP43',
                          'DNAJA4',
                          'SAXO2',
                          'RASGRF1',
                          'HTR1F',
                          'STXBP5L',
                          'LSAMP',
                          'KCNJ6',
                          'MAN1C1',
                          'GPRIN1',
                          'CPLX2',
                          'SNCB',
                          'ADAMTS2',
                          'TTC3',
                          'BRWD1-AS2',
                          'PIGP',
                          'RIPPLY3',
                          'DSCR4',
                          'MIS18A-AS1',
                          'C21orf62-AS1',
                          'PAXBP1-AS1',
                          'KIF1A',
                          'SIL1',
                          'LRRTM2',
                          'WNT4',
                          'CTXN2',
                          'ILDR2',
                          'B3GALT2',
                          'LYSMD2',
                          'SCG3',
                          'GABPB1-AS1',
                          'C15orf61',
                          'APH1B',
                          'C2CD4A',
                          'SNX29',
                          'PPL',
                          'PODXL2',
                          'RAB39B',
                          'KIAA0319',
                          'SCGN',
                          'TIAM1',
                          'LINC01011',
                          'CYYR1-AS1',
                          'SLC17A1',
                          'GPLD1',
                          'CYYR1',
                          'REEP2',
                          'KLHL3',
                          'CCDC188',
                          'DISP2',
                          'PAK6',
                          'MEIS2',
                          'GJD2',
                          'RGS4',
                          'RXRG',
                          'MAP1A',
                          'DLL4',
                          'PDZRN3',
                          'ATP6AP1',
                          'ROBO2',
                          'ROBO1',
                          'SYT17',
                          'CKMT1A',
                          'CATSPER2',
                          'HCFC1-AS1',
                          'LKAAEAR1',
                          'CUTA',
                          'KCNK17',
                          'KCNK16',
                          'DOCK10',
                          'MYT1',
                          'GRM4',
                          'SERPIND1',
                          'SAR1B')
)

# Set up GSVA parameters
params <- ssgseaParam(
  exprData = as.matrix(datasets),
  geneSets = genelist,
  alpha = 0.25, # Weighting parameter for gene ranks
  normalize = FALSE # Don't normalize scores by gene set size
)

# Run GSVA enrichment analysis
gsva_results <- gsva(params)
results_df <- as.data.frame(t(gsva_results))

# Filter to keep only cancer samples
results_df$Group <- gsub(".*_(.*)$", "\\1", rownames(results_df))
results_df$Group <- ifelse(results_df$Group == "HER2+/ER+", "HER2+",
                           ifelse(results_df$Group == "PR+", "TNBC", 
                                  results_df$Group))
results_df$Group <- factor(results_df$Group)

# Visualize the result
results_df %>%
  pivot_longer(cols = -Group, names_to = "Signature", values_to = "Score") %>%
  ggplot(aes(x = Group, y = Score, fill = Group)) +
  geom_boxplot() +
  facet_wrap(~Signature, scales = "free_y") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "ssGSEA Scores by Cell Types", 
       x = "Cell Types", 
       y = "Enrichment Score")

# PCoA analysis 
pcoa_data <- results_df[, sapply(results_df, is.numeric)]
dist_matrix <- vegdist(pcoa_data, method = "bray")

pcoa_result <- tryCatch({
  cmdscale(dist_matrix, eig = TRUE, k = 2)
}, error = function(e) {
  stop("PCoA failed: ", e$message)
})

if(is.null(pcoa_result$points)) {
  stop("PCoA didn't return any points")
}

pcoa_df <- data.frame(
  PC1 = pcoa_result$points[,1],
  PC2 = pcoa_result$points[,2],
  Group = results_df$Group,
  Sample = rownames(results_df),
  stringsAsFactors = FALSE
)

positive_eigs <- pcoa_result$eig[pcoa_result$eig > 0]
if(length(positive_eigs) == 0) {
  stop("All eigenvalues are zero or negative")
}
var_exp <- pcoa_result$eig/sum(positive_eigs) * 100

tryCatch({
  p <- ggplot(pcoa_df, aes(x = PC1, y = PC2, color = Group)) +
    geom_point(size = 3) +
    stat_ellipse(level = 0.95) +
    xlab(paste0("PCo1 (", round(var_exp[1], 1), "%)")) +
    ylab(paste0("PCo2 (", round(var_exp[2], 1), "%)")) +
    ggtitle("PCoA of ssGSEA Scores") +
    theme_minimal(base_size = 12) +
    theme(legend.position = "right")
  
  if(nrow(pcoa_df) <= 20) {
    p <- p + geom_text_repel(aes(label = Sample), size = 3, max.overlaps = Inf)
  }
  
  print(p)
}, error = function(e) {
  message("Plotting failed: ", e$message)
})

# PCA analysis
pca_data <- results_df[, sapply(results_df, is.numeric)]

if(any(is.na(pca_data))) {
  pca_data[is.na(pca_data)] <- colMeans(pca_data, na.rm = TRUE)
}

pca_result <- prcomp(pca_data, scale. = TRUE)
var_exp <- pca_result$sdev^2 / sum(pca_result$sdev^2) * 100

pca_df <- data.frame(
  PC1 = pca_result$x[,1],
  PC2 = pca_result$x[,2],
  Group = results_df$Group,
  Sample = rownames(results_df)
)

ggplot(pca_df, aes(x = PC1, y = PC2, color = Group)) +
  geom_point(size = 3) +
  stat_ellipse(level = 0.95) +
  xlab(paste0("PC1 (", round(var_exp[1], 1), "%)")) +
  ylab(paste0("PC2 (", round(var_exp[2], 1), "%)")) +
  ggtitle("PCA of ssGSEA Scores") +
  theme_minimal()

# Build boxplots
results_df$Sample <- rownames(results_df)

results_df %>%
  pivot_longer(cols = -c(Group, Sample), names_to = "Signature", values_to = "Score") %>%
  ggplot(aes(x = Sample, y = Score, fill = Group, color = Group)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.5) +
  geom_jitter(width = 0.2, size = 1.5, alpha = 0.8) +
  facet_wrap(~Signature, scales = "free_y") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.3, hjust = 1, size = 5),
    strip.text = element_text(size = 5)
  ) +
  labs(
    title = "ssGSEA Scores by Individual Samples", 
    x = "Sample", 
    y = "Enrichment Score"
  )

long_df <- results_df %>%
  pivot_longer(
    cols = -Group,
    names_to = "Tissue",
    values_to = "Value"
  )

tissue_colors <- c(
  "Pancreas Normal_Acinar cell" = "#F781BF",
  "Pancreas Normal_Trunk" = "#B2DF8A",
  "Pancreas Normal_α" = "#A6CEE3"
#  "Ovary Normal_eStromal" = "#984EA3"
)

long_df %>%
  filter(!Group %in% c('Pancreas Normal_Trunk', 'Pancreas Normal_α', 'Pancreas Normal_Acinar cell')) %>%
  ggplot(aes(x = Tissue, y = Value, fill = Tissue)) +
  geom_boxplot(alpha = 0.6, outlier.shape = NA) +
  geom_jitter(width = 0.2, size = 1.5, alpha = 0.8, aes(color = Tissue)) +
  scale_fill_manual(values = tissue_colors, drop = FALSE) +
  scale_color_manual(values = tissue_colors, drop = FALSE) +
  facet_wrap(~ Group, scales = "free_y") +
  theme_minimal(base_size = 10) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1),
    plot.title = element_text(hjust = 0.5),
    strip.background = element_blank(),   
    strip.text = element_text(face = "bold") 
  ) +
  labs(
    title = "ssGSEA Scores by Cancer Types",
    x = "Normal Tissues",
    y = "Enrichment Score"
  )

long_df %>%
  filter(Group == "Cancer") %>%
  ggplot(aes(x = Tissue, y = Value, fill = Tissue)) +
  geom_boxplot(alpha = 0.6, outlier.shape = NA) +
  geom_jitter(width = 0.2, size = 1.5, alpha = 0.8, aes(color = Tissue)) +
  scale_fill_manual(values = tissue_colors, drop = FALSE) +
  scale_color_manual(values = tissue_colors, drop = FALSE) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(hjust = 0.5),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold")
  ) +
  labs(
    title = "ssGSEA Scores in Cancer",
    x = "Tissue",
    y = "Enrichment Score"
  )


