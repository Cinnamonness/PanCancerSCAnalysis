# Import libraries
library(remotes)
#install.packages("GeneNMF") #from CRAN 
#remotes::install_github("carmonalab/GeneNMF")  # or from GitHub

library(GeneNMF)
library(Seurat)
library(ggplot2)
library(UCell)
library(patchwork)
library(Matrix)
library(RcppML)
library(viridis)
library(SeuratDisk)
library(singlet)
library(SingleCellExperiment)
library(dplyr)
library(tidyr)
library(RColorBrewer)

# Load loom file and convert it to Seurat object
loom <- Connect("data/combined_15000_genes.loom", mode = "r")
seurat_obj <- as.Seurat(
  loom,
  features = "var_names", # Use gene names as features    
  cells = "obs_names",        
  normalized = NULL,        
  scaled = NULL             
)

# Close the loom file connection
loom$close()

# Check the conversion result
head(rownames(seurat_obj))
head(colnames(seurat_obj)) 

# Data quality check
print(paste("Число генов:", nrow(seurat_obj)))
print(paste("Число сэмплов:", ncol(seurat_obj)))

# Filter out cells with unwanted annotations
seurat_obj <- subset(
  seurat_obj,
  subset = !(final_cell_annotation %in% c('THCA', 'Thyroid Normal_thyroid follicular cell'))
)

# Run Non-negative Matrix Factorization (NMF)
set.seed(123)
seurat_obj <- singlet::RunNMF(
    object = seurat_obj,
    assay = "RNA",
    slot = "data",
    ranks = 2:100, # Test factorization ranks from 2 to 100
    n_replicates = 3, # Number of replicates per rank
    tol = 1e-6, # Convergence tolerance
    maxit = 200, # Maximum iterations
    verbose = TRUE, # Print progress
    threads = 4 # Use 4 CPU threads
  )

# Analyze NMF results
cv_data <- seurat_obj@reductions$nmf@misc$cv_data
library(tidyr)
error_stats <- cv_data %>%
  group_by(k) %>%
  summarise(mean_error = mean(test_error)) %>%
  complete(k = 1:100, fill = list(mean_error = NA))

# Plot NMF cross-validation results
ggplot() +
  geom_point(
    data = cv_data,
    aes(x = k, y = test_error, color = as.factor(rep)),
    alpha = 0.6,
    size = 2,
    position = position_jitter(width = 0.3)
  ) +
  geom_line(
    data = na.omit(error_stats),
    aes(x = k, y = mean_error),
    color = "black",
    linewidth = 1
  ) +
  geom_vline(
    xintercept = 12,
    linetype = "dashed",
    color = "red"
  ) +
  annotate(
    "text",
    x = 12 + 5,
    y = max(cv_data$test_error, na.rm = TRUE),
    label = paste("Best k =", 12),
    color = "black",
    vjust = 1
  ) +
  labs(
    x = "Factorization Rank (k)",
    y = "Relative Test Set Error",
    color = "Replicate"
  ) +
  scale_x_continuous(limits = c(2, 20), breaks = seq(0, 100, 10)) +
  theme_minimal()

# Visualize NMF results by cell annotation
singlet::MetadataPlot(seurat_obj, "final_cell_annotation", reduction = "nmf")

# Gene Set Enrichment Analysis (GSEA) using GO Biological Processes
gsea_results <- singlet::RunGSEA(seurat_obj, 
                                 category = "C5",
                                 subcategory = "GO:BP",
                                 verbose = FALSE)

gsea_plot <- singlet::GSEAHeatmap(
  gsea_results,
  reduction = "nmf",
  max.terms.per.factor = 3 # Show top 3 terms per factor
)

# Apply custom color scale
custom_coolwarm <- colorRampPalette(c("#3B4CC0", "white", "#B40426"))(100)
gsea_plot + scale_fill_gradientn(colors = custom_coolwarm)
# Prepare data for modified heatmap
plot_data <- gsea_plot$data

# Step 1: Apply -log10 transformation to p-values
plot_data$value <- -log10(plot_data$value)

# Step 2: Bin values for discrete coloring
plot_data$bin <- NA
plot_data$bin[plot_data$value > 1]  <- -1  # Red for significant
plot_data$bin[plot_data$value < 1]  <- 1 # Blue for non-significant
plot_data$bin[is.na(plot_data$value)] <- NA  # Handle missing values

# Step 3: Create binned heatmap
ggplot(plot_data, aes(x = Var2, y = Var1, fill = factor(bin))) +
  geom_tile(color = "white") +
  scale_fill_manual(
    values = c("-1" = "#3B4CC0", "1" = "#B40426"),
    na.value = "white"  # White for NA values
  ) +
  theme_minimal() +
  labs(
    x = "Factor", 
    y = "GO term", 
    fill = "Binned Values"  
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = element_text(size = 6) 
  )

# Cluster cells based on NMF components
seurat_obj <- FindNeighbors(seurat_obj, dims = 1:ncol(seurat_obj@reductions$nmf), reduction = "nmf") %>%
  FindClusters(resolution = 0.8, verbose = FALSE) %>%
  RunUMAP(reduction = "nmf", dims = 1:ncol(seurat_obj@reductions$nmf), verbose = FALSE)

# Visualize specific NMF components
FeaturePlot(seurat_obj, features = paste0("NMF_", 7:12))

# Compare clustering with original annotations
df <- data.frame(
  "nmf_clusters" = seurat_obj@meta.data$seurat_clusters,
  "cell_annotation_clusters" = seurat_obj@meta.data$final_cell_annotation)

# Remove NA annotations
df <- df[!is.na(df$cell_annotation_clusters), ]

# Calculate cluster composition
df <- df %>% 
  group_by(nmf_clusters) %>% 
  count(cell_annotation_clusters) %>% 
  mutate(freq = n / sum(n))

# Create dot plot of cluster composition
ggplot(df, aes(nmf_clusters, cell_annotation_clusters, size = freq, color = n)) + 
  geom_point() + 
  theme_bw() + 
  labs(x = "NMF cluster", 
       y = "Cell annotation", 
       size = "proportion\nof cluster", 
       color = "cells in\nNMF cluster") + 
  scale_color_viridis_c(option = "D")

# Rename clusters by dominant cell type
cluster_names <- df %>% 
  slice(which.max(n)) %>% 
  pull(cell_annotation_clusters)

levels(seurat_obj@meta.data$seurat_clusters) <- make.unique(as.vector(cluster_names))

# Visualize UMAP with various parameters
DimPlot(seurat_obj, 
        reduction = "umap", 
        label = TRUE, 
        group.by = "seurat_clusters", 
        pt.size = 0.5) + NoLegend()

DimPlot(seurat_obj, 
        reduction = "umap", 
        label = FALSE, 
        group.by = "seurat_clusters", 
        pt.size = 0.5)  

DimPlot(seurat_obj, 
        reduction = "umap", 
        label = TRUE, 
        label.size = 4,  
        group.by = "seurat_clusters", 
        pt.size = 0.5) + NoLegend()


DimPlot(seurat_obj, 
        reduction = "umap", 
        label = FALSE, 
        group.by = "final_cell_annotation", 
        pt.size = 0.5)  

DimPlot(seurat_obj, 
        reduction = "umap", 
        label = TRUE, 
        label.size = 2,  
        group.by = "final_cell_annotation", 
        pt.size = 0.5) + NoLegend()

# Extract and save top genes per NMF component
nmf_loadings <- Loadings(seurat_obj[["nmf"]])
dim(nmf_loadings) 
top_genes_per_component <- apply(nmf_loadings, 2, function(x) {
  names(sort(x, decreasing = TRUE))[1:50]
})

# 1. Get NMF components matrix
nmf_matrix <- as.data.frame(seurat_obj@reductions$nmf@cell.embeddings)
nmf_matrix$cell_id <- rownames(nmf_matrix)

# Normalize components so they sum to 1 per cell
nmf_matrix_norm <- nmf_matrix %>%
  select(-cell_id) %>%
  apply(1, function(x) x / sum(x)) %>%
  t() %>%
  as.data.frame()
nmf_matrix_norm$cell_id <- rownames(nmf_matrix)

# 2. Add cell metadata
meta_data <- seurat_obj@meta.data
meta_data$cell_id <- rownames(meta_data)

# Reshape data for plotting
nmf_long <- nmf_matrix_norm %>%
  inner_join(meta_data, by = "cell_id") %>%
  pivot_longer(
    cols = starts_with("NMF"),
    names_to = "component",
    values_to = "value"
  )

# 3. Define order for cell types
ordered_labels <- c(
  'Cholangiocytes', 'CHOL',
  'Skin Normal_Melanocyte', 'MEL', 'UVM',
  'Kidney Normal_kidney loop of Henle thin ascending limb epithelial cell', 'RCC',
  'Colon Normal_TA', 'CRC',
  'Hepatocytes', 'HCC',
  'Fallopian Tube Normal_secretory cell', 'OV',
  'Head and Neck Normal_respiratory basal cell', 'HNSC',
  'Pancreas Normal_Trunk', 'PAAD',
  'Lumsec-basal (Breast)', 'TNBC',
  'LummHR-major (Breast)', 'ER+', 'HER2+'
)

nmf_long$final_cell_annotation <- factor(
  nmf_long$final_cell_annotation,
  levels = ordered_labels
)

# 4. Create stacked bar plot
n_components <- length(unique(nmf_long$component))
pastel_colors <- brewer.pal(min(12, n_components), "Set3")
if (n_components > 12) {
  pastel_colors <- colorRampPalette(brewer.pal(12, "Set3"))(n_components)
}

ggplot(nmf_long, aes(x = final_cell_annotation, y = value, fill = component)) +
  geom_bar(stat = "summary", fun = "mean", position = "stack", width = 0.4) +  
  scale_fill_manual(values = pastel_colors) +
  labs(
    x = "Tissue Type / Cancer",
    y = "Mean NMF Component Fraction",
    fill = "NMF Component"
  ) +
  theme_minimal(base_size = 6) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major.x = element_blank()
  )
