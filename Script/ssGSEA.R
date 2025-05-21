library(GSVA)

datasets <- read.csv("D:/Обучение в ИБ/Диплом/sc-RNA_data/datasets_df_top_10_genes.csv", 
                     row.names = 1, 
                     check.names = FALSE,
                     stringsAsFactors = FALSE)

genelist <- list(
  "LummHR-major (Breast)" = c('DNALI1',
                              'LIN7A',
                              'APH1B',
                              'BCAS4',
                              'FANK1',
                              'IFT22',
                              'EFHC1',
                              'KIF16B',
                              'SMIM14',
                              'TFAP2A'),
  "Lumsec-basal (Breast)" = c('KLK6',
                              'COTL1',
                              'KLK5',
                              'KLK10',
                              'STK25',
                              'KLK7',
                              'MRPL11',
                              'MELTF',
                              'COL22A1',
                              'FUOM'),
  "Lumsec-major (Breast)" = c('PDE4B',
                              'SLC28A3',
                              'CHI3L1',
                              'NCALD',
                              'SLC25A37',
                              'LINC01554',
                              'CNKSR3',
                              'ST8SIA1',
                              'DAPP1',
                              'PLEKHS1'),
  "basal (Breast)" = c('HAS3',
                       'NRP2',
                       'ERG',
                       'TP63',
                       'CD200',
                       'FAM126A',
                       'MSRB3',
                       'GRP',
                       'MYLK',
                       'CCND2')
)

params <- ssgseaParam(
  exprData = as.matrix(datasets),
  geneSets = genelist,
  alpha = 0.25,
  normalize = FALSE 
)

# ssGSEA
gsva_results <- gsva(params)
head(gsva_results)

results_df <- as.data.frame(t(gsva_results))

# Добавление информации о группах опухолей из названий образцов
results_df$Group <- gsub(".*_([^_]+)$", "\\1", rownames(results_df))

# Усреднение результатов по группам
average_scores <- results_df %>%
  group_by(Group) %>%
  summarise(across(everything(), mean, na.rm = TRUE))

library(ggplot2)
library(tidyr)

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

