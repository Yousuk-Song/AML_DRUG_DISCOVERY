library(readxl)
library(dplyr)
library(tidyr)
library(pheatmap)
library(RColorBrewer)
library(Rtsne)
library(ggplot2)
library(factoextra)

# 1. 엑셀 경로
file_path <- "/mnt/S1/sdata/processed_data/scRSEQ_AML/DRUG/blood_bld-2021-011094-mmc2.xlsx"
sheet_names <- excel_sheets(file_path)

# 2. 데이터 읽기 및 전처리
# all_data <- lapply(sheet_names, function(sheet) {
#   df <- read_excel(file_path, sheet = sheet)
#   df %>%
#     select(ID, DRUG_NAME, EC50) %>%
#     mutate(CellLine = sheet)
# }) %>% bind_rows()
# 2. 데이터 읽기 및 전처리
all_data <- lapply(sheet_names, function(sheet) {
  df <- read_excel(file_path, sheet = sheet)
  df %>%
    select(ID, DRUG_NAME, EC50, 'Mechanism/Targets') %>%  # ← 여기 추가
    mutate(CellLine = sheet)
}) %>% bind_rows()


wide_df <- all_data %>%
  filter(EC50 > 0) %>%
  mutate(logEC50 = log10(EC50)) %>%
  select(ID, CellLine, logEC50) %>%
  pivot_wider(names_from = CellLine, values_from = logEC50)

filtered_df <- wide_df %>%
  filter(rowSums(!is.na(.)) >= 2)

drug_mat <- as.data.frame(filtered_df)
rownames(drug_mat) <- drug_mat$ID
drug_mat <- drug_mat %>% select(-ID)
drug_mat <- drug_mat[apply(drug_mat, 1, function(x) sd(x, na.rm = TRUE) > 0), ]

# 3. 상관계수 행렬
cor_matrix <- cor(t(drug_mat), method = "pearson", use = "pairwise.complete.obs")

# 4. 클러스터링 (k=4, elbow는 시각적 판단 기반)
cor_dist <- as.dist(1 - cor_matrix)
optimal_k <- 4
hc_corr <- hclust(cor_dist, method = "complete")
clusters_corr <- cutree(hc_corr, k = optimal_k)

# 5. 클러스터 정보 저장
drug_clusters_corr <- data.frame(
  DrugID = rownames(cor_matrix),
  Cluster = as.factor(clusters_corr)
)
nrow(drug_clusters_corr)

# 6. cor_matrix 정렬: 클러스터 결과 순서 기준으로 정렬
cor_matrix <- cor_matrix[drug_clusters_corr$DrugID, drug_clusters_corr$DrugID]

# 7. 색상 지정 및 히트맵
cluster_colors <- brewer.pal(optimal_k, "Set1")
annotation_colors <- list(Cluster = setNames(cluster_colors, as.character(1:optimal_k)))

pheatmap(
  cor_matrix,
  clustering_distance_rows = cor_dist,
  clustering_distance_cols = cor_dist,
  clustering_method = "complete",
  annotation_row = drug_clusters_corr,
  annotation_colors = annotation_colors,
  show_rownames = FALSE,
  show_colnames = FALSE,
  main = paste("Drug Clustering with", optimal_k, "Clusters (Correlation-Based)")
)

# 8. t-SNE 시각화
tsne_corr <- Rtsne(as.matrix(cor_dist), is_distance = TRUE, perplexity = 30, verbose = TRUE)

tsne_df <- data.frame(
  X = tsne_corr$Y[,1],
  Y = tsne_corr$Y[,2],
  Cluster = drug_clusters_corr$Cluster
)

ggplot(tsne_df, aes(x = X, y = Y, color = Cluster)) +
  geom_point(size = 2, alpha = 0.7) +
  theme_minimal() +
  labs(title = paste("t-SNE of Drug Correlation (Distance-Based),", optimal_k, "Clusters"))

# 9. Enrichment 계산 함수 (DrugName 병합, Cluster 포함)
get_cluster_enrichment <- function(cluster_number, cor_matrix, cluster_df, drug_info) {
  cluster_ids <- cluster_df$DrugID[cluster_df$Cluster == cluster_number]
  submatrix <- cor_matrix[cluster_ids, cluster_ids]
  
  avg_corr <- rowMeans(submatrix, na.rm = TRUE)
  
  df <- data.frame(DrugID = cluster_ids, AvgCorrelation = avg_corr)
  df$Cluster <- cluster_number
  df <- df %>% left_join(drug_info, by = c("DrugID" = "ID"))
  
  list(
    enriched = df %>% arrange(desc(AvgCorrelation)) %>% head(min(10, nrow(df))),
    depleted = df %>% arrange(AvgCorrelation) %>% head(min(10, nrow(df))),
    all = df
  )
}

# 10. DrugID ↔ DrugName 매핑 준비
drug_info <- all_data %>%
  select(ID, DRUG_NAME) %>%
  distinct()

# 11. 클러스터별 Enrichment 결과 계산
cluster_enrichment_results <- lapply(1:optimal_k, function(k) {
  get_cluster_enrichment(k, cor_matrix, drug_clusters_corr, drug_info)
})
names(cluster_enrichment_results) <- paste0("Cluster_", 1:optimal_k)

# 12. 저장 함수 (Top10)
save_top10_enrichment_lists <- function(cluster_results, output_dir = "/mnt/S1/sdata/processed_data/scRSEQ_AML/DRUG/cluster_enrichment_output_top10") {
  dir.create(output_dir, showWarnings = FALSE)
  for (i in seq_along(cluster_results)) {
    cname <- names(cluster_results)[i]
    write.table(cluster_results[[cname]]$enriched, file = file.path(output_dir, paste0(cname, "_enriched_top10.txt")),
                sep = "\t", quote = FALSE, row.names = FALSE)
    write.table(cluster_results[[cname]]$depleted, file = file.path(output_dir, paste0(cname, "_depleted_top10.txt")),
                sep = "\t", quote = FALSE, row.names = FALSE)
  }
}

# 13. 저장 함수 (전체)
save_full_enrichment_lists <- function(cluster_results, output_dir = "/mnt/S1/sdata/processed_data/scRSEQ_AML/DRUG/cluster_enrichment_output_all") {
  dir.create(output_dir, showWarnings = FALSE)
  for (i in seq_along(cluster_results)) {
    cname <- names(cluster_results)[i]
    write.table(cluster_results[[cname]]$all %>% arrange(desc(AvgCorrelation)),
                file = file.path(output_dir, paste0(cname, "_enriched_all.txt")),
                sep = "\t", quote = FALSE, row.names = FALSE)
    write.table(cluster_results[[cname]]$all %>% arrange(AvgCorrelation),
                file = file.path(output_dir, paste0(cname, "_depleted_all.txt")),
                sep = "\t", quote = FALSE, row.names = FALSE)
  }
}

# 14. 파일 저장 실행
save_top10_enrichment_lists(cluster_enrichment_results)
save_full_enrichment_lists(cluster_enrichment_results)





library(jsonlite)
library(dplyr)

# 1. JSON 파일 경로
mechanism_json_path <- "/mnt/S1/sdata/processed_data/scRSEQ_AML/DRUG/final/grouped_drugs.json"

# 2. JSON 로드
raw_mech_data <- fromJSON(mechanism_json_path)



# JSON 구조를 간단한 DrugID ↔ MechanismGroup 형태로 변환
mech_df <- bind_rows(
  lapply(names(raw_mech_data), function(group) {
    mat <- raw_mech_data[[group]]
    if (!is.null(dim(mat)) && ncol(mat) >= 1) {
      data.frame(
        MechanismGroup = group,
        DrugID = mat[, 1],
        stringsAsFactors = FALSE
      )
    } else {
      NULL
    }
  })
)


head(mech_df)

# DrugID를 기준으로 mech_df와 cluster info 병합
drug_clusters_tagged <- drug_clusters_corr %>%
  left_join(mech_df, by = "DrugID")

library(dplyr)
library(tidyr)

mech_cluster_table <- drug_clusters_tagged %>%
  filter(!is.na(MechanismGroup)) %>%
  count(MechanismGroup, Cluster) %>%
  pivot_wider(
    names_from = Cluster,
    values_from = n,
    names_prefix = "Cluster_",
    values_fill = 0
  ) %>%
  arrange(desc(rowSums(across(starts_with("Cluster_")))))

write.csv(
  mech_cluster_table,
  "/mnt/S1/sdata/processed_data/scRSEQ_AML/DRUG/mechanism_vs_cluster_table.csv",
  row.names = FALSE
)







