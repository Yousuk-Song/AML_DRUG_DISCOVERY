library(readxl)
library(dplyr)
library(tidyr)
library(pheatmap)
library(RColorBrewer)
library(ggplot2)

# 1. 엑셀 파일 경로
file_path <- "/mnt/S1/sdata/processed_data/scRSEQ_AML/DRUG/blood_bld-2021-011094-mmc2.xlsx"
sheet_names <- excel_sheets(file_path)

# 2. 데이터 읽기 및 전처리
all_data <- lapply(sheet_names, function(sheet) {
  df <- read_excel(file_path, sheet = sheet)
  df %>%
    select(ID, DRUG_NAME, EC50, 'Mechanism/Targets') %>%
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

# 3. 상관계수 행렬 (Drug 간 Pearson correlation)
cor_matrix <- cor(t(drug_mat), method = "pearson", use = "pairwise.complete.obs")

# 4. pheatmap 실행 (1차) → tree_row 추출을 위한 dendrogram 계산용
pheatmap_result <- pheatmap(
  cor_matrix,
  clustering_method = "complete",
  show_rownames = FALSE,
  show_colnames = FALSE,
  main = "Drug Clustering Preview (correlation-based)"
)

# 5. pheatmap의 dendrogram으로부터 클러스터 추출
optimal_k <- 4
tree_row <- pheatmap_result$tree_row
clusters <- cutree(tree_row, k = optimal_k)

# 6. 클러스터 정보 만들기
drug_clusters_corr <- data.frame(
  DrugID = names(clusters),
  Cluster = as.factor(clusters)
)

# 7. 색상 지정
cluster_colors <- brewer.pal(optimal_k, "Set1")
annotation_colors <- list(Cluster = setNames(cluster_colors, as.character(1:optimal_k)))

# 8. pheatmap 재실행 (annotation 반영)
# 정렬 순서에 맞게 correlation matrix 재정렬
cor_matrix <- cor_matrix[drug_clusters_corr$DrugID, drug_clusters_corr$DrugID]

pheatmap(
  cor_matrix,
  clustering_method = "complete",
  annotation_row = drug_clusters_corr,
  annotation_colors = annotation_colors,
  show_rownames = FALSE,
  show_colnames = FALSE,
  main = paste("Drug Clustering with", optimal_k, "Clusters (cutree from pheatmap)")
)




run_correlation_clustering <- function(cor_matrix,
                                       clustering_method = "complete",
                                       k = 4,
                                       show_names = FALSE,
                                       palette = "Set1") {
  require(pheatmap)
  require(RColorBrewer)
  
  # 1. 초기 pheatmap 실행 (dendrogram 생성용)
  pheatmap_result <- pheatmap(
    cor_matrix,
    clustering_method = clustering_method,
    show_rownames = show_names,
    show_colnames = show_names,
    silent = TRUE  # 결과만 받고 시각화는 나중에
  )
  
  # 2. cutree로 클러스터 추출
  tree <- pheatmap_result$tree_row
  clusters <- cutree(tree, k = k)
  
  # 3. 클러스터 정보 정리
  cluster_df <- data.frame(
    Sample = names(clusters),
    Cluster = as.factor(clusters)
  )
  
  # 4. 색상 설정
  cluster_colors <- brewer.pal(k, palette)
  annotation_colors <- list(Cluster = setNames(cluster_colors, as.character(1:k)))
  
  # 5. cor_matrix 순서 정렬
  cor_matrix_ordered <- cor_matrix[cluster_df$Sample, cluster_df$Sample]
  
  # 6. 최종 pheatmap 실행 (annotation 반영)
  pheatmap(
    cor_matrix_ordered,
    clustering_method = clustering_method,
    annotation_row = cluster_df,
    annotation_colors = annotation_colors,
    show_rownames = show_names,
    show_colnames = show_names,
    main = paste("Correlation-Based Clustering (method =", clustering_method, ", k =", k, ")")
  )
  
  # 7. 클러스터 결과 리턴
  return(cluster_df)
}


# i) complete = Complete linkage: 두 클러스터 사이의 가장 먼 거리 사용 (최대 거리) → 뭉친 클러스터 강조
result_df_complete <- run_correlation_clustering(
  cor_matrix = cor_matrix,
  clustering_method = "complete",
  k = 4
)

# ii) average = Average linkage (UPGMA): 두 클러스터 간 평균 거리 사용 → 균형 잡힌 결과
result_df_average <- run_correlation_clustering(
  cor_matrix = cor_matrix,
  clustering_method = "average",
  k = 4
)

# iii) single = Single linkage: 두 클러스터 사이의 가장 가까운 거리 사용 (최소 거리) → chaining effect 발생 가능
result_df_single <- run_correlation_clustering(
  cor_matrix = cor_matrix,
  clustering_method = "single",
  k = 4
)
# iv) ward.D = Ward's method: 분산 증가 최소화. 유클리디안 거리 전용
result_df_ward.D <- run_correlation_clustering(
  cor_matrix = cor_matrix,
  clustering_method = "ward.D",
  k = 4
)
# v) ward.D2 = Ward의 개선 버전. 역시 유클리디안 거리만 의미 있음
result_df_ward.D2 <- run_correlation_clustering(
  cor_matrix = cor_matrix,
  clustering_method = "ward.D2",
  k = 4
)

# vi) centroid = 클러스터 중심 간 거리 → 극단적인 경우 잘 안 씀
result_df_centroid <- run_correlation_clustering(
  cor_matrix = cor_matrix,
  clustering_method = "centroid",
  k = 4
)

# vii) median = 중간 중심값 기반 (rarely used)
result_df_median <- run_correlation_clustering(
  cor_matrix = cor_matrix,
  clustering_method = "median",
  k = 4
)
