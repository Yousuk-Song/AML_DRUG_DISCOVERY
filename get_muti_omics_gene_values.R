
# Required Libraries
library(Seurat)
library(dplyr)
library(tidyr)
library(readr)
library(biomaRt)
library(cmapR)

#============================#
# 1. scRNA-seq Expression
#============================#

# Load GCT file
expr_matrix <- mat(parse_gctx("/data/processed_data/scRSEQ_AML/DRUG/CCLE/Sequencing/CCLE_RNAseq_genes_counts_20180929.gct"))

# Create Seurat Object
seurat_obj <- CreateSeuratObject(counts = expr_matrix)

# Map Ensembl ID to HGNC Symbol
ensembl_ids <- sub("\\..*$", "", rownames(expr_matrix))
mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
map_tbl <- getBM(
  attributes = c("ensembl_gene_id", "hgnc_symbol"),
  filters = "ensembl_gene_id",
  values = ensembl_ids,
  mart = mart
) %>% 
  filter(hgnc_symbol != "" & !is.na(hgnc_symbol)) %>% 
  distinct(ensembl_gene_id, .keep_all = TRUE)

# Apply mapping
rownames(expr_matrix) <- map_tbl$hgnc_symbol[match(ensembl_ids, map_tbl$ensembl_gene_id)]
expr_matrix <- expr_matrix[!is.na(rownames(expr_matrix)) & !duplicated(rownames(expr_matrix)), ]
seurat_obj[["RNA"]] <- CreateAssayObject(counts = expr_matrix)

# 1. Normalize (LogNormalize)
seurat_obj <- NormalizeData(
  seurat_obj,
  normalization.method = "LogNormalize",
  scale.factor = 10000
)

# 2. Identify variable features
seurat_obj <- FindVariableFeatures(
  seurat_obj,
  selection.method = "vst",
  nfeatures = 2000
)

# 3. Scale data (z-score normalization across cells)
seurat_obj <- ScaleData(seurat_obj)

# 4. PCA
seurat_obj <- RunPCA(seurat_obj, features = VariableFeatures(seurat_obj))

# 5. Neighborhood graph
seurat_obj <- FindNeighbors(seurat_obj, dims = 1:20)

# 6. Clustering
seurat_obj <- FindClusters(seurat_obj, resolution = 0.5)

# 7. UMAP
seurat_obj <- RunUMAP(seurat_obj, dims = 1:20)

#============================#
# 2. Protein: MS & PPRA
#============================#

# Load mapping table
model_map <- read.csv("/data/processed_data/scRSEQ_AML/DRUG/CCLE/Proteomics/Model.csv", stringsAsFactors = FALSE)

# MS Proteomics
MS <- read_csv("/data/processed_data/scRSEQ_AML/DRUG/CCLE/Proteomics/harmonized_MS_CCLE_Gygi.csv")
MS <- as.data.frame(MS)
rownames(MS) <- MS$...1
MS <- MS[, -which(colnames(MS) == "...1")]
MS_mat <- as.matrix(MS)
rownames(MS_mat)[!is.na(match(rownames(MS_mat), model_map$ModelID))] <- 
  model_map$CCLEName[match(rownames(MS_mat), model_map$ModelID)]

# PPRA Proteomics
PPRA <- read.csv("/data/processed_data/scRSEQ_AML/DRUG/CCLE/Proteomics/harmonized_RPPA_CCLE.csv")
PPRA <- as.data.frame(PPRA)
rownames(PPRA) <- PPRA$X
PPRA <- PPRA[, -which(colnames(PPRA) == "X")]
PPRA_mat <- as.matrix(PPRA)
rownames(PPRA_mat)[!is.na(match(rownames(PPRA_mat), model_map$ModelID))] <- 
  model_map$CCLEName[match(rownames(PPRA_mat), model_map$ModelID)]

# Store in Seurat
#seurat_obj@misc$MS_PROTEOMICS   <- MS_mat[intersect(rownames(MS_mat), rownames(seurat_obj@meta.data)), , drop = FALSE]
#seurat_obj@misc$PPRA_PROTEOMICS <- PPRA_mat[intersect(rownames(PPRA_mat), rownames(seurat_obj@meta.data)), , drop = FALSE]
seurat_obj@misc$MS_PROTEOMICS   <- MS_mat
seurat_obj@misc$PPRA_PROTEOMICS <- PPRA_mat

#============================#
# 3. Genetic Perturbation: RNAi & CRISPR
#============================#

# Load RNAi and CRISPR
rnai_matrix <- mat(parse_gctx("/data/processed_data/scRSEQ_AML/DRUG/CCLE/RNAi_Screens/Achilles_QC_v2.4.3.rnai.Gs.gct"))
crispr_matrix <- mat(parse_gctx("/data/processed_data/scRSEQ_AML/DRUG/CCLE/CRISPR_Screens/Achilles_v3.3.8.Gs.gct"))

# Common cell lines
#common_cell_lines <- Reduce(intersect, list(
#  colnames(seurat_obj), colnames(rnai_matrix), colnames(crispr_matrix)))

# Subset
rnai_df <- t(rnai_matrix)
crispr_df <- t(crispr_matrix)

# Extract gene names
rnai_genes <- sapply(strsplit(colnames(rnai_df), "_"), `[`, 1)
crispr_genes <- sapply(strsplit(colnames(crispr_df), "_"), `[`, 1)

# Clean & average
average_by_gene <- function(mat, genes) {
  df <- as.data.frame(mat)
  colnames(df) <- genes
  sapply(unique(genes), function(g) {
    rowMeans(df[, which(genes == g), drop = FALSE], na.rm = TRUE)
  }) %>% as.data.frame() %>% `rownames<-`(rownames(df))
}

rnai_df_avg   <- average_by_gene(rnai_df, rnai_genes)
crispr_df_avg <- average_by_gene(crispr_df, crispr_genes)

# CRISPR: Seurat와 cell line 겹치는 것만 필터링
#common_crispr_cells <- intersect(rownames(seurat_obj@meta.data), rownames(crispr_df_avg))
#seurat_obj@misc$CRISPR_SCORES <- as.matrix(crispr_df_avg[common_crispr_cells, , drop = FALSE])
seurat_obj@misc$CRISPR_SCORES <- as.matrix(crispr_df_avg)

# RNAi: Seurat와 cell line 겹치는 것만 필터링
#common_rnai_cells <- intersect(rownames(seurat_obj@meta.data), rownames(rnai_df_avg))
#seurat_obj@misc$RNAi_SCORES   <- as.matrix(rnai_df_avg[common_rnai_cells, , drop = FALSE])
seurat_obj@misc$RNAi_SCORES   <- as.matrix(rnai_df_avg)


#============================#
# 4. UniProt Mapping & Query Function
#============================#

# Gene ↔ UniProt mapping
# biomaRt 연결
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# HGNC → UniProt 매핑
gene_uniprot_map <- getBM(
  attributes = c("hgnc_symbol", "uniprotswissprot"),
  mart = ensembl
) %>%
  filter(
    hgnc_symbol != "",
    uniprotswissprot != "",
    !is.na(hgnc_symbol),
    !is.na(uniprotswissprot)
  ) %>%
  distinct(hgnc_symbol, .keep_all = TRUE)


# 먼저 gene symbol → UniProt 벡터 매핑 함수
map_genes_to_uniprot <- function(gene_symbols, map_table) {
  filtered <- map_table[map_table$hgnc_symbol %in% gene_symbols, ]
  filtered <- filtered[!duplicated(filtered$hgnc_symbol), ]
  result <- setNames(filtered$uniprotswissprot, filtered$hgnc_symbol)
  result <- result[gene_symbols]
  names(result) <- gene_symbols
  return(result)
}

# 통합된 멀티오믹스 조회 함수
get_multi_omics_gene_values <- function(seurat_obj, cell_line, gene_symbol, map_table = gene_uniprot_map) {
  # 1. gene symbol → UniProt ID
  uniprot_id <- map_genes_to_uniprot(gene_symbol, map_table)[[1]]
  if (is.na(uniprot_id) || is.null(uniprot_id)) uniprot_id <- NULL
  
  # 2. RNA raw counts
  rna_raw <- tryCatch({
    seurat_obj@assays$RNA@counts[gene_symbol, cell_line]
  }, error = function(e) NA)
  
  # 3. RNA log-normalized
  rna_norm <- tryCatch({
    seurat_obj@assays$RNA@data[gene_symbol, cell_line]
  }, error = function(e) NA)
  
  # 4. CRISPR score (matrix)
  crispr_val <- tryCatch({
    if (!is.null(seurat_obj@misc$CRISPR_SCORES) && 
        gene_symbol %in% colnames(seurat_obj@misc$CRISPR_SCORES) &&
        cell_line %in% rownames(seurat_obj@misc$CRISPR_SCORES)) {
      seurat_obj@misc$CRISPR_SCORES[cell_line, gene_symbol]
    } else NA
  }, error = function(e) NA)
  
  # 5. RNAi score (matrix)
  rnai_val <- tryCatch({
    if (!is.null(seurat_obj@misc$RNAi_SCORES) && 
        gene_symbol %in% colnames(seurat_obj@misc$RNAi_SCORES) &&
        cell_line %in% rownames(seurat_obj@misc$RNAi_SCORES)) {
      seurat_obj@misc$RNAi_SCORES[cell_line, gene_symbol]
    } else NA
  }, error = function(e) NA)
  
  # 6. MS proteomics
  ms_val <- tryCatch({
    if (!is.null(uniprot_id) && 
        uniprot_id %in% colnames(seurat_obj@misc$MS_PROTEOMICS) &&
        cell_line %in% rownames(seurat_obj@misc$MS_PROTEOMICS)) {
      seurat_obj@misc$MS_PROTEOMICS[cell_line, uniprot_id]
    } else NA
  }, error = function(e) NA)
  
  # 7. PPRA proteomics
  ppra_val <- tryCatch({
    if (!is.null(uniprot_id) && 
        uniprot_id %in% colnames(seurat_obj@misc$PPRA_PROTEOMICS) &&
        cell_line %in% rownames(seurat_obj@misc$PPRA_PROTEOMICS)) {
      seurat_obj@misc$PPRA_PROTEOMICS[cell_line, uniprot_id]
    } else NA
  }, error = function(e) NA)
  
  # 8. Output
  return(list(
    RNA_raw_counts     = rna_raw,
    RNA_logNormalized  = rna_norm,
    CRISPR_score       = crispr_val,
    RNAi_score         = rnai_val,
    MS_protein         = ms_val,
    PPRA_phospho       = ppra_val,
    UniProt_ID         = uniprot_id
  ))
}

# Example usage:
get_multi_omics_gene_values(seurat_obj, "OCIAML5_HAEMATOPOIETIC_AND_LYMPHOID_TISSUE", "NRAS")
get_multi_omics_gene_values(seurat_obj, "CAL120_BREAST", "TP53")


#saveRDS(seurat_obj, "/data/processed_data/scRSEQ_AML/DRUG/CCLE/CCLE_multi_omics_data.rds")


