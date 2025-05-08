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
seurat_obj@misc$MS_PROTEOMICS   <- MS_mat[intersect(rownames(MS_mat), rownames(seurat_obj@meta.data)), , drop = FALSE]
seurat_obj@misc$PPRA_PROTEOMICS <- PPRA_mat[intersect(rownames(PPRA_mat), rownames(seurat_obj@meta.data)), , drop = FALSE]


#============================#
# 3. Genetic Perturbation: RNAi & CRISPR
#============================#

# Load RNAi and CRISPR
rnai_matrix <- mat(parse_gctx("/data/processed_data/scRSEQ_AML/DRUG/CCLE/RNAi_Screens/Achilles_QC_v2.4.3.rnai.Gs.gct"))
crispr_matrix <- mat(parse_gctx("/data/processed_data/scRSEQ_AML/DRUG/CCLE/CRISPR_Screens/Achilles_v3.3.8.Gs.gct"))

# Common cell lines
common_cell_lines <- Reduce(intersect, list(
  colnames(seurat_obj), colnames(rnai_matrix), colnames(crispr_matrix)))

# Subset
rnai_df <- t(rnai_matrix[, common_cell_lines])
crispr_df <- t(crispr_matrix[, common_cell_lines])

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

# Common genes
common_genes <- Reduce(intersect, list(
  rownames(seurat_obj@assays$RNA@counts),
  colnames(rnai_df_avg),
  colnames(crispr_df_avg)
))

# Save into Seurat
seurat_obj@misc$CRISPR_SCORES <- crispr_df_avg[rownames(seurat_obj@meta.data), common_genes]
seurat_obj@misc$RNAi_SCORES   <- rnai_df_avg[rownames(seurat_obj@meta.data), common_genes]


#============================#
# 4. UniProt Mapping & Query Function
#============================#

# Gene ↔ UniProt mapping
map_genes_to_uniprot <- function(gene_symbols, map_table) {
  filtered <- map_table[map_table$hgnc_symbol %in% gene_symbols, ] %>% 
    filter(!duplicated(hgnc_symbol))
  result <- setNames(filtered$uniprotswissprot, filtered$hgnc_symbol)
  result[gene_symbols]
}

# Query function
get_multi_omics_gene_values <- function(seurat_obj, cell_line, gene_symbol, map_table = gene_uniprot_map) {
  uniprot_id <- map_genes_to_uniprot(gene_symbol, map_table)[[1]]
  if (is.na(uniprot_id)) uniprot_id <- NULL

  list(
    RNA_count     = tryCatch(seurat_obj@assays$RNA@counts[gene_symbol, cell_line], error = function(e) NA),
    CRISPR_score  = tryCatch(seurat_obj@misc$CRISPR_SCORES[cell_line, gene_symbol], error = function(e) NA),
    RNAi_score    = tryCatch(seurat_obj@misc$RNAi_SCORES[cell_line, gene_symbol], error = function(e) NA),
    MS_protein    = tryCatch(if (!is.null(uniprot_id)) seurat_obj@misc$MS_PROTEOMICS[cell_line, uniprot_id] else NA, error = function(e) NA),
    PPRA_phospho  = tryCatch(if (!is.null(uniprot_id)) seurat_obj@misc$PPRA_PROTEOMICS[cell_line, uniprot_id] else NA, error = function(e) NA),
    UniProt_ID    = uniprot_id
  )
}

# Example usage:
# get_multi_omics_gene_values(seurat_obj, "CAL51_BREAST", "TP53")
