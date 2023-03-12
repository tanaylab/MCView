# Generate PBMC metadata
library(anndata)
library(Matrix)
# library(MCView)
devtools::load_all()


# Assume we have the PBMC full.h5ad and cells.h5ad (see metacells vignette)

adata_full <- anndata::read_h5ad("raw/full.h5ad")
mito_genes <- grep("MT-.*", colnames(adata_full$X), value = TRUE)
mito_score <- rowSums(adata_full$X[, mito_genes]) %>% tibble::enframe("cell_id", "mito_score")
cell_size <- rowSums(adata_full$X) %>% tibble::enframe("cell_id", "cell_size")

metadata <- cell_metadata_to_metacell_from_h5ad("raw/cells.h5ad", c("batch"), categorical = "batch")

adata <- anndata::read_h5ad("raw/cells.h5ad")

cell_to_metacell <- adata$obs %>%    
    rownames_to_column("cell_id") %>%
    select(cell_id, metacell = metacell_name)


metadata1 <- cell_metadata_to_metacell(cell_size %>% left_join(mito_score), cell_to_metacell, func=sum)

metadata <- metadata1 %>% left_join(metadata)

fwrite(metadata, "raw/metadata.csv")

