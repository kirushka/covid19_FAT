## Load libraries --------------------------------------------------------------
library(tidyverse)
library(ggtext)
# devtools::install_version("rvcheck", version = "0.1.8") # makes clusterprofiler v3.18.1 work with rvcheck v0.1.8
library(clusterProfiler)
source("revigo_UPD.R")


## Load data with DEGs ---------------------------------------------------------
## Melms dataset was manually fixed 
degs_files <- list.files("SCovid/DEGs", full.names = TRUE)

degs <- map_dfr(degs_files, read_csv)

### Lipid-related GOBP categories
lipid_terms <- c("lipid", "fat", "triglyceride", "cholesterol")
nonlipid_terms <- c("fate", "sulfat")


## Convert SYMBOL to ENTREZ IDs ------------------------------------------------
degs_lng <- degs %>% 
  mutate(UpDown = if_else(log2FC > 0, "Upregulated", "Downregulated")) %>% 
  dplyr::select(Dataset, Tissue, CellType = `Cell Type`, SYMBOL = GENE, UpDown) 

degs_entrez <- bitr(degs_lng$SYMBOL, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db") %>% 
  filter(!duplicated(ENTREZID), !duplicated(SYMBOL)) 

degs_lng <- degs_lng %>% 
  left_join(degs_entrez, by = "SYMBOL")

# write_tsv(degs_lng, "SCovid/degs_combined_table.tsv")

datasets <- degs_lng %>% 
  dplyr::select(Dataset, Tissue) %>% 
  distinct()

## GO BP Enrichment by celltype ------------------------------------------------
go_by_celltype <- compareCluster(
  ENTREZID ~ UpDown + CellType + Dataset, 
  data = degs_lng, 
  fun="enrichGO", 
  OrgDb = "org.Hs.eg.db", 
  ont = "BP", 
  readable = TRUE)

saveRDS(go_by_celltype, "all_go_bp.rds")

go_by_celltype_res <- as_tibble(go_by_celltype@compareClusterResult) %>%
  # Add Tissue info
  left_join(datasets, by = "Dataset")

go_by_celltype_lipid <- go_by_celltype_res %>% 
  # Keep only GO BP terms related to lipids metabolism
  dplyr::filter(str_detect(Description, paste(lipid_terms, collapse = "|"))) %>%
  # Remove some GO BP terms
  dplyr::filter(!str_detect(Description, paste(nonlipid_terms, collapse = "|")))

go_by_celltype_lipid_revigo_res <- go_by_celltype_lipid %>% 
  dplyr::select(ID, qvalue) %>%
  arrange(qvalue) %>% distinct(ID, .keep_all = TRUE) %>% 
  revigo(cutoff = "small") %>% 
  dplyr::filter(str_detect(Eliminated, "False"))


go_by_celltype_lipid_revigo <- go_by_celltype_lipid %>% 
  dplyr::filter(ID %in% go_by_celltype_lipid_revigo_res$`Term ID`)


