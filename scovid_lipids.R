## Load libraries --------------------------------------------------------------
library(tidyverse)
library(ggtext)
# devtools::install_version("rvcheck", version = "0.1.8") # makes clusterprofiler v3.18.1 work with rvcheck v0.1.8
library(clusterProfiler)
source("revigo.R")


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

## Make nested column for dataset:tissue:UpDown groupping ----------------------
degs_by_celltype_tissue <- degs_lng %>% 
  group_by(CellType, Dataset, Tissue, UpDown) %>% 
  nest() %>% 
  ungroup() %>% 
  mutate(CellType_Tissue_Dataset_UpDown = str_c(CellType, Tissue, Dataset, UpDown, sep = "//"))


## Enrichment functions  -------------------------------------------------------
### GO BP Overrepresentation Analysis
run_enrichGO <- function(data){
  data$go_bp = vector(mode = "list", length = nrow(data))
  for (i in 1:nrow(data)){
    data$go_bp[[i]] = enrichGO(gene = data$data[[i]]$ENTREZID, OrgDb = "org.Hs.eg.db", ont = "BP")
    data$go_bp[[i]] = setReadable(data$go_bp[[i]], OrgDb = "org.Hs.eg.db")
  }
  return(data)
}


### KEGG Overrepresentation Analysis
run_enrichKEGG <- function(data){
  data$kegg = vector(mode = "list", length = nrow(data))
  for (i in 1:nrow(data)){
    data$kegg[[i]] = enrichKEGG(gene = data$data[[i]]$ENTREZID, organism = "hsa")
    data$kegg[[i]] = setReadable(data$kegg[[i]], OrgDb = "org.Hs.eg.db", keyType = "ENTREZID")
  }
  return(data)
}



## GO BP Enrichment by celltype ------------------------------------------------
go_by_celltype <- run_enrichGO(degs_by_celltype_tissue)

saveRDS(go_by_celltype, "SCovid/go_by_celltype.rds")
go_by_celltype <- readRDS("SCovid/go_by_celltype.rds")

go_by_celltype_res <- go_by_celltype$go_bp %>% 
  map("result") %>% 
  map2_dfr(degs_by_celltype_tissue$CellType_Tissue_Dataset_UpDown, ., ~ add_column(.y, CellType_Tissue_Dataset_UpDown = .x)) %>% 
  separate(CellType_Tissue_Dataset_UpDown, into = c("CellType", "Tissue", "Dataset", "UpDown"), sep = "//")

go_by_celltype_res_flt <- go_by_celltype_res %>% 
  # Keep only significant GO BP terms
  dplyr::filter(qvalue < 0.05) 

### Filter GO BP categories ----------------------------------------------------
go_lipid <- go_by_celltype_res_flt %>%
  # Keep only GO BP terms related to lipids metabolism
  dplyr::filter(str_detect(Description, paste(lipid_terms, collapse = "|"))) %>%
  # Remove some GO BP terms
  dplyr::filter(!str_detect(Description, paste(nonlipid_terms, collapse = "|")))

go_lipid_revigo_res <- go_lipid %>% 
  dplyr::select(ID, qvalue) %>% 
  arrange(qvalue) %>% distinct(ID, .keep_all = TRUE) %>% 
  revigo(cutoff = "small") %>% 
  dplyr::filter(str_detect(Eliminated, "False"))

go_lipid_revigo <- go_lipid %>% 
  dplyr::filter(ID %in% go_lipid_revigo_res$`Term ID`)

write_tsv(go_lipid_revigo, "SCovid/go_lipid_revigo.tsv")


## GO BP - Half-circles plot ---------------------------------------------------
tissues <- c("Lung", "Airway", "Liver", "Blood", "Lymph node", "Brain", "Heart", "Intestine", "Pancreas", "Kidney")

go_lipid_revigo_2plot <- go_lipid_revigo %>% 
  mutate(
    label = if_else(UpDown == "Upregulated", "\u25D7", "\u25D6"),
    color = if_else(UpDown == "Upregulated", "#f4a261", "#0a9396")) %>% 
  add_count(ID) %>% 
  mutate(
    Tissue = factor(Tissue, levels = tissues),
    Description = factor(Description),
    Description = fct_reorder(Description, n))

min_size <- min(go_lipid_revigo_2plot$Count)
max_size <- max(go_lipid_revigo_2plot$Count)

### By Tissue
ggplot(go_lipid_revigo_2plot, 
       aes(x = Tissue, 
           y = Description)) +
  geom_text(aes(color = color, label = label, size = Count*2), family = "Arial Unicode MS", key_glyph = "point") +
  scale_y_discrete(expand = c(0, 0, 0.03, 0)) +
  scale_color_identity() +
  scale_size_identity(
    "Count", 
    breaks = c(min_size:max_size)[c(FALSE, TRUE)],
    labels = c(min_size:max_size)[c(FALSE, TRUE)],
    guide = "legend") +
  # guides(size = guide_legend(override.aes = list(color = c(rep(c("#f4a261", "#0a9396"), 2), "#f4a261")))) +
  guides(size = guide_legend(override.aes = list(color = "#0a9396"))) +
  labs(
    x = "", y = "", 
    title = "<b style='color:#0a9396'>DOWNREGULATED</b> and <b style='color:#f4a261'>UPREGULATED</b> categories") +
  theme_bw() +
  theme(
    # legend.position = "top",
    legend.justification = c("right"),
    plot.title = element_markdown(lineheight = 1.1, hjust = 1),
    axis.text = element_text(color = "black", size = rel(0.9)),
    axis.text.x = element_text(angle = 35, hjust = 1))  

ggsave("SCovid/GOBP_by_tissue_lipid_keywords.png", dpi = 300, width = 10, height = 10, units = "cm", scale = 2)

### By Dataset - only Lung and Airway
datasets_order <- c("Delorey TM. (Lung)", "Melms JC. (Lung)", "Xu G. (Lung)", "Liao M. (BALF)", "Delorey TM. (Airway)", "Delorey TM. (Trachea)", "Misharin A. (Bronchia)")

go_lipid_revigo_2plot %>% 
  dplyr::filter(Tissue %in% c("Lung", "Airway")) %>% 
  mutate(Dataset = factor(Dataset, levels = datasets_order)) %>% 
  ggplot(aes(x = Dataset, y = Description)) +
  geom_text(aes(color = color, label = label, size = Count*2), family = "Arial Unicode MS", key_glyph = "point") +
  scale_color_identity() +
  scale_size_identity(
    "Count", 
    breaks = c(min_size:max_size)[c(FALSE, TRUE)],
    labels = c(min_size:max_size)[c(FALSE, TRUE)],
    guide = "legend") +
  # guides(size = guide_legend(override.aes = list(color = c(rep(c("#f4a261", "#0a9396"), 2), "#f4a261")))) +
  guides(size = guide_legend(override.aes = list(color = "#0a9396"))) +
  labs(
    x = "", y = "", 
    title = "<b style='color:#0a9396'>DOWNREGULATED</b> and <b style='color:#f4a261'>UPREGULATED</b> categories") +
  theme_bw() +
  theme(
    # legend.position = "top",
    legend.justification = c("right"),
    plot.title = element_markdown(lineheight = 1.1, hjust = 1),
    axis.text = element_text(color = "black", size = rel(0.9)),
    axis.text.x = element_text(angle = 35, hjust = 1))  

ggsave("SCovid/GOBP_by_dataset_lipid_keywords.png", dpi = 300, width = 8, height = 9, units = "cm", scale = 2)

### By Cell Type - only Lung and Airway
go_lipid_revigo_2plot %>% 
  dplyr::filter(Tissue %in% c("Lung", "Airway")) %>% 
  # mutate(Dataset = factor(Dataset, levels = datasets_order)) %>% 
  ggplot(aes(x = CellType, y = Description)) +
  geom_text(aes(color = color, label = label, size = Count*2), family = "Arial Unicode MS", key_glyph = "point") +
  scale_color_identity() +
  scale_size_identity(
    "Count", 
    breaks = c(min_size:max_size)[c(FALSE, TRUE)],
    labels = c(min_size:max_size)[c(FALSE, TRUE)],
    guide = "legend") +
  # guides(size = guide_legend(override.aes = list(color = c(rep(c("#f4a261", "#0a9396"), 2), "#f4a261")))) +
  guides(size = guide_legend(override.aes = list(color = "#0a9396"))) +
  labs(
    x = "", y = "", 
    title = "<b style='color:#0a9396'>DOWNREGULATED</b> and <b style='color:#f4a261'>UPREGULATED</b> categories") +
  theme_bw() +
  theme(
    # legend.position = "top",
    legend.justification = c("right"),
    plot.title = element_markdown(lineheight = 1.1, hjust = 1),
    axis.text = element_text(color = "black", size = rel(0.9)),
    axis.text.x = element_text(angle = 35, hjust = 1))  

ggsave("SCovid/GOBP_by_celltype_lipid_keywords.png", dpi = 300, width = 11, height = 9, units = "cm", scale = 2)


## GO BP - Genes table ---------------------------------------------------------
go_lipid_revigo_2table <- go_lipid_revigo %>% 
  dplyr::select(ID, Description, geneID, UpDown, Tissue, CellType, Dataset) %>% 
  arrange(factor(Tissue, levels = tissues), UpDown, CellType) %>% 
  separate_rows(geneID)

write_tsv(go_lipid_revigo_2table, "SCovid/GOBP_lipid_genes.tsv")


## KEGG Enrichment by tissue ---------------------------------------------------
kegg_by_celltype <- run_enrichKEGG(degs_by_celltype_tissue)

saveRDS(kegg_by_celltype, "SCovid/kegg_by_celltype.rds")
# kegg_by_celltype <- readRDS("SCovid/kegg_by_celltype.rds")

kegg_by_celltype_res <- kegg_by_celltype$kegg %>% 
  map("result") %>% 
  map2_dfr(degs_by_celltype_tissue$CellType_Tissue_Dataset_UpDown, ., ~ add_column(.y, CellType_Tissue_Dataset_UpDown = .x)) %>% 
  separate(CellType_Tissue_Dataset_UpDown, into = c("CellType", "Tissue", "Dataset", "UpDown"), sep = "//")

kegg_by_celltype_res_flt <- kegg_by_celltype_res %>% 
  # Keep only significant GO BP terms
  dplyr::filter(qvalue < 0.05) 

kegg_by_celltype_lipid <- kegg_by_celltype_res_flt %>%
  # Keep only GO BP terms related to lipids metabolism
  dplyr::filter(str_detect(Description, paste(lipid_terms, collapse = "|"))) %>%
  # Remove some GO BP terms
  dplyr::filter(!str_detect(Description, paste(nonlipid_terms, collapse = "|")))

write_tsv(kegg_by_celltype_res_flt, "SCovid/kegg_by_celltype_res_flt.tsv")
write_tsv(kegg_by_celltype_lipid, "SCovid/kegg_by_celltype_lipid.tsv")


## KEGG - Half-circles plot ----------------------------------------------------
kegg_by_celltype_lipid_2plot <- kegg_by_celltype_lipid %>% 
  mutate(
    label = if_else(UpDown == "Upregulated", "\u25D7", "\u25D6"),
    color = if_else(UpDown == "Upregulated", "#f4a261", "#0a9396"))

min_size <- min(kegg_by_celltype_lipid_2plot$Count)
max_size <- max(kegg_by_celltype_lipid_2plot$Count)
tissues <- c("Lung", "Airway", "Liver", "Blood", "Lymph node", "Brain", "Heart", "Intestine", "Pancreas", "Kidney")

ggplot(kegg_by_celltype_lipid_2plot, 
       aes(x = factor(Tissue, levels = tissues), 
           y = Description)) +
  geom_text(aes(color = color, label = label, size = Count*2), family = "Arial Unicode MS", key_glyph = "point") +
  scale_color_identity() +
  scale_size_identity(
    "Count", 
    breaks = c(min_size:max_size)[c(FALSE, TRUE)],
    labels = c(min_size:max_size)[c(FALSE, TRUE)],
    guide = "legend") +
  # guides(size = guide_legend(override.aes = list(color = c(rep(c("#f4a261", "#0a9396"), 2), "#f4a261")))) +
  guides(size = guide_legend(override.aes = list(color = "#0a9396"))) +
  labs(
    x = "", y = "", 
    title = "<b style='color:#0a9396'>DOWNREGULATED</b> and <b style='color:#f4a261'>UPREGULATED</b> KEGG pathways") +
  theme_bw() +
  theme(
    # legend.position = "top",
    legend.justification = c("right"),
    plot.title = element_markdown(lineheight = 1.1, hjust = 1),
    axis.text = element_text(color = "black", size = rel(0.9)),
    axis.text.x = element_text(angle = 35, hjust = 1))  

ggsave("SCovid/KEGG_by_tissue_lipid_keywords.png", dpi = 300, width = 8, height = 4, units = "cm", scale = 2)


## KEGG - Genes table ----------------------------------------------------------
### Genes associated with lipid metabolism
kegg_by_celltype_lipid_2table <- kegg_by_celltype_lipid %>% 
  dplyr::select(ID, Description, geneID, UpDown, Tissue) %>% 
  separate_rows(geneID)

write_tsv(kegg_by_celltype_lipid_2table, "SCovid/kegg_by_tissue_lipid_genes.tsv")

### All enriched KEGG pathways
kegg_by_celltype_res_flt_2table <- kegg_by_celltype_res_flt %>% 
  dplyr::select(ID, Description, Count, geneID, UpDown, Tissue, CellType)

write_tsv(kegg_by_celltype_res_flt_2table, "SCovid/kegg_by_tissue_pathways.tsv")

