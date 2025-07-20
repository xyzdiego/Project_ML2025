library(Seurat)
library(tidyverse)
library(scDblFinder)
library(SingleCellExperiment)
library(SingleR)

txt_files <- list.files(
    path        = ".", 
    pattern     = "\\.txt$",
    full.names  = TRUE,
    ignore.case = TRUE
)
lista_txt <- lapply(txt_files, read.delim, header = TRUE, stringsAsFactors = FALSE)
names(lista_txt) <- basename(txt_files)
merged_df <- Reduce(function(x, y) merge(x, y, by = "gene", all = TRUE),
                    lista_txt)
merged_df[,-1] <- lapply(merged_df[,-1], function(x) ifelse(is.na(x), 0, x))
longitudes <- c()
for(i in 1:length(lista_txt)){
    longitudes = append(longitudes, length(lista_txt[[i]]))
}

grupo <- list(
    "No lesionado - LP" = c("GSM8458830", "GSM8458832", "GSM8458834", "GSM8458836"),
    "Lesionado - LP"    = c("GSM8458831", "GSM8458833", "GSM8458835", "GSM8458837")
)
metadata <- stack(grupo)
names(metadata) <- c("Sample", "Group")
metadata <- metadata[, c("Group", "Sample")]

posiciones <- match(gsub("(GSE274837_|.txt)", "", names(lista_txt)), 
                    metadata$Sample)
info_grupo <- setNames(rep(metadata$Group[posiciones], longitudes), 
                       names(merged_df))

objeto1 <- column_to_rownames(merged_df, var = "gene")

matriz_counts <- as.matrix(objeto1)
seu_objeto <- CreateSeuratObject(counts = matriz_counts, project = "MiProyecto")
raw_countsB <- GetAssayData(object = seu_objeto, layer = "counts")

# Establecer filtro con genes en minimo 10 celulas y celulas con minimo 80 genes
sampleB <- CreateSeuratObject(counts = raw_countsB, project = "sampleA", 
                              min.cells = 3, min.features = 20)
sampleB$batch <- "SampleB"

ref <- celldex::HumanPrimaryCellAtlasData()

## Convert Seurat object to SingleCellExperiment object and perform cell type annotation
sce <- as.SingleCellExperiment(sampleB)
# Add logcounts to SingleCellExperiment
sce <- scuttle::logNormCounts(sce)

## Perform cell type annotation with SingleR using the loaded reference
annotations <- SingleR(test = sce, ref = ref, labels = ref$label.main)

## Add the cell type annotations to the Seurat object
sampleB$cell_type <- annotations$labels
print(table(sampleB$cell_type))

sampleB$batch <- as.character(info_grupo)
sampleB$batch <- factor(sampleB$batch)

A <- data.frame(
    id = 1:dim(sampleB)[2],
    cell_type = sampleB@meta.data$cell_type,
    nCount_RNA = sampleB@meta.data$nCount_RNA,
    nFeature_RNA = sampleB@meta.data$nFeature_RNA,
    cell_type = sampleB@meta.data$cell_type,
    batch = sampleB@meta.data$batch)%>% 
    group_by(cell_type) %>% 
    summarise(conteo = n()) %>% 
    mutate(conteo = ifelse(conteo <= 50, NA, conteo)) %>% 
    drop_na()

which(sampleB@meta.data$cell_type %in% A$cell_type)
sampleB@meta.data[which(sampleB@meta.data$cell_type %in% A$cell_type), ]
objeto1[, which(sampleB@meta.data$cell_type %in% A$cell_type)]

sampleB$tipo_muestra <- paste(sampleB@meta.data$cell_type, 
                              gsub(" - LP", "", sampleB@meta.data$batch),
                              sep = "-")

write.csv(sampleB@meta.data[which(sampleB@meta.data$cell_type %in% A$cell_type), ], 
          file = "Datos_Anotados.csv", row.names = FALSE)

write.csv(objeto1[, which(sampleB@meta.data$cell_type %in% A$cell_type)], 
          file = "Datos_Modelo.csv")

write.csv(t(objeto1[, which(sampleB@meta.data$cell_type %in% A$cell_type)]),
          file = "Datos_Transpuestos.csv")
