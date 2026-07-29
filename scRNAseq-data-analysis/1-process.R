library(tidyr)
library(Seurat)
library(harmony)
library(ggplot2)
library(cowplot)
library(mascarade)
library(dplyr)
library(stringr)
set.seed(1234)
library(RColorBrewer)
brewer_palette <- brewer.pal(n = 12, name = "Paired")
brewer_palette2 <- brewer.pal(n = 5, name = "PuOr")
col <- c(brewer_palette,brewer_palette2,brewer_palette)
library(cowplot)
library(patchwork)
library(ggpubr)
library(org.Mm.eg.db)
library(clusterProfiler)
library(enrichplot)
library(gridExtra)
library(ComplexHeatmap)


# 1.data process
Seurat_26w_A<- readRDS("Intestinal-epithelium-26W-TrnA-G5081A-5199-77.rds")
Seurat_26w_WT <- readRDS("Intestinal-epithelium-26W-WT-5194-0.rds")
Seurat_54w_A <-  readRDS("Intestinal-epithelium-54W-TrnA-G5081A-4993-73.rds")
Seurat_26w_A_2<- readRDS("Intestinal-epithelium-26W-TrnA-2.rds")
Seurat_26w_WT_2<- readRDS("Intestinal-epithelium-26W-WT-2.rds")

Seurat_26w_A@meta.data$group <- "26w_A"
Seurat_26w_WT@meta.data$group <- "26w_WT"
Seurat_26w_A_2@meta.data$group <- "26w_A"
Seurat_26w_WT_2@meta.data$group <- "26w_WT"
Seurat_54w_A@meta.data$group <- "54w_A"
head(Seurat_26w_A_2@meta.data)
Seurat_26w_A$orig.ident <- "A_26w"
Seurat_26w_WT$orig.ident <- "WT_26w"
Seurat_54w_A$orig.ident <- "A_54w"
Seurat_26w_A_2$orig.ident <- "A_26w_2"
Seurat_26w_WT_2$orig.ident <- "WT_26w_2"
head(Seurat_26w_A_2$orig.ident)
Idents(Seurat_26w_A) <- Seurat_26w_A$orig.ident
Idents(Seurat_26w_WT) <- Seurat_26w_WT$orig.ident
Idents(Seurat_26w_A_2) <- Seurat_26w_A_2$orig.ident
Idents(Seurat_26w_WT_2) <- Seurat_26w_WT_2$orig.ident
Idents(Seurat_54w_A) <- Seurat_54w_A$orig.ident
head(Idents(Seurat_54w_A))
Seurat_26w_A <- RenameCells(Seurat_26w_A, add.cell.id  = "A26w")
Seurat_26w_WT <- RenameCells(Seurat_26w_WT, add.cell.id  = "WT26w")
Seurat_54w_A <- RenameCells(Seurat_54w_A, add.cell.id  = "A54w")
Seurat_26w_A_2 <- RenameCells(Seurat_26w_A_2, add.cell.id  = "A26w_2")
Seurat_26w_WT_2 <- RenameCells(Seurat_26w_WT_2, add.cell.id  = "WT26w_2")

combined_sce <- merge(Seurat_26w_A,y=c(Seurat_26w_WT,Seurat_54w_A,Seurat_26w_A_2,Seurat_26w_WT_2))
p3 <- VlnPlot(combined_sce, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo"), ncol = 4, pt.size = 0)
p3
combined_sce_filter <- subset(combined_sce, subset = nFeature_RNA > 500 & nFeature_RNA < 8000 & percent.mt < 10)

combined_sce[["percent.mt"]] <- PercentageFeatureSet(combined_sce, pattern = "^mt-")
n_features <- combined_sce@meta.data$nFeature_RNA
mito_ratio <- combined_sce@meta.data$percent.mt
pdf("quality.pdf",height = 8,width = 8)
p1 <- VlnPlot(combined_sce, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,cols =brewer.pal(n = 5, name = "Set2")) 
p2 <-VlnPlot(combined_sce_filter, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,cols =brewer.pal(n = 5, name = "Set2"))
plot_grid(p1,p2, nrow = 2)

p3 <- ggplot(data = data.frame(n_features), aes(x = n_features)) +
  geom_histogram(binwidth = 50,fill="grey60") +
  labs(title = "Distribution of nFeature_RNA", x = "nFeature_RNA",  y = "count") +
  geom_vline(aes(xintercept = 500), linetype = "dashed", color = "red", linewidth = 0.5)+
  theme_classic()+ 
  scale_y_continuous (expand = c (0,0))+scale_x_continuous (expand = c (0,0))
p3p <- p3   + theme(legend.position = 'none') +xlab("") + ylab("") +scale_x_continuous (expand = c (0,0),limits = c(0, 2500))+ 
  theme(plot.title = element_blank()) + theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5))
p3 + inset_element(p3p, 0.3, 0.4, 0.8, 0.7, on_top = TRUE) 

p4 <- ggplot(data = data.frame(mito_ratio), aes(x = mito_ratio)) +
  geom_histogram(binwidth =0.5,fill="grey60") +
  labs(title = "Distribution of mito ratio",x = "mito ratio", y = "count") +
  theme_classic() +scale_y_continuous (expand = c (0,0))+scale_x_continuous (expand = c (0,0))+
  geom_vline(aes(xintercept = 10), linetype = "dashed", color = "red", linewidth = 0.5)
p4p <- p4 + theme(legend.position = 'none') +xlab("") + ylab("") +scale_x_continuous (expand = c (0,0),limits = c(0, 20))+ 
  theme(plot.title = element_blank()) + theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5))

p4 + inset_element(p4p, 0.3, 0.4, 0.8, 0.7, on_top = TRUE) 

dev.off()

combined_sce_filter<- NormalizeData(combined_sce_filter)
combined_sce_filter <- FindVariableFeatures(combined_sce_filter, selection.method = "vst", nfeatures = 2000)
combined_sce_filter <- ScaleData(combined_sce_filter, vars.to.regress = "nCount_RNA")
combined_sce_filter <- RunPCA(combined_sce_filter, features = VariableFeatures(object = combined_sce_filter))
combined_sce_harmony <- RunHarmony(object = combined_sce_filter, group.by.vars = "orig.ident")
combined_sce_harmony <- FindNeighbors(combined_sce_harmony, dims = 1:30, reduction = "harmony")
combined_sce_harmony <- FindClusters(combined_sce_harmony, resolution = c(0.1,0.2,0.3,0.4, 0.6, 0.8, 1.0, 1.0), reduction = "harmony")
combined_sce_harmony <- RunUMAP(combined_sce_harmony, reduction = "harmony",dims = 1:30)
combined_sce_harmony <- RunTSNE(combined_sce_harmony, reduction = "harmony",dims = 1:30)

p1<-DimPlot(combined_sce_harmony, reduction = "umap", group.by = "RNA_snn_res.1.0",label = T,label.size = 6,repel =T )
p2<-DimPlot(combined_sce_harmony, reduction = "tsne", group.by = "RNA_snn_res.1.0",label = T,label.size = 6,repel =T )
plot_grid(p1,p2)

DimPlot(combined_sce_harmony, reduction = "umap", group.by = "orig.ident")

combined_sce_filter <- RunUMAP(combined_sce_filter, reduction = "pca", dims = 1:40)
pdf("batch_correction.pdf",width = 12)
p1<- DimPlot(combined_sce_filter, reduction = "umap", group.by = "orig.ident",cols =brewer.pal(n = 5, name = "Set2"))+
  ggtitle("Pre-batch correction")+ 
  theme(axis.title = element_text(size=12, face="bold"), panel.grid = element_blank(), 
        legend.position = "none", legend.title = element_blank()  )
p2 <- DimPlot(combined_sce_harmony, reduction = "umap", group.by = "orig.ident",cols =brewer.pal(n = 5, name = "Set2"))+
  ggtitle("After-batch correction")+ 
  theme(axis.title = element_text(size=12, face="bold"), panel.grid = element_blank(), 
        axis.text.y = element_blank(),axis.title.y =element_blank(),axis.line.y = element_blank(),axis.ticks.y=element_blank(),
        legend.position = "right", legend.title = element_blank() )
ggarrange(p1,p2,common.legend = T,legend="right")
dev.off()

saveRDS(combined_sce_harmony,file="combined_sce_harmony.rds")