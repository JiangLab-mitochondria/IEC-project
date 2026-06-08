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

# 2.annotation
top_genes <- c("Lgr5","Slc12a2","Axin2","Olfm4","Gkn3","Mki67","Cdk4","Mcm5","Mcm6","Pcna","Rbp7","Ube2c","Cdc20","Birc5","Cenpa",
               "Alpi","Apoa1","Atoh1","Ccl25","Ldha","Reg3g","Dll1","Dclk1", "Trpm5","Gfi1b","Il25","Chga","Chgb","Neurod1","Cck","Nts", 
               "Muc2","Tff3","Agr2","Lyz1","Defa17","Defa22","Defa24","Ang4", "Ptprc","Cd3e")
pdf("markerplot_1.pdf")
for (i in 1:41) {
  df <- FetchData(object = combined_sce_harmony, vars = c(top_genes[i], "umap_1","umap_2","RNA_snn_res.1.0"),slots = c("data","umap","meta.data") )
  df$value <- df[,1]
  cluster_centers <- df %>% group_by(RNA_snn_res.1.0) %>% summarize(UMAP1 = mean(umap_1), UMAP2 = mean(umap_2), .groups = 'drop')
  masktable<-generateMask( dims=df[,2:3],cluster=df$RNA_snn_res.1, minDensity=1.0,smoothSigma=0.02)
  
  p1 <-ggplot(df, aes(x = umap_1, y = umap_2,color = value)) +
    geom_point(size = 0.1,alpha=1) +
    labs(x = "UMAP 1", y = "UMAP 2", color = "log2(TPM+1)") +
    theme_bw() +
    scale_color_gradientn(colors = c("grey80", "#A1CDE1", "#fff143", "#EC9274", "red"),
                          values = c(0, 0.225, 0.5, 0.75, 1)) +
    geom_path(data=masktable,aes(group=group),linewidth=0.4,linetype=2,col="black") +
    geom_text(data = cluster_centers, aes(x = UMAP1, y = UMAP2, label =RNA_snn_res.1.0),
              size = 5, vjust = 0, hjust = 0.5,color="black") +
    ggtitle(top_genes[i])+
    theme(plot.title = element_text(face = "bold",hjust=0.5) ,axis.title = element_text(size=10, face="bold"), panel.grid = element_blank(), 
          axis.text= element_blank(),axis.ticks=element_blank()) +coord_fixed(ratio = 0.7)
  
  print(p1)
  
}
dev.off()


combined_sce_harmony@meta.data$type <- ifelse(combined_sce_harmony@meta.data$RNA_snn_res.1.4 %in% c(0,27,28,15),"Stem cell",
                                              ifelse(combined_sce_harmony@meta.data$RNA_snn_res.1.4 %in% c(5,7,8,12,13,14,20),"Enterocyte",
                                                     ifelse(combined_sce_harmony@meta.data$RNA_snn_res.1.4 %in% c(9,17),"Enterocyte progenitor",
                                                            #ifelse(combined_sce_harmony@meta.data$RNA_snn_res.1.4 %in% c(10),"Enterocyte progenitor 2",
                                                            # ifelse(combined_sce_harmony@meta.data$RNA_snn_res.1.4 %in% c(20),"Enterocyte progenitor 3",
                                                            ifelse(combined_sce_harmony@meta.data$RNA_snn_res.1.4 %in% c(1,4,16,10),"TA",
                                                                   ifelse(combined_sce_harmony@meta.data$RNA_snn_res.1.4 %in% 11,"Paneth",
                                                                          ifelse(combined_sce_harmony@meta.data$RNA_snn_res.1.4 %in% c(2,3,22),"Goblet",
                                                                                 ifelse(combined_sce_harmony@meta.data$RNA_snn_res.1.4 %in% c(21,6,18),"Secretory progenitor",
                                                                                        ifelse(combined_sce_harmony@meta.data$RNA_snn_res.1.4 %in% c(24,25,19),"EEC",
                                                                                               ifelse(combined_sce_harmony@meta.data$RNA_snn_res.1.4 %in% 23,"Tuft","T cell"))))))))
                                              
)
saveRDS(combined_sce_harmony,file="cluster.rds")

df_A54w <- combined_sce_harmony@meta.data[combined_sce_harmony@meta.data$orig.ident == "A_54w",c(6,16)]
df_A54w <- scRNA@meta.data[scRNA@meta.data$orig.ident == "A_54w",c(6,16)]
df_A54w$barcode <- str_split(rownames(df_A54w),"_",simplify = T)[,2]
write.table(df_A54w,"barcode_54w.txt",row.names = F,quote = F,sep="\t")

combined_sce_harmony_filt <- subset(combined_sce_harmony,type !="T cell") # 去除T细胞污染
saveRDS(combined_sce_harmony_filt,file="combined_sce_harmony_filt.rds")

# 3.DEG and pathway
scRNA <- combined_sce_harmony_filt
Idents(scRNA) <- "group"
scRNA_26w <- subset(scRNA,idents = c("26w_A","26w_WT"))
scRNA_26w<- JoinLayers(scRNA_26w, assay = "RNA")
Idents(scRNA_26w ) <- "type"
cell_types <- unique(Idents(scRNA_26w))
results <- list()
for (cell_type in cell_types) {
  cell_type_data <- subset(scRNA_26w, idents = cell_type)
  Idents(cell_type_data) <- "group"
  diff_genes <- FindMarkers(cell_type_data, ident.1 = "26w_A", ident.2 = "26w_WT",only.pos = F,min.pct = 0.25)
  results[[cell_type]] <- diff_genes
}

ID<-read.csv("/storage/jiangminLab/dengxiaoling/zhangqian/data6/ID.csv",stringsAsFactors = F,row.names = 1)

enrichment_up <- list()
enrichment_down <- list()
enrichment_up_s <- list()
enrichment_down_s <- list()


for (cell_type in cell_types) {
  allDEG <- results[[cell_type]]
  gene_up <- allDEG[allDEG$avg_log2FC>=0.25 & allDEG$p_val < 0.05,]
  gene_down <- allDEG[allDEG$avg_log2FC <= (-0.25) & allDEG$p_val  < 0.05,]
  ego_up <- enrichGO(gene =ID[which(ID$external_gene_name %in% rownames(gene_up)),3] ,  OrgDb = org.Mm.eg.db,keyType = "ENTREZID", ont = "BP",readable = T) 
  ego_up_s <- simplify(ego_up,cutoff = 0.7)
  enrichment_up[[cell_type]] <- ego_up
  enrichment_up_s[[cell_type]] <- ego_up_s
  ego_down <- enrichGO(gene = ID[which(ID$external_gene_name %in% rownames(gene_down)),3], OrgDb = org.Mm.eg.db,keyType = "ENTREZID", ont = "BP",readable = T) 
  ego_down_s <- simplify(ego_down,cutoff = 0.7)
  enrichment_down[[cell_type]] <- ego_down
  enrichment_down_s[[cell_type]] <- ego_down_s
}

for (cell_type in cell_types) {
  write.csv(enrichment_up[[cell_type]]@result, file = paste0(cell_type, "_enrich_up.csv"), row.names = FALSE,quote = F)
  write.csv(enrichment_down[[cell_type]]@result, file = paste0(cell_type, "_enrich_down.csv"), row.names = FALSE,quote = F)
  p1<-dotplot(enrichment_up[[cell_type]],showCategory=20,,font.size=8) + ggtitle("GO_up") + theme(plot.title = element_text(hjust = 0.5))
  p2<-dotplot(enrichment_down[[cell_type]],showCategory=20,,font.size=8) + ggtitle("GO_down") + theme(plot.title = element_text(hjust = 0.5))
  pdf(paste0(cell_type, "_enrich.pdf"))
  if(nrow(p1$data)!=0){print(p1)}
  if(nrow(p2$data)!=0){print(p2)}
  dev.off()
}

### top 10
tmp <- lapply(enrichment_up_s, function(x) {
  df <- as.data.frame(x) 
  top_10_significant <-  df %>% arrange(qvalue) %>% head(5)
  top_10_significant$type <- names(x)  
  return(top_10_significant)
})
merged_df <- do.call(rbind, tmp)
top_up <- merged_df$ID

tmp <- lapply(enrichment_down_s, function(x) {
  df <- as.data.frame(x) 
  top_10_significant <-  df %>% arrange(qvalue) %>% head(5)
  top_10_significant$type <- names(x)  
  return(top_10_significant)
})
merged_df <- do.call(rbind, tmp)
top_down <- merged_df$ID


all_up <- data.frame()
for (i in 1:length(cell_types)) {
  df <- enrichment_up[[i]]@result
  df$type <- names(enrichment_up)[i]
  all_up <- rbind(all_up,df[df$ID %in% top_up,])
}

all_down <- data.frame()
for (i in 1:length(cell_types)) {
  df <- enrichment_down[[i]]@result
  df$type <- names(enrichment_down)[i]
  all_down <- rbind(all_down,df[df$ID %in% top_down,])
}

res_up <- data.frame()
up_term <- unique(all_up$Description )
for (i in 1:length(up_term)) {
  df <- all_up[ all_up$Description %in% up_term[i] ,c(2,10,13)]
  df <- df %>% group_by(Description) %>% pivot_wider(names_from = type,values_from = qvalue )
  res_up <- rbind(res_up,df)
}
res_up  <- na.omit(res_up)

res_down <- data.frame()
down_term <- unique(all_down$Description )
for (i in 1:length(down_term)) {
  df <- all_down[ all_down$Description %in% down_term[i] ,c(2,10,13)]
  df <- df %>% group_by(Description) %>% pivot_wider(names_from = type,values_from = qvalue )
  res_down <- rbind(res_down,df)
}
res_down  <- na.omit(res_down)

### plot
res_up <-as.data.frame(res_up)
rownames(res_up) <- res_up[,1]
res_up <- res_up[,-1]
res_up <- -log10(res_up)
head(res_up)
res_down <-as.data.frame(res_down)
rownames(res_down) <- res_down[,1]
res_down <- res_down[,-1]
res_down <- -log10(res_down)
library(circlize)

exp <- apply(res_up, 1, scale)
rownames(exp) <- colnames(res_up)
exp <- t(exp)

p1 <- Heatmap(exp, col = colorRamp2(c(-2, 0, 2), c("#366799", "white", "#C43F33")),heatmap_legend_param = list(title = "-log10(qvalue)"),,row_names_max_width = unit(10, 'cm'),
              show_column_dend = F,show_row_dend = F,column_title ="Upregulation pathway",column_title_gp = gpar(fontsize = 20,fontface = "bold"),
              column_names_rot = 45,column_names_side = "top",column_names_gp = gpar(fontsize = 12),column_names_centered = F,row_names_gp = gpar(fontsize = 10),
              column_order = c("Stem cell","TA","Enterocyte progenitor","Enterocyte","Secretory progenitor","Tuft","EEC","Goblet","Paneth"))
p1

exp <- apply(res_down , 1, scale)
rownames(exp) <- colnames(res_down)
exp <- t(exp)

p2 <- Heatmap(exp, col = colorRamp2(c(-2, 0, 2), c("#366799", "white", "#C43F33")),heatmap_legend_param = list(title = "-log10(qvalue)"),row_names_max_width = unit(10, 'cm'),
              show_column_dend = F,show_row_dend = F,column_title ="Downregulation pathway",column_title_gp = gpar(fontsize = 20,fontface = "bold"),
              column_names_rot = 45,column_names_side = "top",column_names_gp = gpar(fontsize = 12),column_names_centered = F,row_names_gp = gpar(fontsize = 10),
              column_order = c("Stem cell","TA","Enterocyte progenitor","Enterocyte","Secretory progenitor","Tuft","EEC","Goblet","Paneth"))
p2

pdf("enrich_heatmap_5_simplify.pdf",width = 14,height = 12)
p1
p2
dev.off()

# 4.marker plot
top_genes <- c("Lgr5","Slc12a2","Axin2","Olfm4","Gkn3","Mki67","Cdk4","Mcm5","Mcm6","Pcna","Rbp7","Ube2c","Cdc20","Birc5","Cenpa",
               "Alpi","Apoa1","Atoh1","Ccl25","Ldha","Reg3g","Dll1","Dclk1", "Trpm5","Gfi1b","Il25","Chga","Chgb","Neurod1","Cck","Nts", 
               "Muc2","Tff3","Agr2","Lyz1","Defa17","Defa22","Defa24","Ang4", "Ptprc","Cd3e")


Idents(scRNA) <- "group"
scRNA_26w <- subset(scRNA,idents = c("26w_A","26w_WT"))
Idents(scRNA_26w ) <- "type"
Idents(scRNA_26w ) <- factor(Idents(scRNA_26w ),levels=rev(c("Stem cell","TA" ,"Enterocyte progenitor","Enterocyte","Secretory progenitor","Tuft","EEC","Goblet","Paneth","T cell")))



pdf("marker.pdf",width = 20)
DotPlot(scRNA_26w, features = top_genes, dot.scale = 10,scale = TRUE) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, size = 12), axis.text.y = element_text(face = "italic", size = 15),
        legend.text = element_text(size = 9), plot.title = element_text(hjust = 0.5, face = "bold"),
        axis.title = element_blank()) +
  scale_color_gradientn(colors = c("#A1CDE1", "white", "red"), values = c(0,  0.5,  1)) +
  
  labs(title = "26w scRNA-seq Data Marker Expression")
dev.off()
