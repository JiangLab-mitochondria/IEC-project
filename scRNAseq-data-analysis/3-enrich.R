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


# 3.DEG and pathway
scRNA <- readRDS("combined_sce_harmony_filt.rds")
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

ID<-read.csv("ID.csv",stringsAsFactors = F,row.names = 1)

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