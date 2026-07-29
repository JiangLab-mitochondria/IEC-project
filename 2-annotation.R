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

combined_sce_harmony <- readRDS("combined_sce_harmony.rds")
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