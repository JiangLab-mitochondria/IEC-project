library(ArchR)
addArchRThreads(threads = 6)
library(dplyr)
library(Seurat)
library(stringr)
library(mascarade)
library(RColorBrewer)
library(BSgenome.Mmusculus.UCSC.mm10)
library(openxlsx)
brewer_palette <- brewer.pal(n = 12, name = "Paired")
brewer_palette2 <- brewer.pal(n = 5, name = "PuOr")
col <- c(brewer_palette,brewer_palette2,brewer_palette)
set.seed(1234)
library(openxlsx)
library(patchwork)

#2. annotation
scRNA <- readRDS("res/scrna.rds")
dna<- loadArchRProject(path="data/Arch_cluster_1.0")

#### 2.1 1st
dna_anno <- addGeneIntegrationMatrix(ArchRProj = dna, seRNA = scRNA,groupRNA = "type",addToArrow = FALSE,threads = 1,
                                     useMatrix = "GeneScoreMatrix",matrixName = "GeneIntegrationMatrix",reducedDims = "Harmony")
dna.cell.info <- as.data.frame(getCellColData(dna ))
dna.cell.info$barcode <- rownames(dna.cell.info)
dna.umap <- getEmbedding(ArchRProj = dna , embedding = "UMAPHarmony", returnDF = T)
colnames(dna.umap) <- c("UMAP_1", "UMAP_2")
dna.umap$barcode <- rownames(dna.umap)
df <- merge(dna.cell.info ,dna.umap,by="barcode")                      
head(df)
cluster_centers <- df %>% group_by(Clusters) %>% summarize(UMAP1 = mean(UMAP_1), UMAP2 = mean(UMAP_2), .groups = 'drop')
df$Clusters <- factor(df$Clusters,levels=paste0("C",1:21))
#masktable<-generateMask( dims=df[,4:5],cluster=df$Clusters, minDensity=1.0,smoothSigma=0.02)

p1 <- ggplot(df, aes(x = UMAP_1, y = UMAP_2, color = Clusters)) +
  geom_point(size = 0.5,alpha=1) +
  labs(x = "UMAP 1", y = "UMAP 2", color = "Cluster") +
  theme_minimal() +
  geom_text(data = cluster_centers, aes(x = UMAP1, y = UMAP2, label = Clusters),
            size = 5, vjust = 1, hjust = 0.5,color="black",fontface="bold")+
  guides(color = guide_legend(override.aes = list(size = 3)))+
  scale_color_manual(values = col) + 
  coord_fixed(ratio=1)

p1
cluster_centers <- df %>% group_by(predictedGroup) %>% summarize(UMAP1 = mean(UMAP_1)-0.5, UMAP2 = median(UMAP_2), .groups = 'drop')
axis <- ggh4x::guide_axis_truncated(trunc_lower = unit(0, "npc"),trunc_upper = unit(3, "cm"))
#masktable<-generateMask( dims=df[,4:5],cluster=df$ kendall_pred, minDensity=1.0,smoothSigma=0.02)
p2 <- ggplot(df, aes(x = UMAP_1, y = UMAP_2, color =  predictedGroup)) +
  geom_point(size = 1,alpha=0.8) +
  labs(x = "UMAP 1", y = "UMAP 2", color = "Annote_nature") +
  theme_minimal() +
  geom_text(data = cluster_centers, aes(x = UMAP1, y = UMAP2, label = predictedGroup),
            size = 4, vjust = 0 ,hjust = 0.5,color="black",fontface="bold")+
  guides(color = guide_legend(override.aes = list(size = 3)))+
  scale_color_manual(values = col)  +coord_fixed(ratio=1) +
  theme(aspect.ratio = 1, panel.background = element_blank(), panel.grid = element_blank(),
        axis.line = element_line(arrow = arrow(type = "closed",length = unit(0.2,"cm"))),axis.title = element_text(hjust = 0.05, face = "italic")) +
  scale_x_continuous(breaks = NULL)+
  scale_y_continuous(breaks = NULL) +
  guides(color = FALSE, x = axis, y = axis)
p2



pdf("archr_single_new.pdf",width = 20)
grid.arrange(p1, p2, ncol = 2)
dev.off()
#### 2.2 2nd
dna <- dna_anno
cM <- as.matrix(confusionMatrix(dna$Clusters, dna$predictedGroup)) 
preClust <- colnames(cM)[apply(cM, 1 , which.max)]

rnaGoblet <- as.character(colnames(scRNA)[grep("Goblet", scRNA$type )])
atacGoblet <- as.character(dna$cellNames[dna$Clusters %in% c("C13","C14")])
rnaEnterocyte <- as.character(colnames(scRNA)[grep("Enterocyte", scRNA$type )])
atacEnterocyte <- as.character(dna$cellNames[dna$Clusters %in% c("C1","C2","C3","C4","C5","C6","C7","C8")])
rnaPaneth<- as.character(colnames(scRNA)[grep("Paneth", scRNA$type )])
atacPaneth <- as.character(dna$cellNames[dna$Clusters %in% "C10"])
rnaEEC<- as.character(colnames(scRNA)[grep("EEC", scRNA$type )])
atacEEC <- as.character(dna$cellNames[dna$Clusters %in% "C11"])

rna_other <- as.character(colnames(scRNA)[ !colnames(scRNA) %in% c(rnaPaneth,rnaEEC) ])
atac_other <- as.character(dna$cellNames[ !dna$cellNames %in% c(atacGoblet,atacPaneth,atacEEC,atacEnterocyte)])

groupList <- SimpleList(
  Goblet = SimpleList(ATAC =atacGoblet ,RNA =rnaGoblet ),
  Enterocyte = SimpleList(ATAC = atacEnterocyte,RNA = rnaEnterocyte),
  Paneth = SimpleList(ATAC = atacPaneth,RNA = rnaPaneth),
  EEC = SimpleList(ATAC = atacEEC,RNA =rnaEEC ),
  other = SimpleList(ATAC = atac_other,RNA =rna_other )
)

dna_anno <- addGeneIntegrationMatrix(ArchRProj = dna, seRNA = scRNA,groupRNA = "type",addToArrow = FALSE,threads = 1,groupList = groupList,
                                     useMatrix = "GeneScoreMatrix",matrixName = "GeneIntegrationMatrix",reducedDims = "Harmony",
                                     nameCell = "predictedCell_Co",nameGroup = "predictedGroup_Co",nameScore = "predictedScore_Co")

saveArchRProject(ArchRProj = dna_anno ,  load = TRUE,dropCells = TRUE,
                 outputDirectory =  paste0("data/", "Arch_cluster1_single"))

dna.cell.info <- as.data.frame(getCellColData(dna_anno ))
dna.cell.info$barcode <- rownames(dna.cell.info)
dna.umap <- getEmbedding(ArchRProj = dna_anno , embedding = "UMAPHarmony", returnDF = T)
colnames(dna.umap) <- c("UMAP_1", "UMAP_2")
dna.umap$barcode <- rownames(dna.umap)
df <- merge(dna.cell.info ,dna.umap,by="barcode") 

cluster_centers <- df %>% group_by(predictedGroup_Co) %>% summarize(UMAP1 = mean(UMAP_1)-0.5, UMAP2 = median(UMAP_2), .groups = 'drop')
axis <- ggh4x::guide_axis_truncated(trunc_lower = unit(0, "npc"),trunc_upper = unit(3, "cm"))
#masktable<-generateMask( dims=df[,4:5],cluster=df$ kendall_pred, minDensity=1.0,smoothSigma=0.02)
p2 <- ggplot(df, aes(x = UMAP_1, y = UMAP_2, color =  predictedGroup_Co)) +
  geom_point(size = 1,alpha=0.8) +
  labs(x = "UMAP 1", y = "UMAP 2", color = "Annote_nature") +
  theme_minimal() +
  geom_text(data = cluster_centers, aes(x = UMAP1, y = UMAP2, label = predictedGroup_Co),
            size = 4, vjust = 0 ,hjust = 0.5,color="black",fontface="bold")+
  guides(color = guide_legend(override.aes = list(size = 3)))+
  scale_color_manual(values = col)  +coord_fixed(ratio=1) +
  theme(aspect.ratio = 1, panel.background = element_blank(), panel.grid = element_blank(),
        axis.line = element_line(arrow = arrow(type = "closed",length = unit(0.2,"cm"))),axis.title = element_text(hjust = 0.05, face = "italic")) +
  scale_x_continuous(breaks = NULL)+
  scale_y_continuous(breaks = NULL) +
  guides(color = FALSE, x = axis, y = axis)
p2

pdf("archr_single_new.pdf",width = 10)
p2
dev.off()

cluster_centers <- df %>% group_by(Clusters) %>% summarize(UMAP1 = mean(UMAP_1)-0.5, UMAP2 = median(UMAP_2), .groups = 'drop')
axis <- ggh4x::guide_axis_truncated(trunc_lower = unit(0, "npc"),trunc_upper = unit(3, "cm"))
p3 <- ggplot(df, aes(x = UMAP_1, y = UMAP_2, color = Clusters)) +
  geom_point(size = 1,alpha=0.8) +
  labs(x = "UMAP 1", y = "UMAP 2", color = "Annote_nature") +
  theme_minimal() +
  geom_text(data = cluster_centers, aes(x = UMAP1, y = UMAP2, label = Clusters),
            size = 4, vjust = 0 ,hjust = 0.5,color="black",fontface="bold")+
  guides(color = guide_legend(override.aes = list(size = 3)))+
  scale_color_manual(values = col)  +coord_fixed(ratio=1) +
  theme(aspect.ratio = 1, panel.background = element_blank(), panel.grid = element_blank(),
        axis.line = element_line(arrow = arrow(type = "closed",length = unit(0.2,"cm"))),axis.title = element_text(hjust = 0.05, face = "italic")) +
  scale_x_continuous(breaks = NULL)+
  scale_y_continuous(breaks = NULL) +
  guides(color = FALSE, x = axis, y = axis)
p3

pdf("archr_single_cluster_new.pdf",width = 10)
p3
dev.off()

#### 2.3 54w multi
dna.cell.info <- as.data.frame(getCellColData(dna_anno))
dna.cell.info$barcode <- rownames(dna.cell.info)
dna.cell.info <- dna.cell.info %>% separate(barcode, c("sample2", "barcode2"), remove = F, sep = "#")
dna_trna <- dna.cell.info[grep("Intestinal-epithelium-54W-TrnA",dna.cell.info$barcode),c(ncol(dna.cell.info),21)]
head(dna_trna)
#dna_trna$barcode <- str_split(dna_trna$barcode,"#",simplify = T)[,2]

data <- read.table("barcode_54w.txt", sep="\t",header = TRUE, stringsAsFactors = FALSE)
cell_type <- data[,2:3]
head(cell_type)
anno_data <- merge(cell_type,dna_trna,by.x="barcode",by.y="barcode2")

dna.umap <- getEmbedding(ArchRProj = dna, embedding = "UMAPHarmony", returnDF = T)
colnames(dna.umap) <- c("UMAP_1", "UMAP_2")
dna.umap$barcode <- rownames(dna.umap)
dna.umap_trna <- dna.umap[grep("Intestinal-epithelium-54W-TrnA",dna.umap$barcode),]
dna.umap_trna$barcode <- str_split(dna.umap_trna$barcode,"#",simplify = T)[,2]
df <- merge(anno_data,dna.umap_trna,by="barcode")
#df <- merge(dna.cell.info,dna.umap,by="barcode")
write.csv(df,"multi_pro.csv")

## plot

cluster_centers <- df %>% group_by(predictedGroup_Co) %>% summarize(UMAP1 = mean(UMAP_1), UMAP2 = mean(UMAP_2), .groups = 'drop') #与tidyr包冲突
p1 <- ggplot(df, aes(x = UMAP_1, y = UMAP_2, color =  predictedGroup_Co)) +
  geom_point(size = 1,alpha=0.8) +
  labs(x = "UMAP 1", y = "UMAP 2", color = "Cell type") +
  theme_minimal() +
  geom_text(data = cluster_centers, aes(x = UMAP1, y = UMAP2, label = predictedGroup_Co),
            size = 5, vjust = 1.6 ,hjust = 0.5,color="black",fontface="bold")+
  guides(color = guide_legend(override.aes = list(size = 3)))+
  scale_color_manual(values = col)  +
  theme(aspect.ratio = 0.9, panel.background = element_blank(), panel.grid = element_blank(),axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank()) +
  geom_segment(aes(x = min(df$UMAP_1) , y = min(df$UMAP_2) ,  xend = min(df$UMAP_1) +3, yend = min(df$UMAP_2) ),colour = "black", size=0.5,arrow = arrow(length = unit(0.3,"cm")))+ 
  geom_segment(aes(x = min(df$UMAP_1)  , y = min(df$UMAP_2) , xend = min(df$UMAP_1) , yend = min(df$UMAP_2) + 3),colour = "black", size=0.5,arrow = arrow(length = unit(0.3,"cm"))) +
  annotate("text", x = min(df$UMAP_1) +1.5, y = min(df$UMAP_2) -0.7, label = "UMAP1", color="black",size = 4) + 
  annotate("text", x = min(df$UMAP_1) -0.7, y = min(df$UMAP_2) + 1.5, label = "UMAP2", color="black",size =4 ,angle=90) 


p1


cluster_centers <- df %>% group_by(type) %>% summarize(UMAP1 = mean(UMAP_1), UMAP2 = mean(UMAP_2), .groups = 'drop') #与tidyr包冲突
p2 <- ggplot(df, aes(x = UMAP_1, y = UMAP_2, color =  type)) +
  geom_point(size = 1,alpha=0.8) +
  labs(x = "UMAP 1", y = "UMAP 2", color = "Cell type") +
  theme_minimal() +
  geom_text(data = cluster_centers, aes(x = UMAP1, y = UMAP2, label = type),
            size = 5, vjust = 1.6 ,hjust = 0.5,color="black",fontface="bold")+
  guides(color = guide_legend(override.aes = list(size = 3)))+
  scale_color_manual(values = col)  +
  theme(aspect.ratio = 0.9, panel.background = element_blank(), panel.grid = element_blank(),axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank()) +
  geom_segment(aes(x = min(df$UMAP_1) , y = min(df$UMAP_2) ,  xend = min(df$UMAP_1) +3, yend = min(df$UMAP_2) ),colour = "black", size=0.5,arrow = arrow(length = unit(0.3,"cm")))+ 
  geom_segment(aes(x = min(df$UMAP_1)  , y = min(df$UMAP_2) , xend = min(df$UMAP_1) , yend = min(df$UMAP_2) + 3),colour = "black", size=0.5,arrow = arrow(length = unit(0.3,"cm"))) +
  annotate("text", x = min(df$UMAP_1) +1.5, y = min(df$UMAP_2) -0.7, label = "UMAP1", color="black",size = 4) + 
  annotate("text", x = min(df$UMAP_1) -0.7, y = min(df$UMAP_2) + 1.5, label = "UMAP2", color="black",size =4 ,angle=90) 

p2
pdf(file = "cluster_multi.pdf", width = 16, height = 10)
grid.arrange(p1, p2, ncol = 2)
dev.off()