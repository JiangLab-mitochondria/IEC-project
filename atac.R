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

# 1.data process
archr_clustering <- function(sample) {
  sample <- as.character(sample)
  print(sample)
  print("==================start=======================")
  # inputfile <- paste0("/mnt/transposon1/zhangyanxiaoLab/chaiguoshi/projects/jiangmin-lab-project/work/zhang-qian-project/resutls/mtscATAC/cellranger-atac-results/",
  #                     sample, "/outs/possorted_bam_rmdup_DNA_fragments_sorted.bed.gz")
  inputfile <- paste0("/mnt/transposon1/zhangyanxiaoLab/chaiguoshi/HPC_projects/zhang-qian-projects/results/mtscATAC/cellranger-atac-results/",
                      sample, "/outs/possorted_bam_rmdup_DNA_fragments_sorted.bed.gz")
  sample.name <- sample
  ArrowFiles <- createArrowFiles(
    inputFiles = inputfile,
    sampleNames = sample.name,
    force = T
  )
  
  
  ArrowFiles <- c(paste0(sample.name, ".arrow"))
  doubScores <- addDoubletScores(
    input = ArrowFiles,
    k = 10, #Refers to how many cells near a "pseudo-doublet" to count.
    knnMethod = "UMAP", #Refers to the embedding to use for nearest neighbor search with doublet projection.
    LSIMethod = 1
  )
  
  
  projMB1 <- ArchRProject(
    ArrowFiles = ArrowFiles,
    copyArrows = F)
  projMB2 <- filterDoublets(projMB1)
  
  projMB2.IterativeLSI <- addIterativeLSI(
    ArchRProj = projMB2,
    useMatrix = "TileMatrix",
    name = "IterativeLSI",
    iterations = 8,
    clusterParams = list(
      resolution = c(0.2),
      sampleCells = 10000,
      n.start = 10
    ),
    varFeatures = 25000,
    dimsToUse = 1:30,
    force = T
  )
  projMB2.IterativeLSI.addclusters <- addClusters(
    input = projMB2.IterativeLSI,
    reducedDims = "IterativeLSI",
    method = "Seurat",
    name = "Clusters",
    resolution = 0.8,
    force = T
  )
  projMB2.IterativeLSI.addclusters.addumap <- addUMAP(
    ArchRProj = projMB2.IterativeLSI.addclusters,
    reducedDims = "IterativeLSI",
    name = "UMAP",
    nNeighbors = 30,
    minDist = 0.5,
    metric = "cosine",
    force = T
  )
  projMB2.IterativeLSI.addclusters.addumap.addtsne <- addTSNE(
    ArchRProj = projMB2.IterativeLSI.addclusters.addumap,
    reducedDims = "IterativeLSI",
    name = "TSNE",
    perplexity = 30,
    force = T
  )
  saveArchRProject(ArchRProj = projMB2.IterativeLSI.addclusters.addumap.addtsne,
                   outputDirectory = paste0("clustering/", sample.name),
                   load = TRUE,
                   dropCells = TRUE
  )
  
  rm(ArrowFiles, doubScores, projMB1, projMB2, projMB2.IterativeLSI,
     projMB2.IterativeLSI.addclusters, projMB2.IterativeLSI.addclusters.addumap,
     projMB2.IterativeLSI.addclusters.addumap.addtsne
  )
  
  print(sample)
  print("==================end=======================")
  print("")
  return()
}
result <- archr_clustering(sample = "Intestinal-epithelium-100W-ND5-G12918A-2135-70")
saveArchRProject(ArchRProj = dna ,  load = TRUE, dropCells = TRUE,outputDirectory = paste0("data/", "Intestinal-epithelium-100W-ND5-G12918A-2135-70"))

archr_clustering_addHarmony <- function(sample_vector, sample_name) {
  sample_name <- as.character(sample_name)
  ArrowFiles <- sample_vector
  projMB1 <- ArchRProject(
    ArrowFiles = ArrowFiles, 
    copyArrows = F)
  
  projMB1.IterativeLSI <- addIterativeLSI(
    ArchRProj = projMB1,
    useMatrix = "TileMatrix", 
    name = "IterativeLSI", 
    iterations = 8, 
    clusterParams = list( 
      resolution = c(0.2), 
      sampleCells = 10000, 
      n.start = 10
    ), 
    varFeatures = 25000, 
    dimsToUse = 1:30,
    force = T
  )
  
  projMB1.IterativeLSI.addHarmony <- addHarmony(
    ArchRProj = projMB1.IterativeLSI,
    reducedDims = "IterativeLSI",
    name = "Harmony",
    groupBy = "Sample"
  )
  
  projMB1.IterativeLSI.addHarmony.addclusters <- addClusters(
    input = projMB1.IterativeLSI.addHarmony,
    reducedDims = "Harmony",
    method = "Seurat",
    name = "Clusters",
    resolution = 0.8,
    force = T
  )
  projMB1.IterativeLSI.addHarmony.addclusters.addumap <- addUMAP(
    ArchRProj = projMB1.IterativeLSI.addHarmony.addclusters, 
    reducedDims = "Harmony", 
    name = "UMAPHarmony", 
    nNeighbors = 30, 
    minDist = 0.5, 
    metric = "cosine",
    force = T
  )
  projMB1.IterativeLSI.addHarmony.addclusters.addumap.addtsne <- addTSNE(
    ArchRProj = projMB1.IterativeLSI.addHarmony.addclusters.addumap, 
    reducedDims = "Harmony", 
    name = "TSNEHarmony", 
    perplexity = 30,
    force = T
  )
  saveArchRProject(ArchRProj = projMB1.IterativeLSI.addHarmony.addclusters.addumap.addtsne, 
                   outputDirectory = paste0("res/", sample_name), 
                   load = TRUE,
                   dropCells = TRUE
  )
  return()
}
result <- archr_clustering_addHarmony(sample_vector = c("Intestinal-epithelium-11-samples_merged_harmony/ArrowFiles/Intestinal-epithelium-3W-ND5-G12918A-6-70.arrow",
                                                        "Intestinal-epithelium-11-samples_merged_harmony/ArrowFiles/Intestinal-epithelium-3W-TrnA-G5081A-5292-79.arrow",
                                                        "Intestinal-epithelium-11-samples_merged_harmony/ArrowFiles/Intestinal-epithelium-25W-TrnA-G5081A-4508-78.arrow",
                                                        "Intestinal-epithelium-11-samples_merged_harmony/ArrowFiles/Intestinal-epithelium-25w-TrnA-G5081A-5470-80.arrow",
                                                        "Intestinal-epithelium-11-samples_merged_harmony/ArrowFiles/Intestinal-epithelium-25w-WT-5476-0.arrow",
                                                        "Intestinal-epithelium-11-samples_merged_harmony/ArrowFiles/Intestinal-epithelium-54W-TrnA-G5081A-4993-73.arrow",
                                                        "Intestinal-epithelium-11-samples_merged_harmony/ArrowFiles/Intestinal-epithelium-55W-TrnA-G5081A-4841-82.arrow",
                                                        "Intestinal-epithelium-11-samples_merged_harmony/ArrowFiles/Intestinal-epithelium-100W-ND5-G12918A-2135-70.arrow",
                                                        "Intestinal-epithelium-100W-TrnA-G5081A-26-84-masked/ArrowFiles/Intestinal-epithelium-100W-TrnA-G5081A-26-84-masked.arrow"
                                                        
),sample_name = "9samples_merged")

dna_1 <-  addClusters( input = result ,reducedDims = "Harmony",method = "Seurat",name = "Clusters",resolution = 1.0,force = T)
plotEmbedding(ArchRProj = dna_1, colorBy = "cellColData", name = "Clusters", embedding = "UMAPHarmony")
dna <- dna_1
saveArchRProject(ArchRProj = dna ,  load = TRUE, dropCells = TRUE,outputDirectory = paste0("data/", "Arch_cluster_1.0"))
                 
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

## 绘图

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


#3. TF analysis
dna_peak <- addGroupCoverages(dna_anno,groupBy = "predictedGroup_Co")

dna_peak <- addReproduciblePeakSet(ArchRProj = dna_peak,groupBy = "predictedGroup_Co",
                                   pathToMacs2 = "/home/jiangminLab/dengxiaoling/miniconda3/envs/seq/bin/macs2")
dna_peak_matrix  <- addPeakMatrix(dna_peak)

saveArchRProject(ArchRProj = dna_peak_matrix , load = FALSE, dropCells = TRUE,
                 outputDirectory = paste0("data/", "Arch_cluster1_single_peak"))
### 3.1 mutation load merge
dna <- dna_peak_matrix
dna.cell.info <- as.data.frame(getCellColData(dna))
dna.cell.info$barcode <- rownames(dna.cell.info)
data <-read.table("mutation_load_cell_type.txt", header = TRUE, sep = "\t", comment.char = "", stringsAsFactors = FALSE)
df <- data[,c(18,21:23)]
df$mutation_load <- round(df$A.depth / df$coverage * 100,2)
df$norm_mutation_load <- round(df$mutation_load /df$tail_mutation_load,2)
df <- merge(df,dna.cell.info,by="barcode")
df <- df[df$coverage >=10,] 
df$sample2 <- paste(str_split(df$barcode,"-",simplify = T)[,3],str_split(df$barcode,"-",simplify = T)[,4],sep ="_")
df$cell <- df$predictedGroup_Co
df_all <- df 

### 3.2 TF analysis

type <- unique(dna_anno$predictedGroup_Co)
sample <- unique(df_all$sample2)

for (sample_name in sample) {
  res_high <- df_all %>% group_by(sample2, predictedGroup_Co)  %>% arrange(mutation_load) %>% mutate( n = n(),rank = row_number()) %>%
    filter(rank >= 0.75 * n , sample2 %in%  sample_name) %>% mutate(condition="high")
  res_low <- df_all %>% group_by(sample2, predictedGroup_Co) %>% arrange(mutation_load) %>% mutate( n = n(),rank = row_number()) %>%
    filter(rank <= (0.25 * n+1), sample2 %in% sample_name) %>% mutate(condition="low")
  res <- rbind(res_high,res_low)
  idxSample <- BiocGenerics::which(rownames(dna.cell.info) %in% res$barcode)
  cellsSample <- dna$cellNames[idxSample]
  dna_subset  <- dna[cellsSample, ]
  res <- res[match(rownames(getCellColData(dna_subset)), res$barcode), ] 
  identical(rownames(getCellColData(dna_subset)),res$barcode)
  
  dna_subset$mutation_load <- res$mutation_load
  dna_subset$condition <- res$condition
  dna_subset$norm_mutation_load <- res$norm_mutation_load
  for (i in type) {
    idxSample <- BiocGenerics::which( dna_subset$predictedGroup_Co %in% i)
    cellsSample <- dna_subset$cellNames[idxSample]
    dna_subset_en  <- dna_subset[cellsSample, ]
    markersPeaks <- getMarkerFeatures(ArchRProj = dna_subset_en, useMatrix = "PeakMatrix",  groupBy = "condition", bias = c("TSSEnrichment", "log10(nFrags)"),
                                      testMethod = "wilcoxon",useGroups = "high",bgdGroups = "low")
    dna_motif  <- addMotifAnnotations(ArchRProj = dna_subset_en, motifSet = "cisbp", name = "Motif")
    motif_up<- peakAnnoEnrichment(seMarker = markersPeaks,ArchRProj = dna_motif,peakAnnotation = "Motif",cutOff = "FDR <= 0.1 & Log2FC >= 0.25")
    motif_down<- peakAnnoEnrichment(seMarker = markersPeaks,ArchRProj = dna_motif,peakAnnotation = "Motif",cutOff = "FDR <= 0.1 & Log2FC <=( -0.25)")
    
    data_up <- data.frame(Padj=10^(-motif_up@assays@data$mlog10Padj[,1]),P =10^(-motif_up@assays@data$mlog10p[,1]),Enrichment=motif_up@assays@data$Enrichment[,1],row.names = rownames(motif_up))
    data_down <- data.frame(Padj=10^(-motif_down@assays@data$mlog10Padj[,1]),P =10^(-motif_down@assays@data$mlog10p[,1]),Enrichment=motif_down@assays@data$Enrichment[,1],row.names = rownames(motif_down))
    markerList <- getMarkers(markersPeaks, cutOff = "FDR <= 0.1 & abs(Log2FC) >= 0.25")
    data_peak <- as.data.frame(markerList$high)
    data <- readRDS("data/Arch_cluster1_single_peak/Annotations/Motif-Matches-In-Peaks.rds")
    gene <- c(rownames(data_up)[data_up$Padj < 0.05],rownames(data_down)[data_down$Padj < 0.05])
    data_res <- lapply(gene, function(x){
      tmp <- data[data@assays@data$matches[,x]==1]
      data_tf <- unname(tmp@rowRanges)
      data_tf <- as.data.frame(data_tf)
      data_tf <- data_tf[data_tf$start %in% data_peak$start, ]
      data_tf$TF <- x
      return(data_tf)
    })
    data_tf  <- do.call(rbind, data_res)
    
    if(nrow(data_peak!=0)){
      wb <- createWorkbook()
      addWorksheet(wb, "tf_up")
      writeData(wb, "tf_up", data_up, rowNames = TRUE )
      addWorksheet(wb, "tf_down")
      writeData(wb, "tf_down", data_down, rowNames = TRUE)
      addWorksheet(wb, "tf")
      writeData(wb, "tf", data_tf, rowNames = FALSE)
      addWorksheet(wb, "peak")
      writeData(wb, "peak", data_peak, rowNames = FALSE)
      saveWorkbook(wb, paste(sample_name,i,"tf.xlsx",sep = "_"), overwrite = TRUE)
      
      pdf(paste(sample_name,i,"tf.pdf",sep = "_"))
      p1 <- plotMarkers(seMarker = markersPeaks,  cutOff = "FDR <= 0.1 & abs(Log2FC) >= 0.25", plotAs = "Volcano",name="high")
      print(p1)
      df <- data.frame(TF = rownames(motif_up), mlog10Padj = assay(motif_up)[,1])
      df <- df[order(df$mlog10Padj, decreasing = TRUE),]
      df$rank <- seq_len(nrow(df))
      p2 <- ggplot(df, aes(rank, mlog10Padj, color = mlog10Padj)) + 
        geom_point(size = 1) +
        ggrepel::geom_label_repel(data = df[rev(seq_len(30)), ], aes(x = rank, y = mlog10Padj, label = TF), size = 1.5,nudge_x = 2,color = "black",max.overlaps = 30) + 
        ylab("-log10(P-adj) Motif Enrichment") + xlab("Rank Sorted TFs Enriched") + ggtitle("motif up")+
        theme_ArchR() +  scale_color_gradientn(colors = paletteContinuous(set = "captain"))+
        geom_hline(yintercept = -log10(0.1), color = "black", linetype = "dashed",linewidth=0.1)
      print(p2)
      df <- data.frame(TF = rownames(motif_down), mlog10Padj = assay(motif_down)[,1])
      df <- df[order(df$mlog10Padj, decreasing = TRUE),]
      df$rank <- seq_len(nrow(df))
      p3 <- ggplot(df, aes(rank, mlog10Padj, color = mlog10Padj)) + 
        geom_point(size = 1) +
        ggrepel::geom_label_repel(data = df[rev(seq_len(30)), ], aes(x = rank, y = mlog10Padj, label = TF), size = 1.5,nudge_x = 2,color = "black",max.overlaps = 30) + 
        ylab("-log10(P-adj) Motif Enrichment") + xlab("Rank Sorted TFs Enriched") + ggtitle("motif down")+
        theme_ArchR() +  scale_color_gradientn(colors = paletteContinuous(set = "captain"))+
        geom_hline(yintercept = -log10(0.1), color = "black", linetype = "dashed",linewidth=0.1)
      print(p3)
      dev.off()
    }
  }
}

#4. plot
### 4.1 quality control
dna.cell.info <- as.data.frame(getCellColData(dna))
dna.cell.info$sample2 <- paste(str_split(dna.cell.info$Sample,"-",simplify = T)[,3],str_split(dna.cell.info$Sample,"-",simplify = T)[,4],sep ="_")
dna.cell.info$sample2 <- factor(dna.cell.info$sample2 ,levels = c("3W_ND5","3W_TrnA","25W_TrnA_2","25w_TrnA" , "25w_WT","54W_TrnA","55W_TrnA" ,"100W_ND5","100W_TrnA"))
df <- dna.cell.info[dna.cell.info$sample2 !="54W_TrnA" ,]
library(ggplot2)
library(ggbeeswarm)
library(ggpubr)
color<- brewer.pal( name = "Set3",n=9)

df_subset <- df[df$sample2 %in% c("3W_TrnA","25W_TrnA_2","55W_TrnA","100W_TrnA") ,]

pdf("quality_A.pdf")
ggplot(data=df_subset,aes(y=TSSEnrichment,x=sample2)) +
  geom_violin(linewidth =0.5,trim = F,alpha=0.6,aes(fill=sample2))+
  geom_boxplot(width=0.2,aes(fill=sample2),alpha=0.9) +
  xlab(label ="") +ylab(label = "Score") + theme_classic() +ggtitle("TSS Enrichment Score")+
  theme (legend.position ="none",
         axis.line = element_line(color ="black", size = 0.5),plot.title = element_text(hjust = 0.5, face = "bold",size=20),
         axis.ticks.length = unit(0.3,"cm"),axis.ticks = element_line(size = 0.5), 
         axis.title = element_text(size = 15,,face ="bold"),axis.text = element_text(size = 10,face ="bold") ) +
  scale_fill_manual(values =color )

ggplot(data=df_subset,aes(y=nFrags,x=sample2)) +
  geom_violin(linewidth =0.5,trim = F,alpha=0.6,aes(fill=sample2))+
  geom_boxplot(width=0.1,aes(fill=sample2),alpha=0.9) +
  xlab(label ="") +ylab(label = "Fragments per cell") + theme_classic() +ggtitle("High-quality Fragments Per Cell")+
  theme (legend.position ="none",
         axis.line = element_line(color ="black", size = 0.5),plot.title = element_text(hjust = 0.5, face = "bold",size=20),
         axis.ticks.length= unit(0.3,"cm"),axis.ticks = element_line(size = 0.5), 
         axis.title = element_text(size = 15,,face ="bold"),axis.text = element_text(size = 10,face ="bold") ) +
  scale_fill_manual(values =color )
dev.off()

df_subset <- df[df$sample2 %in% c("3W_ND5","100W_ND5") ,]
pdf("quality_nd5.pdf")
ggplot(data=df_subset,aes(y=TSSEnrichment,x=sample2)) +
  geom_violin(linewidth =0.5,trim = F,alpha=0.6,aes(fill=sample2))+
  geom_boxplot(width=0.2,aes(fill=sample2),alpha=0.9) +
  xlab(label ="") +ylab(label = "Score") + theme_classic() +ggtitle("TSS Enrichment Score")+
  theme (legend.position ="none",
         axis.line = element_line(color ="black", size = 0.5),plot.title = element_text(hjust = 0.5, face = "bold",size=20),
         axis.ticks.length = unit(0.3,"cm"),axis.ticks = element_line(size = 0.5), 
         axis.title = element_text(size = 15,,face ="bold"),axis.text = element_text(size = 10,face ="bold") ) +
  scale_fill_manual(values =color[5:6] )

ggplot(data=df_subset,aes(y=nFrags,x=sample2)) +
  geom_violin(linewidth =0.5,trim = F,alpha=0.6,aes(fill=sample2))+
  geom_boxplot(width=0.1,aes(fill=sample2),alpha=0.9) +
  xlab(label ="") +ylab(label = "Fragments per cell") + theme_classic() +ggtitle("High-quality Fragments Per Cell")+
  theme (legend.position ="none",
         axis.line = element_line(color ="black", size = 0.5),plot.title = element_text(hjust = 0.5, face = "bold",size=20),
         axis.ticks.length= unit(0.3,"cm"),axis.ticks = element_line(size = 0.5), 
         axis.title = element_text(size = 15,,face ="bold"),axis.text = element_text(size = 10,face ="bold") ) +
  scale_fill_manual(values =color[5:6] )
dev.off()

## 合并coverage
library(dplyr)
library(stringr)
data <- read.table("coverage.txt",header = T)
data_subset <- data[data$sample %in% c("3W_ND5","3W_TrnA","25W_TrnA_2","25W_TrnA" , "25W_WT","55W_TrnA" ,"100W_ND5","100W_TrnA") ,]
data_subset$sample <- factor(data_subset$sample ,levels = c("3W_ND5","3W_TrnA","25W_TrnA_2","25w_TrnA" , "25w_WT","55W_TrnA" ,"100W_ND5","100W_TrnA"))

df <- data_subset[data_subset$sample %in% c("3W_TrnA","25W_TrnA_2","55W_TrnA","100W_TrnA") ,]
pdf("coverage_A.pdf",width = 5,height = 2.5)
ggplot(df, aes(x = position, y = mean_coverage,group = sample,colour =sample)) +
  geom_line(alpha=0.8,linewidth=0.4) + theme_bw()+
  labs(x = "", y = "mean coverage per cell") +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),plot.title = element_text(hjust = 0.5, size = 15, face = "bold"),
        axis.title = element_text(size = 10,,face ="bold"),axis.text = element_text(size = 6,face ="bold") )+
  scale_color_manual(values =color )
dev.off()

df <- data_subset[data_subset$sample %in% c("3W_ND5","100W_ND5") ,]
pdf("coverage_nd5.pdf",width = 5,height = 2.5)
ggplot(df, aes(x = position, y = mean_coverage,group = sample,colour =sample)) +
  geom_line(alpha=0.8,linewidth=0.4) + theme_bw()+
  labs(x = "", y = "mean coverage per cell") +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),plot.title = element_text(hjust = 0.5, size = 15, face = "bold"),
        axis.title = element_text(size = 10,,face ="bold"),axis.text = element_text(size = 6,face ="bold") )+
  scale_color_manual(values =color[5:6] )
dev.off()


### 4.2 peak
dna$group <- factor(dna$predictedGroup_Co,levels = c("Enterocyte","Enterocyte progenitor","Stem cell","TA","Goblet","Secretory progenitor","EEC","Paneth","Tuft","T cell"))

cluster_col <- c("#A6CEE3","#1F78B4","#B2DF8A","#33A02C","#FB9A99","#E31A1C","#FDBF6F","#FF7F00","#CAB2D6", "#6A3D9A")

dna$group[dna$group=="Stem cell"]<-"AStem cell"
dna$group[dna$group=="TA"]<-"BTA"
dna$group[dna$group=="Enterocyte progenitor"]<-"CEnterocyte progenitor"
dna$group[dna$group=="Enterocyte"]<-"DEnterocyte"
dna$group[dna$group=="Secretory progenitor"]<-"ESecretory progenitor"
dna$group[dna$group=="Tuft"]<-"FTuft"
dna$group[dna$group=="EEC"]<-"GEEC"
dna$group[dna$group=="Goblet"]<-"HGoblet"
dna$group[dna$group=="Paneth"]<-"IPaneth"
dna$group[dna$group=="T cell"]<-"JT cell"

cluster_col <- c("#FDBF6F","#CAB2D6","#B2DF8A","#1F78B4","#E31A1C", "#6A3D9A","#A6CEE3","#33A02C","#FB9A99","#FF7F00")


pdf("peak.pdf")
gene <- "Lgr5"
p <- plotBrowserTrack(ArchRProj = dna, groupBy = "group",geneSymbol = gene,title = gene,pal=cluster_col,
                      upstream = 150000, downstream = 10000,plotSummary = c("bulkTrack", "geneTrack")
)
grid::grid.newpage()
grid::grid.draw(p[[gene]])


gene <- "Olfm4"
p <- plotBrowserTrack(ArchRProj = dna, groupBy = "group",geneSymbol = gene,title = gene,pal=cluster_col, 
                      highlight = GRanges(seqnames = Rle("chr14", 1),ranges = IRanges(start = 79999000, end = 80000900)),
                      upstream = 10000, downstream = 30000,plotSummary = c("bulkTrack", "geneTrack")
)
grid::grid.newpage()
grid::grid.draw(p[[gene]])  


gene <- "Mki67"
p <- plotBrowserTrack(ArchRProj = dna, groupBy = "group",geneSymbol = gene,title = gene,pal=cluster_col,
                      upstream = 30000, downstream = 3000,plotSummary = c("bulkTrack", "geneTrack")
)
grid::grid.newpage()
grid::grid.draw(p[[gene]])

gene <- "Pcna"
p <- plotBrowserTrack(ArchRProj = dna, groupBy = "group",geneSymbol = gene,title = gene,pal=cluster_col,
                      upstream = 5000, downstream = 1000,plotSummary = c("bulkTrack", "geneTrack")
)
grid::grid.draw(p[[gene]])


gene <- "Rbp7"
p <- plotBrowserTrack(ArchRProj = dna, groupBy = "group",geneSymbol = gene,title = gene,pal=cluster_col,
                      upstream = 7000, downstream = 1000,plotSummary = c("bulkTrack", "geneTrack")
)
grid::grid.newpage()
grid::grid.draw(p[[gene]])


gene <- "Cenpa"
p <- plotBrowserTrack(ArchRProj = dna, groupBy = "group",geneSymbol = gene,title = gene,pal=cluster_col,
                      upstream = 500, downstream = 8000,plotSummary = c("bulkTrack", "geneTrack")
)
grid::grid.newpage()
grid::grid.draw(p[[gene]])


gene <- "Alpi"
p <- plotBrowserTrack(ArchRProj = dna, groupBy = "group",geneSymbol = gene,title = gene,pal=cluster_col,
                      upstream = 5000, downstream = 500,plotSummary = c("bulkTrack", "geneTrack")
)
grid::grid.newpage()
grid::grid.draw(p[[gene]])

gene <- "Apoa1"
p <- plotBrowserTrack(ArchRProj = dna, groupBy = "group",geneSymbol = gene,title = gene,pal=cluster_col,
                      upstream = 1000, downstream = 3000,plotSummary = c("bulkTrack", "geneTrack")
)
grid::grid.newpage()
grid::grid.draw(p[[gene]])


gene <- "Atoh1"
p <- plotBrowserTrack(ArchRProj = dna, groupBy = "group",geneSymbol = gene,title = gene,pal=cluster_col,
                      upstream = 1000, downstream = 3000,plotSummary = c("bulkTrack", "geneTrack")
)
grid::grid.newpage()
grid::grid.draw(p[[gene]])

gene <- "Dll1"
p <- plotBrowserTrack(ArchRProj = dna, groupBy = "group",geneSymbol = gene,title = gene,pal=cluster_col,
                      upstream = 10000, downstream = 500,plotSummary = c("bulkTrack", "geneTrack")
)
grid::grid.newpage()
grid::grid.draw(p[[gene]])

gene <- "Dclk1"
p <- plotBrowserTrack(ArchRProj = dna, groupBy = "group",geneSymbol = gene,title = gene,pal=cluster_col,
                      upstream = 10000, downstream = 300000,plotSummary = c("bulkTrack", "geneTrack")
)
grid::grid.newpage()
grid::grid.draw(p[[gene]])

gene <- "Chga"
p <- plotBrowserTrack(ArchRProj = dna, groupBy = "group",geneSymbol = gene,title = gene,pal=cluster_col,
                      upstream = 1000, downstream = 20000,plotSummary = c("bulkTrack", "geneTrack")
)
grid::grid.newpage()
grid::grid.draw(p[[gene]])


gene <- "Muc2"
p <- plotBrowserTrack(ArchRProj = dna, groupBy = "group",geneSymbol = gene,title = gene,pal=cluster_col,
                      upstream = 5000, downstream = 80000,plotSummary = c("bulkTrack", "geneTrack")
)
grid::grid.newpage()
grid::grid.draw(p[[gene]])

gene <- "Tff3"
p <- plotBrowserTrack(ArchRProj = dna, groupBy = "group",geneSymbol = gene,title = gene,pal=cluster_col,
                      upstream = 5000, downstream = 1000,plotSummary = c("bulkTrack", "geneTrack")
)
grid::grid.newpage()
grid::grid.draw(p[[gene]])

gene <- "Lyz1"
p <- plotBrowserTrack(ArchRProj = dna, groupBy = "group",geneSymbol = gene,title = gene,pal=cluster_col,
                      upstream = 8000, downstream = 1000,plotSummary = c("bulkTrack", "geneTrack")
)
grid::grid.newpage()
grid::grid.draw(p[[gene]])

gene <- "Defa17"
p <- plotBrowserTrack(ArchRProj = dna, groupBy = "group",geneSymbol = gene,title = gene,pal=cluster_col,
                      upstream = 1000, downstream = 3000,plotSummary = c("bulkTrack", "geneTrack")
)
grid::grid.newpage()
grid::grid.draw(p[[gene]])

gene <- "Ptprc"
p <- plotBrowserTrack(ArchRProj = dna, groupBy = "group",geneSymbol = gene,title = gene,pal=cluster_col,
                      upstream = 200000, downstream = 5000,plotSummary = c("bulkTrack", "geneTrack")
)
grid::grid.newpage()
grid::grid.draw(p[[gene]])

gene <- "Cd3e"
p <- plotBrowserTrack(ArchRProj = dna, groupBy = "group",geneSymbol = gene,title = gene,pal=cluster_col,
                      upstream = 20000, downstream = 3000,plotSummary = c("bulkTrack", "geneTrack")
)
grid::grid.newpage()
grid::grid.draw(p[[gene]])
dev.off()


pdf("peak_final.pdf")

gene <- "Lgr5"
p <- plotBrowserTrack(ArchRProj = dna, groupBy = "group",geneSymbol = gene,title = gene,pal=cluster_col, highlightFill ="#B3B3B3",
                      region = GRanges(seqnames = Rle("chr10", 1),ranges = IRanges(start = 115470000, end = 115530000)),
                      highlight = GRanges(seqnames = Rle("chr10", 1),ranges = IRanges(start = 115483000, end = 115486000)),
                      plotSummary = c("bulkTrack", "geneTrack")
)
grid::grid.newpage()
grid::grid.draw(p)  

gene <- "Mki67"
p <- plotBrowserTrack(ArchRProj = dna, groupBy = "group",geneSymbol = gene,title = gene,pal=cluster_col, highlightFill ="#B3B3B3",
                      region = GRanges(seqnames = Rle("chr7", 1),ranges = IRanges(start = 135670000, end = 135730000)),
                      highlight = GRanges(seqnames = Rle("chr7", 1),ranges = IRanges(start = 135702000, end = 135705000)),
                      plotSummary = c("bulkTrack", "geneTrack")
)
grid::grid.newpage()
grid::grid.draw(p)  


gene <- "Rbp7"
p <- plotBrowserTrack(ArchRProj = dna, groupBy = "group",geneSymbol = gene,title = gene,pal=cluster_col, highlightFill ="#B3B3B3",
                      region = GRanges(seqnames = Rle("chr4", 1),ranges = IRanges(start = 149430000, end = 149490000)),
                      highlight = GRanges(seqnames = Rle("chr4", 1),ranges = IRanges(start = 149453000, end = 149456000)),
                      plotSummary = c("bulkTrack", "geneTrack")
)
grid::grid.newpage()
grid::grid.draw(p) 

gene <- "Alpi"
p <- plotBrowserTrack(ArchRProj = dna, groupBy = "group",geneSymbol = gene,title = gene,pal=cluster_col, highlightFill ="#B3B3B3",
                      region = GRanges(seqnames = Rle("chr1", 1),ranges = IRanges(start = 87070000, end = 87130000)),
                      highlight = GRanges(seqnames = Rle("chr1", 1),ranges = IRanges(start = 87100000, end = 87103000)),
                      plotSummary = c("bulkTrack", "geneTrack")
)
grid::grid.newpage()
grid::grid.draw(p)

pdf("peak_Dll1.pdf")
gene <- "Dll1"
p <- plotBrowserTrack(ArchRProj = dna, groupBy = "group",geneSymbol = gene,title = gene,pal=cluster_col, highlightFill ="#B3B3B3",
                      region = GRanges(seqnames = Rle("chr17", 1),ranges = IRanges(start = 15360000, end = 15420000)),
                      highlight = GRanges(seqnames = Rle("chr17", 1),ranges = IRanges(start = 15375000, end = 15378000)),
                      plotSummary = c("bulkTrack", "geneTrack")
)
grid::grid.newpage()
grid::grid.draw(p)
dev.off()
gene <- "Dclk1"
p <- plotBrowserTrack(ArchRProj = dna, groupBy = "group",geneSymbol = gene,title = gene,pal=cluster_col, highlightFill ="#B3B3B3",
                      region = GRanges(seqnames = Rle("chr3", 1),ranges = IRanges(start = 55220000, end =55280000)),
                      highlight = GRanges(seqnames = Rle("chr3", 1),ranges = IRanges(start = 55241000, end = 55244000)),
                      plotSummary = c("bulkTrack", "geneTrack")
)
grid::grid.newpage()
grid::grid.draw(p)

gene <- "Chga"
p <- plotBrowserTrack(ArchRProj = dna, groupBy = "group",geneSymbol = gene,title = gene,pal=cluster_col, highlightFill ="#B3B3B3",
                      region = GRanges(seqnames = Rle("chr12", 1),ranges = IRanges(start = 102530000, end =102590000)),
                      highlight = GRanges(seqnames = Rle("chr12", 1),ranges = IRanges(start = 102553500, end = 102556500)),
                      plotSummary = c("bulkTrack", "geneTrack")
)
grid::grid.newpage()
grid::grid.draw(p)

gene <- "Muc2"
p <- plotBrowserTrack(ArchRProj = dna, groupBy = "group",geneSymbol = gene,title = gene,pal=cluster_col, highlightFill ="#B3B3B3",
                      region = GRanges(seqnames = Rle("chr7", 1),ranges = IRanges(start = 141670000, end =141730000)),
                      highlight = GRanges(seqnames = Rle("chr7", 1),ranges = IRanges(start = 141702000, end = 141705000)),
                      plotSummary = c("bulkTrack", "geneTrack")
)
grid::grid.newpage()
grid::grid.draw(p)

gene <- "Lyz1"
p <- plotBrowserTrack(ArchRProj = dna, groupBy = "group",geneSymbol = gene,title = gene,pal=cluster_col, highlightFill ="#B3B3B3",
                      region = GRanges(seqnames = Rle("chr10", 1),ranges = IRanges(start = 117260000, end =117320000)),
                      highlight = GRanges(seqnames = Rle("chr10", 1),ranges = IRanges(start = 117286000, end = 117289000)),
                      plotSummary = c("bulkTrack", "geneTrack")
)
grid::grid.newpage()
grid::grid.draw(p)

gene <- "Cd3e"
p <- plotBrowserTrack(ArchRProj = dna, groupBy = "group",geneSymbol = gene,title = gene,pal=cluster_col, highlightFill ="#B3B3B3",
                      region = GRanges(seqnames = Rle("chr9", 1),ranges = IRanges(start = 44970000, end =45030000)),
                      highlight = GRanges(seqnames = Rle("chr9", 1),ranges = IRanges(start = 45008000, end = 45011000)),
                      plotSummary = c("bulkTrack", "geneTrack")
)
grid::grid.newpage()
grid::grid.draw(p)

dev.off()

### 4.2 propotion
freq_df <- res %>% group_by(sample, group) %>%  summarise(Number = n() )
result <- freq_df %>%  group_by(sample) %>% reframe(type = group, number=Number, percentage = (Number/ sum(Number)) * 100 )
write.csv(result,"freq_all.csv",row.names = F,quote = F)
result <- result[result$sample %in% c("Intestinal-epithelium-25w-WT-5476-0","Intestinal-epithelium-25w-TrnA-G5081A-5470-80","Intestinal-epithelium-25W-TrnA-G5081A-4508-78"),]
result$group <- paste(str_split(result$sample,"-",simplify = T)[,3],str_split(result$sample,"-",simplify = T)[,4],sep="_")


result$group <- factor(result$group,levels = c("25w_WT","25w_TrnA","25W_TrnA_2"))
result$type <- factor(result$type,levels = c("Stem cell","TA","Enterocyte progenitor","Enterocyte","Secretory progenitor","Tuft","EEC","Goblet","Paneth","T cell"))
col_2 <- c("#FDBF6F","#CAB2D6","#B2DF8A","#1F78B4","#E31A1C","#6A3D9A","#A6CEE3","#33A02C","#FB9A99","#FF7F00")
write.csv(result ,"atac_freq_all_25wA.csv",row.names = F)
pdf("propotion_all_25wA.pdf",width = 4)
ggplot(result, aes( x = group,y=percentage,fill = type))+
  geom_bar(position = "stack", stat = "identity", width = 0.65) +
  theme_classic()+   
  scale_fill_manual(values = col_2 )+  
  scale_y_continuous(expand = c(0,0))+# 调整y轴属性，使柱子与X轴坐标接触
  labs(x=NULL,y="Proportion") +
  guides(fill=guide_legend(title=NULL))+
  theme(axis.text.x = element_text(size = 8),axis.title.y = element_text(size = 14))
dev.off() 


