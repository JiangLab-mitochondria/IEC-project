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