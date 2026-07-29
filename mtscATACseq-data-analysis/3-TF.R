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


dna_anno<- loadArchRProject(path="data/Arch_cluster1_single")
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