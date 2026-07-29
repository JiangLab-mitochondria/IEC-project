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

dna <- loadArchRProject(path="data/Arch_cluster1_single_peak")
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
  scale_y_continuous(expand = c(0,0))+
  labs(x=NULL,y="Proportion") +
  guides(fill=guide_legend(title=NULL))+
  theme(axis.text.x = element_text(size = 8),axis.title.y = element_text(size = 14))
dev.off() 