##TRANSCRIPTOMIC ANALYSIS REVEALS SPECIFIC METABOLIC PATHWAYS OF ENTEROHEMORRHAGIC E. COLI O157:H7 IN BOVINE DIGESTIVE CONTENTS
##RNAseq data analysis
#December 2, 2016

#see https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.pdf

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 1. Installing packages and loading libraries
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#install.packages("RSvgDevice")
library(RSvgDevice)
#source("http://bioconductor.org/biocLite.R")
#biocLite("DESeq2")
library(DESeq2)
library(ggplot2)
#install.packages("pheatmap")
library("pheatmap")
library('RColorBrewer')
#install.packages('VennDiagram')
library(VennDiagram)
#install.packages("readODS")
library(readODS)
#see https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.pdf


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 2. Set working directory and input files
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
<<<<<<< HEAD
working_dir="/home/pauline/Bureau/TransEHEC" 						#set working directory
=======
working_dir="~/transEHEC" 						#set working directory
>>>>>>> 43bfabe0a6aa3158b62794f55b5affbde36d5c45
setwd(working_dir)
counts<-read.table("TransEHEC_featureCounts_R.txt",row.names=1,head=T) ; head(counts) 	#count table coming from featureCounts
info<-read.table("sampleSheet.txt",head=T,row.names=1,sep=";") ; head(info)		#sample sheet with sample informations
rnagenes<-read.table("RNAgenes.txt",head=T) ; head(rnagenes)				#list of EDL933 RNA genes
quickC<-read.table("quick_correspondence_281216.txt",h=T,sep=";") ; dim(quickC)		#gene names
annot <- read.csv("Annotations_trimmed_281216.txt", h=T, sep=";") ; head(annot)		#gene annotations

#Select only EDL933 columns in count table and info table
counts_EDL<-counts[1:5652,1:21]
info_EDL<-info[1:21,]

#Remove rRNA, tRNA and ncRNA from count table and quickC table
messagers<-setdiff(as.character(rownames(counts_EDL)), as.character(rnagenes$Locus_tag)) 
counts_EDL <- counts_EDL[which(rownames(counts_EDL)%in%messagers),] ; dim(counts_EDL)
messagers2<-setdiff(as.character(quickC[,1]), as.character(rnagenes$Locus_tag)) 
quickC <- quickC[which(quickC[,1]%in%messagers2),] ; dim(quickC)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 3. Construct DESeq2 object
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#ST variable is a combination of the two factors of interest, DC and growth phases (see sampleSheet).
dds_EDL <- DESeqDataSetFromMatrix(countData = counts_EDL, colData = info_EDL, design = ~ ST) ; dds_EDL
dds_EDL <- dds_EDL[ rowSums(counts(dds_EDL)) > 1, ] ; dds_EDL

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 4. Plot PCA
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
rld_EDL <- rlog(dds_EDL, blind=FALSE)
pdf("PCA.pdf")
data <- plotPCA(rld_EDL, intgroup=c("Source","Time"), returnData=TRUE)#, ntop=50000)
percentVar <- round(100*attr(data, "percentVar"))
ggplot(data, aes(PC1, PC2, color=Source, shape=Time)) +
  geom_point(size=1) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  coord_fixed() + ggtitle("Principal Component Analysis plot") +
  geom_hline(aes(yintercept=0)) +
  geom_vline(aes(xintercept=0)) +
  #scale_shape_manual(values=c(21,24)) +
  geom_text(aes(label=colnames(dds_EDL)),hjust=0.5, vjust=1.8, size=1.3)
dev.off()


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 5. Differential analysis
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
dds_EDL <- DESeq(dds_EDL)

#Results table for the 8 comparisons 
res_EDL_1 <- results(dds_EDL,contrast=c("ST","BICEM","M9EM"))
res_EDL_2 <- results(dds_EDL,contrast=c("ST","RAEM","M9EM"))
res_EDL_3 <- results(dds_EDL,contrast=c("ST","BICEE","M9EE"))
res_EDL_4 <- results(dds_EDL,contrast=c("ST","RAEE","M9EE"))
res_EDL_5 <- results(dds_EDL,contrast=c("ST","RCEE","M9EE"))
res_EDL_6 <- results(dds_EDL,contrast=c("ST","BICEM","BICEE"))
res_EDL_7 <- results(dds_EDL,contrast=c("ST","RAEM","RAEE"))
res_EDL_8 <- results(dds_EDL,contrast=c("ST","M9EM","M9EE"))
res_EDL_9 <- results(dds_EDL,contrast=c("ST","BICEE","RAEE"))
res_EDL_10 <- results(dds_EDL,contrast=c("ST","BICEE","RCEE"))
res_EDL_11 <- results(dds_EDL,contrast=c("ST","BICEM","RAEM"))
res_EDL_12 <- results(dds_EDL,contrast=c("ST","RAEE","RCEE"))


#q-value and fold-change cutoff
res_EDL_1<-res_EDL_1[which(res_EDL_1$padj<=0.05),] ; res_EDL_1<-res_EDL_1[which(abs(res_EDL_1$log2FoldChange) >= 2),] 
res_EDL_2<-res_EDL_2[which(res_EDL_2$padj<=0.05),] ; res_EDL_2<-res_EDL_2[which(abs(res_EDL_2$log2FoldChange) >= 2),] 
res_EDL_3<-res_EDL_3[which(res_EDL_3$padj<=0.05),] ; res_EDL_3<-res_EDL_3[which(abs(res_EDL_3$log2FoldChange) >= 2),] 
res_EDL_4<-res_EDL_4[which(res_EDL_4$padj<=0.05),] ; res_EDL_4<-res_EDL_4[which(abs(res_EDL_4$log2FoldChange) >= 2),] 
res_EDL_5<-res_EDL_5[which(res_EDL_5$padj<=0.05),] ; res_EDL_5<-res_EDL_5[which(abs(res_EDL_5$log2FoldChange) >= 2),] 
res_EDL_6<-res_EDL_6[which(res_EDL_6$padj<=0.05),] ; res_EDL_6<-res_EDL_6[which(abs(res_EDL_6$log2FoldChange) >= 2),] 
res_EDL_7<-res_EDL_7[which(res_EDL_7$padj<=0.05),] ; res_EDL_7<-res_EDL_7[which(abs(res_EDL_7$log2FoldChange) >= 2),] 
res_EDL_8<-res_EDL_8[which(res_EDL_8$padj<=0.05),] ; res_EDL_8<-res_EDL_8[which(abs(res_EDL_8$log2FoldChange) >= 2),] 
res_EDL_9<-res_EDL_9[which(res_EDL_9$padj<=0.05),] ; res_EDL_9<-res_EDL_9[which(abs(res_EDL_9$log2FoldChange) >= 2),] 
res_EDL_10<-res_EDL_10[which(res_EDL_10$padj<=0.05),] ; res_EDL_10<-res_EDL_10[which(abs(res_EDL_10$log2FoldChange) >= 2),] 
res_EDL_11<-res_EDL_11[which(res_EDL_11$padj<=0.05),] ; res_EDL_11<-res_EDL_11[which(abs(res_EDL_11$log2FoldChange) >= 2),] 
res_EDL_12<-res_EDL_12[which(res_EDL_12$padj<=0.05),] ; res_EDL_12<-res_EDL_12[which(abs(res_EDL_12$log2FoldChange) >= 2),] 


#Export DE genes lists
write.table(res_EDL_1,"res_EDL_BICEM_vs_M9EM", sep="\t", dec=",")
write.table(res_EDL_2,"res_EDL_RAEM_vs_M9EM", sep="\t", dec=",")
write.table(res_EDL_3,"res_EDL_BICEE_vs_M9EE", sep="\t", dec=",")
write.table(res_EDL_4,"res_EDL_RAEE_vs_M9EE", sep="\t", dec=",")
write.table(res_EDL_5,"res_EDL_RCEE_vs_M9EE", sep="\t", dec=",")
write.table(res_EDL_6,"res_EDL_BICEM_vs_BICEE", sep="\t", dec=",")
write.table(res_EDL_7,"res_EDL_RAEM_vs_RAEE", sep="\t", dec=",")
write.table(res_EDL_8,"res_EDL_M9EM_vs_M9EE", sep="\t", dec=",")
write.table(res_EDL_9,"res_EDL_BICEE_vs_RAEE", sep="\t", dec=",")
write.table(res_EDL_10,"res_EDL_BICEE_vs_RCEE", sep="\t", dec=",")
write.table(res_EDL_11,"res_EDL_BICEM_vs_RAEM", sep="\t", dec=",")
write.table(res_EDL_12,"res_EDL_RAEE_vs_RCEE", sep="\t", dec=",")



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 6. Get norm counts
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
counts_norm <- counts(dds_EDL, normalized=TRUE) ; head(counts_norm)
write.table(counts_norm,"counts_norm.txt", sep="\t", dec=",")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 7. MA plots
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#defining a new function
plotMA.DESeqResults <- function(object, alpha, main="", xlab="mean of normalized counts", ylim, MLE=FALSE, ...) {
  if (missing(alpha)) {
    alpha <- if (is.null(metadata(object)$alpha)) {
      0.1
    } else {
      metadata(object)$alpha
    }
  }
  df <- if (MLE) {
    # test if MLE is there
    if (is.null(object$lfcMLE)) {
      stop("lfcMLE column is not present: you should first run results() with addMLE=TRUE")
    }
    data.frame(mean = object$baseMean,
               lfc = object$lfcMLE,
               isDE = ifelse(is.na(object$padj), FALSE, object$padj < alpha))
  } else {
    data.frame(mean = object$baseMean,
               lfc = object$log2FoldChange,
               #isDE = ifelse(is.na(object$padj), FALSE, object$padj < alpha))
               isDE = ifelse(object$padj < alpha & abs(object$log2FoldChange) >= 2, TRUE, FALSE))
    
  }
  print(dim(df[df$isDE==TRUE,]))
  if (missing(ylim)) {
    plotMA(df, main=main, xlab=xlab, ...)
  } else {
    plotMA(df, main=main, xlab=xlab, ylim=ylim, ...)
  }  
}

par(mfrow=c(6,2))
pdf("MAplots.pdf")
plotMA(res_EDL_1, main="BICEM vs M9EM", ylim=c(-10,10),alpha=0.05)
plotMA(res_EDL_2, main="RAEM vs M9EM", ylim=c(-10,10),alpha=0.05)
plotMA(res_EDL_3, main="BICEE vs M9EE", ylim=c(-10,10),alpha=0.05)
plotMA(res_EDL_4, main="RAEE vs M9EE", ylim=c(-10,10),alpha=0.05)
plotMA(res_EDL_5, main="RCEE vs M9EE", ylim=c(-10,10),alpha=0.05)
plotMA(res_EDL_6, main="BICEM vs BICEE", ylim=c(-10,10),alpha=0.05)
plotMA(res_EDL_7, main="RAEM vs RAEE", ylim=c(-10,10),alpha=0.05)
plotMA(res_EDL_8, main="M9EM vs M9EE", ylim=c(-10,10),alpha=0.05)
plotMA(res_EDL_9, main="BICEE vs RAEE", ylim=c(-10,10),alpha=0.05)
plotMA(res_EDL_10, main="BICEE vs RCEE", ylim=c(-10,10),alpha=0.05)
plotMA(res_EDL_11, main="BICEM vs RAEM", ylim=c(-10,10),alpha=0.05)
plotMA(res_EDL_12, main="RAEE vs RCEE", ylim=c(-10,10),alpha=0.05)
dev.off()


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 8. Heatmap
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
nCounts_pf <- counts(dds_EDL, normalized=TRUE)

#Rename dds_EDL locus tags with quickC gene names
rownames(quickC)<-quickC[,1]
temp<-c()
for(elem in rownames(dds_EDL)) {
  if(!(elem %in% quickC[,1])) {
    temp <- c(temp,elem)
  }
  else if(is.na(quickC[elem,2])) {
      temp <- c(temp,elem)
  }
  else {
      temp <- c(temp,as.character(quickC[elem,2]))
    }
}
head(rownames(dds_EDL))
head(temp)
rownames(dds_EDL) <- temp

#Heatmap
nCounts <- counts(dds_EDL, normalized=TRUE)

#Method_1
par(mfrow=c(1,1))
select <- nCounts[order(rowMeans(nCounts),decreasing=TRUE)[1:100],]
heatmap(as.matrix(select), Rowv = NA, col = hmcol, mar = c(8,2))
dev.off()

#Method_2
pdf("heatmap.pdf")
select <- order(rowMeans(nCounts),decreasing=TRUE)[1:75]
nt <- normTransform(dds_EDL)
log2.norm.counts <- assay(nt)[select,]
df <- as.data.frame(colData(dds_EDL)[,c("Source","Time")])
dist1 <- "euclidean"
clust <- "average"
pheatmap(log2.norm.counts, clustering_distance_cols = dist1, clustering_method = clust, cluster_rows=FALSE,cluster_cols=TRUE, annotation_col=df,show_rownames=TRUE,fontsize_row=5,fontsize_col=8,main=paste("Clustered heatmap of 75 most abundant genes\n",dist1," distance with ",clust, " clustering method",sep=""),fontsize=7)
dev.off()

#Method_3 : only DE genes
n = 50
resOrdered1 <- res_EDL_1[order(res_EDL_1$padj),]
topResults <- rbind(resOrdered1[ resOrdered1[,'log2FoldChange'] > 0, ][1:n,], resOrdered1[ resOrdered1[,'log2FoldChange'] < 0, ][n:1,] )
topResults[c(1:5,(2*n-4):(2*n)), c('baseMean','log2FoldChange','padj')]
hmcol <- brewer.pal(11,'RdBu')
heatmap(as.matrix(nCounts_pf[row.names(topResults), ]), Rowv = NA, col = hmcol, mar = c(8,2))


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 9. Venn Diagrams
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

BICEM_M9EM<-rownames(res_EDL_1)
RAEM_M9EM<-rownames(res_EDL_2)
BICEE_M9EE<-rownames(res_EDL_3)
RAEE_M9EE<-rownames(res_EDL_4)
RCEE_M9EE<-rownames(res_EDL_5)
BICEM_BICEE<-rownames(res_EDL_6)
RAEM_RAEE<-rownames(res_EDL_7)
M9EM_M9EE<-rownames(res_EDL_8)
BICEE_RAEE<-rownames(res_EDL_9)
BICEE_RCEE<-rownames(res_EDL_10)
BICEM_RAEM<-rownames(res_EDL_11)
RAEE_RCEE<-rownames(res_EDL_12)

#Venn Diagram 1
pdf("Venn_diagram_1.pdf")
inter_BICEM_M9EM_RAEM_M9EM <- intersect(BICEM_M9EM, RAEM_M9EM) ; length(inter_BICEM_M9EM_RAEM_M9EM)
draw.pairwise.venn(length(BICEM_M9EM),length(RAEM_M9EM),length(inter_BICEM_M9EM_RAEM_M9EM),category=c("BICEM_vs_M9EM","RAEM_vs_M9EM"),fill=c("coral","deepskyblue"), cat.col=c("coral","deepskyblue"), cat.cex=0.8, cat.pos=c(3,6), cat.dist=c(0.02,0.02),cat.fontfamily=c("sans","sans"))
dev.off()


#Venn Diagram 2
pdf("Venn_diagram_EE.pdf")
temp1 <- intersect(BICEE_M9EE, RCEE_M9EE) ; length(temp1)
temp2 <- intersect(BICEE_M9EE, RAEE_M9EE) ; length(temp2)
temp3 <- intersect(RCEE_M9EE, RAEE_M9EE) ; length(temp3)
temp123 <- intersect(temp1,temp2); length(temp123)
draw.triple.venn(length(BICEE_M9EE),length(RCEE_M9EE),length(RAEE_M9EE),length(temp1),length(temp3),length(temp2),length(temp123), 
                 category=c("BICEE_vs_M9EE","RCEE_vs_M9EE","RAEE_vs_M9EE"),fill=c("cadetblue","chocolate1","chartreuse4"), 
                 cat.col=c("cadetblue","chocolate1","chartreuse4"), cat.fontfamily=c("sans","sans","sans"),margin=0.05,cat.dist=c(0.07,0.07,0.05))
dev.off()


#Venn Diagram 3
pdf("Venn_diagram_EM_vs_EE.pdf")
temp1 <- intersect(BICEM_BICEE, RAEM_RAEE) ; length(temp1)
temp2 <- intersect(RAEM_RAEE, M9EM_M9EE) ; length(temp2)
temp3 <- intersect(BICEM_BICEE, M9EM_M9EE) ; length(temp3)
temp123 <- intersect(temp1,temp2); length(temp123)
draw.triple.venn(length(BICEM_BICEE),length(RAEM_RAEE),length(M9EM_M9EE),length(temp1),length(temp2),length(temp3),length(temp123), 
                 category=c("BICEM_vs_BICEE","RAEM_vs_RAEE","M9EM_vs_M9EE"),fill=c("coral","deepskyblue","chartreuse3"), 
                 cat.col=c("coral","deepskyblue","chartreuse3"), cat.fontfamily=c("sans","sans","sans"),margin=0.05,cat.dist=c(0.07,0.07,0.05),main="End of exponential phase (EE)\n versus middle of exponental phase (EM)")
dev.off()

