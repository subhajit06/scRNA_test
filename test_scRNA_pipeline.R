##### testing of scRNA code
library(DropletUtils)
library(biomaRt)
library(dplyr)
library(scater)

setwd("~/Desktop/ToDo/CPG/JD/scRNASeq_test/")

path1 = "~/Desktop/ToDo/CPG/JD/scRNASeq_test/"
path2 = paste0(path1, "matrices_mex/hg19/")
sce   = read10xCounts(path2, col.names=TRUE)

sce

anno.file = "~/Desktop/ToDo/CPG/JD/scRNASeq_test/gene.annoation.rds"
if(file.exists(anno.file)){
  gene.annotation = readRDS(anno.file)
}else{
  ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl")
  
  attr.string = c('ensembl_gene_id', 'hgnc_symbol', 'chromosome_name')
  attr.string = c(attr.string, 'start_position', 'end_position', 'strand')
  attr.string = c(attr.string, 'description', 'percentage_gene_gc_content')
  attr.string = c(attr.string, 'gene_biotype')
  
  rowData(sce)[1:2,]
  gene.annotation = getBM(attributes=attr.string, 
                          filters =  'ensembl_gene_id', 
                          values = rowData(sce)$ID, 
                          mart = ensembl)
}

dim(gene.annotation)

gene.annotation[1:2,]

t1 = table(gene.annotation$ensembl_gene_id)
t2 = t1[t1 > 1]
t2 

gene.annotation[which(gene.annotation$ensembl_gene_id %in% names(t2)),]
gene.annotation = distinct(gene.annotation, ensembl_gene_id, 
                           .keep_all = TRUE)

dim(gene.annotation)

#####
gene.annotation[1:2,]

table(gene.annotation$chromosome_name)
table(gene.annotation$gene_biotype)

## some genes do not have annotation because their ids are retired
gene.missing = setdiff(rowData(sce)$ID, gene.annotation$ensembl_gene_id)
length(gene.missing)

gene.missing[1:6]

w2kp = match(gene.annotation$ensembl_gene_id, rowData(sce)$ID)
sce  = sce[w2kp,]
dim(sce)

table(gene.annotation$ensembl_gene_id == rowData(sce)$ID)
rowData(sce)  = gene.annotation
rownames(sce) = uniquifyFeatureNames(rowData(sce)$ensembl_gene_id, 
                                     rowData(sce)$hgnc_symbol)

######
bcrank = barcodeRanks(counts(sce))

# Only showing unique points for plotting speed.
uniq = !duplicated(bcrank$rank)

pdf("./plots/inflection_knee.pdf")
par(mar=c(5,4,2,1), bty="n")
plot(bcrank$rank[uniq], bcrank$total[uniq], log="xy", 
     xlab="Rank", ylab="Total UMI count", cex=0.5, cex.lab=1.2)

abline(h=bcrank$inflection, col="darkgreen", lty=2)
abline(h=bcrank$knee, col="dodgerblue", lty=2)

legend("left", legend=c("Inflection", "Knee"), bty="n", 
       col=c("darkgreen", "dodgerblue"), lty=2, cex=1.2)

dev.off()

bcrank$inflection
bcrank$knee

summary(bcrank$total)
table(bcrank$total >= bcrank$knee)
table(bcrank$total >= bcrank$inflection)

set.seed(100)
e.out = emptyDrops(counts(sce))
e.out

is.cell = (e.out$FDR <= 0.01)

pdf("./plots/UMI_and_prob.pdf")
par(mar=c(5,4,1,1), mfrow=c(1,2), bty="n")
plot(e.out$Total, -e.out$LogProb, col=ifelse(is.cell, "red", "black"),
     xlab="Total UMI count", ylab="-Log Probability", cex=0.2)
abline(v = bcrank$inflection, col="darkgreen")
abline(v = bcrank$knee, col="dodgerblue")
legend("bottomright", legend=c("Inflection", "Knee"), bty="n", 
       col=c("darkgreen", "dodgerblue"), lty=1, cex=1.2)

plot(e.out$Total, -e.out$LogProb, col=ifelse(is.cell, "red", "black"),
     xlab="Total UMI count", ylab="-Log Probability", cex=0.2, xlim=c(0,2000), ylim=c(0,2000))
abline(v = bcrank$inflection, col="darkgreen")
abline(v = bcrank$knee, col="dodgerblue")
dev.off()

table(colnames(sce) == rownames(e.out))
table(e.out$FDR <= 0.01, useNA="ifany")

table(is.cell, e.out$Total >= bcrank$inflection)

w2kp = which(is.cell & e.out$Total >= bcrank$inflection)
sce = sce[,w2kp]
dim(sce)

library(data.table)
ribo.file = "~/Desktop/ToDo/CPG/JD/scRNASeq_test/ribosome_genes.txt"
ribo = fread(ribo.file)
dim(ribo)

ribo[1:2,]
is.mito = which(rowData(sce)$chromosome_name == "MT")
is.ribo = which(rowData(sce)$hgnc_symbol %in% ribo$`Approved Symbol`)
length(is.mito)
length(is.ribo)

sce = calculateQCMetrics(sce, feature_controls=list(Mt=is.mito, Ri=is.ribo))
colnames(colData(sce))

pdf("./plots/ribo_MT.pdf")
par(mfrow=c(2,2), mar=c(5, 4, 1, 1), bty="n")
hist(log10(colData(sce)$total_counts), xlab="log10(Library sizes)", main="", 
     breaks=20, col="grey80", ylab="Number of cells")

hist(log10(colData(sce)$total_features), xlab="log10(# of expressed genes)", 
     main="", breaks=20, col="grey80", ylab="Number of cells")

hist(colData(sce)$pct_counts_Ri, xlab="Ribosome prop. (%)",
     ylab="Number of cells", breaks=40, main="", col="grey80")

hist(colData(sce)$pct_counts_Mt, xlab="Mitochondrial prop. (%)", 
     ylab="Number of cells", breaks=80, main="", col="grey80")

dev.off()

par(mfrow=c(1,3), mar=c(5,4,1,1))
hist(log10(rowData(sce)$mean_counts+1e-6), col="grey80",  main="", 
     breaks=40, xlab="log10(ave # of UMI + 1e-6)")
hist(log10(rowData(sce)$n_cells_counts+1), col="grey80", main="", 
     breaks=40, xlab="log10(# of expressed cells + 1)")
smoothScatter(log10(rowData(sce)$mean_counts+1e-6), 
              log10(rowData(sce)$n_cells_counts + 1), 
              xlab="log10(ave # of UMI + 1e-6)", 
              ylab="log10(# of expressed cells + 1)")

