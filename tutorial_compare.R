#### Tutorial comparison ----

#Going to integrate the key elements from both the Chen et al. and Love et al. tutorials to see how they differ in their gene comparisons.

#Chen, Yunshun, A. T. L. Lun, and G. K. Smyth. 2016. “From Reads to Genes to Pathways: Differential Expression Analysis of RNA-Seq Experiments Using Rsubread and the edgeR Quasi-Likelihood Pipeline.” F1000Research 5:1438.

#Love, Michael, I., S. Anders, V. Kim & W. Huber. 2019. “RNA-seq workflow: gene-level exploratory analysis and differential expression.” F1000Research 5:1438.
#https://www.bioconductor.org/packages/devel/workflows/vignettes/rnaseqGene/inst/doc/rnaseqGene.html

#### Acquire Data ----
setwd("/Users/lisa/Documents/Bioinformatics/6210/A5/scripts_data/a5_data/")

#This package is for downloading the requisite file (I think), do not need to load if using saved file.
#BiocManager::install("RnaSeqGeneEdgeRQL")
library(RnaSeqGeneEdgeRQL)

#I'll follow the same method as the Chen et al. tutorial to download the count data from my experiment

#Count_URL <- paste("http://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE114702", "format=file", "file=GSE114702_HTSeq_gene_counts.txt.gz", sep="&")
#download.file(Count_URL, "GSE114702_HTSeq_gene_counts.txt.gz")

#Read in file; first column is the gene identifiers, which will be used for row names
countdata_mice <- read.delim("GSE114702_HTSeq_gene_counts.txt.gz", row.names = 1)
head(countdata_mice, 3)
class(countdata_mice)
colnames(countdata_mice)
head(rownames(countdata_mice))

#To start, I plan on comparing the naive CON and GF mice. If there is time, I will add the social interaction CON and GF mice.

col_P_CON_GF <- grep("^CON.P|^GF.P", colnames(countdata_mice))
counts_P_CON_GF <- countdata_mice[,col_P_CON_GF]
colnames(counts_P_CON_GF)

rm(rows_P_CON_GF)


#### Formatting counts for analysis ----

#First, I will specify the number of significant figures to be displayed on the screen for all of the outputs. This is to ensure readability of results.
options(digits = 3)

#Specify groups for samples at factor level

library(stringr)
groups <- str_remove(colnames(counts_P_CON_GF), pattern = ".P[0-9]+")
groups <- factor(groups)
groups
table(groups)

#Next Love et al. requires a dataframe of column information for the construction of the DESeqDataSet object
coldata_mice <- data.frame(Group = groups)
class(coldata_mice)

#Similarly, Chen et al. requires a dataframe of gene information for an annotation component to the DGEList object. Since the example data in the tutorial had length information to use, they started with just length, and added annotations later, but I will extract some annotations first to set up the dataframe. I'm going to add Entrez gene IDs for each ENSEMBL gene ID, as well as gene symbols

#BiocManager::install("org.Mm.eg.db")
library(org.Mm.eg.db)

entrez_id_anot <-mapIds(org.Mm.eg.db, rownames(counts_P_CON_GF), keytype = "ENSEMBL", column = "ENTREZID")

symbols_anot <-mapIds(org.Mm.eg.db, rownames(counts_P_CON_GF), keytype = "ENSEMBL", column = "SYMBOL")

#for assignment submission, especially if I want to add any other mappings, consider creating a function or otherwise simplifying these repeated lines

gene_mice <- data.frame(Entrez_ID = entrez_id_anot, Symbol = symbols_anot)
head(gene_mice, 3)
class(gene_mice)

#The Chen et al. tutorial filters out any genes that do not have a corresponding Symbol, but I will leave it for now.

#Chen et al.'s tutorial uses EdgeR, and creates a DGEList object:
library(edgeR)

DGE_Counts <- DGEList(counts_P_CON_GF, group = groups, genes = gene_mice)
class(DGE_Counts)
DGE_Counts$samples
head(DGE_Counts$counts)
head(DGE_Counts$genes)


#Love et al.'s tutorial uses DESeq2, and creates a DESeqDataSet object here. For design, the tutorial has a few suggestions, but I will try ~0 + groups because that is the formula used in creating the design matrix in the Chen et al. tutorial
library(DESeq2)
DDS_Counts <- DESeqDataSetFromMatrix(countData = counts_P_CON_GF, colData = coldata_mice, design = ~ 0 + Group)
class(DDS_Counts)

#The design matrix for Chen et al. is held in a separate variable.
design_m <- model.matrix(~0+groups)
colnames(design_m) <- levels(groups)
design_m


#### Data Filtering  ----

#The EdgeR tutorial uses the following default parameters for filtering: min.count=10 (minimum numeric count for at least some samples), min.total.count = 15, also, min.prop=0.7. These are variables I can consider adjusting to see how they impact the results

#The result is a logical vector specifying true or false for if each gene meets the parameters considred.
filt_data_dge <- filterByExpr(DGE_Counts, design_m)
table(filt_data_dge) #TRUE=14811
head(filt_data_dge)

DGE_Filt <- DGE_Counts[filt_data_dge, , keep.lib.sizes=FALSE]

#The DESeq2 tutorial filters out genes that have 1 or fewer reads total. The authors describe this only as a pre-filtering step to eliminate empty and nearly empty row, with further filtering occuring downstream.
keep_rows_dds <- rowSums(counts(DDS_Counts)) > 1
DDS_Filt <- DDS_Counts[keep_rows_dds,]
nrow(DDS_Filt) #keeps 23762

#I could try changing the count to >= 15, which would be similar to the min.total.count argument above

keep_rows_dds_15 <- rowSums(counts(DDS_Counts)) >= 15
DDS_Filt_15 <- DDS_Counts[keep_rows_dds_15,]
nrow(DDS_Filt_15) #Keeps 18526

#This is much more similar. The difference is likely coming from the min.count argument, or the min.prop argument, which are not specified in this format for filtering.

#Alternately, I could try changing the parameters for the filter on the EdgeR filtering
#Since the Stilling et al. paper that is the source of the data set did not specify their filtering parameters, but did mention using DESeq2 with default parameters for their analysis, I will filter the EdgeR set by min count of 2 as well.

filt_data_dge_2 <- filterByExpr(DGE_Counts, design_m, min.count = 2, min.total.count = 2, min.prop = 0)
table(filt_data_dge_2) #TRUE = 18300

#This is still much less than the comparative filter in DESeq2. This could be because the function for filtering in EdgeR takes into account library sizes and the design of the experiment provided.
#For the purposes of this exploration, I could use the same filtering format as used for DESeq, but with the EdgeR data
keep_rows_dge <- rowSums(DGE_Counts$counts) > 1
DGE_Filt_2 <- DGE_Counts[keep_rows_dge,]
nrow(DGE_Filt_2) #23762

#This achieves the same result.

#The Chen et al. tutorial for EdgeR also suggests filtering by log counts per million as an alternative in order to take into account library sizes when filtering.

AveLogCPM_DGE <- aveLogCPM(DGE_Counts)
hist(AveLogCPM_DGE)

#recommendations about filtering found at: https://support.bioconductor.org/p/108798/
keep_test_dds <- rowSums(cpm(DDS_Counts) > 0.5) >= 2
keep_test_dge <- rowSums(cpm(DGE_Counts) > 0.5) >= 2
length(keep_test_dds)
length(keep_test_dge)

#I think I will skip this for now, and keep the matching sized sets of:
#DDS_Filt for the DESeq2 data set; and
#DGE_Filt_2 for the EdgeR data set

#Both tutorials then do normalization on the library sizes.

#In EdgeR, this is done by TMM (trimmed mean of M) normalization, and the information is updated directly into the DGEList object.
DGE_Filt_N <- calcNormFactors(DGE_Filt_2)
DGE_Filt_N$samples

#DESeq uses "median ratio method" to estimate size factors
DDS_Filt <- estimateSizeFactors(DDS_Filt)
DDS_Filt$sizeFactor

#The approximate values and relative values of the factors are comparable, even if the different methods produce different results.
#A discussion on the different techniques for normalization, and why they aren't directly comparable:
#https://support.bioconductor.org/p/51479/
#http://www.nathalievialaneix.eu/doc/html/TP1_normalization.html


#### Data Visualization ----

#The Chen et al. tutorial for EdgeR  does plots for MDS for all samples together, MD for each sample, BCV for dispersion estimates, and QL dispersion on fitted GLM coefficients

#The Love et al. tutorial for DESeq2 does heatmap on data that has been transformed with variance stabilization, PCA plot on the transformed data, and MDS plot on transformed data

#For the direct comparison, let's start with the MDS plot as done in each tutorial

#EdgeR:
groups
pch_groups <- c(16,17)
colors_groups <-c("pink", "blue")
plotMDS(DGE_Filt_N, col=colors_groups[groups], pch=pch_groups[groups])
legend("bottomright", legend = levels(groups), pch=pch_groups, col=colors_groups, ncol = 2)

#DESeq2
library(ggplot2)
#The Love et al. DESeq2 tutorial recommends either variance stabilizing transformation, or rlog, with the authors favouring rlog for smaller data sets. Since I'm only working with 8 samples in this set, let's try rlog:
DDS_Filt_T <- rlog(DDS_Filt, blind = FALSE)
head(assay(DDS_Filt_T), 3)
sampleDists_vst <- dist(t(assay(DDS_Filt_T)))
sampleDists_vst
sampleDistMatrix_vst <- as.matrix(sampleDists_vst)
DDS_MDS <- as.data.frame(colData(DDS_Filt_T)) %>%
  cbind(cmdscale(sampleDistMatrix_vst))

ggplot(DDS_MDS, aes(x = `1`, y = `2`, color = groups, shape = groups)) +
  geom_point(size = 3) + coord_fixed() + ggtitle("MDS with Rlog data")

#Both have the GF mice grouped closer together than the CON mice, but the scales and appearances are quite different.
#This is worth looking into more, to explain what the graphs are showing in each case

#Come back to this section later to get into some of the other visualizations that show the structure of the data


#### Differential Expression Analysis ----

#EdgeR
EdgeR.CONvsGF <- makeContrasts(CON-GF, levels = design_m)
class(EdgeR.CONvsGF)

DGE_Disp <- estimateDisp(DGE_Filt_N, design_m, robust=TRUE)

fit_dgeglm <- glmQLFit(DGE_Disp, design_m, robust=TRUE)
res_ER.CONvsGF <- glmQLFTest(fit_dgeglm, contrast = EdgeR.CONvsGF)
class(res_ER.CONvsGF)
topTags(res_ER.CONvsGF)


#DESeq
DESeq.CONvsGF <- DESeq(DDS_Filt)

res_DS.CONvsGF <- results(DESeq.CONvsGF)
summary(res_DS.CONvsGF)

resSig_mice <- subset(res_DS.CONvsGF, padj < 0.1)
head(resSig_mice[ order(resSig_mice$log2FoldChange), ])
head(resSig_mice[ order(resSig_mice$log2FoldChange, decreasing = TRUE), ]) 
