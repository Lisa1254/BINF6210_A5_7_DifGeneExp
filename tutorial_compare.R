#### Tutorial comparison ----

#Going to integrate the key elements from both the Chen et al. and Love et al. tutorials to see how they differ in their gene comparisons.

#Chen, Yunshun, A. T. L. Lun, and G. K. Smyth. 2016. “From Reads to Genes to Pathways: Differential Expression Analysis of RNA-Seq Experiments Using Rsubread and the edgeR Quasi-Likelihood Pipeline.” F1000Research 5:1438.

#Love, Michael, I., S. Anders, V. Kim & W. Huber. 2019. “RNA-seq workflow: gene-level exploratory analysis and differential expression.” F1000Research 5:1438.
#https://www.bioconductor.org/packages/devel/workflows/vignettes/rnaseqGene/inst/doc/rnaseqGene.html

#### Acquire Data ----
setwd("/Users/lisa/Documents/Bioinformatics/6210/A5/scripts_data/a5_data/")

#This package was loaded before getting the metadata from a system file. Unclear from tutorial if needed for anything else. I'll try not loading it.
#BiocManager::install("RnaSeqGeneEdgeRQL")
#library(RnaSeqGeneEdgeRQL)

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

rm(col_P_CON_GF)


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

#The result is a logical vector specifying true or false for if each gene meets the parameters considered.
filt_data_dge <- filterByExpr(DGE_Counts, design_m)
table(filt_data_dge) #TRUE=14811
head(filt_data_dge)

DGE_Filt <- DGE_Counts[filt_data_dge, , keep.lib.sizes=FALSE]

#The DESeq2 tutorial filters out genes that have 1 or fewer reads total. The authors describe this only as a pre-filtering step to eliminate empty and nearly empty row, with further filtering occurring downstream.
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
DGE_Filt_1 <- DGE_Counts[keep_rows_dge,]
nrow(DGE_Filt_1) #23762

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
#DGE_Filt_1 for the EdgeR data set

#Both tutorials then do normalization on the library sizes.

#In EdgeR, this is done by TMM (trimmed mean of M) normalization, and the information is updated directly into the DGEList object.
DGE_Filt_N <- calcNormFactors(DGE_Filt_1)
DGE_Filt_N$samples

#DESeq uses "median ratio method" to estimate size factors
DDS_Filt_N <- estimateSizeFactors(DDS_Filt)
DDS_Filt_N$sizeFactor

#The approximate values and relative values of the factors are comparable, even if the different methods produce different results.
#A discussion on the different techniques for normalization, and why they aren't directly comparable:
#https://support.bioconductor.org/p/51479/
#http://www.nathalievialaneix.eu/doc/html/TP1_normalization.html


#### Data Visualization ----

#The Chen et al. tutorial for EdgeR  does plots for MDS for all samples together, MD for each sample, BCV for dispersion estimates, and QL dispersion on fitted GLM coefficients

#The Love et al. tutorial for DESeq2 does heatmap on data that has been transformed with variance stabilization, PCA plot on the transformed data, and MDS plot on transformed data

#For the direct comparison, let's start with the MDS plot as done in each tutorial
#Refering to http://www.stat.columbia.edu/~tzheng/files/Rcolor.pdf for colour options

#EdgeR:
groups
pch_groups <- c(16,17)
colors_groups <-c("orchid2", "blue")
MDS_DGE.out <- plotMDS(DGE_Filt_N, col=colors_groups[groups], pch=pch_groups[groups])
legend("bottomright", legend = levels(groups), pch=pch_groups, col=colors_groups, ncol = 2)
MDS_DGE.out$cmdscale.out

#I can use the DDS object for the same plotMDS command to make a comparable plot
MDS_DDS.out <- plotMDS(DDS_Filt_N, col=colors_groups[groups], pch=pch_groups[groups])
legend("bottomright", legend = levels(groups), pch=pch_groups, col=colors_groups, ncol = 2)
MDS_DDS.out$cmdscale.out

#They look quite similar, with the differences likely being attributed to the difference in library size normalization.
#From the documentation on plotMDS https://www.rdocumentation.org/packages/edgeR/versions/3.14.0/topics/plotMDS.DGEList
#The default method of using logFC can be changed to BCV values. The log2FC values are calculated using counts per million and Euclidean distance
#The default use of top 500 genes can also be changed. A recommended number is to be roughly the number of genes expected to have a large fold change.

#DESeq2

#The Love et al. DESeq2 tutorial recommends either variance stabilizing transformation or rlog prior to making their MDS plot with a distance matrix and ggplot. The authors suggest rlog for smaller data sets. Let's try a few different options, and see how it impacts the plot:
library(ggplot2)
#install.packages("gridExtra")
library(gridExtra)

#Start with rlog transformation
DDS_Filt_Tr <- rlog(DDS_Filt_N, blind = FALSE)

dists_dds_rlog <- dist(t(assay(DDS_Filt_Tr)), diag = TRUE, upper = TRUE)
#Note: assay() takes "a SummarizedExperiment object" (the DESeqDataSet), and extracts the count information;
#t() transposes the matrix;
#dist() "computes the distances between the rows of a data matrix". The default method is euclidian
#In the tutorial, they followed this up with an as.matrix command to use a separate object to store the symmetrical matrix, but adding the diag = TRUE & upper = TRUE arguments does the same thing to the original object. The tutorial had separate objects because the distances without completing the matrix was used in a heatmap visualization

DDS_MDS_rlog <- as.data.frame(colData(DDS_Filt_Tr)) %>%
  cbind(cmdscale(dists_dds_rlog))
#This dataframe combines the group and size factor information with the cmdscale product of the distance matrix.
#cmdscale() is "classical multidimensional scaling of a data matrix, also known as principal coordinates analysis"
#A bit of research suggests that the main difference between principle components analysis (PCA) and MDS is that PCA preserves the spread of the data, while MDS preserves the pairwise distances

plot_DDSRlog <- ggplot(DDS_MDS_rlog, aes(x = `1`, y = `2`, color = groups, shape = groups)) +
  geom_point(size = 3) + coord_fixed() + ggtitle("MDS with Rlog data")
plot_DDSRlog

#Compare to the vst 
DDS_Filt_Tv <- vst(DDS_Filt_N, blind = FALSE)
head(assay(DDS_Filt_Tv), 3)

dists_dds_vst <- dist(t(assay(DDS_Filt_Tv)), diag = TRUE, upper = TRUE)
DDS_MDS_vst <- as.data.frame(colData(DDS_Filt_Tv)) %>%
  cbind(cmdscale(dists_dds_vst))

plot_DDSvst <- ggplot(DDS_MDS_vst, aes(x = `1`, y = `2`, color = groups, shape = groups)) +
  geom_point(size = 3) + coord_fixed() + ggtitle("MDS with VST data")
plot_DDSvst


#The DGE object does not fit into the rlog or vst functions of the DESeq2 package, but the plot could be viewed with untransformed data

dists_dge <- dist(t(DGE_Filt_N$counts), diag = TRUE, upper = TRUE)
DGE_CMD <- as.data.frame(cmdscale(dists_dge))

plot_DGE <- ggplot(DGE_CMD, aes(x = V1, y = V2, color = groups, shape = groups)) +
  geom_point(size = 3) + coord_fixed() + ggtitle("MDS with Untransformed data (DGE)")
plot_DGE

#The untransformed DDS for comparison:
dists_dds <- dist(t(assay(DDS_Filt_N)), diag = TRUE, upper = TRUE)
DDS_CMD <- as.data.frame(cmdscale(dists_dds))

plot_DDS <- ggplot(DDS_CMD, aes(x = V1, y = V2, color = groups, shape = groups)) +
  geom_point(size = 3) + coord_fixed() + ggtitle("MDS with Untransformed data (DDS)")
plot_DGE

grid.arrange(plot_DDS, plot_DGE, plot_DDSRlog, plot_DDSvst, nrow = 2, top = "MDS plots with transformations on DESeqDataSet and DGEList")

#Although the transformation changes the scales and specific values, it is clear that the general structure of the two sets is the same here, with the GF mice being much more clustered than the CON mice, and the GF mice more to the left of the plot. 
#This is worth looking into more, to explain what the graphs are showing in each case


#Might be worth checking out the preliminary heatmap, as used in the DESeq2 tutorial

#Come back to this section later to get into some of the other visualizations that show the structure of the data


#### Differential Expression Analysis ----

#EdgeR
#First constructs a contrast matrix. This step seems most relevant when there are multiple possible contrasts in the complete data set, but as I have only retained naive mice for CON and GF types, it is basically just reformatting the existing design
EdgeR.CONvsGF <- makeContrasts(CON-GF, levels = design_m)
class(EdgeR.CONvsGF)
EdgeR.CONvsGF

#Estimates dispersion by maximizing negative binomial likelihood
DGE_Disp <- estimateDisp(DGE_Filt_N, design_m, robust=TRUE)
#The dispersion can be used in a plot
plotBCV(DGE_Disp)

#EdgeR uses the dispersion to estimate the GLM for each gene
fit_dgeglm <- glmQLFit(DGE_Disp, design_m, robust=TRUE)
#This can also be plotted
plotQLDisp(fit_dgeglm)

#The results are then calculated using the contrast matrix and the GLM fit
res_ER.CONvsGF <- glmQLFTest(fit_dgeglm, contrast = EdgeR.CONvsGF)
class(res_ER.CONvsGF)
head(res_ER.CONvsGF$table)
topTags(res_ER.CONvsGF)
#The top tags are listed by order of lowest p-value

#Identify the top 10 genes by p-value for comparison with the DESeq method
top_p_ER <- row.names(topTags(res_ER.CONvsGF))
top_p_ER


#DESeq
#The first step does the estimation of size factors, dispersion estimate, and GLM fitting all internally, and thus uses the non-normalized dataset as the input. The reason the tutorial had done data normalization and transformations previously was to be able to do preliminary investigations and visualizations
DESeq.CONvsGF <- DESeq(DDS_Filt)

#The results function extracts a results table from the DESeq object
res_DS.CONvsGF <- results(DESeq.CONvsGF)
class(res_DS.CONvsGF)
head(res_DS.CONvsGF[order(res_DS.CONvsGF$pvalue), ], 10)

#Identify the top 10 genes by p-value for comparison with the EdgeR method
top_p_DS <- row.names(head(res_DS.CONvsGF[order(res_DS.CONvsGF$pvalue), ], 10))
top_p_DS

#The DESeq tutorial uses the padj for filtering, so lets see if the same genes are identified if ordering by adjusted p-value
top_padj_DS <- row.names(head(res_DS.CONvsGF[order(res_DS.CONvsGF$padj), ], 10))
top_padj_DS

top_p_padj <- intersect(top_p_DS, top_padj_DS)
top_p_padj
#Same top 10

top_p_both <- intersect(top_p_DS, top_p_ER)
#This shows that 5/10 genes with the most significant p-values are the same in both methods.

#The DESeq tutorial filters out any genes with a padj of 0.1 or greater to represent a FDR of 10%
DS_p10 <- subset(res_DS.CONvsGF, padj < 0.1)

#Then the top genes upregulated and downregulated are considered by the highest logfold change:
top_down_DS <- row.names(head(DS_p10[order(DS_p10$log2FoldChange), ], 30))
top_up_DS <- row.names(head(DS_p10[order(DS_p10$log2FoldChange, decreasing = TRUE), ], 30))

#Try the same with the EdgeR data:
ER_p10 <-subset(res_ER.CONvsGF$table, PValue < 0.1)
top_down_ER <- row.names(head(ER_p10[order(ER_p10$logFC), ], 30))
top_up_ER <- row.names(head(ER_p10[order(ER_p10$logFC, decreasing = TRUE), ], 30))

top_down_both <- intersect(top_up_DS, top_up_ER)
top_down_both #No overlap

top_up_both <- intersect(top_down_DS, top_down_ER)
top_up_both #No overlap

#This is interesting that the 30 most differentially expressed genes up and down according to the logfold change are unique for each method.

#Both methods use a heat map visualization on the results data: 
#EdgeR tutorial does CPM on dispersion data of the top 30 genes according to p-value, and uses the results in the heatmap
#DESeq2 tutorial only did a heat map with the top 20 genes in terms of variance across the samples, using the vst data, not the results. Following along with the EdgeR tutorial, let's try to do the same