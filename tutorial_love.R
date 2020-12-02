#### A5 Tutorial by Love et al. ----

#Starting by practicing skills and packages, following along with tutorial:
#Love, Michael, I., S. Anders, V. Kim & W. Huber. 2019. “RNA-seq workflow: gene-level exploratory analysis and differential expression.” F1000Research 5:1438.
#https://www.bioconductor.org/packages/devel/workflows/vignettes/rnaseqGene/inst/doc/rnaseqGene.html


#Using dataset from "Social interaction-induced activation of RNA splicing in the amygdala of microbiome-deficient mice" from 2018 paper by Stilling et al.

#### Obtain & Format the Data ---- 

setwd("/Users/lisa/Documents/Bioinformatics/6210/A5/scripts_data/a5_data/")

#I'll follow the same method as the Chen et al. tutorial to download the count data from my experiment

#Count_URL <- paste("http://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE114702", "format=file", "file=GSE114702_HTSeq_gene_counts.txt.gz", sep="&")
#download.file(Count_URL, "GSE114702_HTSeq_gene_counts.txt.gz")

#Read in file; first column is the gene identifiers, which will be used for row names
countdata_mice <- read.delim("GSE114702_HTSeq_gene_counts.txt.gz", row.names = 1)
head(countdata_mice)
class(countdata_mice)
colnames(countdata_mice)

#Like in the Chen et al. tutorial, I will go down to 2 bioloigcal replicates per condition. I'm subsetting by column numbers for the first two of each type, but like mentioned in the other tutorial, I would like to troubleshoot using a sampling method to randomly select.

countdata_mice2 <- countdata_mice[,c(1,2,5,6,13,14,17,18,21,22,29,30)]

#Section 3.2 of this tutorial details starting from count matrices, so the preceding sections about formatting from raw data or other formats can be skipped.
#If I want to try following along the processing, I can download sequencing data from:  
#https://trace.ncbi.nlm.nih.gov/Traces/study/?acc=SRP148566&o=acc_s%3Aa

library(DESeq2)

#First need to format count matrix to a DESeqDataSet object. This will require the count data defined above, as well as arguments "colData" and "design"
#From function's documentation: "colData: for matrix input: a DataFrame or data.frame with at least a single column. Rows of colData correspond to columns of countData"
#In Section 2.5, colData(gse) is a dataframe with one sample per row, and columns describing the donor and treatment condition. I can construct a dataframe to describe which group each sample belongs to:

#From the function's documentation, design can be a formula or a matrix, and "the formula expresses how the counts for each gene depend on the variables in colData." The tutorial uses the forumla "~ cell + dex". Variables must be columns in colData. Looking into the tutorial, dex is the experimental condition. it is a column that is factor level, and lists "untreated" and "dexamethasone" for treatment type. The cell column is the donor. Each cell donor has an untreated and a treated. So looks like they are modeling if treatment has an impact on the cell. 
#From the tutorial: "The simplest design formula for differential expression would be ~ condition , where condition is a column in colData(dds) that specifies which of two (or more groups) the samples belong to. For the airway experiment, we will specify ~ cell + dex meaning that we want to test for the effect of dexamethasone ( dex ) controlling for the effect of different cell line ( cell )"
#So, I will try making the coldata dataframe to have separate columns for social condition (naive, social interaction), mouse type (conventional, germ-free, and ex-germ-free), and for group. I will try using the design as ~ Type + Interaction to see how the different mice types have response to social interaction.
#Alternatively could just do ~ Group 

library(stringr)
samples_mice <-colnames(countdata_mice2)
social_mice <- str_extract(string = samples_mice, pattern = "[SP]")
social_mice <- str_replace(social_mice, "S", "social")
social_mice <- str_replace(social_mice, "P", "naive")
social_mice <- factor(social_mice)
type_mice <- factor(str_extract(string = samples_mice, pattern = "^[A-z]+"))
groups_mice <- factor(str_extract(string = samples_mice, pattern = "^[A-z]+.[SP]"))
coldata_mice <- data.frame(Sample = samples_mice, Interaction = social_mice, Type = type_mice, Group = groups_mice)
class(coldata_mice)
coldata_mice
#lot of repetition in the way I coded this. use piping or other simplifying features if putting in assignment

dds_mice <- DESeqDataSetFromMatrix(countData = countdata_mice2, colData = coldata_mice, design = ~ Type + Interaction)


#### Data Exploration and Visualization ----

#First, the tutorial removes rows with 1 or fewer fragments across samples
nrow(dds_mice)
keep_rows <- rowSums(counts(dds_mice)) > 1
dds_mice <- dds_mice[keep_rows,]
nrow(dds_mice)
#They retained 54% of their samples, mine retains 56%

#The tutorial did variance stabilizing next, and recommends VST for samples >30, and rlog for samples less than 30. Since this data set has n=40, I will use VST

vsd_mice <- vst(dds_mice, blind = FALSE)
head(assay(vsd_mice), 3)
colData(vsd_mice)

#From the tutorial: "For differential testing we recommend the DESeq function applied to raw counts, as described later in this workflow, which also takes into account the dependence of the variance of counts on the mean value during the dispersion estimation step."
#Since I am using the counts for differential testing, I may not need this transformation.

#The tutorial explored overall distance between samples next:
#They did the sample distances on the VST object, because it was normalized. Let's see what the difference is for using each

sampleDists_dds <- dist(t(assay(dds_mice)))
sampleDists_dds

sampleDists_vst <- dist(t(assay(vsd_mice)))
sampleDists_vst

#The ratio difference between the top two numbers is quite large. There must be a large difference in library size between the two. A bit more research on the package and different normalization methods is probably a good idea. 

#They next visualize the differences with a heat map.

library("pheatmap")
library("RColorBrewer")

# They use the transformed data and its corresponding distances calculated above, so I will try with both transformed and not to see
#I'm using my Interaction information (Social/Naive) in place of the "dex" treatment information, and Type (CON/exGF/GF) in place of cell type
sampleDistMatrix_dds <- as.matrix(sampleDists_dds)
rownames(sampleDistMatrix_dds) <- paste(dds_mice$Interaction, dds_mice$Type, sep = " - " )
colnames(sampleDistMatrix_dds) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix_dds,
         clustering_distance_rows = sampleDists_dds,
         clustering_distance_cols = sampleDists_dds,
         col = colors)

sampleDistMatrix_vst <- as.matrix(sampleDists_vst)
rownames(sampleDistMatrix_vst) <- paste(vsd_mice$Interaction, vsd_mice$Type, sep = " - " )
colnames(sampleDistMatrix_vst) <- NULL
pheatmap(sampleDistMatrix_vst,
         clustering_distance_rows = sampleDists_vst,
         clustering_distance_cols = sampleDists_vst,
         col = colors)

#A lot lighter colour with the transformed one!
#mostly the naive mice are grouped together and the social together, with the exception of the one social GF mouse.
#mostly biological replicates are grouped together as well, but there are some differences.

#I'm going to skip the Poisson distances, but if it seems relevant, I will research the statistics a bit more and com back to it later

#I see a lot of PCA in papers, so I will try that visualization as well.
#I did not work on the dds object, only on the VST transformed one. 

plotPCA(vsd_mice, intgroup = c("Interaction", "Type"))
pcaData_mice <- plotPCA(vsd_mice, intgroup = c( "Interaction", "Type"), returnData = TRUE)
pcaData_mice
#Should have saved images from other tutorial so I could compare.
#Also, use the formatting of group markers from Chen et al. tutorial so I can easier compare differenct conditions


percentVar <- round(100 * attr(pcaData_mice, "percentVar"))

library(ggplot2)
ggplot(pcaData_mice, aes(x = PC1, y = PC2, color = Interaction, shape = Type)) +
  geom_point(size =3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  coord_fixed() +
  ggtitle("PCA with VST data")

#This looks much better as a visualization

## MDS plot visualization
mds_mice <- as.data.frame(colData(vsd_mice)) %>%
  cbind(cmdscale(sampleDistMatrix_vst))

ggplot(mds_mice, aes(x = `1`, y = `2`, color = Interaction, shape = Type)) +
  geom_point(size = 3) + coord_fixed() + ggtitle("MDS with VST data")


#### Differential Expression Analysis ----

#Even though doing some of the visualisations for quality control and preliminary analysis on the transformed data, the next step of creating the DESeq object is done on the DESeqDataSet object, not the transformed vsd. Since the DESeq command does the appropriate transformations within its algorithm, that earlier step was just for preliminary visualizations, and is no longer needed. The DESeq command does not even work with the vsd object.

des_mice <- DESeq(dds_mice)


#Building a results table

res_mice <- results(des_mice)
res_mice

#this results table models social vs naive mice, because the original design was listed as ~Type + Interaction. The tutorial also describes the command below as doing the same as not specifying contrast
#res <- results(dds, contrast=c("dex","trt","untrt")) 
#which shows that I could choose a different contrast to model.

#This compares GF mice to CON mice
res_con_gf <- results(des_mice, contrast = c("Type", "GF", "CON"))
res_con_gf

#To view the metadata of the columns, can use mcol
mcols(res_mice, use.names = TRUE)

#Can also check summary
summary(res_mice)
summary(res_con_gf)



#Let's see how it goes if I do the initial design by Group instead

dds_group <- DESeqDataSetFromMatrix(countData = countdata_mice2, colData = coldata_mice, design = ~ Group)

keep_row_g <- rowSums(counts(dds_group)) > 1
dds_group <- dds_group[keep_row_g,]

vsd_group <- vst(dds_group, blind = FALSE)
head(assay(vsd_group), 3)
colData(vsd_group)

sampleDists_group <- dist(t(assay(vsd_group)))
sampleDists_group

sampleDistMatrix_group <- as.matrix(sampleDists_group)
rownames(sampleDistMatrix_group) <- paste(vsd_group$Interaction, vsd_group$Type, sep = " - " )
colnames(sampleDistMatrix_group) <- NULL
pheatmap(sampleDistMatrix_group,
         clustering_distance_rows = sampleDists_group,
         clustering_distance_cols = sampleDists_group,
         col = colors)

pcaData_group <- plotPCA(vsd_group, intgroup = c( "Interaction", "Type"), returnData = TRUE)
pcaData_group
percentVar_f <- round(100 * attr(pcaData_group, "percentVar"))
library(ggplot2)
ggplot(pcaData_group, aes(x = PC1, y = PC2, color = Interaction, shape = Type)) +
  geom_point(size =3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  coord_fixed() +
  ggtitle("PCA with VST data")


#Looks similar?

mds_group <- as.data.frame(colData(vsd_group)) %>%
  cbind(cmdscale(sampleDistMatrix_group))

ggplot(mds_group, aes(x = `1`, y = `2`, color = Interaction, shape = Type)) +
  geom_point(size = 3) + coord_fixed() + ggtitle("MDS with VST data")

##Now for results

des_group <- DESeq(dds_group)


#Building a results table

res_group <- results(des_group)
res_group

#This compares GF mice to CON mice within the social interaction category
res_con_gf_g <- results(des_group, contrast = c("Group", "GF.S", "CON.S"))
res_con_gf_g

#The Stilling et al. article says they used 10% FDR adjusted p-values to determine significant genes:

sum(res_mice$padj < 0.1, na.rm=TRUE)
sum(res_group$padj < 0.1, na.rm=TRUE)
sum(res_con_gf$padj < 0.1, na.rm=TRUE)
sum(res_con_gf_g$padj < 0.1, na.rm=TRUE)

#Different modeling and comparisons show different amounts of genes considered significant. 

#Let's see which genes are most significantly down and up regulated for each model above

# ~Interaction + Type: total
resSig_mice <- subset(res_mice, padj < 0.1)
head(resSig_mice[ order(resSig_mice$log2FoldChange), ])
head(resSig_mice[ order(resSig_mice$log2FoldChange, decreasing = TRUE), ]) 

# ~Interaction + Type: con vs gf
resSig_con_gf <- subset(res_con_gf, padj < 0.1)
head(resSig_con_gf[ order(resSig_con_gf$log2FoldChange), ])
head(resSig_con_gf[ order(resSig_con_gf$log2FoldChange, decreasing = TRUE), ]) 

# ~Group: total
resSig_group <- subset(res_group, padj < 0.1)
head(resSig_group[ order(resSig_group$log2FoldChange), ])
head(resSig_group[ order(resSig_group$log2FoldChange, decreasing = TRUE), ]) 

# ~Group: social con vs gf
resSig_con_gf_g <- subset(res_con_gf_g, padj < 0.1)
head(resSig_con_gf_g[ order(resSig_con_gf_g$log2FoldChange), ])
head(resSig_con_gf_g[ order(resSig_con_gf_g$log2FoldChange, decreasing = TRUE), ]) 


#### Plot results ----

topGene_mice <- rownames(res_mice)[which.min(res_mice$padj)]
plotCounts(dds_mice, gene = topGene_mice, intgroup=c("Interaction", "Type"))
#I wonder if I can modify this plot to colour. Might need to find ggplot alternative to do so, but I can choose to use one or more than one variable in the intgroup argument.

#The chart above ordered the genes by log2FC, provided that they met a particular p-value threshold, but this chart is plotting the expression counts for the gene with the lowest p-value. So, if I understand correctly, the other function showed me the genes that were most different in expression for a particular comparison, but there is a chance some are false-discoveries, and this shows the gene with the most significant p-value?

topGene_con_gf_g <- rownames(res_con_gf_g)[which.min(res_con_gf_g$padj)]
plotCounts(dds_group, gene = topGene_con_gf_g, intgroup=c("Group"))

#Let's look at a heatmap

BiocManager::install("genefilter")
library("genefilter")

topVarGenes_mice <- head(order(rowVars(assay(vsd_mice)), decreasing = TRUE), 20)
mat_mice <- assay(vsd_mice)[ topVarGenes_mice, ]
mat_mice <- mat_mice - rowMeans(mat_mice)
anno_mice <- as.data.frame(colData(vsd_mice)[, c("Type","Interaction")])
pheatmap(mat_mice, annotation_col = anno_mice)

topVarGenes_group <- head(order(rowVars(assay(vsd_group)), decreasing = TRUE), 20)
mat_group <- assay(vsd_group)[ topVarGenes_group, ]
mat_group <- mat_group - rowMeans(mat_group)
anno_group <- as.data.frame(colData(vsd_group)[, c("Type","Interaction")])
pheatmap(mat_group, annotation_col = anno_group)

