#### A5 Draft work ----

#Starting by practicing skills and packages, following along with tutorial:
#Chen, Yunshun, A. T. L. Lun, and G. K. Smyth. 2016. “From Reads to Genes to Pathways: Differential Expression Analysis of RNA-Seq Experiments Using Rsubread and the edgeR Quasi-Likelihood Pipeline.” F1000Research 5:1438.

#Using dataset from "Social interaction-induced activation of RNA splicing in the amygdala of microbiome-deficient mice" from 2018 paper by Stilling et al.

#### Obtaining & Format the Data ---- 

setwd("/Users/lisa/Documents/Bioinformatics/6210/A5/scripts_data/a5_data/")

#BiocManager::install("RnaSeqGeneEdgeRQL")
library(RnaSeqGeneEdgeRQL)

#First thing the tutorial did was extract a file of metadata from the package to use for annotation of cell types. I'll have to create the groups from the sequence data file, as there does not appear to be a preset file of this information.
#Could possibly figure out how to download metadata or accession list from 
#https://www.ncbi.nlm.nih.gov/Traces/study/?acc=PRJNA472338&o=acc_s%3Aa 

#Next, the tutorial downloads the matrix of sequence counts
#Count_URL <- paste("http://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE114702", "format=file", "file=GSE114702_HTSeq_gene_counts.txt.gz", sep="&")
#download.file(Count_URL, "GSE114702_HTSeq_gene_counts.txt.gz")

#Read in file; first column is the gene identifiers, which will be used for row names
GenewiseCounts <- read.delim("GSE114702_HTSeq_gene_counts.txt.gz", row.names = 1)

dim(GenewiseCounts)
#There are 43629 genes and 40 samples. I might want to subset this data for the project to practice with a more manageable size. Since this data has naive (labeled with "p") and social interaction (labeled with "si") mice for 3 mouse types (CON, ex-GF, GF), this is a similar number of conditions to the tutorial. But the tutorial used 2 biological replicates for each condition, and this data has 4-12. I could subset for 2 replicates per condition.

head(GenewiseCounts)
#The data looks like it loaded in correctly. The second gene on the list (ENSMUSG00000000003) has 0 count for all of the samples. This shows that there will likely be some preprocessing of the counts needed before analysis.

#Now that I have the data, I can set the groups like the tutorial had at the beginning
#The following resource was helpful in understanding what the factor command and object type is:
#https://www.datamentor.io/r-programming/factor/
library(stringr)
groups <-colnames(GenewiseCounts)
class(groups)
groups <- str_extract(string = groups, pattern = "^[A-z]+.[SP]")
groups <- factor(groups)
groups
table(groups)

#Next the counts are formatted into a DGEList object for use in edgeR. It says that the gene argument is optional to create a dataframe for annotations. I like how they added gene symbols into this dataframe later, so I will use the argument, but since I don't have gene length information in my dataset, I will use the gene names as a placeholder to create the dataframe.
library(edgeR)

DGE_Counts <- DGEList(GenewiseCounts, group = groups, genes = row.names(GenewiseCounts))
options(digits = 3)
DGE_Counts$samples
head(DGE_Counts$counts)
head(DGE_Counts$genes)

#Add correct annotations to the "gene" dataframe
#BiocManager::install("org.Mm.eg.db")
library(org.Mm.eg.db)

#Add EntrezID and gene symbol to the annotations frame, removing placeholder column
DGE_Counts$genes$EntrezID <-mapIds(org.Mm.eg.db, rownames(DGE_Counts), keytype = "ENSEMBL", column = "ENTREZID")
DGE_Counts$genes$Symbol <-mapIds(org.Mm.eg.db, rownames(DGE_Counts), keytype = "ENSEMBL", column = "SYMBOL")
#Maybe also add KEGG path link, which could be useful downstream for analysis?
DGE_Counts$genes$Path <-mapIds(org.Mm.eg.db, rownames(DGE_Counts), keytype = "ENSEMBL", column = "PATH")
#This repetitious bit might be a good candidate for an apply function?
DGE_Counts$genes$genes <- NULL
head(DGE_Counts$genes)

#The tutorial dropped genes that don't have a symbol from analysis, but did not explain why. Perhaps just to keep the data manageable for instruction purposes? Instead, here I will subset the first 2 biological replicates per group to be comparable size to the tutorial set.
#If time, troubleshoot randomly sampling from group list

DGE_Counts$samples$group
DGE_Sub_Counts <- DGE_Counts[,c(1,2,5,6,13,14,17,18,21,22,29,30)]
DGE_Sub_Counts$samples

## Design matrix
#This is a design matrix, which encodes as a binary condition which group each sample belongs to. The tutorial uses the group variable, but as I removed some of the samples, I need to redefine group to reflect the current data set.

groups <- factor(DGE_Sub_Counts$samples$group)
groups

#Making the design matrix
#This is how they did it in the tutorial. Not sure what the ~0+ is for. Let's experiment with different numbers. When I did ~1+ it replaced the CON.P column with (Intercept) and all values were equal to 1
design_m <- model.matrix(~0+groups)
colnames(design_m) <- levels(groups)
design_m


## Filtering to remove low counts
#This step is a good candidate for playing with the filters. Documentation doesn't list other defaults when the DGEList class is specified, but the "Default S3 Method" has a few. Default arguments not considered in the tutorial are: min.count=10 (minimum numeric count for at least some samples), min.total.count = 15, large.n = 10 (number of samples per group that is considered to be "large")-- I wonder if this would have been where I could filter my extra replicates. Is that what is meant by samples per group? Also, min.prop=0.7 (minimum proportion of samples in the smallest group that express the gene)

#Start changing min.count=10
filt_data_def <- filterByExpr(DGE_Sub_Counts, design_m)
table(filt_data_def)

filt_data_mincount5 <- filterByExpr(DGE_Sub_Counts, design_m, min.count = 5)
table(filt_data_mincount5)

filt_data_mincount20 <- filterByExpr(DGE_Sub_Counts, design_m, min.count = 20)
table(filt_data_mincount20)

#Bigger numbers in min.count make more "False" for the table, i.e. bigger number means more rows (genes) are being discarded

#Now let's change min.total.count=15
table(filt_data_def)

filt_data_mintotal5 <- filterByExpr(DGE_Sub_Counts, design_m, min.total.count = 5)
table(filt_data_mintotal5)

filt_data_mintotal25 <- filterByExpr(DGE_Sub_Counts, design_m, min.total.count = 25)
table(filt_data_mintotal25)
#Same for all. guess this one is not applicable to this data set

#Let's change large.n=10
table(filt_data_def)

filt_data_largen5 <- filterByExpr(DGE_Sub_Counts, design_m, large.n = 5)
table(filt_data_largen5)

filt_data_largen1 <- filterByExpr(DGE_Sub_Counts, design_m, large.n = 1)
table(filt_data_largen1)

filt_data_largen15 <- filterByExpr(DGE_Sub_Counts, design_m, large.n = 15)
table(filt_data_largen15)
#All the same

#Other possible check to determine filter
AveLogCPM <- aveLogCPM(DGE_Sub_Counts)
hist(AveLogCPM)

#Could simply choos something like if the read count is more than 50
filt_data_50 <- rowSums(DGE_Sub_Counts$counts) > 50
table(filt_data_50)
#Pretty similar, but keeps a few more
#Maybe if time research when I would look at the average log CPM or the read count to determine the filter cutoff

#Let's continue with the default setting then.

DGE_Filt <- DGE_Sub_Counts[filt_data_def, , keep.lib.sizes=FALSE]

#This is interesting because it is only keeping about a third of the genes. Is this normal? The tutorial kept ~60%, so I got around half of what they did.


## Library size normalization

DGE_Filt_TMM <- calcNormFactors(DGE_Filt)
DGE_Filt_TMM$samples
#The really small number at CON.P1 means that a few genes are being very heavily expressed relative to the others in the set
#The tutorial suggests that although TMM normalization is good most of the time, if zero counts are especially large, try TMMwsp normalization. How do I know if my zero counts are large?
DGE_Filt_TMMwsp <- calcNormFactors(DGE_Filt, method = "TMMwsp")
DGE_Filt_TMMwsp$samples
#It's the same.


## Exploring differences between libraries

#MDS plot provides a "type of unsupervised clustering of the samples"
groups
pch_groups <- c(0,15,1,16,2,17)
colors_groups <-c("darkgreen", "darkgreen", "red", "red", "blue", "blue")
plotMDS(DGE_Filt_TMM, col=colors_groups[groups], pch=pch_groups[groups])
legend("bottom", legend = levels(groups), pch=pch_groups, col=colors_groups, ncol = 2)

#Biological replicates mostly appear closer to each other than to the other groups, except for CON.P. The GF S&P seem pretty mixed up together as well. It also seems like the naive mice are to the right side, and the social interaction mice are to the left side, mostly. There is a greater difference between CON.P and CON.S than exGF.P and exGF.S, with GF all being quite central. 

#mean-difference plot visualizes the library-size adjusted log-fold change between two libraries against the average log-expression across those libraries" (difference to mean). It is comparing the sample of your choice to "an artificial reference library constructed from the average of all the other samples"
plotMD(DGE_Filt_TMM, column=1)
abline(h=0, col="red", lty=2, lwd=2)
#This one had the normalization factor of 0.556, which was low, and here we can see that a number of the genes are highly upregulated. There is a definite positive skew.

plotMD(DGE_Filt_TMM, column=2)
abline(h=0, col="red", lty=2, lwd=2)

plotMD(DGE_Filt_TMM, column=5)
abline(h=0, col="red", lty=2, lwd=2)
#This is one of the larger numbers, 1.100, and it seems reasonably balanced, with less upregulated genes

plotMD(DGE_Filt_TMM, column=7)
abline(h=0, col="red", lty=2, lwd=2)
plotMD(DGE_Filt_TMM, column=11)
abline(h=0, col="red", lty=2, lwd=2)


## Dispersion estimation

DGE_Disp <- estimateDisp(DGE_Filt_TMM, design_m, robust=TRUE)
plotBCV(DGE_Disp)
#Tutorial suggests that the comon line (asymptotic value) is between 0.05 and 0.2 for genetically identical mice, and above 0.3 for humans. The value in this chart of ~0.1 seems reasonable, given that the mice are not mentioned to be genetically identical, but are offspring of Swiss Webster breeding pairs. Quick search does not seem to indicate if that means that they are genetically identical or not. Maybe look into more if it seems relevant.

fit_dgeglm <- glmQLFit(DGE_Disp, design_m, robust=TRUE)
head(fit_dgeglm$coefficients)
#output is a DGEGLM object with estimated GLM coefficients for each gene.
#maybe research further into GLM if using this metric in assignment

#DGEGLM object can also plot a QL dispersion 
plotQLDisp(fit_dgeglm)
summary(fit_dgeglm$df.prior)


#### Differential Expression Analysis ----

#All that previous stuff was preliminary processing and analysis to ensure quality of data. Now the experimental evaluations can begin
#Demonstration of analysis in tutorial was between basal cells of lactating mice and basal cells of pregnant mice.
#I will choose a comparison of the social interaction mice that were conventionally raised, and those that were germ free. In the results spreadsheet downloaded, it is sheet 5.

S.CONvsGF <- makeContrasts(CON.S-GF.S, levels = design_m)
class(S.CONvsGF)
res_S.CONvsGF <- glmQLFTest(fit_dgeglm, contrast = S.CONvsGF)
class(res_S.CONvsGF)
topTags(res_S.CONvsGF)
#The first one on the list here is #232 in the results spreadsheet. The second is at 110. Wait though, I don't know what order the spreadsheet is in. But still, comparing the logFC values, the results are quite different. I wonder if this has to do with the subsampling of only 2 biological replicates, or something else about the methods.

#Next the tutorial estimated the false discovery rate. From the documentation, "Identify which genes are significantly differentially expressed from an edgeR fit object containing p-values and test statistics."
de_S.CONvsGF <- decideTestsDGE(res_S.CONvsGF)
class(de_S.CONvsGF)
summary(de_S.CONvsGF)

#use the plotMD function on the DGELRT object to visualize the magnitude of DE
plotMD(res_S.CONvsGF, status = de_S.CONvsGF)

#test whether DE is above a certain threshold.
#tutorial uses log2(1.5) as the cutoff, but this could be a good point to try different values.

tr_S.CONvsGF <-glmTreat(fit_dgeglm, contrast = S.CONvsGF, lfc = log2(1.5))
topTags(tr_S.CONvsGF)

#Now the tutorial does the FDR check again from the new object with the log cutoff applied
is.de.log <- decideTestsDGE(tr_S.CONvsGF)
summary(is.de.log)
#only 11 down and 15 up for DE genes now in my set (compared to 687 down and 368 up in theirs). Does this mean that there are less differences in my samples?

#visualize, as before:
plotMD(tr_S.CONvsGF, status = is.de.log)

#A different way to visualize and express DE genes is through heatmaps. To do so, first the data must be converted to log2-CPM. The tutorial uses the most recent DGEList object for this purpose:
logCPM <- cpm(DGE_Disp, prior.count = 2, log = TRUE)
rownames(logCPM) <- DGE_Disp$genes$Symbol
colnames(logCPM) <- paste(DGE_Disp$samples$group, 1:2, sep = "-")

#The tutorial chooses the top 30 DE genes from the TREAT test above, then makes the heatmap of those genes:
ordered_genes <- order(tr_S.CONvsGF$table$PValue)
logCPM <- logCPM[ordered_genes[1:30],]

coolmap(logCPM, margins = c(7,7), lhei = c(1,6), lwid = c(1,3))
#I'm getting some graphics errors that I'm not sure what they mean, and playing around with the numbers in the arguments provided does not seem to help. But I did get an image that looks like theirs, except for the Z-score legend, which is missing in mine. My biological replicates did cluster together, with the most interesting differences being between CON.S and GF.S, which makes sense, since I chose the top genes in that comparison.

## Analysis of deviance

#Now try comparing 3 at a time, I'll do the three naive mice groups
con_3P <- makeContrasts(
          P.CONvsGF = CON.P - GF.P, 
          P.CONvsEXGF = CON.P - exGF.P,
          P.EXGFvsGF = exGF.P - GF.P, levels = design_m)

res_3P <- glmQLFTest(fit_dgeglm, contrast = con_3P)
topTags(res_3P)


##Complicated contrasts
#The tutorial looks at the interaction effect between mouse status (pregnant, lactating) and cell type, which will describe if the expression change between pregnant and lactating mice is the same for both cell types.
#I will look at if the change in expression between naive and social interaction mice is the same for exGF and CON mice.

con_P.S_CON.exGF <- makeContrasts(
                  (CON.P-CON.S)-(exGF.P-exGF.S),
                  levels = design_m)

res_P.S_CON.exGF <- glmQLFTest(fit_dgeglm, contrast = con_P.S_CON.exGF)
topTags(res_P.S_CON.exGF)

#Tutorial does not describe the interpretation too much, like does this mean within the CON mice from P and S, the change in ENSMUSG00000030103 is significantly different than its change in exGF mice in the P and S conditions?


#### Pathway analysis & Read Alignment----

#The tutorial next goes through pathway analysis, but I will skip it because it is beyond the scope of this assignment

#The tutorial also goes over read alignment and quantification from raw sequence files, which seems useful, but it is probably also beyond the scope of this assignment. I might return to this if there is time.



#### Other Stuff I might use ----

#Was trying to subset 2 biological replicates for each group, but having trouble. Return to this if time
sampled_groups <- vector()
for (i in 1:length(levels(groups))){
  cnt <- 0
  while (cnt<2){
    for (j in 1:length(DGE_Counts$samples$group)){
      if (levels(groups)[i] == DGE_Counts$samples$group[j]){
        sampled_groups <- c(sampled_groups, j)
        cnt <- cnt + 1
      }
    }
  }
}

sampled_groups

#This is the link to download their results, so I can compare what I get to their final product
#Results_URL <- paste("http://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE114702", "format=file", "file=GSE114702_Processed_data_file_2_-_DESeq2_results.xlsx", sep="&")
#download.file(Results_URL, "GSE114702_DESeq2_results.xlsx")

#install.packages("readxl")
library(readxl)

sheets <- excel_sheets("GSE114702_DESeq2_results.xlsx")
#Can make a function, or lapply to read in desired sheets into a list of the dataframes/tibbles
#https://stackoverflow.com/questions/12945687/read-all-worksheets-in-an-excel-workbook-into-an-r-list-with-data-frames
Results_data_1 <- read_excel("GSE114702_DESeq2_results.xlsx", sheet = sheets[1])
Results_data_2 <- read_excel("GSE114702_DESeq2_results.xlsx", sheet = sheets[2])
#also, maybe make more efficient by using filename <- "GSE114702_DESeq2_results.xlsx" or similar
