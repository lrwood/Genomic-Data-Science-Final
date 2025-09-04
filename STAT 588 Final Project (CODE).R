## ----setup, include=FALSE---------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)


## ---------------------------------------------------------------------------
#Install needed packages
library(BiocManager)
library(GEOquery)
library(affy)
library(oligo)
library(umap)
library(limma)
library(siggenes)
library(pd.hg.u133a)
library(data.table)
library(genefilter)
library(edgeR)
library(Glimma)
library(gplots)
library(org.Mm.eg.db)
library(RColorBrewer)


## ---------------------------------------------------------------------------
#Import data
rawdata<-ReadAffy()


## ---------------------------------------------------------------------------
#Import information about data
gset <- getGEO("GSE994", GSEMatrix =TRUE, getGPL=FALSE)
if (length(gset) > 1) idx <- grep("GPL96", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]


## ---------------------------------------------------------------------------
#View distribution of raw data
par(mfrow=c(1,2))
boxplot(rawdata)
hist(rawdata)


## ---------------------------------------------------------------------------
#Normalize the data
normdata<-mas5(rawdata)


## ---------------------------------------------------------------------------
#PLot normalized data
par(mfrow=c(1,2))
boxplot(normdata)
hist(normdata)


## ---------------------------------------------------------------------------
#normdata.2<-affy::rma(rawdata)


## ---------------------------------------------------------------------------
#boxplot(normdata.2)
#hist(normdata.2)


## ---------------------------------------------------------------------------
dat<-gset@phenoData@data


## ---------------------------------------------------------------------------
#Extract expression data
expression<-data.frame(exprs(normdata))
newnames<-gsub(".CEL.gz", "", colnames(expression))
setnames(expression, new=newnames, old=colnames(expression))


## ---------------------------------------------------------------------------
#Ordering the data matrix in order to begin eBayes analysis on Smoker vs Never vs Past smoker
dat<-dat[order(as.numeric(row.names(dat))),]
dat[1:34,1]<-"Current Smoker"
dat[35:58,1]<-"Never Smoker"
dat[59:75,1]<-"Former Smoker"
dat[,1]<-as.factor(dat$title)
sample_info<-data.frame(
  sample=dat$geo_accession,
  condtion=dat$title
)
condition<-sample_info$condtion
design <- model.matrix(~ condition, data = sample_info)


## ---------------------------------------------------------------------------
#Get those differentially expressed genes!
fit<-lmFit(normdata, design)
fit<-eBayes(fit)
results <- topTable(fit, coef = condition, adjust.method = "fdr", number = 75)
results


## ---------------------------------------------------------------------------
#Store the top genes!
top_genes <- results[results$adj.P.Val < 0.05, ]
top_gene_matrix <- expression[rownames(top_genes), ]


## ---------------------------------------------------------------------------
# Create the heatmap with hierarchical clustering
dist_matrix <- dist(top_gene_matrix)
hclust_result <- hclust(dist_matrix, method = "complete")

ordered_data <- top_gene_matrix[hclust_result$order, ]
ordered_data<-as.matrix(ordered_data)

#jpeg("FinalHeatMap.jpg")
heatmap.2(ordered_data,
          scale = "row",    # Scale rows (genes)
          col = colorRampPalette(c("green", "black", "red"))(100),  # Color palette
          key = TRUE,       # Show color key
          keysize = 1.5,     # Size of the color key
          symkey = FALSE,    # Include negative values in the color key
          density.info = "none",  # Remove density plot
          trace = "none",    # Remove row and column dendrograms
          margins = c(5, 10)) # Adjust margins

#dev.off()


## ---------------------------------------------------------------------------
#T-test between current and former smokers
new_matrix<-top_gene_matrix[,1:34]
current_averages<-rowMeans(new_matrix)
never_matrix<-top_gene_matrix[,35:58]
never_averages<-rowMeans(never_matrix)
current.v.never<-data.frame(current_averages,never_averages)
test<-t.test(current_averages,never_averages,alternative="two.sided", mu = 0, paired = FALSE, var.equal = FALSE, conf.level = 0.95)
test$p.value


## ---------------------------------------------------------------------------
#Finding pathway analysis 

library(hgu133plus2.db)
library(hgu133a.db)
library(clusterProfiler)

affyUniverse = featureNames(normdata)
uniId = hgu133aENTREZID[affyUniverse]
entrezUniverse = unique(as.character(uniId))
entrezUniverse<-c(entrezUniverse)

enrich_result <- enrichGO(
  gene = entrezUniverse,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "BP",  # Biological Process. Other options: "MF" (Molecular Function), "CC" (Cellular Component)
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.2
)
print(enrich_result)
dotplot(enrich_result, showCategory = 15)
barplot(enrich_result, showCategory = 15)
enrich_mat <- as.matrix(enrich_result[, c("ID", "Description", "GeneRatio", "qvalue")])



## ---------------------------------------------------------------------------
#Begin analysis about the effect of sex
sample_info2<-data.frame(
  sample=dat$geo_accession,
  sex=dat$characteristics_ch1.1
)
sex<-sample_info2$sex
design2 <- model.matrix(~ sex, data = sample_info2)
sex<-na.omit(sex)


## ---------------------------------------------------------------------------
fit_sex<-lmFit(normdata, design2)
fit_sex<-eBayes(fit_sex)
results_sex <- topTable(fit_sex, coef = condition, adjust.method = "fdr", number = Inf)
results_sex$GeneName<-rownames(results_sex)
TOPGENE<-"203445_s_at"
row_of_interest <- results_sex[results_sex$GeneName == TOPGENE, ]


## ---------------------------------------------------------------------------
#Begin analysis about the effect of race
sample_info3<-data.frame(
  sample=dat$geo_accession,
  race=dat$characteristics_ch1.2
)
race<-sample_info3$race
design3 <- model.matrix(~ race, data = sample_info3)
race<-na.omit(race)


## ---------------------------------------------------------------------------
fit_race<-lmFit(normdata, design3)
fit_race<-eBayes(fit_race)
results_race <- topTable(fit_race, coef = condition, adjust.method = "fdr", number = Inf)
results_race$GeneName<-rownames(results_race)
row_of_interest3 <- results_race[results_race$GeneName == TOPGENE, ]



## ---------------------------------------------------------------------------
#Begin analysis about the effect of age
ages<-data.frame(dat$geo_accession)
ages$age<-dat$characteristics_ch1.3
age_test<-top_gene_matrix[1:10,]
age_test$"32"<-age_test$GSM15693
age_test$"23"<-age_test$GSM15694
age_test$"28"<-age_test$GSM15722
age_test$"31"<-age_test$GSM15723
age_test$"28"<-age_test$GSM15724
age_test$"27"<-age_test$GSM15726
age_test$"68"<-age_test$GSM15727
age_test$"25"<-age_test$GSM15730
age_test$"45"<-age_test$GSM15748
test_age<-age_test[,76:83]



## ---------------------------------------------------------------------------
library(tidyr)
your_data_long <- gather(test_age, key = "Variable", value = "Value")
plot(your_data_long$Variable, your_data_long$Value, xlab = "Age", ylab = "Gene Expression")


