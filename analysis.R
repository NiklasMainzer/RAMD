library(BiocGenerics)
library(Biobase)
library(limma)

##################
### Referenzen ###
##################

# Folgende Online Referenzen, wurden bei der Programmierung zur Hilfe genommen:

# https://ucdavis-bioinformatics-training.github.io/2018-June-RNA-Seq-Workshop/thursday/DE.html
# https://jdblischak.github.io/dc-bioc-limma/vdx.html
# https://learn.gencore.bio.nyu.edu/rna-seq-analysis/gene-set-enrichment-analysis/
# https://academic.oup.com/nar/article/43/7/e47/2414268
# Official Bioconductor documentation

# Alle Links wurden zuletzt am 1.04.2022 um 12:00 geöffnet

####################
### Prepare Data ###
####################

data <- read.table(file = "./Project8_crohns.txt", header = TRUE)

data_control <- read.table(file = "./Project8_control.txt", header = TRUE)

colnames(data) <- c('hgnc', 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
colnames(data_control) <- c('hgnc', 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1)

total_data <- merge(data, data_control)

# 0 for MC patients and 1 for control group
colnames(total_data)= c('hgnc', 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1)

###########################
### Explorative Analyse ###
###########################

# define important genes by calculating mean/median btw sick/not sick and chose those where the difference is significant
# clustering might identify which genes belong to one group (pathway analysis and stuff?)

total_data_matrix <- data.matrix(total_data)
total_data_matrix <- total_data_matrix[,-c(1)]
plotDensities(total_data_matrix)

# MA Plot
# Points further from M=0 line indicate genes with significant expression - down/upregulated
limma::plotMA(total_data_matrix, main="MA Plot")

plotMDS(total_data_matrix)

# Hierarchical cluster (all genes would be too much), choose the ones that have the highest difference between the median of ibd patients & control group
ibd_median <- apply(data[,-1], 1, median)
control_median <- apply(data_control[,-1], 1, median)
median_both <- data.frame(ibd = ibd_median, control=control_median)
median_both$hgnc <- data$hgnc
median_both$diff <- median_both$ibd - median_both$control
high_genes <- median_both
# high_genes[high_genes$diff >= 1,]

# Filtered 124 genes to continue clustering
filtered_genes <- high_genes[abs(high_genes$diff) >= 1,]

data_filtered = data.frame()
data_filtered <- total_data[total_data$hgnc %in% filtered_genes$hgnc,]

data_filtered_dist_genes <- dist(data_filtered)
filtered_matrix <- data.matrix(data_filtered)
filtered_matrix <- filtered_matrix[,-c(1)]
data_filtered_dist_samples <- dist(t(filtered_matrix))

gene_hclust <- hclust(data_filtered_dist_genes, method = "complete")
plot(gene_hclust)

# Cut the tree at cluster size k=2
cutree(gene_hclust, k = 2)

# Hierarchical Clustering with 500 random genes
random_total_data <- total_data[sample(nrow(total_data), 500), ]
random_total_data_matrix <- data.matrix(random_total_data)
random_total_data_matrix <- random_total_data_matrix[,-c(1)]
random_total_dist_samples = dist(t(random_total_data_matrix))
random_clust = hclust(random_total_dist_samples, method = "complete")
plot(random_clust)

#############
###  DEA  ###
#############

groups = gsub("_.*", "", colnames(total_data))
colnames(total_data)= c('hgnc', 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1)

groups = colnames(total_data)
groups <- factor(groups, levels = c(0,1))

# Design matrix has all samples as rows and the columns indicate 1 for belonging to a group and 0 to not belonging, while group0 has MC and group1 is the control group
design <- model.matrix( ~ 0 + groups )

# Fit a linear model to every row/gene
fit <- lmFit(total_data, design)
head(coef(fit))
dim(coef(fit))
coef_matrix <- coef(fit)

# Test DE between group0 (MC patients) and group1 (control)
contmatrix = makeContrasts(groups0 - groups1, levels=colnames(design))
fit2 <- contrasts.fit(fit, contmatrix)
fit2

# Empirischer Bayes
# Berechnet t-Statistik, F-Statistik, log-odds und differential expression
fit2 <- eBayes(fit2)

# Multiple Testing: If t-statistic (and p value) is signifincantly neg/pos/not significant
# Compares MC with control 
summary(decideTests(fit2, adjust.method = "BH"))
#       groups0 - groups1
# Down                3554  -> comparing ibd to control these genes are neg. significant 
# NotSig             16655  -> comparing ibd to control these genes are not significant 
# Up                  2267  -> comparing ibd to control these genes are pos. significant 

# Gibt die top 20 differentially expressed Genes an 
top.table <- topTable(fit2, sort.by = "P", n = Inf)
head(top.table, 20) # see topGenes_DEA.txt file for output

top.table_fc <- topTable(fit2, sort.by = "logFC", n = Inf)
head(top.table_fc, 20)

#############
### Plots ###
#############

statistics <- topTable(fit2, number = nrow(fit2), sort.by = "none")
head(top.table, 20)

# Histogram by p value of DE
hist(statistics[, "P.Value"])

# Statistical significance (p value) vs magnitude of change (fold change)
# most significant genes are at top (high p value) and left/right the most up/downregulated genes
volcanoplot(fit2, highlight = 5, names = fit2$genes[, "hgnc"])

#############
### GSEA  ###
#############
library('org.Hs.eg.db') # Genome wide annotation for Hs (homo sapiens)
library(ggplot2)
library(enrichplot)
library(clusterProfiler)
library(EnrichmentBrowser)

columns(org.Hs.eg.db)

test_func = function(x) {return(x[1])}
tmp <- apply(total_data,1,test_func)

# Map hgnc Nomenklatur zur entrez ID 
entrez <- mapIds(org.Hs.eg.db, tmp, 'ENTREZID', 'SYMBOL')

original_gene_list <- top.table$logFC
names(original_gene_list) <- total_data$hgnc
gene_list<-na.omit(original_gene_list)
gene_list = sort(gene_list, decreasing = TRUE)

# Gene Set Enrichment Analysis of Gene Ontology
# CC: Cellular Component, MF: Molecular function, BP: Biological Process
gse <- gseGO(geneList=gene_list, 
             ont ="ALL", 
             keyType = "SYMBOL", 
             minGSSize = 3, 
             maxGSSize = 800, 
             pvalueCutoff = 0.05, 
             verbose = TRUE, 
             OrgDb = 'org.Hs.eg.db', 
             pAdjustMethod = "none")

#############
### Plots ###
#############

require(DOSE)

dotplot(gse, showCategory=10, split=".sign") + facet_grid(.~.sign)

gse2 <- pairwise_termsim(gse, semData="org.Hs.eg.db") 
emapplot(gse2)

cnetplot(gse, categorySize="pvalue", foldChange=gene_list, showCategory = 3)
cnetplot(gse, categorySize="pvalue", foldChange=gene_list) # categorySize="geneNum"

ridgeplot(gse) + labs(x = "enrichment distribution")

# GSEAplot/Barcodeplot, index ist das jeweilige Gen in der GSEA
gseaplot(gse, by = "all", title = gse$Description[2], geneSetID = 2)
# Find the gene set to be plotted by the Name
which(gse_desc %in% "cellular response to cholesterol" )

pmcplot(terms, 2010:2021, proportion=FALSE)

###################################
### Overrepresentation Analysis ###
###################################

# Problematic: just counting ignores p-values and fold-changes
# GO over representation analysis
fit3 <- fit2
enrich_go <- goana(fit3, geneid = entrez, species = "Hs") # Hs for Homo Sapiens
topGO(enrich_go, ontology = "BP")
topGO(enrich_go)

# Convert SYMBOL to ENTREZ
de <- names(gene_list)[abs(gene_list) > 2]
ego <- enrichGO(de, OrgDb = "org.Hs.eg.db", keyType="SYMBOL", ont="BP", readable=FALSE)
goplot(ego)

ggo <- groupGO(gene     = de,
                OrgDb    = org.Hs.eg.db,
                ont      = "BP",
                level    = 3,
                keyType="SYMBOL",
               readable = FALSE)

# GSEA KEGG gseKegg
ids<-bitr(names(original_gene_list), fromType = "SYMBOL", toType = "ENTREZID", OrgDb="org.Hs.eg.db")
dedup_ids = ids[!duplicated(ids[c("SYMBOL")]),]
df2 = top.table[total_data$hgnc %in% dedup_ids$SYMBOL,]
df2$Y = dedup_ids$ENTREZID
kegg_gene_list <- df2$logFC
names(kegg_gene_list) <- df2$Y
kegg_gene_list<-na.omit(kegg_gene_list)
kegg_gene_list = sort(kegg_gene_list, decreasing = TRUE)
kegg_organism = "hsa" # homo sapiens

kk2 <- gseKEGG(geneList     = kegg_gene_list,
               organism     = kegg_organism,
               minGSSize    = 3,
               maxGSSize    = 800,
               pvalueCutoff = 0.05,
               pAdjustMethod = "none",
               keyType       = "ncbi-geneid")

kk2@result[["ID"]]

#Plots
dotplot(kk2, showCategory = 5, title = "Enriched Pathways" , split=".sign") + facet_grid(.~.sign)
kk_new <- pairwise_termsim(kk2, semData="org.Hs.eg.db")
emapplot(kk_new)
cnetplot(kk2, categorySize="pvalue", foldChange=gene_list)
ridgeplot(kk2) + labs(x = "enrichment distribution")
out_kk2 <- as.matrix(kk2@result)

#######################
### Pathway Analyse ###
#######################

# Reformat the hgnc-entrez dictionary to a dataframe only containing the corrensponding entrez
library(data.table)
entrez2 <- as.data.frame(t(entrez))
entrez2 <- transpose(entrez2)
colnames(entrez2) <- 'entrez'

# Replace the hgnc genes with the entrez
fit3$genes <- entrez2

# KEGG Pathway analysis (kegga)
entrez_kegga = fit3$genes[, "entrez"]
enrich_kegg <- kegga(fit3, geneid = entrez_kegga, species = "Hs")
topKEGG(enrich_kegg) # kegga_pathway.txt

# Plot Pathways in graphics
library(pathview)
# take pathways form kegga_pathway.txt
dme <- pathview(gene.data=kegg_gene_list, pathway.id="hsa02010", species = kegg_organism)
