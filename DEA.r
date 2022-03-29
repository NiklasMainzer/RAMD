library(Biobase)
library(limma)

data <- read.table(file = "/Users/deborahhoeltje/Desktop/TUBS4/medDaten/Prüfungsprojekte/Project8_crohns.txt", header = TRUE)

data_control <- read.table(file = "/Users/deborahhoeltje/Desktop/TUBS4/medDaten/Prüfungsprojekte/Project8_control.txt", header = TRUE)

total_data <- data.frame(data, data_control)

#############
### DEA  ###
#############

groups = gsub("_.*", "", colnames(total_data))

# 0 for MD patients and 1 for control group
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

bayes_fit <- eBayes(fit)
head(bayes_fit)
head(bayes_fit$t, 3)
results <- decideTests(bayes_fit[, 1])

# Test DE between group0 (MD patients) and group1 (control)
contmatrix = makeContrasts(groups0 - groups1, levels=colnames(design))
fit2 <- contrasts.fit(fit, contmatrix)
fit2

# empirischer Bayes wird verwendet, wenn man wenige samples hat. Hier wird die globale Varianz über alle Gene geschätzt in in die Genweite Varianz miteinbezogen (Gene mit hoher Varianz werden nach unten angepasst und Gene mit einer niedrigen Varianz nach oben)
# Berechnet t-Statistik, F-Statistik, log-odds und differential expression
fit2 <- eBayes(fit2)

# Multiple Testing: If t-statistic (and p value) is signifincantly neg/pos/not significant
# Compares ibd with control 
summary(decideTests(fit2))
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
stats <- topTable(fit2, number = nrow(fit2), sort.by = "none")
head(top.table, 20)

# random deviates of the uniform distribution 
hist(runif(10000))

# Histogram by p value of DE
hist(stats[, "P.Value"])

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
# Only gene Set size in [minGSSize, maxGSSize] will be tested
# GSEA use permutation test -> nPerm ---> also tried without perm (see plots annotated with *_noPerm.pdf)
gse <- gseGO(geneList=gene_list, 
             ont ="ALL", 
             keyType = "SYMBOL", 
             #nPerm = 10000, 
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

gse2 <- pairwise_termsim(gse, semData="org.Hs.eg.db") # Falsche DB -> mit neues DB emap_neu.pdf
emapplot(gse2)

cnetplot(gse, categorySize="pvalue", foldChange=gene_list, showCategory = 3)
cnetplot(gse, categorySize="pvalue", foldChange=gene_list)

ridgeplot(gse) + labs(x = "enrichment distribution")

# http://www.gsea-msigdb.org/gsea/doc/GSEAUserGuideFrame.html?_Interpreting_GSEA_Results
# GSEAplot, index ist das jeweilige Gen in der GESA
gseaplot(gse, by = "all", title = gse$Description[2], geneSetID = 2)

pmcplot(terms, 2010:2021, proportion=FALSE)


#######################
### Pathway Analyse ###
#######################

# Reformat the hgnc-entrez dictionary to a dataframe only containing the corrensponding entrez
library(data.table)
entrez2 <- as.data.frame(t(entrez))
entrez2 <- transpose(entrez2)
colnames(entrez2) <- 'entrez'

# Replace the hgnc genes with the entrez
fit3 <- fit2
fit3$genes <- entrez2

# KEGG Pathway analysis (kegga)
entrez_kegga = fit3$genes[, "entrez"]
enrich_kegg <- kegga(fit3, geneid = entrez_kegga, species = "Hs")
topKEGG(enrich_kegg)

# https://jdblischak.github.io/dc-bioc-limma/vdx.html
# GO over representation analysis
enrich_go <- goana(fit3, geneid = entrez, species = "Hs") # Hs for Homo Sapiens
topGO(enrich_go, ontology = "BP")
topGO(enrich_go)

# Convert SYMBOL to ENTREZ
#names(gene_list) <- entrez2$entrez
de <- names(gene_list)[abs(gene_list) > 2]
ego <- enrichGO(de, OrgDb = "org.Hs.eg.db", keyType="SYMBOL", ont="BP", readable=FALSE)
goplot(ego)

ggo <- groupGO(gene     = de,
                OrgDb    = org.Hs.eg.db,
                ont      = "BP",
                level    = 3,
                keyType="SYMBOL",
               readable = FALSE)

# KEGG gseKegg
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

# Plot Pathways in graphics
library(pathview)
# take pathways form kegga_pathway.txt
dme <- pathview(gene.data=kegg_gene_list, pathway.id="hsa02010", species = kegg_organism)

# TODO: 
# - Fishers t test to obtain p-values
# - overrepresentation analysis: https://yulab-smu.top/biomedical-knowledge-mining-book/enrichment-overview.html
# - gsea top table bekommen

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

#Hierarchical cluster (all genes would be too much), choose the ones that have the highest difference between the median of ibd patients & control group
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

cutree(gene_hclust, k = 2)

