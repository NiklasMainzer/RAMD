library(Biobase)
library(limma)

data <- read.table(file = "/Users/deborahhoeltje/Desktop/TUBS4/medDaten/Prüfungsprojekte/Project8_crohns.txt", header = TRUE)

data_control <- read.table(file = "/Users/deborahhoeltje/Desktop/TUBS4/medDaten/Prüfungsprojekte/Project8_control.txt", header = TRUE)

#############
### DEA  ###
#############

groups = gsub("_.*", "", colnames( total_data))
colnames(total_data)= c('hgnc', 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1)

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

# Test DE between group0 (MC patients) and group1 (control)
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
library('org.Hs.eg.db')
library(ggplot2)
library(enrichplot)
library(clusterProfiler)

columns(org.Hs.eg.db)

test_func = function(x) {return(x[1])}
tmp <- apply(total_data,1,test_func)

# Map hgnc Nomenklatur zur entrez ID 
entrez <- mapIds(org.Hs.eg.db, tmp, 'ENTREZID', 'SYMBOL')

original_gene_list <- top.table$logFC
names(original_gene_list) <- total_data$hgnc
gene_list<-na.omit(original_gene_list)
gene_list = sort(gene_list, decreasing = TRUE)

gse <- gseGO(geneList=gene_list, 
             ont ="ALL", 
             keyType = "SYMBOL", 
             nPerm = 10000, 
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

gse2 <- pairwise_termsim(gse, semData="org.Dm.eg.db")
emapplot(gse2)

cnetplot(gse, categorySize="pvalue", foldChange=gene_list, showCategory = 3)
cnetplot(gse, categorySize="pvalue", foldChange=gene_list)

ridgeplot(gse) + labs(x = "enrichment distribution")

# http://www.gsea-msigdb.org/gsea/doc/GSEAUserGuideFrame.html?_Interpreting_GSEA_Results
# GSEAplot, index ist das jeweilige Gen in der GESA
gseaplot(gse, by = "all", title = gse$Description[2], geneSetID = 2)

pmcplot(terms, 2010:2022, proportion=FALSE)

