library(Biobase)
library(limma)

data <- read.table(file = "/Users/deborahhoeltje/Desktop/TUBS4/medDaten/Prüfungsprojekte/Project8_crohns.txt", header = TRUE)

data_control <- read.table(file = "/Users/deborahhoeltje/Desktop/TUBS4/medDaten/Prüfungsprojekte/Project8_control.txt", header = TRUE)

groups = gsub("_.*", "", colnames( total_data))
colnames(total_data)= c('hgnc', 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1)

groups = colnames(total_data)
groups <- factor(groups, levels = c(0,1))

design <- model.matrix( ~ 0 + groups )

fit <- lmFit(total_data, design)
head(coef(fit))
dim(coef(fit))
coef_matrix <- coef(fit)
bayes_fit <- eBayes(fit)
head(bayes_fit)
head(bayes_fit$t, 3)

results <- decideTests(bayes_fit[, 1])
contmatrix = makeContrasts(groups0 - groups1, levels=colnames(design))
fit2 <- contrasts.fit(fit, contmatrix)
fit2
fit2 <- eBayes(fit2)
summary(decideTests(fit2))

top.table <- topTable(fit2, sort.by = "P", n = Inf)
head(top.table, 20)


#############
### Plots ###
#############
stats <- topTable(fit2, number = nrow(fit2), sort.by = "none")
head(top.table, 20)
hist(runif(10000))
hist(stats[, "P.Value"])

volcanoplot(fit2, highlight = 5, names = fit2$genes[, "hgnc"])

