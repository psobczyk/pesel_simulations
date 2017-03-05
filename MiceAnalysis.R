
library(FactoMineR)
library(pesel)
mouse = read.table("http://factominer.free.fr/docs/souris.csv", header = T, sep = ";",
                   row.names = 1)
expressions = mouse[, 24: ncol(mouse)]
dim(expressions)
# 40 120

## Pesel analysis
res <- pesel(expressions, npc.min = 0, npc.max = min(ncol(expressions) -2, nrow(expressions)-2))
res # 5
plot(res)
plot(res, posterior = F)

## Pesel with exponential prior
pc_prior <- dgeom(0:min(ncol(expressions) -2, nrow(expressions)-2), 0.5)
res <- pesel(expressions, npc.min = 0, npc.max = min(ncol(expressions) -2, nrow(expressions)-2), 
             prior = pc_prior)
res # 5
plot(res)
plot(res, posterior = F)


## GCV
res.gcv <- estim_ncp(expressions)
plot(res.gcv$crit)
res.gcv$ncp # 12

## PCA plots for mice data
res.pca <- PCA(cbind.data.frame(mouse[,1:2],expressions), quali.sup= 1:2, graph = F)
plot.PCA(res.pca, habillage = 2, invisible = "quali")
plotellipses(res.pca, keepvar = "Regime", axes= c(3,4))
plot.PCA(res.pca, choix =  "var", axes= c(1,2), select = "contrib 20", cex = 0.7)
plot.PCA(res.pca, choix =  "var", axes= c(3,4), select = "contrib 20", cex = 0.7)




