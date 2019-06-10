library(limma)
library(ggplot2)
library(reshape2)
library(annotate)


pheno <- read.csv('gpheno_new.csv', row.names=1)
design.factor = factor(pheno$progression, levels=c("Regressive","Progressive"))
design.matrix = model.matrix(~0+design.factor)
colnames(design.matrix) = levels(design.factor)

exprSet <- read.csv('gdata_new.csv', row.names=1)
fit <- lmFit(exprSet, design.matrix)#线性拟合
cont.matrix <- makeContrasts(Regressive-Progressive, levels=design.matrix)#case和control进行比较
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)
DEG <- topTable(fit2,adjust="BH",sort.by="B",number=100)
exprSet_DEG = exprSet[rownames(DEG),]
write.csv(exprSet_DEG,"exprSet_DEG100.csv",quote=F)
