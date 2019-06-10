library(limma)


tcga.gdata <- read.csv('tcga.gdata.csv', row.names=1)
tcga.gpheno <- read.csv('tcga.gpheno.csv', row.names=1)
design.factor = factor(tcga.gpheno$progression, levels=c('LUSC','Control'))
design.matrix = model.matrix(~0+design.factor)
colnames(design.matrix) = levels(design.factor)

fit <- lmFit(tcga.gdata , design.matrix)#ÏßÐÔÄâºÏ
cont.matrix <- makeContrasts(LUSC-Control, levels=design.matrix)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)
DEG <- topTable(fit2,adjust="fdr",sort.by="B",number=100)
exprSet_DEG = tcga.gdata[rownames(DEG),]
dim(exprSet_DEG)
write.csv(exprSet_DEG,"exprSet_TCGA100.csv",quote=F)
