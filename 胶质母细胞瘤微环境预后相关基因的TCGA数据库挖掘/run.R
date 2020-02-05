#安装R包
#source('https://bioconductor.org/biocLite.R')
#biocLite(pkg_name)
#library(BiocManager)
#BiocManager::install()

#一、读取表达矩阵
suppressPackageStartupMessages(library(CLL))
data(sCLLex)
exprSet=exprs(sCLLex)  
samples=sampleNames(sCLLex)
pdata=pData(sCLLex)
group_list=as.character(pdata[,2])
dim(exprSet)
exprSet[1:5,1:5]
group_list
#二、制作分组矩阵
suppressMessages(library(limma))
design <- model.matrix(~0+factor(group_list))
colnames(design)=levels(factor(group_list))
rownames(design)=colnames(exprSet)
design
#三、制作差异比较矩阵
contrast.matrix<-makeContrasts(paste0(unique(group_list),
　　　　　　　　collapse = "-"),levels = design)
contrast.matrix

#使用limma包来进行差异分析
fit <- lmFit(exprSet,design)

fit2 <- contrasts.fit(fit, contrast.matrix) 
fit2 <- eBayes(fit2)

tempOutput = topTable(fit2, coef=1, n=Inf)
nrDEG = na.omit(tempOutput) 
head(nrDEG)
write.csv(nrDEG,"limma_notrend.results.csv",quote = F)
