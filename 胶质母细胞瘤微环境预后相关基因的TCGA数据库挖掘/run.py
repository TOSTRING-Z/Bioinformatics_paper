import pandas as pd
# pandas和R转换 https://rpy2.github.io/doc/latest/html/pandas.html
# pandas与R对比 https://pandas.pydata.org/pandas-docs/stable/getting_started/comparison/comparison_with_r.html
# 让R与Python共舞 https://www.cnblogs.com/lantingg/p/9600280.html
from rpy2.robjects import r
from rpy2 import robjects
import rpy2.robjects as ro
from rpy2.robjects.packages import importr
from rpy2.robjects import pandas2ri
from rpy2.robjects.conversion import localconverter

if __name__ == "__main__":
    # 设置数据路径
    path = "/install/git/Bioinformatics_paper/胶质母细胞瘤微环境预后相关基因的TCGA数据库挖掘/"
    r.setwd(path)
    # 读取处理好的数据
    sample = pd.read_csv(f"{path}sample.txt", sep="\t")
    sample_Group = sample["Stromal_Group"]
    # 读取处理好的基因表达数据
    HT_HG_U133A_sample = pd.read_csv(f"{path}HT_HG_U133A_sample.txt", sep="\t")

    ################# 方差分析(ANOVA) GeneExp_Subtype #################
    # https://www.bioinfo-scrounger.com/archives/588/
    with localconverter(ro.default_converter + pandas2ri.converter):
        ANOVA_data_R = ro.conversion.py2rpy(sample[["Stromal_score","GeneExp_Subtype"]])
        print(r.summary(r.aov(r("Stromal_score~GeneExp_Subtype"), data = ANOVA_data_R)))

    ################# t检验 IDH1 #################
    r('''suppressMessages(library(MASS))''')
    with localconverter(ro.default_converter + pandas2ri.converter):
        Ttest_data_R = ro.conversion.py2rpy(sample[["Stromal_score","IDH1"]].query("IDH1==1 or IDH1==0"))
        print(r["t.test"](r("Stromal_score~IDH1"), data=Ttest_data_R))

    ################# 生存分析 #################
    # https://www.jianshu.com/p/4ad9ba730719
    # r('''suppressMessages(library(survival))''')
    importr("survival")
    survminer = importr("ggfortify")
    with localconverter(ro.default_converter + pandas2ri.converter):
        # 构建生存对象
        robjects.globalenv["cox_data"] = ro.conversion.py2rpy(sample[["Stromal_score","Stromal_Group","OS.time","OS"]].dropna())
    robjects.globalenv["surv"] = r("Surv(time=cox_data$OS.time, event=cox_data$OS)")
    # Kaplan-Meier生存曲线
    KM_Stromal_fit = r("survfit(surv~cox_data$Stromal_Group)")
    r.ggsave(r.autoplot(KM_Stromal_fit), file="KM_Stromal_fit.pdf")
    # Long-rank检验(对数秩和检验)
    print(r("survdiff(surv~cox_data$Stromal_Group,rho = 0)"))
    # 单因素Cox
    print(r("coxph(surv~cox_data$Stromal_score)"))

    ################# 基因差异表达分析 #################
    ## 读取表达矩阵
    group_list = "-".join(list(sample["Stromal_Group"].astype("str").unique()))
    exprSet_R = r(f'''read.table("{path}HT_HG_U133A_sample.txt",header=TRUE,sep="\t",row.names="sample")''')
    ## 制作分组矩阵
    design = sample[["Stromal_Group_H", "Stromal_Group_L"]].set_index(sample["ID"])
    with localconverter(ro.default_converter + pandas2ri.converter):
        design_R = ro.conversion.py2rpy(design.rename(columns={"Stromal_Group_H": "H", "Stromal_Group_L": "L"}))
    ## 制作差异比较矩阵
    r('''suppressMessages(library(limma))''')
    contrast_matrix = r['makeContrasts'](group_list, levels=design_R)
    ## 使用limma包来进行差异分析
    fit = r["lmFit"](exprSet_R, design_R)
    fit2 = r["contrasts.fit"](fit, contrast_matrix)
    fit2 = r["eBayes"](fit2)
    tempOutput = r["topTable"](fit2, coef=1, number=float('inf'), lfc=0.5849625007211562)
    with localconverter(ro.default_converter + pandas2ri.converter):
        nrDEG = ro.conversion.rpy2py(tempOutput)
        nrDEG.dropna(axis=0, how='any', inplace=True)
    ## 筛选差异显著基因abs(log2FoldChange)>1.5,P.Value<0.05
    Differentially_significant_genes = nrDEG[nrDEG["P.Value"] < 0.05].query(
        '''logFC>0.5849625007211562 or logFC<-0.5849625007211562''')
    Differentially_significant_genes.to_csv(f"{path}Differentially_significant_genes.txt", sep="\t")

    ################# 基因功能富集分析 #################
    # 富集分析 https://www.jianshu.com/p/988d90484f77
    # BiocManager::install("org.Hs.eg.db")
    r("suppressMessages(library(clusterProfiler))")
    r("suppressMessages(library(org.Hs.eg.db))")
    symbol_id = robjects.StrVector(list(Differentially_significant_genes.index))
    # symbol_id 转 entre_id
    entre_id = r.unlist(r.bitr(symbol_id, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db").rx2('ENTREZID'))
    # GO富集分析
    #ALL_R = r.enrichGO(entre_id, "org.Hs.eg.db", keyType="ENTREZID", ont='ALL', pvalueCutoff=0.05, pAdjustMethod="BH",
    #                   qvalueCutoff=0.1, readable=True)  # 一步到位
    BP_R = r.enrichGO(entre_id, "org.Hs.eg.db", keyType="ENTREZID", ont="BP", pvalueCutoff=0.05, pAdjustMethod="BH",
                      qvalueCutoff=0.1, readable=True)  # 3种分开进行富集
    MF_R = r.enrichGO(entre_id, "org.Hs.eg.db", keyType="ENTREZID", ont="MF", pvalueCutoff=0.05, pAdjustMethod="BH",
                      qvalueCutoff=0.1, readable=True)
    CC_R = r.enrichGO(entre_id, "org.Hs.eg.db", keyType="ENTREZID", ont="CC", pvalueCutoff=0.05, pAdjustMethod="BH",
                      qvalueCutoff=0.1, readable=True)
    # KEGG富集分析
    KEGG_R = r.enrichKEGG(gene=entre_id, organism='hsa', qvalueCutoff=0.05)
    with localconverter(ro.default_converter + pandas2ri.converter):
        BP = ro.conversion.rpy2py(r["as.data.frame"](r.slot(BP_R, "result")))
        MF = ro.conversion.rpy2py(r["as.data.frame"](r.slot(MF_R, "result")))
        CC = ro.conversion.rpy2py(r["as.data.frame"](r.slot(CC_R, "result")))
        KEGG = ro.conversion.rpy2py(r["as.data.frame"](r.slot(KEGG_R, "result")))
        #GO_ALL = ro.conversion.rpy2py(r["as.data.frame"](r.slot(KEGG_R, "result")))
        BP.query("pvalue<0.05").to_csv(f"{path}BP.txt", sep="\t")
        MF.query("pvalue<0.05").to_csv(f"{path}MF.txt", sep="\t")
        CC.query("pvalue<0.05").to_csv(f"{path}CC.txt", sep="\t")
        KEGG.query("pvalue<0.05").to_csv(f"{path}KEGG.txt", sep="\t")


