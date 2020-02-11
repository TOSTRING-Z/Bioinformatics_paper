#文件说明

[https://xenabrowser.net/datapages/][UCSC]

GBM_clinicalMatrix `分子亚型数据` 

GBM_mc3_gene_level `基因突变数据` 

GBM_survival `预后数据` 

HT_HG-U133A `基因表达谱数据` 

[https://bioinformatics.mdanderson.org/estimate/disease.html?glioblastoma%20multiforme_Affymetrix][ESTIMATE] 

glioblastoma_multiforme_Affymetrix `ESTIMATE得分数据`

[http://www.cgga.org.cn/download.jsp][CGGA]

[https://string-db.org/][STRING] `基因互作信息`

##软件使用

GraphPad `画图，统计分析`
SangerBox `数据分析`

[UCSC]: https://xenabrowser.net/datapages/

[ESTIMATE]: https://bioinformatics.mdanderson.org/estimate/disease.html?glioblastoma%20multiforme_Affymetrix

[CGGA]: http://www.cgga.org.cn/download.jsp

[STRING]: https://string-db.org/

##分析流程

1.数据准备

2.使用ESTIMATE计算 基质细胞 和 免疫细胞 得分

3.方差分析 亚型 是否对 得分 有统计学意义 --> 有

4.T检验 IDH1 是否对 得分 有统计学意义 --> 无

5.将 得分数据 排序分以为 高低2组，绘制KM生存曲线，使用Log-rank算法计算Pvalue，看 得分高低 是否对 病人生存时间 具有统计学意义 --> 无

6.将 得分数据 排序分以为 高低2组，进行 基因差异表达分析 筛选具有 显著差异的基因 --> 差异基因占比较集大部分

7.将 上调的 差异表达基因 单个单个 排序分以为 高低2组，绘制KM生存曲线，使用Log-rank算法计算Pvalue，看 基因表达高低 是否对 病人生存时间 具有统计学意义，筛选具有统计学意义基因

8.使用 差异表达基因 绘制 蛋白-蛋白 互作网络，分为几个组，选择其中 度 最多的节点

9.使用CGGA已经实验验证的数据对 当前结果进行验证 --> 验证包含了结果基因 并且 发现了实验未证实的基因
