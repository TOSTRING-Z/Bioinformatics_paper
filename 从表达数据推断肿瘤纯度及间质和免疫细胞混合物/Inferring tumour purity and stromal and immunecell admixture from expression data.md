#Inferring tumour purity and stromal and immune cell admixture from expression data
从表达数据推断肿瘤[ˈtuːmər]纯度及基质(stromal)细胞和免疫[ɪˈmjuːn]细胞混合物[ədˈmɪkstʃər]

Infiltrating stromal and immune cells form **the major fraction** of normal cells in tumour tissue \
浸润性（infiltrating）基质细胞和免疫细胞 **主要 部分**[ˈfrækʃn]来自在肿瘤组织的正常细胞

and not only perturb the tumour signal in molecular studies but also **have an important role** in cancer biology. \
不仅 干扰[pərˈtɜːrb]在 分子生物学[məˈlekjələr]的肿瘤信号 而且在肿瘤生物学中 **有重要作用** 

Here we describe ‘Estimation of STromal and Immune cells in MAlignant Tumours using Expression data’ \
这里我们描述了 `使用表达数据 估计 在恶性[məˈlɪɡnənt]肿瘤中的 基质和肿瘤细胞`

(ESTIMATE)—a method that uses gene expression signatures to infer the fraction of stromal and immune cells in tumour samples. \
(ESTIMATE)—使用 基因表达信号 去推断 在肿瘤样本中的 基质和免疫细胞的部分

ESTIMATE scores **correlate with** DNA copy number-based tumour purity across samples from 11 different tumour types, \
ESTIMATE分数 **和** 来自11个不同肿瘤类型的样本的 基于DNA拷贝数的 肿瘤纯度 **有关**

profiled on Agilent, Affymetrix platforms or based on RNA sequencing and available through The Cancer Genome Atlas. \
图谱[ˈproʊfaɪlz]在Agilent,Affymetrix平台 或者基于RNA测序 并可通过癌症基因组 图谱(Atlas) 获得

The prediction accuracy is further corroborated using 3,809 transcriptional profiles available elsewhere in the public domain. \
预测准确度 使用 在公共领域其它可用的 3809个转录图谱 被 进一步[ˈfɜːrðər]证实[kəˈrɑːbəreɪtɪd]

The ESTIMATE method allows consideration of **tumour-associated** normal cells in genomic and transcriptomic studies. \ 
ESTIMATE方法允许 在基因组(genomic)和转录组(transcriptomic)研究中 考虑 **肿瘤相关的** 正常细胞

An R-library **is available on** https://sourceforge.net/projects/estimateproject/. \
**在**https://sourceforge.net/projects/estimateproject/**上提供了** 一个R库

## ARTICLE
Malignant solid tumour tissues consist of not only tumour cells but also tumour-associated normal epithelial and stromal cells,immune cells and vascular cells.
恶性(Malignant)实体瘤 组织 组成不仅是肿瘤细胞 而且有和肿瘤相关的正常 上皮[ɛpɪˈθɛljəl]细胞 和基质细胞,免疫细胞和血管[ˈvæskjələr]细胞

Stromal cells are thought to have important roles in tumour growth,disease progression 1,2 and drug resistance 3. 
基质细胞被 认为[θɔːt]对在肿瘤生长,疾病发展1和药物抗性3中有重要影响

Infiltrating immune cells act in a context-dependent manner,
浸润性免疫细胞 以一种依赖上下文方式 移动

and whereas antitumor effects of infiltrating T-lymphocytes have been observed in ovarian cancer 4–6 , 
鉴于(whereas) 浸润性T淋巴细胞 在卵巢[oʊˈveriən]癌中 抗肿瘤作用 已经 被发现

associations with tumour growth, invasion and metastasis were described in colorectal cancer 7,8 .
在 直结肠[ˌkoʊloʊˈrɛktəl]癌7,8中 和肿瘤生长, 侵袭[ɪnˈveɪʒn]和 转移[məˈtæstəsɪs]有关

The comprehensive understanding of tumour-associated normal cells in tumour tissues may provide important insights into tumour biology and **aid in** the development of robust prognostic and predictive models.
对肿瘤组织中与 肿瘤相关的正常细胞的 全面的理解 可以 对肿瘤生物学提供重要的 见解[ˈɪnˌsaɪts]和 **在**可靠[roʊˈbʌst]预后和预测模型 的建立**方面提供了帮助**

Gene expression profiling of cancer has **resulted in** the identification of molecular subtypes and the development of models for predication prognostic and has enriched our knowledge of the molecular pathways of tumorigenesis.
癌症基因表达谱 已经**导致**分子亚型的识别 和预测预后模型的发展 并丰富了我们对 肿瘤发生分子通路 的认识

Increasing evidence suggests that the infiltration of tumour-associated normal cells influences the analysis of clinical tumour samples by genomic approaches,
越来越多的证据表明 肿瘤相关的正常细胞 的浸润 影响了 基因组方法 分析 临床肿瘤样本

such as gene expression profiles or copy number data,and biological interpretation of the results requires considerable **attention to** sample heterogeneity.
例如基因表达谱或者拷贝数数据和生物学解释结果 需要 相当[kənˈsɪdərəbl]**注意** 样本多样性[ˌhɛtərədʒɪˈniəti]

Several methods have been proposed to estimate the fraction of tumour cells in clinical tumour samples by using DNA copy number array data or by using next-generation sequencing data.
多种方法 已经 被提出[prəˈpoʊzd] 通过使用拷贝数阵列数据或者下一代基因测序数据 去估计 临床肿瘤样本的 肿瘤细胞 比例[ˈfrækʃn]

DNA copy number-based estimation of tumour purity is rapidly gaining traction in predicting the purity of tumour samples;
基于DNA拷贝数 估计肿瘤纯度 在肿瘤样本纯度预测中 迅速得到 重视[ˈtrækʃn]

however,such methods **are limited to** samples with available copy number profiles.
但是,这种方法 **仅限于**具有可用的拷贝数谱的 样本

Previous studies have attempted to deconvolve gene expression data into gene expression profiles from their constituent cellular fractions,
先前的研究尝试 从它们的组成[kənˈstɪtʃuənt]细胞部分 去分解(deconvolve)基因表达数据到基因表达谱

whereas have focused on deconvolution of microarray data obtained from normal tissue into cell-type-specific profiles,by calculating enrichment scores.
然而,重点是 通过 计算[ˈkælkjuleɪtɪŋ]富集分数 将 来自正常组织获得的 微阵列数据分解为 细胞类型特异性图谱

These methods **take advantage of** the differences in transcriptome properties of distinct cell types.
这些方法 **利用** 不同细胞类型的 转录组 特性[ˈprɑpərtiz]差异

Here we present a new algorithm that advantage of the unique properties of the transcriptional profiles of cancer samples to infer tumour cellularity **as well as** the different infiltrating normal cells,
这里我们提出了一种新算法 用癌症样品的转录谱的 独特 特性[ˈprɑpərtiz]的优势 去推测肿瘤细胞 **以及** 不同浸润的正常细胞

called ESTIMATE(Estimation of STromal and Immune cells in MAlignant Tumour tissues using Expression data).
被称作 ESTIMSTE(使用在 表达数据 估计 恶性[məˈlɪɡnənt]肿瘤组织的 基质和免疫细胞)

We **focus on** stromal and immune cells that form the major non-tumour constituents of tumour samples and identify specific signatures related to the infiltration of stromal and immune cells in tumour tissues.
我们**关注** 肿瘤样本的 主要的 非肿瘤成分的 基质和免疫细胞 和 识别 在肿瘤组织中的 基质和免疫细胞的 浸润性 相关的 特异性 特征[ˈsɪgnətʃərz]

By performing single-sample gene set-enrichment analysis(ssGSEA),
通过使用 单样本基因集富集分析

we calculate stromal and immune scores to predict the level of infiltrating stromal and immune cells and these **form the basis for** the ESTIMATE score to infer tumour purity in tumour tissue.
我们计算了基质和肿瘤得分 去预测 浸润性基质和免疫细胞的水平 这些**形成了**ESTIMATE分数推断 在肿瘤组织中的肿瘤纯度**的基础**

Finally,we describe the biological characteristics of stromal and immune scores in The Cancer Genome Atlas(TCGA) data sets.
最后,我们 描述了 TCGA数据集中 基质和免疫评分的 生物学特性