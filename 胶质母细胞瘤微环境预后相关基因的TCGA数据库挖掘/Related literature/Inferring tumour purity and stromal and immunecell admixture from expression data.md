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