# ChIP_seq-Pipeline=bowtie2+macs2

## 一、总览

- 去接头前的质控报告：fastqc
- 去接头：trim_galore
- 去接头后的质控报告：fastqc
- 比对：bowtie2
- 去除PCR重复：picard
- 搜峰：macs2

## 二、双末端流程

### 2.1、输入数据

~~~jy
01_raw_data/case_1.fastq
01_raw_data/case_2.fastq
~~~

### 2.2、去接头前的质控报告

~~~jy
fastqc \
    01_raw_data/case_1.fastq 01_raw_data/case_2.fastq \
    -o 02_qc_before_fastqc
~~~

### 2.3、去接头

~~~jy
trim_galore -q 20 --phred33 --stringency 3 --length 20 -e 0.1 --paired --gzip \
    01_raw_data/case_1.fastq 01_raw_data/case_2.fastq \
    -o 03_qc
~~~

### 2.4、去接头后的质控报告

~~~jy
fastqc \
    03_qc/case_1_val_1.fq.gz 03_qc/case_2_val_2.fq.gz \
    -o 04_qc_after_fastqc
~~~

### 2.5、Alignment to genome

将read比对到参考基因组，最消耗时间的一步

~~~jy
bowtie2 -p 10 -k 1 \
    -x bowtie2_index/hg19/hg19 \
    -1 03_qc/case_1_val_1.fq.gz -2 03_qc/case_2_val_2.fq.gz \
    -S 05_alignment_to_genome/case.sam
~~~

---

### 2.6、sam2bam

将sam文件转化为更好存储的bam文件，bam文件是sam文件的二进制文件，可以更好存储数据。同时对bam文件进行排序，然后再构建索引。

#### sam文件转bam文件

~~~jy
samtools view -b -S 05_alignment_to_genome/case.sam > 05_alignment_to_genome/case.bam 
~~~

#### 对bam文件进行排序

~~~jy
samtools sort 05_alignment_to_genome/case.bam > 05_alignment_to_genome/case.sort.bam
~~~

#### 对排序后的bam文件构建索引

~~~jy
samtools index 05_alignment_to_genome/case.sort.bam 05_alignment_to_genome/case.sort.bam.bai
~~~

---

### 2.7、去除PCR重复

~~~jy
picard MarkDuplicates METRICS_FILE= H3K27Ac.markDup.metric REMOVE_DUPLICATES=true ASSUME_SORTED=true CREATE_INDEX=true \
    INPUT=05_alignment_to_genome/case.sort.bam \
    OUTPUT=06_rm_pcr_duplicates/case.final.sort.bam
~~~

### 2.8、Calling peak

将case组和control组，得到相关峰

~~~jy
macs2 callpeak -q 0.05 -f BAMPE -g hs \
-t case.sort.bam -c control.sort.bam \
-n out_name
~~~

-t case组，-c control组，-n 生成文件的名字

## 三、单末端流程

### 3.1、输入数据

~~~jy
01_raw_data/case.fastq
~~~

### 3.2、去接头前的质控报告

~~~jy
fastqc \
    01_raw_data/case.fastq \
    -o 02_qc_before_fastqc
~~~

### 3.3、去接头

~~~jy
trim_galore -q 20 --phred33 --stringency 3 --length 20 -e 0.1 --paired --gzip \
    01_raw_data/case.fastq \
    -o 03_qc
~~~

### 3.4、去接头后的质控报告

~~~jy
fastqc \
    03_qc/case_val_1.fq.gz \
    -o 04_qc_after_fastqc
~~~

### 3.5、Alignment to genome

将read比对到参考基因组，最消耗时间的一步

~~~jy
bowtie2 -p 10 -k 1 \
    -x bowtie2_index/hg19/hg19 \
    -U 03_qc/case_val_1.fq.gz \
    -S 05_alignment_to_genome/case.sam
~~~

---

### 3.6、sam2bam

将sam文件转化为更好存储的bam文件，bam文件是sam文件的二进制文件，可以更好存储数据。同时对bam文件进行排序，然后再构建索引。

#### sam文件转bam文件

~~~jy
samtools view -b -S 05_alignment_to_genome/case.sam > 05_alignment_to_genome/case.bam 
~~~

#### 对bam文件进行排序

~~~jy
samtools sort 05_alignment_to_genome/case.bam > 05_alignment_to_genome/case.sort.bam
~~~

#### 对排序后的bam文件构建索引

~~~jy
samtools index 05_alignment_to_genome/case.sort.bam 05_alignment_to_genome/case.sort.bam.bai
~~~

---

### 3.7、去除PCR重复

~~~jy
picard MarkDuplicates METRICS_FILE= H3K27Ac.markDup.metric REMOVE_DUPLICATES=true ASSUME_SORTED=true CREATE_INDEX=true \
    INPUT=05_alignment_to_genome/case.sort.bam \
    OUTPUT=06_rm_pcr_duplicates/case.final.sort.bam
~~~

### 3.8、Calling peak

将case组和control组，得到相关峰

~~~jy
macs2 callpeak -q 0.05 -f BAM -g hs \
-t 05_alignment_to_genome/case.sort.bam -c 05_alignment_to_genome/control.sort.bam \
-n out_name
~~~

-t case组，-c control组，-n 生成文件的名字