# big-data-analysis
homework code 
Methods: Raw RNA-seq read were wuality assesed using FastQC (v0.12.1) (S. Andrews, 2010) and MultiQC (v1.33) (Ewels et al., 2016). Adapter trimming was performed using Trim Galore (v0.6.10) (Krueger, 2021) with default paired-end parameters. RNA-Seq reads were mapped to the H19 reference genome using HISAT2 (v. 2.2.1) (K. Daehwan et al., 2019) with default parameters. Alignment quality was assesed using samtools (1.19.2R) (Li et al., 2009) and RseQC (v5.0.4) (Wang et al., 2012). Read quantification was performed using feature counts (v. 2.0.6) (Liaoet al., 2014), using ensembl GRCh37.87 annotation with paired-end counting and unstranded library settings. Differential expression analysis was performed using R DESeq2 (v. 1.50.2) package (MI Love et al, 2014). Gene ontology over-representation analysis and gene set enrichment analysis were performed using clusterProfiler (v. 4.18.4) (Yu et al., 2012).

Samples:
GSM4505877 – normal N2 RNA-Seq
GSM4505880 – normal N8 RNA-Seq
GSM4505883 – normal N12 RNA-Seq

GSM4505887 – Tumor T2 RNA-Seq
GSM4505890 – Tumor T8 RNA-Seq
GSM4505893 – Tumor T12 RNA-Seq

Three normal (healthy cells) samples, three esophageal squamous cell carcinoma (tumour cells) samples. Number indicates replica, patient (ex. N2 and T2 has been taken from the same patient).
