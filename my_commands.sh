#fastqc/multiqc
mkdir fastqcoutput
fastqc -t 5 *.fastq.gz -o fastqcoutput/
python -m multiqc /home/genetics/homework1/data/raw/fastqcoutput/ -o /home/genetics/homework1/data/multiqc_report/

#trimming

mkdir trimmed
for sample in SRR11647686 SRR11647689 SRR11647692 SRR11647696 SRR11647699 SRR11647702; do trim_galore --paired --cores 8 -o trimmed/ ${sample}_1.fastq.gz ${sample}_2.fastq.gz; done


#genomo indeksavimas

wget https://ftp.ensembl.org/pub/grch37/current/fasta/homo_sapiens/dna/
gunzip Homo_sapiens.GRCh37.dna.primary_assembly.fa.gz
mkdir indexed_genome
hisat2-build -p 8 Homo_sapiens.GRCh37.dna.primary_assembly.fa indexed_genome/indexed

#mapinimas
mkdir data/raw/sam
for SAMPLE in SRR11647686 SRR11647689 SRR11647692 SRR11647696 SRR11647699 SRR11647702; do
    hisat2 -p 6 -q -x indexed_genome/indexed \
    -1 trimmed/${SAMPLE}_1_val_1.fq.gz \
    -2 trimmed/${SAMPLE}_2_val_2.fq.gz \
    -S sam/${SAMPLE}_output.sam --verbose
done

mv *.sam sam/

mkdir bam
for SAMPLE in SRR11647686 SRR11647689 SRR11647692 SRR11647696 SRR11647699 SRR11647702; do   
   samtools view -bS sam/${SAMPLE}_output.sam > bam/${SAMPLE}_output.bam
done

#markdup + qc
for SHomo_sapiens.GRCh37.dna.primary_assembly.fa.gzAMPLE in SRR11647686 SRR11647689 SRR11647692 SRR11647696 SRR11647699 SRR11647702; do
    samtools collate -@ 4 -O -u bam/${SAMPLE}_output.bam | \
    samtools fixmate -@ 4 -m -u - - | \
    samtools sort -@ 4 -u - | \
    samtools markdup -@ 4 - bam/${SAMPLE}_markdup.bam

    echo "${SAMPLE} duplicates:" >> qc_stats.txt
    samtools view -c -f 1024 bam/${SAMPLE}_markdup.bam >> qc_stats.txt
    echo "${SAMPLE} flagstat:" >> qc_stats.txt
    samtools flagstat bam/${SAMPLE}_markdup.bam >> qc_stats.txt
    echo "---" >> qc_stats.txt
done

#coverage
for SAMPLE in SRR11647686 SRR11647689 SRR11647692 SRR11647696 SRR11647699 SRR11647702; do
    echo "${SAMPLE} coverage:" >> coverage_stats.txt
    samtools coverage bam/${SAMPLE}_markdup.bam >> coverage_stats.txt
    echo "---" >> coverage_stats.txt
done
mkdir stats
mv qc_stats.txt coverage_stats.txt stats/

#gene body coverage

##bam indexing
for SAMPLE in SRR11647686 SRR11647689 SRR11647692 SRR11647696 SRR11647699 SRR11647702; do
    samtools index -@ 8 bam/${SAMPLE}_markdup.bam
done

##bed
wget https://sourceforge.net/projects/rseqc/files/BED/Human_Homo_sapiens/hg19_Ensembl_gene.bed.gz
gunzip hg19_Ensembl_gene.bed.gz
sed 's/^chr//' hg19_Ensembl_gene.bed > hg19_Ensembl_gene_nochr.bed

mkdir genebody
for SAMPLE in SRR11647686 SRR11647689 SRR11647692 SRR11647696 SRR11647699 SRR11647702; do
   geneBody_coverage.py -r hg19_Ensembl_gene_nochr.bed     
   -i bam/${SAMPLE}_markdup.bam     
   -o genebody/genebody_${SAMPLE}; 
done

for SAMPLE in SRR11647686 SRR11647689 SRR11647692 SRR11647696 SRR11647699 SRR11647702; do
   geneBody_coverage.py -r hg19_Ensembl_gene_nochr.bed     
   -i bam/${SAMPLE}_markdup.bam     
   -o genebody/genebody_${SAMPLE};
done

# inner distance
mkdir inner_distance
for SAMPLE in SRR11647686 SRR11647689 SRR11647692 SRR11647696 SRR11647699 SRR11647702; do
    inner_distance.py -r hg19_Ensembl_gene_nochr.bed \
    -i bam/${SAMPLE}_markdup.bam \
    -o inner_distance/${SAMPLE}
done

# clipping profile
mkdir clipping_profile
for SAMPLE in SRR11647686 SRR11647689 SRR11647692 SRR11647696 SRR11647699 SRR11647702; do
    clipping_profile.py \
    -i bam/${SAMPLE}_markdup.bam \
    -o clipping_profile/${SAMPLE} \
    -s PE 
done


# junction anotation
mkdir junctions
for SAMPLE in SRR11647686 SRR11647689 SRR11647692 SRR11647696 SRR11647699 SRR11647702; do
    junction_annotation.py -r hg19_Ensembl_gene_nochr.bed \
    -i bam/${SAMPLE}_markdup.bam \
    -o junctions/${SAMPLE}
done

# ── 12. CORRELATION + PCA ────────────────────────────────────
multiBamSummary bins --bamfiles bam/*_markdup.bam -o out.npz
plotCorrelation -in out.npz -c spearman -p heatmap -o correlation_plot.png
plotPCA -in out.npz -o pca_plot.png