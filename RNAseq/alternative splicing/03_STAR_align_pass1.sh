
source samples.sh


for i in ${SAMPLES[@]};
do
    mkdir results/star_pass1/${i};
    STAR --runThreadN 16 \
	--outSAMtype BAM Unsorted \
	--sjdbOverhang 124 \
	--outFilterType BySJout \
	--outFilterMismatchNoverLmax 0.04 \
	--alignEndsType EndToEnd \
	--alignSJDBoverhangMin 4 \
	--alignIntronMax 1000000 \
	--alignMatesGapMax 1000000 \
	--alignSJoverhangMin 8 \
	--alignIntronMin 20 \
	--quantMode GeneCounts \
	--outSAMstrandField intronMotif \
	--genomeDir star279a_124_gencode_vM26/ \
	--limitBAMsortRAM 10000000000 \
	--outFileNamePrefix results/star_pass1/${i}/${i}. \
	--readFilesCommand zcat \
	--readFilesIn fastq/trim/${i}_R1trim125.fastq.gz fastq/trim/${i}_R2trim125.fastq.gz
done
