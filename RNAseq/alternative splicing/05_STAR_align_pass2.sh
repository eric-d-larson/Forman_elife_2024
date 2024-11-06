
source scripts/samples.sh


for i in ${SAMPLES[@]};
do
    mkdir results/star_pass2/${i}; 
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
	--sjdbFileChrStartEnd results/star_pass1/21_S40/21_S40.SJ.out.exp.tab results/star_pass1/22_S41/22_S41.SJ.out.exp.tab results/star_pass1/23_S42/23_S42.SJ.out.exp.tab results/star_pass1/24_S43/24_S43.SJ.out.exp.tab \
	--outFileNamePrefix results/star_pass2/${i}/${i}. \
	--readFilesCommand zcat \
	--readFilesIn fastq/trim/${i}_R1trim125.fastq.gz fastq/trim/${i}_R2trim125.fastq.gz
done

