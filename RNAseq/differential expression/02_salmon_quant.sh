
source config.sh
source scripts/samples.sh


for i in ${SAMPLES[@]};
do
salmon quant -i $sal_index_mm -l A \
    --validateMappings \
    -1 fastq/trim/${i}_R1trim.fastq.gz \
    -2 fastq/trim/${i}_R2trim.fastq.gz \
    -p 12 -o results/salmon/${i} \
	--numBootstraps 10
done
