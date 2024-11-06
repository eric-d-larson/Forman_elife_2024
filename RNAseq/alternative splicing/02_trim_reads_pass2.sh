
source config.sh
source samples.sh


for i in ${SAMPLES[$(($LSB_JOBINDEX - 1))]};
do 
    $bbduk -Xmx7g t=4 \
    in1=fastq/trim/${i}_R1trim.fastq.gz \
    in2=fastq/trim/${i}_R2trim.fastq.gz \
    out1=fastq/trim/${i}_R1trim125.fastq.gz \
    out2=fastq/trim/${i}_R2trim125.fastq.gz \
    ftr=124 \
    minlength=125
done

