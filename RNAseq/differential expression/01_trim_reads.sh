
source config.sh
source samples.sh


for i in ${SAMPLES[$(($LSB_JOBINDEX - 1))]};
do 
    $bbduk -Xmx7g t=4 \
    in1=fastq/${i}_L003_R1_001.fastq.gz \
    in2=fastq/${i}_L003_R2_001.fastq.gz \
    out1=fastq/trim/${i}_R1trim.fastq.gz \
    out2=fastq/trim/${i}_R2trim.fastq.gz \
    ref=$bbadapters \
    ktrim=r \
    k=23 \
    mink=11 \
    hdist=1 \
    tpe=t \
    tbo=t \
    qtrim=rl \
    trimq=5 \
    minlength=37
done

