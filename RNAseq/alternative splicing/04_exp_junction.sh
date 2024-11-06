
source scripts/samples.sh


for i in ${SAMPLES[@]};
do
    cat results/star_pass1/${i}/${i}.SJ.out.tab | awk '$7 >0' > results/star_pass1/${i}/${i}.SJ.out.exp.tab
done

