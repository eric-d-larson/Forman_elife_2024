source samples.sh
source config.sh

python $rmats --b1 wt_125.txt \
	--b2 trt_125.txt \
	--gtf GRCm38_gencode_vM26.gtf \
	--od results/star_pass2/rmats -t paired \
	--nthread 8 \
	--readLength 125 \
	--cstat 0.0001 \
	--libType fr-secondstrand
