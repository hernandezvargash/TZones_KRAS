
###############################
# TZ Zones project
# run on @nanotower
# 130321
###############################

# run 354
# KRAS Chronic wound 4 weeks (Carcinoma)
# fastq files downloaded on 090321
# 16,000 cells have been processed

###############################

export PATH=/home/hh/software/cellranger-3.1.0:$PATH

cd /home/hh/scP/10x_GG

ref=/home/hh/software/cellranger_refs/refdata-cellranger-mm10-3.0.0
fastqdir=/home/hh/scP/10x_GG/202010029_R354_GGG/fastq/

# align:

cellranger count --id=wound_2 \
--fastqs=$fastqdir \
--sample=1_KRAS_chronic_wound_ARNm \
--transcriptome=$ref \
--expect-cells=10000 \
--localcores=12

# run velocity

conda activate velocyto

velocyto run10x wound_2/ /home/hh/software/cellranger_refs/refdata-cellranger-mm10-3.0.0/genes/genes.gtf

conda deactivate

###############################

# run 355
# normal TZ

###############################

export PATH=/home/hh/software/cellranger-3.1.0:$PATH

cd /home/hh/scP/10x_GG/

ref=/home/hh/software/cellranger_refs/refdata-cellranger-mm10-3.0.0
fastqdir=/home/hh/scP/10x_GG/normal_2/fastq/

# align:

cellranger count --id=normal_2b \
--fastqs=$fastqdir \
--sample=2_GFP-Cre4weeksTAM_ARNm \
--transcriptome=$ref \
--expect-cells=10000 \
--localcores=8

# run velocity

conda activate velocyto

velocyto run10x normal_2b/ /home/hh/software/cellranger_refs/refdata-cellranger-mm10-3.0.0/genes/genes.gtf

conda deactivate

###############################





