
###############################
# TZ Zones project
# run on @nanotower
# 250821
###############################

# run 380
# 13/07/2021
# Normal epithelial anorectal TZ with stroma

###############################

export PATH=/home/hh/software/cellranger-6.0.1:$PATH

cd /home/hh/scP/10x_GG

ref=/home/hh/software/cellranger_refs/refdata-cellranger-mm10-3.0.0
fastqdir=/home/hh/scP/10x_GG/Run_380/fastq/

# align:

cellranger count --id=normal_130721 \
--fastqs=$fastqdir \
--sample=1_NormalTZandStroma_ARNm \
--transcriptome=$ref \
--expect-cells=5000 \
--localcores=12

# run velocity

conda activate velocyto

velocyto run10x normal_130721/ /home/hh/software/cellranger_refs/refdata-cellranger-mm10-3.0.0/genes/genes.gtf

conda deactivate

###############################

# run 381
# 26/07/2021
# Hyperplasia epithelial anorectal TZ with stroma (KRAS expression alone in K17pos cells)
# aka "Hyperplasia 1"

###############################

export PATH=/home/hh/software/cellranger-6.0.1:$PATH

cd /home/hh/scP/10x_GG/

ref=/home/hh/software/cellranger_refs/refdata-cellranger-mm10-3.0.0
fastqdir=/home/hh/scP/10x_GG/Run_381/fastq/

# align:

cellranger count --id=KRAS_260721 \
--fastqs=$fastqdir \
--sample=2_Stroma_ARNm \
--transcriptome=$ref \
--expect-cells=5000 \
--localcores=12

# run velocity

conda activate velocyto

velocyto run10x KRAS_260721/ /home/hh/software/cellranger_refs/refdata-cellranger-mm10-3.0.0/genes/genes.gtf

conda deactivate

###############################





