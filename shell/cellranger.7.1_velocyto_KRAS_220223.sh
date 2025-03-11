
###############################
# TZ Zones project
# run on @nanotower
# 080323
###############################

# run 202301005_R447_GGG
# 21/02/2023
# KRAS Hyperplasia (3rd try)

###############################

#export PATH=/home/hh/software/cellranger-6.1.2:$PATH
export PATH=/mnt/sdb/software/cellranger-7.1.0:$PATH

cd '/mnt/sdb/projects/single_cell/TZones/KRAS_model/202301005_R447_GGG/' 

#ref=/home/hh/software/cellranger_refs/refdata-gex-mm10-2020-A
ref=/mnt/sdb/refs/cellranger_refs/refdata-gex-mm10-2020-A
fastqdir=/mnt/sdb/projects/single_cell/TZones/KRAS_model/202301005_R447_GGG/cellranger_mkfastq/Run_447_NS500-357_16-02-2023_GGG/S005546

# align:

cellranger count --id=hyperplasia_080323 \
--fastqs=$fastqdir \
--sample=1_Kras_hyperplasie_ARNm \
--transcriptome=$ref \
--expect-cells=10000 \
--localcores=12

# run velocity

conda activate velocyto

velocyto run10x hyperplasia_080323/ /mnt/sdb/refs/cellranger_refs/refdata-gex-mm10-2020-A/genes/genes.gtf

conda deactivate

###############################




