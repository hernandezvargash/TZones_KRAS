
###############################
# TZ Zones project
# run on @nanotower
# 220223
###############################

# run 202301005_R447_GGG
# 21/02/2023
# KRAS Hyperplasia (3rd try)

###############################

export PATH=/home/hh/software/cellranger-6.1.2:$PATH

cd '/mnt/sdb/projects/single_cell/TZones/KRAS_model/202301005_R447_GGG/' 

ref=/home/hh/software/cellranger_refs/refdata-gex-mm10-2020-A
fastqdir=/mnt/sdb/projects/single_cell/TZones/KRAS_model/202301005_R447_GGG/cellranger_mkfastq/Run_447_NS500-357_16-02-2023_GGG/S005546

# align:

cellranger count --id=hyperplasia_220223 \
--fastqs=$fastqdir \
--sample=1_Kras_hyperplasie_ARNm \
--transcriptome=$ref \
--expect-cells=10000 \
--localcores=12

# run velocity


conda activate velocyto

# MT genes are absent in the last runs of velocyto
# to fix this, modified two lines of the counter.py code as explained here: https://github.com/velocyto-team/velocyto.py/issues/318 , 
# and here https://github.com/velocyto-team/velocyto.py/blob/0963dd2d/velocyto/counter.py#L282
velocyto run10x hyperplasia_220223/ /home/hh/software/cellranger_refs/refdata-gex-mm10-2020-A/genes/genes.gtf
#velocyto run10x hyperplasia_220223/ /mnt/sdb/refs/mm10.refGene.gtf.gz

conda deactivate

###############################



