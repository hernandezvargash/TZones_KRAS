
###############################
# TZ Zones project
# run on @nanotower
###############################

# run 202205025_R425_GGG
# 10/08/2022
# KRAS Hyperplasia
# Run 425

###############################

export PATH=/home/hh/software/cellranger-6.1.2:$PATH

cd /home/hh/scP/TZones

ref=/home/hh/software/cellranger_refs/refdata-gex-mm10-2020-A
fastqdir=/mnt/sdb/projects/single_cell/TZones/KRAS_model/202205025_R425_GGG/fastq_path/

# align:

cellranger count --id=hyperplasia_100822 \
--fastqs=$fastqdir \
--sample=1_KRAS_Hyperplasie_ \
--transcriptome=$ref \
--expect-cells=10000 \
--localcores=12

# run velocity

conda activate velocyto

velocyto run10x hyperplasia_100822/ /home/hh/software/cellranger_refs/refdata-gex-mm10-2020-A/genes/genes.gtf

conda deactivate

###############################



