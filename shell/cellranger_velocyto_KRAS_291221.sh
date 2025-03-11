
###############################
# TZ Zones project
# run on @nanotower
# 291221
# KRAS Dysplasia (last run of 4 for the Guasch part of the TZ Zones project)
###############################

# run 202110043_R400_GGG
# 28/12/2021
# KRAS Dysplasia
# Run 400

###############################

export PATH=/home/hh/software/cellranger-6.1.2:$PATH

cd /home/hh/scP/10x_GG

ref=/home/hh/software/cellranger_refs/refdata-cellranger-mm10-3.0.0
fastqdir=/mnt/sdb/projects/single_cell/TZones/202110043_R400_GGG/cellranger_mkfastq/fastq_path/Run_400_NS500-308_16-12-2021_GGG/S004453/

# align:

cellranger count --id=dysplasia_281221 \
--fastqs=$fastqdir \
--sample=1_KRAS_Dysplasie_bas_grade \
--transcriptome=$ref \
--expect-cells=10000 \
--localcores=12

# run velocity

conda activate velocyto

velocyto run10x dysplasia_281221/ /home/hh/software/cellranger_refs/refdata-cellranger-mm10-3.0.0/genes/genes.gtf

conda deactivate

###############################



