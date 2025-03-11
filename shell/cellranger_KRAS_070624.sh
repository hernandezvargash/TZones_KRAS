
###############################
# TZ Zones project
# run on @nanotower
# 070624
# KRAS Dysplasia
# test with force-cells to recover neutrophils
###############################

# run 202110043_R400_GGG
# 28/12/2021
# KRAS Dysplasia
# Run 400

###############################

export PATH='/home/hh/Software/cellranger-7.2.0':$PATH

cd /home/hh/Documents

ref=/home/hh/Software/cellranger_refs/refdata-gex-mm10-2020-A
fastqdir='/media/hh/MyBook/single_cell/TZones/KRAS/Dysplasia/202110043_R400_GGG/cellranger_mkfastq/Run_400_NS500-308_16-12-2021_GGG/S004453'

cellranger count --id=dysplasia_281221 \
--fastqs=$fastqdir \
--sample=1_KRAS_Dysplasie_bas_grade \
--transcriptome=$ref \
--force-cells=10000 \
--include-introns=true \
--localcores=12

cellranger count --id=dysplasia_281221b \
--fastqs=$fastqdir \
--sample=1_KRAS_Dysplasie_bas_grade \
--transcriptome=$ref \
--force-cells=15000 \
--include-introns=true \
--localcores=12


###############################



