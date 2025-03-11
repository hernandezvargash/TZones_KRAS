
###############################
# TZ Zones project
# run on @nanotower
# 260723
# align all KRAS samples to ZsGreen reporter gene
###############################


###############################
# Custom Reference
###############################

cd '/mnt/sdb/refs/cellranger_refs'

# custom reference genome (combining mm10 and GFP) was created as in:
# https://github.com/igordot/genomics/blob/master/workflows/ref-genome-gfp.md

## genome.ZsGreen.fa
# kept few bases before and after coding sequence
https://www.addgene.org/browse/sequence_vdb/4762/
>pZsGreen1-N1

# remove hidden characters introduced by libre office
dos2unix genome.ZsGreen.fa

# to get sequence length do:
cat genome.ZsGreen.fa | grep -v "^>" | tr -d "\n" | wc -c
# 994 in the "UTR" version
# 696 in the selected version

## genes.ZsGreen.gtf
echo -e 'ZsGreen\tunknown\texon\t1\t696\t.\t+\t.\tgene_id "ZsGreen"; transcript_id "ZsGreen"; gene_name "ZsGreen"; gene_biotype "protein_coding";' > genes.ZsGreen.gtf

# old version of gtf file:
#ZsGreen	unknown	gene	1	994	.	+	.	gene_id "ZsGreen"; gene_name "ZsGreen"; gene_biotype "protein_coding";
#ZsGreen	unknown	transcript	1	994	.	+	.	gene_id "ZsGreen"; transcript_id "ZsGreen"; gene_name "ZsGreen"; gene_biotype "protein_coding";
#ZsGreen	unknown	exon	1	994	.	+	.	gene_id "ZsGreen"; transcript_id "ZsGreen"; gene_name "ZsGreen"; gene_biotype "protein_coding";

# join genes and genome files
cat ./refdata-gex-mm10-2020-A/fasta/genome.fa genome.ZsGreen.fa > genome.mm10.ZsGreen.fa
cat ./refdata-gex-mm10-2020-A/genes/genes.gtf genes.ZsGreen.gtf > genes.mm10.ZsGreen.gtf

# check concatenation
grep ">" genome.mm10.ZsGreen.fa
tail genes.mm10.ZsGreen.gtf

# create custom reference
export PATH=/home/hh/software/cellranger-6.1.2:$PATH
cellranger mkref --genome=refdata-cellranger-mm10.ZsGreen --fasta=genome.mm10.ZsGreen.fa --genes=genes.mm10.ZsGreen.gtf


###############################

# run 380
# 13/07/2021
# Normal epithelial anorectal TZ with stroma

###############################

export PATH=/home/hh/software/cellranger-6.1.2:$PATH

cd '/mnt/sdb/projects/single_cell/TZones/KRAS_model'

ref=/mnt/sdb/refs/cellranger_refs/refdata-cellranger-mm10.ZsGreen

fastqdir='/media/hh/MyBook/single_cell/TZones/KRAS/Normal/Run_380/cellranger_mkfastq/Run_380_NS500-288_13-07-2021_GGG/S004286'

# align:

cellranger count --id=normal_130721 \
--fastqs=$fastqdir \
--sample=1_NormalTZandStroma_ARNm \
--transcriptome=$ref \
--expect-cells=5000 \
--localcores=12


###############################

# run 381
# 26/07/2021
# Hyperplasia epithelial anorectal TZ with stroma (KRAS expression alone in K17pos cells)
# aka "Hyperplasia 1"

###############################

export PATH=/home/hh/software/cellranger-6.1.2:$PATH

cd '/mnt/sdb/projects/single_cell/TZones/KRAS_model'

ref=/mnt/sdb/refs/cellranger_refs/refdata-cellranger-mm10.ZsGreen

fastqdir='/media/hh/MyBook/single_cell/TZones/KRAS/Hyperplasia/Run_381/cellranger_mkfastq/Run_381_NS500-289_26-07-2021_GGG/S004306'

# align:

cellranger count --id=KRAS_260721 \
--fastqs=$fastqdir \
--sample=2_Stroma_ARNm \
--transcriptome=$ref \
--expect-cells=5000 \
--localcores=12


###############################

# run 202110043_R400_GGG
# 28/12/2021
# KRAS Dysplasia
# Run 400

###############################


export PATH=/home/hh/software/cellranger-6.1.2:$PATH

cd '/mnt/sdb/projects/single_cell/TZones/KRAS_model'

ref=/mnt/sdb/refs/cellranger_refs/refdata-cellranger-mm10.ZsGreen

fastqdir='/media/hh/MyBook/single_cell/TZones/KRAS/Dysplasia/202110043_R400_GGG/cellranger_mkfastq/Run_400_NS500-308_16-12-2021_GGG/S004453'

# align:

cellranger count --id=dysplasia_281221 \
--fastqs=$fastqdir \
--sample=1_KRAS_Dysplasie_bas_grade \
--transcriptome=$ref \
--expect-cells=10000 \
--localcores=12



###############################

# run 202301005_R447_GGG
# 21/02/2023
# KRAS Hyperplasia (3rd try)
# Run 447

###############################

export PATH=/home/hh/software/cellranger-6.1.2:$PATH

cd '/mnt/sdb/projects/single_cell/TZones/KRAS_model'

ref=/mnt/sdb/refs/cellranger_refs/refdata-cellranger-mm10.ZsGreen

fastqdir='/media/hh/MyBook/single_cell/TZones/KRAS/Hyperplasia/202301005_R447_GGG/cellranger_mkfastq/Run_447_NS500-357_16-02-2023_GGG/S005546'

# align:

cellranger count --id=hyperplasia_220223 \
--fastqs=$fastqdir \
--sample=1_Kras_hyperplasie_ARNm \
--transcriptome=$ref \
--expect-cells=10000 \
--localcores=12



###############################

# run 354
# KRAS Chronic wound 4 weeks (Carcinoma)
# fastq files downloaded on 090321
# 16,000 cells have been processed

###############################

export PATH=/home/hh/software/cellranger-6.1.2:$PATH

cd '/mnt/sdb/projects/single_cell/TZones/KRAS_model'

ref=/mnt/sdb/refs/cellranger_refs/refdata-cellranger-mm10.ZsGreen

fastqdir='/media/hh/MyBook/single_cell/TZones/KRAS/Carcinoma/202010029_R354_GGG/fastq'

# align:

cellranger count --id=carcinoma_090321 \
--fastqs=$fastqdir \
--sample=1_KRAS_chronic_wound_ARNm \
--transcriptome=$ref \
--expect-cells=10000 \
--localcores=12


###############################

# velocity did not work with the ZsGreen gtf file

# add exon_number to gtf file for velocity

cd '/mnt/sdb/refs/cellranger_refs'

echo -e 'ZsGreen\tunknown\texon\t1\t696\t.\t+\t.\tgene_id "ZsGreen"; exon_number 1; transcript_id "ZsGreen"; gene_name "ZsGreen"; gene_biotype "protein_coding";' > genes.ZsGreen.gtf
cat ./refdata-gex-mm10-2020-A/genes/genes.gtf genes.ZsGreen.gtf > genes.mm10.ZsGreen.velo.gtf
tail genes.mm10.ZsGreen.velo.gtf


# run velocity

cd '/mnt/sdb/projects/single_cell/TZones/KRAS_model'

conda activate velocyto

#ref=/mnt/sdb/refs/cellranger_refs/genes.mm10.ZsGreen.velo.gtf
ref=/mnt/sdb/refs/cellranger_refs/refdata-gex-mm10-2020-A/genes/genes.gtf

velocyto run10x normal_130721/ $ref
velocyto run10x KRAS_260721/ $ref
velocyto run10x dysplasia_281221/ $ref
velocyto run10x hyperplasia_220223/ $ref
velocyto run10x carcinoma_090321/ $ref


conda deactivate



###############################








