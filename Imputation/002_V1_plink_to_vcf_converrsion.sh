## A protocol is available here: https://www.protocols.io/view/genotype-imputation-workflow-v3-0-nmndc5e?version_warning=no&step=1
## I am doing following work to prepare file to be imputed using sanger imputation service https://imputation.sanger.ac.uk/ 


# Files are bacth wise
# Merging all plink file


plink2=/rds/project/rds-csoP2nj6Y6Y/msr52/SOFTWARES/Apps/plink

# V1 chip non-imputed plink files: 
/rds/project/rds-csoP2nj6Y6Y/msr52/Imputation_required/v1_chip_extra


for i in {1..10}
do
echo /rds/project/rds-csoP2nj6Y6Y/msr52/Imputation_required/v1_chip_extra/batch${i}.nbr30_ids >> /rds/project/rds-csoP2nj6Y6Y/msr52/Imputation_required/v1_chip_extra/mergelist.txt
done

batch${i}.nbr30_ids

/rds/project/rds-csoP2nj6Y6Y/msr52/SOFTWARES/Apps/plink1.9 \
--merge-list /rds/project/rds-csoP2nj6Y6Y/msr52/Imputation_required/v1_chip_extra/mergelist.txt \
--make-bed \
--snps-only just-acgt \
--out /rds/project/rds-csoP2nj6Y6Y/msr52/Imputation_required/v1_chip_extra/v1_chip_extra_merged


# convert plink file to VCF

/rds/project/rds-csoP2nj6Y6Y/msr52/SOFTWARES/Apps/plink2 \
--bfile /rds/project/rds-csoP2nj6Y6Y/msr52/Imputation_required/v1_chip_extra/v1_chip_extra_merged \
--export vcf bgz id-paste=iid \
--out /rds/project/rds-csoP2nj6Y6Y/msr52/Imputation_required/v1_chip_extra/v1_chip_extra_merged



#Creating index file for new VCF
=
bcftools index v1_chip_extra_merged.vcf.gz

#Cheking merged VCF file
bcftools view v1_chip_extra_merged.vcf.gz | less
bcftools stats v1_chip_extra_merged.vcf.gz | awk '/^SN/'

bcftools query -l v1_chip_extra_merged.vcf.gz # sample column 

bcftools view v1_chip_extra_merged.vcf.gz | less

# filtering chromosome column: https://samtools.github.io/bcftools/howtos/query.html
query -f '%CHROM\n' v1_chip_extra_merged.vcf.gz > chr.txt

# There was chromosome x, y and MT


# Auto detect sample sex and fix ploidy in VCF


./configure --prefix=/rds/project/rds-csoP2nj6Y6Y/msr52/SOFTWARES/Apps
make
make install


export BCFTOOLS_PLUGINS=/rds/project/rds-csoP2nj6Y6Y/msr52/SOFTWARES/bcftools/plugins
/rds/project/rds-csoP2nj6Y6Y/msr52/SOFTWARES/Apps/bcftools/bcftools +guess-ploidy -g b37 /rds/project/rds-csoP2nj6Y6Y/msr52/Imputation_required/v1_chip_extra/v1_chip_extra_merged.vcf.gz > samples.txt



/rds/project/rds-csoP2nj6Y6Y/msr52/SOFTWARES/Apps/bcftools/bcftools +fixref /rds/project/rds-csoP2nj6Y6Y/msr52/Imputation_required/v1_chip_extra/v1_chip_extra_merged.vcf.gz -- -f /rds/project/rds-csoP2nj6Y6Y/msr52/Imputation_required/GRCh37.p13.genome.fa



#https://www.ncbi.nlm.nih.gov/variation/docs/human_variation_vcf/
#wget ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606_b146_GRCh37p13/VCF/All_20151104.vcf.gz

wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/human_g1k_v37.fasta.gz



# checking refrence alleles


## Chromosome name in the vcf file and fasta file need to be matched while using fixref plugin and fasta files. 

# chnage chromosome names
bcftools annotate -Oz --rename-chrs chr_change.txt v1_chip_extra_merged.vcf.gz > v1_renamed_merged.vcf.gz
# chr_change.txt - created manually to change chr names from one version to other
bcftools index v1_renamed_merged.vcf.gz
bcftools stats v1_renamed_merged.vcf.gz | awk '/^SN/'

#v1_chip_extra_merged.vcf.gz -- -f human_g1k_v37.fasta.gz

# First environments for plugins need to be set using export command: here is the path 'BCFTOOLS_PLUGINS=/rds/project/rds-csoP2nj6Y6Y/msr52/SOFTWARES/bcftools/plugins"

bgzip -c human_g1k_v37.fasta.gz > human_g1k_v37.fasta.gz # this file does not work, taken following link from sanger imputaion server
samtools faidx hg19.fa
samtools faidx GRCh37.p13.genome.fa


export BCFTOOLS_PLUGINS=/rds/project/rds-csoP2nj6Y6Y/msr52/SOFTWARES/bcftools/plugins

/rds/project/rds-csoP2nj6Y6Y/msr52/SOFTWARES/Apps/bcftools/bcftools +fixref /rds/project/rds-csoP2nj6Y6Y/msr52/Imputation_required/v1_chip_extra/v1_renamed_merged.vcf.gz -- -f /rds/project/rds-csoP2nj6Y6Y/msr52/Imputation_required/hg19.fa
# NS, Number of sites:
NS	total        	731205
NS	ref match    	619272	86.8%
NS	ref mismatch 	94318	13.2%
NS	skipped      	17615
NS	non-ACGT     	0
NS	non-SNP      	17615
NS	non-biallelic	0


wget 'ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606_b146_GRCh37p13/VCF/All_20151104.vcf.gz'


# Get the dbSNP annotation file. Make sure the correct reference build is used (e.g. b37)
#   https://www.ncbi.nlm.nih.gov/variation/docs/human_variation_vcf/
#   wget ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606_b146_GRCh37p13/VCF/All_20151104.vcf.gz  # NOT PRESENT


wget "https://ftp.ncbi.nih.gov/snp/organisms/human_9606_b150_GRCh37p13/VCF/00-All.vcf.gz"
wget "https://ftp.ncbi.nih.gov/snp/organisms/human_9606_b150_GRCh37p13/VCF/00-All.vcf.gz.tbi"


# again chromosome name was differnt than the lates vcf file used, so need to update chromosome file.
# as hg19 file has chr in it, dbsnp file has to be same

# Swap the alleles

# checking how chromosome names were in the file
/rds/project/rds-csoP2nj6Y6Y/msr52/SOFTWARES/Apps/bcftools/bcftools  query -f '%CHROM\n' 00-All.vcf.gz > 00_All_chr.txt
# .text file was checked with R

# changing chromosome names in dbsnp
module add bcftools-1.9-gcc-5.4.0-b2hdt5n
bcftools annotate -Oz --rename-chrs chr_renamed_numtocharec.txt /rds/project/rds-csoP2nj6Y6Y/msr52/Imputation_required/00-All.vcf.gz > 00_All_renamed_dbsnpchr.vcf.gz


# creating .tbi file 
module add tabix-2013-12-16-gcc-5.4.0-xn3xiv7

bcftools view 00_All_renamed_dbsnpchr.vcf.gz -Oz -o 00_All_renamed_dbsnpchr.vcf.gz
bcftools index -f 00_All_renamed_dbsnpchr.vcf.gz



# creatingg file with 22 autosomes only as our vcf does not have X,Y and MT varriants.

bcftools view 00_All_renamed_dbsnpchr.vcf.gz --regions chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 > 00_All_renamed_dbsnpchr1_22.vcf.gz
bcftools view 00_All_renamed_dbsnpchr1_22.vcf.gz -Oz -o 00_All_renamed_dbsnpchr1_22.vcf.gz
bcftools index 00_All_renamed_dbsnpchr1_22.vcf.gz


# Swap the alleles

# Results are not different 
/rds/project/rds-csoP2nj6Y6Y/msr52/SOFTWARES/Apps/bcftools/bcftools +fixref v1_renamed_merged.vcf.gz -Oz -o v1_h19_aligned_merged.vcf.gz -- -d -f /rds/project/rds-csoP2nj6Y6Y/msr52/Imputation_required/hg19.fa -i 00_All_renamed_dbsnpchr1_22.vcf.gz
# this file shows truncated
bcftools view v1_h19_aligned_merged.vcf.gz-Oz -o v1_h19_aligned_merged.vcf.gz
bcftools index v1_h19_aligned_merged.vcf.gz
# using version one file: v1_h19_aligned_merged.vcf.gz
/rds/project/rds-csoP2nj6Y6Y/msr52/SOFTWARES/Apps/bcftools/bcftools +fixref /rds/project/rds-csoP2nj6Y6Y/msr52/Imputation_required/v1_chip_extra/v1_h19_aligned_merged.vcf.gz -- -f /rds/project/rds-csoP2nj6Y6Y/msr52/Imputation_required/hg19.fa

# SC, guessed strand convention
SC	TOP-compatible	0
SC	BOT-compatible	0
# NS, Number of sites:
NS	total        	626814
NS	ref match    	626814	100.0%
NS	ref mismatch 	0	0.0%
NS	skipped      	0
NS	non-ACGT     	0
NS	non-SNP      	0
NS	non-biallelic	0

/rds/project/rds-csoP2nj6Y6Y/msr52/SOFTWARES/Apps/bcftools/bcftools +fixref v1_renamed_merged.vcf.gz -Oz -o Version2_v1_h19_aligned_merged.vcf.gz -- -d -f /rds/project/rds-csoP2nj6Y6Y/msr52/Imputation_required/hg19.fa -i 00_All_renamed_dbsnpchr.vcf.gz
bcftools view Version2_v1_h19_aligned_merged.vcf.gz -Oz -o Version2_v1_h19_aligned_merged.vcf.gz
bcftools index Version2_v1_h19_aligned_merged.vcf.gz
# Results again same as previous one 

# using version 2 file: Version2_v1_h19_aligned_merged.vcf.gz
/rds/project/rds-csoP2nj6Y6Y/msr52/SOFTWARES/Apps/bcftools/bcftools +fixref /rds/project/rds-csoP2nj6Y6Y/msr52/Imputation_required/v1_chip_extra/Version2_v1_h19_aligned_merged.vcf.gz -- -f /rds/project/rds-csoP2nj6Y6Y/msr52/Imputation_required/hg19.fa
#There is something wrong with this file


# This is very extreme: just filippedd all other mismatched SNPs 
/rds/project/rds-csoP2nj6Y6Y/msr52/SOFTWARES/Apps/bcftools/bcftools +fixref v1_renamed_merged.vcf.gz -Oz -o vresion3_v1_h19_aligned_merged.vcf.gz -- -f /rds/project/rds-csoP2nj6Y6Y/msr52/Imputation_required/hg19.fa -m flip -d
/rds/project/rds-csoP2nj6Y6Y/msr52/SOFTWARES/Apps/bcftools/bcftools index  vresion3_v1_h19_aligned_merged.vcf.gz

# SC, guessed strand convention
SC	TOP-compatible	0
SC	BOT-compatible	0
# NS, Number of sites:
NS	total        	731205
NS	ref match    	619272	86.8%
NS	ref mismatch 	94318	13.2%
NS	flipped      	6	0.0%
NS	swapped      	89313	12.5%
NS	flip+swap    	1	0.0%
NS	unresolved   	45417	6.4%
NS	fixed pos    	0	0.0%
NS	skipped      	17615
NS	non-ACGT     	0
NS	non-SNP      	17615
NS	non-biallelic	0

# using this final version: version3_v1_h19_aligned_merged.vcf.gz
/rds/project/rds-csoP2nj6Y6Y/msr52/SOFTWARES/Apps/bcftools/bcftools +fixref /rds/project/rds-csoP2nj6Y6Y/msr52/Imputation_required/v1_chip_extra/vresion3_v1_h19_aligned_merged.vcf.gz -- -f /rds/project/rds-csoP2nj6Y6Y/msr52/Imputation_required/hg19.fa
# SC, guessed strand convention
SC	TOP-compatible	0
SC	BOT-compatible	0
# NS, Number of sites:
NS	total        	668173
NS	ref match    	668173	100.0%
NS	ref mismatch 	0	0.0%
NS	skipped      	0
NS	non-ACGT     	0
NS	non-SNP      	0
NS	non-biallelic	0




wget -O 1000g_af.vcf.gz "ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.wgs.phase3_shapeit2_mvncall_integrated_v5c.20130502.sites.vcf.gz"

/rds/project/rds-csoP2nj6Y6Y/msr52/SOFTWARES/Apps/bcftools/bcftools 1000g_af.vcf.gz

/rds/project/rds-csoP2nj6Y6Y/msr52/SOFTWARES/Apps/bcftools/bcftools annotate -c INFO/AF -a 1000g_af.vcf.gz vresion3_v1_h19_aligned_merged.vcf.gz | /rds/project/rds-csoP2nj6Y6Y/msr52/SOFTWARES/Apps/bcftools/bcftools +af-dist | grep ^PROB > data.dist.txt



wget -O hg19_dbsnp_freq.vcf.gz "https://ftp.ncbi.nih.gov/snp/latest_release/VCF/GCF_000001405.25.gz"
wget -O hg19_dbsnp_freq.vcf.gz.tbi "https://ftp.ncbi.nih.gov/snp/latest_release/VCF/GCF_000001405.25.gz.tbi"

wget "https://www.well.ox.ac.uk/~wrayner/tools/1000GP_Phase3_combined.legend.gz"

# Details available here: 
https://www.well.ox.ac.uk/~wrayner/tools/ 

## converting VCF file to plink file to be checked with HRC-1000G-check-bim.pl program:
 
# First need to convert chr names from ChrXX to XX
module load bcftools-1.9-gcc-5.4.0-b2hdt5n
bcftools annotate -Oz --rename-chrs chr_rename_chractonum_1_22.txt vresion3_v1_h19_aligned_merged.vcf.gz > vresion3_v1_h19_aligned_merged_renamed.vcf.gz


# VCF to Plink conversion
mocule load plink/2.00-alpha

plink \
--vcf vresion3_v1_h19_aligned_merged_renamed.vcf.gz \
--make-bed \
--out vresion3_v1_h19_aligned_merged_renamed_temp

# creating .freq file required for frequency checking
plink \
--bfile vresion3_v1_h19_aligned_merged_renamed_temp \
--freq \
--out vresion3_v1_h19_aligned_merged_renamed_temp 

plink \
--bfile v1_chip_extra_merged \
--freq \
--out v1_chip_extra_merged 


# Allele frequency checking 

# 1000 genome as reference panel: .bim file used here created from VCF filr after QC
perl /rds/project/rds-csoP2nj6Y6Y/msr52/SOFTWARES/Apps/HRC-1000G-check-bim.pl \
-b vresion3_v1_h19_aligned_merged_renamed_temp.bim \
-f vresion3_v1_h19_aligned_merged_renamed_temp.frq \
-r /rds/project/rds-csoP2nj6Y6Y/msr52/Imputation_required/1000GP_Phase3_combined.legend.gz \
-g \
-p EUR 

# HRC as reference panel: .bim file used here created from VCF filr after QC
perl /rds/project/rds-csoP2nj6Y6Y/msr52/SOFTWARES/Apps/HRC-1000G-check-bim.pl \
-b vresion3_v1_h19_aligned_merged_renamed_temp.bim \
-f vresion3_v1_h19_aligned_merged_renamed_temp.frq \
-r /rds/project/rds-csoP2nj6Y6Y/msr52/Imputation_required/HRC.r1-1.GRCh37.wgs.mac5.sites.tab \
-h \
-p EUR 

# 1000 genome as reference panel: .bim file used here comes from original files when merged together, no QC was performed 
perl /rds/project/rds-csoP2nj6Y6Y/msr52/SOFTWARES/Apps/HRC-1000G-check-bim.pl \
-b v1_chip_extra_merged.bim \
-f v1_chip_extra_merged.frq \
-r /rds/project/rds-csoP2nj6Y6Y/msr52/Imputation_required/1000GP_Phase3_combined.legend.gz \
-g \
-p EUR 

# HRC as reference panel: .bim file used here comes from original files when merged together, no QC was performed 
perl /rds/project/rds-csoP2nj6Y6Y/msr52/SOFTWARES/Apps/HRC-1000G-check-bim.pl \
-b v1_chip_extra_merged.bim \
-f v1_chip_extra_merged.frq \
-r /rds/project/rds-csoP2nj6Y6Y/msr52/Imputation_required/HRC.r1-1.GRCh37.wgs.mac5.sites.tab \
-h \
-p EUR 











