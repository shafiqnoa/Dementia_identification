
#Directory orginal data:
/rds/project/rds-csoP2nj6Y6Y/msr52

#hpc work directory
/home/msr52/rds/hpc-work

#Imputation of NBR30 batch 16 and 17
#Directory
/rds/project/rds-csoP2nj6Y6Y/msr52/Imputation_required/V2.1_B16_17
#Batch 16
v2chip_b16_V3_calls_1.vcf.gz
#Batch 17
v2chip_b17_V3_calls_1.vcf.gz


module add bcftools-1.9-gcc-5.4.0-b2hdt5n


#Combine two batches

# Making list of files to merge
ls -d *.gz >mergelist.txt

# merging using bcftools 

bcftools merge -O z -m both -l mergelist.txt --force-samples -o merged_b16_17.vcf.gz
# or
bcftools merge -Oz -l mergelist.txt --force-samples -o merged_b16_17.vcf.gz

## This merging attempt shows allleic mismatch thus followed up each file

# Some bacics using bcftools
bcftools view v2chip_b16_V3_calls_1.vcf.gz | less
bcftools view v2chip_b16_V3_calls_1.vcf.gz -s SP00300081923W | less # sampling one sample

# Filtering by position and sample
bcftools view v2chip_b16_V3_calls_1.vcf.gz -r 1:25405592 -s SP00300081923W | less # This is the position where mismatch occured, its a CNV

# Only bi-allelic SNPs are kept which also passed filter:
bcftools view v2chip_b16_V3_calls_1.vcf.gz -i'FILTER=="PASS"' -m2 -M2 -v snps -Oz > v2chip_b16_V3_SNPs.vcf.gz
bcftools view v2chip_b17_V3_calls_1.vcf.gz -i'FILTER=="PASS"' -m2 -M2 -v snps -Oz > v2chip_b17_V3_SNPs.vcf.gz

# How many samples and variants in each of those?
bcftools stats v2chip_b16_V3_calls_1.vcf.gz | awk '/^SN/'
bcftools stats v2chip_b16_V3_SNPs.vcf.gz | awk '/^SN/' # only bi-allelic SNPs are avialble in the file

bcftools stats v2chip_b17_V3_calls_1.vcf.gz | awk '/^SN/'
bcftools stats v2chip_b17_V3_SNPs.vcf.gz | awk '/^SN/' # only bi-allelic SNPs file

#For more similar check look #https://eriqande.github.io/eca-bioinf-handbook/basic-handling-of-vcf-files.html
#https://github.com/dantaki/videos/tree/master/bcftools#view

# Remove variants with any missing genotypes ("." or "./." or ".|." or "./0")
# The ~ is a matches operator, different from ==

bcftools view v2chip_b16_V3_SNPs.vcf.gz -e'GT~"\."' | less

bcftools view v2chip_b17_V3_SNPs.vcf.gz -e'GT~"\."' | less

# Indexing new file
bcftools index v2chip_b16_V3_SNPs.vcf.gz
bcftools index v2chip_b17_V3_SNPs.vcf.gz

# Mergging newly filtered files
ls -d *SNPs.vcf.gz > mergelist.txt

bcftools merge -Oz -l mergelist.txt -o merged_b16_17.vcf.gz

#Cheking merged VCF file
bcftools view merged_b16_17.vcf.gz | less
bcftools stats merged_b16_17.vcf.gz | awk '/^SN/'
# Number of remaining multi-allelic SNP sites 10


# Lift over from hg38 to hg19 (not needed for UKBBv1_0 chip, which comes referenced to hg19)
# Files are downloaded according to https://gitlab.haem.cam.ac.uk/cjp64/affy_impute 
# Downloading Chain file

wget --timestamping 'ftp://hgdownload.cse.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHg19.over.chain.gz' -O hg38ToHg19.over.chain.gz



wget --timestamping 'ftp://hgdownload.cse.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHg19.over.chain.gz' -O hg38ToHg19.over.chain.gz




# Latest fasta file was dowloded from http://hgdownload.cse.ucsc.edu/goldenpath/hg19/bigZips/ 
wget http://hgdownload.cse.ucsc.edu/goldenpath/hg19/bigZips/latest/hg19.fa.gz 
gunzip hg19.fa.gz 

# The following script were used by Chirs Penket from NIHR Bioresource 


##### Lift over from hg38 to hg19 (not needed for UKBBv1_0 chip, which comes referenced to hg19)

module load samtools/1.9 > /dev/null 2>&1
module load CrossMap/0.2.9 > /dev/null 2>&1


module add python/3.5 
module add samtools-1.9-gcc-5.4.0-vf6vvem
module load ceuadmin/crossmap/0.4 



#Installed Picard following https://github.com/broadinstitute/picard/
java -jar /rds/project/rds-csoP2nj6Y6Y/msr52/SOFTWARES/picard/build/libs/picard.jar -h 
#picard.jar file copied to Apps folder
#Enviornment was set following https://linuxize.com/post/how-to-set-and-list-environment-variables-in-linux/
PICARD='/rds/project/rds-csoP2nj6Y6Y/msr52/SOFTWARES/Apps/'
# Now the picard is accessible with following command 
java -jar $PICARD/picard.jar

#https://gatk.broadinstitute.org/hc/en-us/articles/360037060932-LiftoverVcf-Picard-

# Follow this for error
#https://gatk.broadinstitute.org/hc/en-us/community/posts/360062428371-LiftoverVcf-fails-to-lift-over-a-GVCF-file-Error-java-lang-ArrayIndexOutOfBoundsException-1

#java -jar $PICARD/picard.jar LiftoverVcf \\

#hg19.fa file require a .dic file to work on. the following command makes it. 
java -jar /rds/project/rds-csoP2nj6Y6Y/msr52/SOFTWARES/picard/build/libs/picard.jar CreateSequenceDictionary \
R=/rds/project/rds-csoP2nj6Y6Y/msr52/Imputation_required/hg19.fa \
O=/rds/project/rds-csoP2nj6Y6Y/msr52/Imputation_required/hg19.dict 


java -jar /rds/project/rds-csoP2nj6Y6Y/msr52/SOFTWARES/picard/build/libs/picard.jar CreateSequenceDictionary \
R=/rds/project/rds-csoP2nj6Y6Y/msr52/Imputation_required/GRCh37.p13.genome.fa \
O=/rds/project/rds-csoP2nj6Y6Y/msr52/Imputation_required/GRCh37.p13.genome.dict


java -jar /rds/project/rds-csoP2nj6Y6Y/msr52/SOFTWARES/picard/build/libs/picard.jar LiftoverVcf \
I=/rds/project/rds-csoP2nj6Y6Y/msr52/Imputation_required/V2.1_B16_17/merged_b16_17.vcf.gz \
O=/rds/project/rds-csoP2nj6Y6Y/msr52/Imputation_required/V2.1_B16_17/lifted_over_merged_b16_17.vcf.gz \
CHAIN=/rds/project/rds-csoP2nj6Y6Y/msr52/Imputation_required/hg38ToHg19.over.chain.gz \
REJECT=/rds/project/rds-csoP2nj6Y6Y/msr52/Imputation_required/V2.1_B16_17/rejected_variants_for_b16_17.vcf \
R=/rds/project/rds-csoP2nj6Y6Y/msr52/Imputation_required/hg19.fa


java -jar /rds/project/rds-csoP2nj6Y6Y/msr52/SOFTWARES/picard/build/libs/picard.jar LiftoverVcf \
I=/rds/project/rds-csoP2nj6Y6Y/msr52/Imputation_required/V2.1_B16_17/merged_b16_17.vcf.gz \
O=/rds/project/rds-csoP2nj6Y6Y/msr52/Imputation_required/V2.1_B16_17/lifted_over_merged_b16_17.vcf.gz \
CHAIN=/rds/project/rds-csoP2nj6Y6Y/msr52/Imputation_required/hg38ToHg19.over.chain.gz \
R=/rds/project/rds-csoP2nj6Y6Y/msr52/Imputation_required/hg19.fa


java -jar /rds/project/rds-csoP2nj6Y6Y/msr52/SOFTWARES/picard/build/libs/picard.jar LiftoverVcf \
I=/rds/project/rds-csoP2nj6Y6Y/msr52/Imputation_required/V2.1_B16_17/merged_b16_17.vcf.gz \
O=/rds/project/rds-csoP2nj6Y6Y/msr52/Imputation_required/V2.1_B16_17/lifted_over_merged_b16_17.vcf.gz \
CHAIN=/rds/project/rds-csoP2nj6Y6Y/msr52/Imputation_required/hg38ToHg19.over.chain.gz \
REJECT=/rds/project/rds-csoP2nj6Y6Y/msr52/Imputation_required/V2.1_B16_17/rejected_variants_for_b16_17.vcf \
R=/rds/project/rds-csoP2nj6Y6Y/msr52/Imputation_required/GRCh37.p13.genome.fa


/rds/project/rds-csoP2nj6Y6Y/msr52/SOFTWARES/Apps/bcftools/bcftools query -f '%CHROM\n' merged_b16_17.vcf.gz > chr.txt


> unique(chr$V1)
                  1 1_KI270706v1_random                  10                  11 
              61097                   1               36055               37739 
                 12                  13                  14                  15 
              35518               24050               24956               24893 
                 16                  17   17_KI270909v1_alt                  18 
              28482               27912                   4               20699 
                 19   19_KI270938v1_alt                   2                  20 
              26595                   7               58802               18809 
                 21                  22   22_KI270879v1_alt   22_KI270928v1_alt 
              11129               12976                   4                   1 
                  3                   4                   5                   6 
              49742               44622               42723               53727 
                  7    7_KI270803v1_alt                   8    8_KI270821v1_alt 
              39802                  17               37326                   7 
                  9                  MT                   X                   Y 
              33588                 305               18937                 545 


##contig=<ID=Y,species="Homo sapiens">
##contig=<ID=MT,species="Homo sapiens">
##contig=<ID=17_KI270909v1_alt,species="Homo sapiens">
##contig=<ID=19_KI270938v1_alt,species="Homo sapiens">
##contig=<ID=1_KI270706v1_random,species="Homo sapiens">
##contig=<ID=22_KI270879v1_alt,species="Homo sapiens">
##contig=<ID=22_KI270928v1_alt,species="Homo sapiens">
##contig=<ID=7_KI270803v1_alt,species="Homo sapiens">
##contig=<ID=8_KI270821v1_alt,species="Homo sapiens">
##INFO=<ID=AC,Number=A,Type=Integer,Description="Allele count in genotypes">
##INFO=<ID=AN,Number=1,Type=Integer,Description="Total number of alleles in called genotypes">
##bcftools_viewVersion=1.9+htslib-1.9
##bcftools_viewCommand=view -S /home/cb864/nbr30/v2_chip/b16_eligible1.csv -Oz -o /home/cb864/nbr30/v2_chip/v2_chip_b16_V3_calls_1.vcf.gz /rds/project/rds-CqN4aZ7hKNY/data/batches/UKBBv2_1/V3/b16
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  SP00300107677P  SP00300088656P  SP00300087554R  SP00300081923W  SP00300068747Z  SP00300085899S  SP00300087943P  SP003000709
1       792461  AX-37361813     G       A       .       PASS    CR=97;ConversionType=PolyHighResolution;RSID=rs116587930;PSID=AX-37361813;ASID=Affx-35298040;AC=90;AN=1774      GT      0/0     0/0
1       794252  AX-32137419     C       T       .       PASS    CR=100;ConversionType=NoMinorHom;RSID=rs116720794;PSID=AX-32137419;ASID=Affx-13637449;AC=78;AN=1850     GT      0/0     0/1     0/0
1       817341  AX-13191280     A       G       .       PASS    CR=100;ConversionType=PolyHighResolution;RSID=rs3131972;PSID=AX-13191280;ASID=Affx-13945728;AC=1520;AN=1852     GT      0/1     0/1
1       818725  AX-11194291     C       T       .       PASS    CR=100;ConversionType=NoMinorHom;RSID=rs12184325;PSID=AX-11194291;ASID=Affx-13963217;AC=76;AN=1852      GT      0/0     0/1     0/0
1       821224  AX-32225497     A       G       .       PASS    CR=100;ConversionType=PolyHighResolution;RSID=rs3131962;PSID=AX-32225497;ASID=Affx-13995532;AC=1578;AN=1850     GT      0/1     1/1
1       823656  AX-32233025     G       A       .       PASS    CR=95;ConversionType=NoMinorHom;RSID=rs114525117;PSID=AX-32233025;ASID=Affx-14027812;AC=93;AN=1744      GT      0/0     0/1     0/0
1       825767  AX-40202607     T       C       .       PASS    CR=97;ConversionType=PolyHighResolution;RSID=rs3115850;PSID=AX-40202607;ASID=Affx-14055733;AC=1548;AN=1824      GT      0/1     1/1
1       833068  AX-11214939     G       A       .       PASS    CR=100;ConversionType=PolyHighResolution;RSID=rs12562034;PSID=AX-11214939;ASID=Affx-14143477;AC=189;AN=1852     GT      0/1     0/0
1       837547  AX-32276299     C       T       .       PASS   



