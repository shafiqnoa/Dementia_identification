plink --bfile /rds/project/rds-csoP2nj6Y6Y/msr52/Imputation_required/v1_chip_extra/v1_chip_extra_merged --exclude /rds/project/rds-csoP2nj6Y6Y/msr52/Imputation_required/v1_chip_extra/Allele_freq_check/Exclude-v1_chip_extra_merged-HRC.txt --make-bed --out /rds/project/rds-csoP2nj6Y6Y/msr52/Imputation_required/v1_chip_extra/TEMP1
plink --bfile /rds/project/rds-csoP2nj6Y6Y/msr52/Imputation_required/v1_chip_extra/TEMP1 --update-map /rds/project/rds-csoP2nj6Y6Y/msr52/Imputation_required/v1_chip_extra/Allele_freq_check/Chromosome-v1_chip_extra_merged-HRC.txt --update-chr --make-bed --out /rds/project/rds-csoP2nj6Y6Y/msr52/Imputation_required/v1_chip_extra/TEMP2
plink --bfile /rds/project/rds-csoP2nj6Y6Y/msr52/Imputation_required/v1_chip_extra/TEMP2 --update-map /rds/project/rds-csoP2nj6Y6Y/msr52/Imputation_required/v1_chip_extra/Allele_freq_check/Position-v1_chip_extra_merged-HRC.txt --make-bed --out /rds/project/rds-csoP2nj6Y6Y/msr52/Imputation_required/v1_chip_extra/TEMP3
plink --bfile /rds/project/rds-csoP2nj6Y6Y/msr52/Imputation_required/v1_chip_extra/TEMP3 --flip /rds/project/rds-csoP2nj6Y6Y/msr52/Imputation_required/v1_chip_extra/Allele_freq_check/Strand-Flip-v1_chip_extra_merged-HRC.txt --make-bed --out /rds/project/rds-csoP2nj6Y6Y/msr52/Imputation_required/v1_chip_extra/TEMP4
plink --bfile /rds/project/rds-csoP2nj6Y6Y/msr52/Imputation_required/v1_chip_extra/TEMP4 --a2-allele /rds/project/rds-csoP2nj6Y6Y/msr52/Imputation_required/v1_chip_extra/Allele_freq_check/Force-Allele1-v1_chip_extra_merged-HRC.txt --make-bed --out /rds/project/rds-csoP2nj6Y6Y/msr52/Imputation_required/v1_chip_extra/Allele_freq_check/v1_chip_extra_merged-updated
plink --bfile /rds/project/rds-csoP2nj6Y6Y/msr52/Imputation_required/v1_chip_extra/Allele_freq_check/v1_chip_extra_merged-updated --real-ref-alleles --make-bed --chr 1 --out /rds/project/rds-csoP2nj6Y6Y/msr52/Imputation_required/v1_chip_extra/Allele_freq_check/v1_chip_extra_merged-updated-chr1
plink --bfile /rds/project/rds-csoP2nj6Y6Y/msr52/Imputation_required/v1_chip_extra/Allele_freq_check/v1_chip_extra_merged-updated --real-ref-alleles --recode vcf --chr 1 --out /rds/project/rds-csoP2nj6Y6Y/msr52/Imputation_required/v1_chip_extra/Allele_freq_check/v1_chip_extra_merged-updated-chr1
plink --bfile /rds/project/rds-csoP2nj6Y6Y/msr52/Imputation_required/v1_chip_extra/Allele_freq_check/v1_chip_extra_merged-updated --real-ref-alleles --make-bed --chr 2 --out /rds/project/rds-csoP2nj6Y6Y/msr52/Imputation_required/v1_chip_extra/Allele_freq_check/v1_chip_extra_merged-updated-chr2
plink --bfile /rds/project/rds-csoP2nj6Y6Y/msr52/Imputation_required/v1_chip_extra/Allele_freq_check/v1_chip_extra_merged-updated --real-ref-alleles --recode vcf --chr 2 --out /rds/project/rds-csoP2nj6Y6Y/msr52/Imputation_required/v1_chip_extra/Allele_freq_check/v1_chip_extra_merged-updated-chr2
plink --bfile /rds/project/rds-csoP2nj6Y6Y/msr52/Imputation_required/v1_chip_extra/Allele_freq_check/v1_chip_extra_merged-updated --real-ref-alleles --make-bed --chr 3 --out /rds/project/rds-csoP2nj6Y6Y/msr52/Imputation_required/v1_chip_extra/Allele_freq_check/v1_chip_extra_merged-updated-chr3
plink --bfile /rds/project/rds-csoP2nj6Y6Y/msr52/Imputation_required/v1_chip_extra/Allele_freq_check/v1_chip_extra_merged-updated --real-ref-alleles --recode vcf --chr 3 --out /rds/project/rds-csoP2nj6Y6Y/msr52/Imputation_required/v1_chip_extra/Allele_freq_check/v1_chip_extra_merged-updated-chr3
plink --bfile /rds/project/rds-csoP2nj6Y6Y/msr52/Imputation_required/v1_chip_extra/Allele_freq_check/v1_chip_extra_merged-updated --real-ref-alleles --make-bed --chr 4 --out /rds/project/rds-csoP2nj6Y6Y/msr52/Imputation_required/v1_chip_extra/Allele_freq_check/v1_chip_extra_merged-updated-chr4
plink --bfile /rds/project/rds-csoP2nj6Y6Y/msr52/Imputation_required/v1_chip_extra/Allele_freq_check/v1_chip_extra_merged-updated --real-ref-alleles --recode vcf --chr 4 --out /rds/project/rds-csoP2nj6Y6Y/msr52/Imputation_required/v1_chip_extra/Allele_freq_check/v1_chip_extra_merged-updated-chr4
plink --bfile /rds/project/rds-csoP2nj6Y6Y/msr52/Imputation_required/v1_chip_extra/Allele_freq_check/v1_chip_extra_merged-updated --real-ref-alleles --make-bed --chr 5 --out /rds/project/rds-csoP2nj6Y6Y/msr52/Imputation_required/v1_chip_extra/Allele_freq_check/v1_chip_extra_merged-updated-chr5
plink --bfile /rds/project/rds-csoP2nj6Y6Y/msr52/Imputation_required/v1_chip_extra/Allele_freq_check/v1_chip_extra_merged-updated --real-ref-alleles --recode vcf --chr 5 --out /rds/project/rds-csoP2nj6Y6Y/msr52/Imputation_required/v1_chip_extra/Allele_freq_check/v1_chip_extra_merged-updated-chr5
plink --bfile /rds/project/rds-csoP2nj6Y6Y/msr52/Imputation_required/v1_chip_extra/Allele_freq_check/v1_chip_extra_merged-updated --real-ref-alleles --make-bed --chr 6 --out /rds/project/rds-csoP2nj6Y6Y/msr52/Imputation_required/v1_chip_extra/Allele_freq_check/v1_chip_extra_merged-updated-chr6
plink --bfile /rds/project/rds-csoP2nj6Y6Y/msr52/Imputation_required/v1_chip_extra/Allele_freq_check/v1_chip_extra_merged-updated --real-ref-alleles --recode vcf --chr 6 --out /rds/project/rds-csoP2nj6Y6Y/msr52/Imputation_required/v1_chip_extra/Allele_freq_check/v1_chip_extra_merged-updated-chr6
plink --bfile /rds/project/rds-csoP2nj6Y6Y/msr52/Imputation_required/v1_chip_extra/Allele_freq_check/v1_chip_extra_merged-updated --real-ref-alleles --make-bed --chr 7 --out /rds/project/rds-csoP2nj6Y6Y/msr52/Imputation_required/v1_chip_extra/Allele_freq_check/v1_chip_extra_merged-updated-chr7
plink --bfile /rds/project/rds-csoP2nj6Y6Y/msr52/Imputation_required/v1_chip_extra/Allele_freq_check/v1_chip_extra_merged-updated --real-ref-alleles --recode vcf --chr 7 --out /rds/project/rds-csoP2nj6Y6Y/msr52/Imputation_required/v1_chip_extra/Allele_freq_check/v1_chip_extra_merged-updated-chr7
plink --bfile /rds/project/rds-csoP2nj6Y6Y/msr52/Imputation_required/v1_chip_extra/Allele_freq_check/v1_chip_extra_merged-updated --real-ref-alleles --make-bed --chr 8 --out /rds/project/rds-csoP2nj6Y6Y/msr52/Imputation_required/v1_chip_extra/Allele_freq_check/v1_chip_extra_merged-updated-chr8
plink --bfile /rds/project/rds-csoP2nj6Y6Y/msr52/Imputation_required/v1_chip_extra/Allele_freq_check/v1_chip_extra_merged-updated --real-ref-alleles --recode vcf --chr 8 --out /rds/project/rds-csoP2nj6Y6Y/msr52/Imputation_required/v1_chip_extra/Allele_freq_check/v1_chip_extra_merged-updated-chr8
plink --bfile /rds/project/rds-csoP2nj6Y6Y/msr52/Imputation_required/v1_chip_extra/Allele_freq_check/v1_chip_extra_merged-updated --real-ref-alleles --make-bed --chr 9 --out /rds/project/rds-csoP2nj6Y6Y/msr52/Imputation_required/v1_chip_extra/Allele_freq_check/v1_chip_extra_merged-updated-chr9
plink --bfile /rds/project/rds-csoP2nj6Y6Y/msr52/Imputation_required/v1_chip_extra/Allele_freq_check/v1_chip_extra_merged-updated --real-ref-alleles --recode vcf --chr 9 --out /rds/project/rds-csoP2nj6Y6Y/msr52/Imputation_required/v1_chip_extra/Allele_freq_check/v1_chip_extra_merged-updated-chr9
plink --bfile /rds/project/rds-csoP2nj6Y6Y/msr52/Imputation_required/v1_chip_extra/Allele_freq_check/v1_chip_extra_merged-updated --real-ref-alleles --make-bed --chr 10 --out /rds/project/rds-csoP2nj6Y6Y/msr52/Imputation_required/v1_chip_extra/Allele_freq_check/v1_chip_extra_merged-updated-chr10
plink --bfile /rds/project/rds-csoP2nj6Y6Y/msr52/Imputation_required/v1_chip_extra/Allele_freq_check/v1_chip_extra_merged-updated --real-ref-alleles --recode vcf --chr 10 --out /rds/project/rds-csoP2nj6Y6Y/msr52/Imputation_required/v1_chip_extra/Allele_freq_check/v1_chip_extra_merged-updated-chr10
plink --bfile /rds/project/rds-csoP2nj6Y6Y/msr52/Imputation_required/v1_chip_extra/Allele_freq_check/v1_chip_extra_merged-updated --real-ref-alleles --make-bed --chr 11 --out /rds/project/rds-csoP2nj6Y6Y/msr52/Imputation_required/v1_chip_extra/Allele_freq_check/v1_chip_extra_merged-updated-chr11
plink --bfile /rds/project/rds-csoP2nj6Y6Y/msr52/Imputation_required/v1_chip_extra/Allele_freq_check/v1_chip_extra_merged-updated --real-ref-alleles --recode vcf --chr 11 --out /rds/project/rds-csoP2nj6Y6Y/msr52/Imputation_required/v1_chip_extra/Allele_freq_check/v1_chip_extra_merged-updated-chr11
plink --bfile /rds/project/rds-csoP2nj6Y6Y/msr52/Imputation_required/v1_chip_extra/Allele_freq_check/v1_chip_extra_merged-updated --real-ref-alleles --make-bed --chr 12 --out /rds/project/rds-csoP2nj6Y6Y/msr52/Imputation_required/v1_chip_extra/Allele_freq_check/v1_chip_extra_merged-updated-chr12
plink --bfile /rds/project/rds-csoP2nj6Y6Y/msr52/Imputation_required/v1_chip_extra/Allele_freq_check/v1_chip_extra_merged-updated --real-ref-alleles --recode vcf --chr 12 --out /rds/project/rds-csoP2nj6Y6Y/msr52/Imputation_required/v1_chip_extra/Allele_freq_check/v1_chip_extra_merged-updated-chr12
plink --bfile /rds/project/rds-csoP2nj6Y6Y/msr52/Imputation_required/v1_chip_extra/Allele_freq_check/v1_chip_extra_merged-updated --real-ref-alleles --make-bed --chr 13 --out /rds/project/rds-csoP2nj6Y6Y/msr52/Imputation_required/v1_chip_extra/Allele_freq_check/v1_chip_extra_merged-updated-chr13
plink --bfile /rds/project/rds-csoP2nj6Y6Y/msr52/Imputation_required/v1_chip_extra/Allele_freq_check/v1_chip_extra_merged-updated --real-ref-alleles --recode vcf --chr 13 --out /rds/project/rds-csoP2nj6Y6Y/msr52/Imputation_required/v1_chip_extra/Allele_freq_check/v1_chip_extra_merged-updated-chr13
plink --bfile /rds/project/rds-csoP2nj6Y6Y/msr52/Imputation_required/v1_chip_extra/Allele_freq_check/v1_chip_extra_merged-updated --real-ref-alleles --make-bed --chr 14 --out /rds/project/rds-csoP2nj6Y6Y/msr52/Imputation_required/v1_chip_extra/Allele_freq_check/v1_chip_extra_merged-updated-chr14
plink --bfile /rds/project/rds-csoP2nj6Y6Y/msr52/Imputation_required/v1_chip_extra/Allele_freq_check/v1_chip_extra_merged-updated --real-ref-alleles --recode vcf --chr 14 --out /rds/project/rds-csoP2nj6Y6Y/msr52/Imputation_required/v1_chip_extra/Allele_freq_check/v1_chip_extra_merged-updated-chr14
plink --bfile /rds/project/rds-csoP2nj6Y6Y/msr52/Imputation_required/v1_chip_extra/Allele_freq_check/v1_chip_extra_merged-updated --real-ref-alleles --make-bed --chr 15 --out /rds/project/rds-csoP2nj6Y6Y/msr52/Imputation_required/v1_chip_extra/Allele_freq_check/v1_chip_extra_merged-updated-chr15
plink --bfile /rds/project/rds-csoP2nj6Y6Y/msr52/Imputation_required/v1_chip_extra/Allele_freq_check/v1_chip_extra_merged-updated --real-ref-alleles --recode vcf --chr 15 --out /rds/project/rds-csoP2nj6Y6Y/msr52/Imputation_required/v1_chip_extra/Allele_freq_check/v1_chip_extra_merged-updated-chr15
plink --bfile /rds/project/rds-csoP2nj6Y6Y/msr52/Imputation_required/v1_chip_extra/Allele_freq_check/v1_chip_extra_merged-updated --real-ref-alleles --make-bed --chr 16 --out /rds/project/rds-csoP2nj6Y6Y/msr52/Imputation_required/v1_chip_extra/Allele_freq_check/v1_chip_extra_merged-updated-chr16
plink --bfile /rds/project/rds-csoP2nj6Y6Y/msr52/Imputation_required/v1_chip_extra/Allele_freq_check/v1_chip_extra_merged-updated --real-ref-alleles --recode vcf --chr 16 --out /rds/project/rds-csoP2nj6Y6Y/msr52/Imputation_required/v1_chip_extra/Allele_freq_check/v1_chip_extra_merged-updated-chr16
plink --bfile /rds/project/rds-csoP2nj6Y6Y/msr52/Imputation_required/v1_chip_extra/Allele_freq_check/v1_chip_extra_merged-updated --real-ref-alleles --make-bed --chr 17 --out /rds/project/rds-csoP2nj6Y6Y/msr52/Imputation_required/v1_chip_extra/Allele_freq_check/v1_chip_extra_merged-updated-chr17
plink --bfile /rds/project/rds-csoP2nj6Y6Y/msr52/Imputation_required/v1_chip_extra/Allele_freq_check/v1_chip_extra_merged-updated --real-ref-alleles --recode vcf --chr 17 --out /rds/project/rds-csoP2nj6Y6Y/msr52/Imputation_required/v1_chip_extra/Allele_freq_check/v1_chip_extra_merged-updated-chr17
plink --bfile /rds/project/rds-csoP2nj6Y6Y/msr52/Imputation_required/v1_chip_extra/Allele_freq_check/v1_chip_extra_merged-updated --real-ref-alleles --make-bed --chr 18 --out /rds/project/rds-csoP2nj6Y6Y/msr52/Imputation_required/v1_chip_extra/Allele_freq_check/v1_chip_extra_merged-updated-chr18
plink --bfile /rds/project/rds-csoP2nj6Y6Y/msr52/Imputation_required/v1_chip_extra/Allele_freq_check/v1_chip_extra_merged-updated --real-ref-alleles --recode vcf --chr 18 --out /rds/project/rds-csoP2nj6Y6Y/msr52/Imputation_required/v1_chip_extra/Allele_freq_check/v1_chip_extra_merged-updated-chr18
plink --bfile /rds/project/rds-csoP2nj6Y6Y/msr52/Imputation_required/v1_chip_extra/Allele_freq_check/v1_chip_extra_merged-updated --real-ref-alleles --make-bed --chr 19 --out /rds/project/rds-csoP2nj6Y6Y/msr52/Imputation_required/v1_chip_extra/Allele_freq_check/v1_chip_extra_merged-updated-chr19
plink --bfile /rds/project/rds-csoP2nj6Y6Y/msr52/Imputation_required/v1_chip_extra/Allele_freq_check/v1_chip_extra_merged-updated --real-ref-alleles --recode vcf --chr 19 --out /rds/project/rds-csoP2nj6Y6Y/msr52/Imputation_required/v1_chip_extra/Allele_freq_check/v1_chip_extra_merged-updated-chr19
plink --bfile /rds/project/rds-csoP2nj6Y6Y/msr52/Imputation_required/v1_chip_extra/Allele_freq_check/v1_chip_extra_merged-updated --real-ref-alleles --make-bed --chr 20 --out /rds/project/rds-csoP2nj6Y6Y/msr52/Imputation_required/v1_chip_extra/Allele_freq_check/v1_chip_extra_merged-updated-chr20
plink --bfile /rds/project/rds-csoP2nj6Y6Y/msr52/Imputation_required/v1_chip_extra/Allele_freq_check/v1_chip_extra_merged-updated --real-ref-alleles --recode vcf --chr 20 --out /rds/project/rds-csoP2nj6Y6Y/msr52/Imputation_required/v1_chip_extra/Allele_freq_check/v1_chip_extra_merged-updated-chr20
plink --bfile /rds/project/rds-csoP2nj6Y6Y/msr52/Imputation_required/v1_chip_extra/Allele_freq_check/v1_chip_extra_merged-updated --real-ref-alleles --make-bed --chr 21 --out /rds/project/rds-csoP2nj6Y6Y/msr52/Imputation_required/v1_chip_extra/Allele_freq_check/v1_chip_extra_merged-updated-chr21
plink --bfile /rds/project/rds-csoP2nj6Y6Y/msr52/Imputation_required/v1_chip_extra/Allele_freq_check/v1_chip_extra_merged-updated --real-ref-alleles --recode vcf --chr 21 --out /rds/project/rds-csoP2nj6Y6Y/msr52/Imputation_required/v1_chip_extra/Allele_freq_check/v1_chip_extra_merged-updated-chr21
plink --bfile /rds/project/rds-csoP2nj6Y6Y/msr52/Imputation_required/v1_chip_extra/Allele_freq_check/v1_chip_extra_merged-updated --real-ref-alleles --make-bed --chr 22 --out /rds/project/rds-csoP2nj6Y6Y/msr52/Imputation_required/v1_chip_extra/Allele_freq_check/v1_chip_extra_merged-updated-chr22
plink --bfile /rds/project/rds-csoP2nj6Y6Y/msr52/Imputation_required/v1_chip_extra/Allele_freq_check/v1_chip_extra_merged-updated --real-ref-alleles --recode vcf --chr 22 --out /rds/project/rds-csoP2nj6Y6Y/msr52/Imputation_required/v1_chip_extra/Allele_freq_check/v1_chip_extra_merged-updated-chr22
rm TEMP*
