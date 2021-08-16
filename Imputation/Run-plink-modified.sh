#https://rpubs.com/maffleur/452627 this blog helped

plink --bfile /rds/project/rds-csoP2nj6Y6Y/msr52/Imputation_required/v1_chip_extra/v1_chip_extra_merged --exclude /rds/project/rds-csoP2nj6Y6Y/msr52/Imputation_required/v1_chip_extra/Allele_freq_check/Exclude-v1_chip_extra_merged-HRC.txt --make-bed --out /rds/project/rds-csoP2nj6Y6Y/msr52/Imputation_required/v1_chip_extra/Allele_freq_check/TEMP1

# output file was created in another directory so had to change input directory in next line
plink --bfile /rds/project/rds-csoP2nj6Y6Y/msr52/Imputation_required/v1_chip_extra/Allele_freq_check/TEMP1 --update-map /rds/project/rds-csoP2nj6Y6Y/msr52/Imputation_required/v1_chip_extra/Allele_freq_check/Chromosome-v1_chip_extra_merged-HRC.txt --update-chr --make-bed --out /rds/project/rds-csoP2nj6Y6Y/msr52/Imputation_required/v1_chip_extra/Allele_freq_check/TEMP2
# updated input directory in all of these codes
plink --bfile /rds/project/rds-csoP2nj6Y6Y/msr52/Imputation_required/v1_chip_extra/Allele_freq_check/TEMP2 --update-map /rds/project/rds-csoP2nj6Y6Y/msr52/Imputation_required/v1_chip_extra/Allele_freq_check/Position-v1_chip_extra_merged-HRC.txt --make-bed --out /rds/project/rds-csoP2nj6Y6Y/msr52/Imputation_required/v1_chip_extra/Allele_freq_check/TEMP3
plink --bfile /rds/project/rds-csoP2nj6Y6Y/msr52/Imputation_required/v1_chip_extra/Allele_freq_check/TEMP3 --flip /rds/project/rds-csoP2nj6Y6Y/msr52/Imputation_required/v1_chip_extra/Allele_freq_check/Strand-Flip-v1_chip_extra_merged-HRC.txt --make-bed --out /rds/project/rds-csoP2nj6Y6Y/msr52/Imputation_required/v1_chip_extra/Allele_freq_check/TEMP4
#plink --bfile /rds/project/rds-csoP2nj6Y6Y/msr52/Imputation_required/v1_chip_extra/Allele_freq_check/TEMP4 --a2-allele /rds/project/rds-csoP2nj6Y6Y/msr52/Imputation_required/v1_chip_extra/Allele_freq_check/Force-Allele1-v1_chip_extra_merged-HRC.txt --make-bed --out /rds/project/rds-csoP2nj6Y6Y/msr52/Imputation_required/v1_chip_extra/Allele_freq_check/v1_chip_extra_merged-updated

plink --bfile /rds/project/rds-csoP2nj6Y6Y/msr52/Imputation_required/v1_chip_extra/Allele_freq_check/TEMP4 --a2-allele /rds/project/rds-csoP2nj6Y6Y/msr52/Imputation_required/v1_chip_extra/Allele_freq_check/Force-Allele1-v1_chip_extra_merged-HRC.txt --autosome --recode vcf-iid bgz --out /rds/project/rds-csoP2nj6Y6Y/msr52/Imputation_required/v1_chip_extra/Allele_freq_check/v1_chip_extra_merged-updated

rm TEMP*

