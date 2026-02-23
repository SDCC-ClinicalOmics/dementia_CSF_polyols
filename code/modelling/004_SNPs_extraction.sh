## we have verified that in ADNIGO/2 only these 3 SNPs were profiled with 003_SNPs_extraction.py
grep -E "rs2229540|rs7239738|rs2847153" ../data/ExternalCohorts/ADNIGO_2_genetics/ADNI_GO_2_Forward_Bin.bim
grep -E "rs12767543|rs3755137|rs164104|rs743818|rs1050828|rs11231825|rs11265548|rs11417418|rs1213768" ../data/ExternalCohorts/ADNIGO_2_genetics/ADNI_GO_2_Forward_Bin.bim

# 1. Create our SNP list
echo -e "rs2229540\nrs2847153\nrs7239738\nrs10508282\nrs10508287\nrs10508288\nrs11231825" > target_snps.txt

# 2. Extract and recode to additive format (0, 1, 2)
./plink2/plink2.exe \
     --bfile ../data/ExternalCohorts/ADNIGO_2_genetics/ADNI_GO_2_Forward_Bin \
      --extract target_snps.txt \
      --recode A \
      --out ../results/ADNI_genetics/adni_polyol_genotypes