#https://github.com/wenchichou/eQTL
# git add *
# git commit -m "updates"
# git push

# 102716
cd /home/unix/wcchou/gsapWenChi/gautvik/results

for CHR in {1..22}; do 
	echo $CHR; 
	cat bagSNP.chr${CHR}.print.txt| awk -F " " -v CHR=${CHR} '{for(i = 1; i < NF; i++){printf "%s:%s ",CHR,$i};printf "%s:%s",CHR,$NF; printf "\n"}' > bagSNP.chr${CHR}.print.txt2
done 

reuse UGER
#ish
qrsh -q interactive -l h_vmem=32g
cd /home/unix/wcchou/gsapWenChi/gautvik/results
use R-3.3

R

Rscript hypergeometric_test_after_bagging.R "/home/unix/wcchou/gsapWenChi/gautvik/results/bagging/" "/home/unix/wcchou/gsapWenChi/gautvik/data/gwas.2.5m/LSBMD_CHRPOS_PVALUE_done.gz" 1e-8 1e-6 1e-1000 test110216.2.txt 





