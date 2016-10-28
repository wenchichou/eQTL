# 102716
cd /home/unix/wcchou/gsapWenChi/gautvik/results

CHR=1; cat bagSNP.chr${CHR}.print.txt| awk -F " " -v CHR=${CHR} '{for(i = 1; i < NF; i++){printf "%s:%s ",CHR,$i};printf "%s:%s",CHR,$NF; printf "\n"}' > bagSNP.chr${CHR}.print.txt2
