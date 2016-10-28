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
ish
cd /home/unix/wcchou/gsapWenChi/gautvik/results
use R-3.3

R

