#use UGER
#qrsh -q interactive -l m_mem_free=10g
#use .python-2.7.8-sqlite3-rtrees
#for CHR in {21..1}; do echo $CHR;
#python ../bin/eQTL/baggingSNPsByLD.py $CHR 0.8 199 ../data/pairwiseR2/ ../data/smallestPvalue.sorted.k1.folder/;
#done

cd /home/unix/wcchou/gsapWenChi/gautvik/results
for CHR in {21..1}; do echo "#! /bin/bash
#$ -cwd
#$ -q long
#$ -P gscid 
#$ -pe smp 2 -R y
#$ -l m_mem_free=10g
#$ -e bagging.${CHR}.sh.err
#$ -o bagging.${CHR}.sh.out
source /broad/software/scripts/useuse;
reuse -q .python-2.7.8-sqlite3-rtrees
python ../bin/eQTL/baggingSNPsByLD.py $CHR 0.8 199 ../data/pairwiseR2/ ../data/smallestPvalue.sorted.k1.folder/
"> bagging.${CHR}.sh
done

for shFile in `find -maxdepth 1 -name "bagging.*.sh"`; do echo $shFile; chmod u+x $shFile; qsub $shFile; done

watch -n5 qstat



