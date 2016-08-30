use UGER
qrsh -q interactive -l m_mem_free=10g
cd /home/unix/wcchou/gsapWenChi/gautvik/results
use .python-2.7.8-sqlite3-rtrees
for CHR in {21..1}; do echo $CHR;
python ../bin/eQTL/baggingSNPsByLD.py $CHR 0.8 199 ../data/pairwiseR2/ ../data/smallestPvalue.sorted.k1.folder/;
done

