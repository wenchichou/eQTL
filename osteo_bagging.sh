#!/bin/bash

PROJECT=chengcancer
echo "PROJECT=${PROJECT}"
#JOBNAME=`(basename $0) |  sed -e 's/\.[^.]*$//'` #tophat_mapping_y
JOBNAME='bag'
echo "JOBNAME=${JOBNAME}"
NCORE_REQ=2 #based on Charlie's benmark, the increased number of cores above 4 didn't play much in performance gain for tophat,user can experiment this more to set the best value
echo "NCORE_REQ=${NCORE_REQ}"
RT_LIMIT="48:00:00"
MEM_LIMIT="125G" # it seems each job do not consume 10G to my observation; so set 125G shall be safe, even after increase the cores.

#FILE_COUNT=`find /restricted/projectnb/chengcancer/Simon/project/T2D/bulk_RNA-seq/merge_170415_NB501164_170506_NB501164/2mapping/QT/osteoblast/smallest_pval_data/*.cis_osteoblast_eQTL_position_pval.txt.gz  | wc -l`
#echo "FILE_COUNT=${FILE_COUNT}"


#for file in /restricted/projectnb/chengcancer/Simon/project/T2D/bulk_RNA-seq/merge_170415_NB501164_170506_NB501164/2mapping/QT/osteoblast/smallest_pval_data/*.cis_osteoblast_eQTL_position_pval.txt.gz ;
#do
 # mkdir  $(basename ${file} _Analysis_chrpos_pval.txt.gz)
for CHR in {1..22}; 
do
   file="/restricted/projectnb/chengcancer/Simon/project/T2D/bulk_RNA-seq/merge_170415_NB501164_170506_NB501164/2mapping/QT/osteoblast/smallest_pval_data/chr${CHR}.cis_osteoblast_eQTL_position_pval.txt.gz"
   echo $CHR
   echo $file
   qsub -P $PROJECT -N $JOBNAME -pe omp $NCORE_REQ -l mem_total=${MEM_LIMIT} -l h_rt=${RT_LIMIT}  -v CHR=${CHR} -v file=${file}  -v outdir=$(basename ${file} _eQTL_position_pval.txt.gz) osteo_bagging.qsub  
done

#done
