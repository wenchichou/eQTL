#!/bin/bash

PROJECT=chengcancer
echo "PROJECT=${PROJECT}"
#JOBNAME=`(basename $0) |  sed -e 's/\.[^.]*$//'` #tophat_mapping_y
JOBNAME='enrich'
echo "JOBNAME=${JOBNAME}"
NCORE_REQ=2 #based on Charlie's benmark, the increased number of cores above 4 didn't play much in performance gain for tophat,user can experiment this more to set the best value
echo "NCORE_REQ=${NCORE_REQ}"
RT_LIMIT="48:00:00"
MEM_LIMIT="125G" # it seems each job do not consume 10G to my observation; so set 125G shall be safe, even after increase the cores.

#FILE_COUNT=`find /restricted/projectnb/chengcancer/Simon/project/T2D/bulk_RNA-seq/merge_170415_NB501164_170506_NB501164/2mapping/QT/GTeX_bagged_eQTL/* -maxdepth 0 -type d  | wc -l`
#echo "FILE_COUNT=${FILE_COUNT}"

# Let each directory of GTeX bagged  eQTL be processed by various p-value cut-off bins
cd /restricted/projectnb/chengcancer/Simon/project/T2D/bulk_RNA-seq/merge_170415_NB501164_170506_NB501164/2mapping/QT/osteoblast/smallest_pval_data/osteo_bagged
#for GTeX in /restricted/projectnb/chengcancer/Simon/project/T2D/bulk_RNA-seq/merge_170415_NB501164_170506_NB501164/2mapping/QT/GTeX_bagged_eQTL/*/ ;

#do
 # mkdir  $(basename ${file} _Analysis_chrpos_pval.txt.gz)
#GTeX=${GTeX%*/}
GTeX="/restricted/projectnb/chengcancer/Simon/project/T2D/bulk_RNA-seq/merge_170415_NB501164_170506_NB501164/2mapping/QT/osteoblast/smallest_pval_data/osteo_bagged"
echo ${GTeX##*/}
  # The GWAS here is actually bone eQTL 
GWAS="/restricted//projectnb//chengcancer//Simon//project//T2D//bulk_RNA-seq//merge_170415_NB501164_170506_NB501164//2mapping//QT//eQTL"
  # make one enrichment analysis with open cut-off  
GTeX_eQTLpvalueUpperBound=1e-2
GTeX_eQTLpvalueLowerBound=1e-100
sigGWASpvalueCutoff=1e-4
echo ${GTeX_eQTLpvalueUpperBound}
echo ${GTeX_eQTLpvalueLowerBound}

outputFilePath="/restricted/projectnb/chengcancer/Simon/project/T2D/bulk_RNA-seq/merge_170415_NB501164_170506_NB501164/2mapping/QT/osteoblast/smallest_pval_data/"
qsub -P $PROJECT -N $JOBNAME -pe omp $NCORE_REQ -l mem_total=${MEM_LIMIT} -l h_rt=${RT_LIMIT}  -v eQTLfilePath=${GTeX} -v GWASfilePath=${GWAS} -v GTeX_eQTLpvalueUpperBound=${GTeX_eQTLpvalueUpperBound}  -v GTeX_eQTLpvalueLowerBound=${GTeX_eQTLpvalueLowerBound}  -v sigGWASpvalueCutoff=${sigGWASpvalueCutoff}  -v outdir=${outputFilePath}  /restricted/projectnb/chengcancer/Simon/project/T2D/bulk_RNA-seq/merge_170415_NB501164_170506_NB501164/2mapping/QT/osteoblast/smallest_pval_data/osteo_enrichment.qsub


for power  in {2..10..2};
do

   pow2=$[${power}+2]
   GTeX_eQTLpvalueUpperBound=1e-${power}
   GTeX_eQTLpvalueLowerBound=1e-${pow2}
   echo ${GTeX_eQTLpvalueUpperBound}
   echo ${GTeX_eQTLpvalueLowerBound}
   qsub -P $PROJECT -N $JOBNAME -pe omp $NCORE_REQ -l mem_total=${MEM_LIMIT} -l h_rt=${RT_LIMIT}  -v eQTLfilePath=${GTeX} -v GWASfilePath=${GWAS} -v GTeX_eQTLpvalueUpperBound=${GTeX_eQTLpvalueUpperBound}  -v GTeX_eQTLpvalueLowerBound=${GTeX_eQTLpvalueLowerBound}  -v sigGWASpvalueCutoff=${sigGWASpvalueCutoff}  -v outdir=${outputFilePath}  /restricted/projectnb/chengcancer/Simon/project/T2D/bulk_RNA-seq/merge_170415_NB501164_170506_NB501164/2mapping/QT/osteoblast/smallest_pval_data/osteo_enrichment.qsub
done




