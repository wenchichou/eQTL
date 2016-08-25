#!/usr/bin/python
#CHR=sys.argv[0]
##args <- commandArgs(TRUE)
##eQTLPvaluePath <- args[1]
##additiveSumDir <- args[2]
##rDataPath <- args[3]
##r2cutoff <- as.numeric(args[4])
#
import sys
import csv
import re
import glob
import gzip
#import numpy as np
#import pandas as pd
import getopt

CHR = str(sys.argv[1])
rsquaredCutOff = float(sys.argv[2])
printOutUnit = int(sys.argv[3])
#def main(argv):
#   CHR = ''
#   rsquaredCutOff = ''
#   printOutUnit = ''
#   try:
#      opts, args = getopt.getopt(argv,"hc:r:p:",["chromosome=","rsquaredCutOff=","printOutUnit="])
#   except getopt.GetoptError:
#      print 'baggingByLD.print.arg.py -c <chromosome> -r <rsquaredCutOff> -p <printOutUnit>'
#      sys.exit(2)
#   for opt, arg in opts:
#      if opt == '-h':
#         print 'baggingByLD.print.arg.py -c <chromosome> -r <rsquaredCutOff> -p <printOutUnit>'
#         sys.exit()
#      elif opt in ("-c", "--chromosome"):
#         CHR = arg
#      elif opt in ("-r", "--rsquaredCutOff"):
#         rsquaredCutOff = arg
#      elif opt in ("-p", "--printOutUnit"):
#         printOutUnit = arg
#   print 'Running Chromosome ', CHR
#   print 'with rsquaredCutOff ', rsquaredCutOff
#   print 'with printOutUnit ', printOutUnit


## 1. load eQTL, oneFile and twofile merged rsquare table data
eqtlPvalue_hash = {}
with open('/Volumes/Seagate3TB/orchestraBackup.081716/gautvik.eQTL/results/111314/filesForeQTL.043015/completedExtraction/completedExtraction2/eQTL.cis/smallestPvalue.sorted.k1.folder/CHR%s.cis.gz.smallestPvalue.sorted.k1' % CHR) as f:
    for line in f:
	(key, val, tmp) = line.split()
	key = re.sub(r"^.*:", "", key)
	eqtlPvalue_hash[str(key)] = float(val)

#smallestValue = min(eqtlPvalue_hash.values())
#keyWithSmallestValue = min(eqtlPvalue_hash, key=eqtlPvalue_hash.get)

#for x in eqtlPvalue_hash:
#    print (x)
#    print (eqtlPvalue_hash[x])

#rsquared_hash = {}
position1 = []
position2 = []
pairwiseR2oneFiles = glob.glob('/Volumes/Seagate3TB/orchestraBackup.081716/gautvik.eQTL/results/111314/pairwiseR2/pairwiseR2.oneFile/chr%s.*.impute2.r2.gz.gt0.4.gz' % CHR)
#print pairwiseR2oneFiles
for file in pairwiseR2oneFiles:
	print file
	with gzip.open(file,'r') as fin:
		for line in fin:
			(pos1, pos2, r2) = line.split()
			if r2 > rsquaredCutOff:
				#positions = (pos1, pos2)
				#rsquared_hash[" ".join(positions)] = float(r2)
				position1.append(pos1)
				position2.append(pos2)

pairwiseR2twoFiles = glob.glob('/Volumes/Seagate3TB/orchestraBackup.081716/gautvik.eQTL/results/111314/pairwiseR2/pairwiseR2.twoFiles/chr%s/chr%s.*.chr%s.*.r2.gz' % (CHR, CHR, CHR))
#print pairwiseR2oneFiles
for file in pairwiseR2twoFiles:
	print file
	with gzip.open(file,'r') as fin:
		for line in fin:
			(pos1, pos2, r2) = line.split()
			if r2 > 0.8:
				#positions = (pos1, pos2)
				#rsquared_hash[" ".join(positions)] = float(r2)
				position1.append(pos1)
				position2.append(pos2)
#	iter_csv = pd.read_csv(file, iterator=True, chunksize=1000, compression='gzip',  header=0, sep=' ')
#	df_one = pd.concat([chunk[chunk.ix[:,2] > 0.8] for chunk in iter_csv])
#	df_all = pd.concat([df_all, df_one])

#print "eQTL Length : %d" % len (eqtlPvalue_hash)
#print "position1 Length : %d" % len (position1)
#print "position2 Length : %d" % len (position2)
#print df_all.shape

#for x in eqtlPvalue_hash:
#    print (x)
#    print (eqtlPvalue_hash[x])


#bag = [[] for i in range(len(eqtlPvalue_hash))]
bag = [[] for i in range(0)]
pvalue = []
bagi = 0
unitBagi = 0
while len(eqtlPvalue_hash) > 1:
#while bagi < 10:
    if bagi >= len(bag):
            bag.append([])
    #bag[bagi] = []
    rsquaredPositions = []
    smallestValue = min(eqtlPvalue_hash.values())
    positionWithSmallestValue = min(eqtlPvalue_hash, key=eqtlPvalue_hash.get) # the position with smallest eQTL pvalue
    print ("\nbagi: %s" % (bagi+unitBagi*(printOutUnit+1)))
    print "eQTLpvalue Length : ", len(eqtlPvalue_hash)
    print "position1 Length  : ", len(position1)
    print "position2 Length  : ", len(position2)
    print ("pos: %s; pvalue: %s" % (positionWithSmallestValue, smallestValue))
    indices_position1 = [i for i, x in enumerate(position1) if x == positionWithSmallestValue]
    indices_position2 = [i for i, x in enumerate(position2) if x == positionWithSmallestValue]
    if len(indices_position1) > 0:
        rsquaredPositions.extend([ position2[index] for index in indices_position1 ])
#        print [ position1[index] for index in indices_position1 ]
#        print [ position2[index] for index in indices_position1 ]
    if len(indices_position2) > 0:
        rsquaredPositions.extend([ position1[index] for index in indices_position2 ])
#        print [ position2[index] for index in indices_position2 ]
#        print [ position1[index] for index in indices_position2 ]
#    rsquaredPositions.append(positionWithSmallestValue)
    rsquaredPositions.insert(0, positionWithSmallestValue)
    pvalue.append(smallestValue)
    allIndices = indices_position1 + indices_position2
    for index in sorted(allIndices, reverse=True):
        del position1[index]
        del position2[index]
    for k in rsquaredPositions:
        if k in eqtlPvalue_hash:
            bag[bagi].append(k)
            eqtlPvalue_hash.pop(k, None)
    print "Length of rsquaredPositions: ", len(rsquaredPositions)
    print "Length of this bag: ", len(bag[bagi])
    if (bagi == printOutUnit):
        bagOutFile = open('./bagSNP.chr%s.print.txt' % CHR, 'a')
        for item in bag:
            bagOutFile.write("%s\n" % ' '.join(item))
        pvalueOutFile = open('./bagPvalue.chr%s.print.txt' % CHR, 'a')
        for item in pvalue:
            pvalueOutFile.write("%s\n" % item)
        bag = [[] for i in range(0)]
        pvalue = []
        bagi = 0
        unitBagi = unitBagi + 1
    else:
        bagi = bagi + 1

if (bagi < printOutUnit) and (bagi > 0):
    bagOutFile = open('./bagSNP.chr%s.print.txt' % CHR, 'a')
    for item in bag:
        bagOutFile.write("%s\n" % ' '.join(item))
    pvalueOutFile = open('./bagPvalue.chr%s.print.txt' % CHR, 'a')
    for item in pvalue:
        pvalueOutFile.write("%s\n" % item)

