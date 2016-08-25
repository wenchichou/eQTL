#!/usr/bin/python
import sys
import csv
import re
import glob
import gzip
import getopt
import os
#use .python-2.7.8-sqlite3-rtrees
#python ../bin/eQTL/baggingSNPsByLD.py 22 0.8 199 ../data/pairwiseR2/ ../data/smallestPvalue.sorted.k1.folder/

CHR = str(sys.argv[1])
rsquaredCutOff = float(sys.argv[2])
printOutUnit = int(sys.argv[3])
allRsquaredFolder = str(sys.argv[4]) #/home/unix/wcchou/gsapWenChi/gautvik/data/pairwiseR2
allSmallesteQTLPvalueFolder = str(sys.argv[5]) #/home/unix/wcchou/gsapWenChi/gautvik/data/smallestPvalue.sorted.k1.folder
#CHR = ''
#rsquaredCutOff = ''
#printOutUnit = ''
#allRsquaredFolder = ''
#allSmallesteQTLPvalueFolder = ''
#try:
#    opts, args = getopt.getopt(sys.argv,"hc:r:p:dr:de:",["chromosome=","rsquaredCutOff=","printOutUnit=","allRsquaredFolder=","allSmallesteQTLPvalueFolder="])
#except getopt.GetoptError:
#    print 'baggingByLD.print.arg.py -c <chromosome> -r <rsquaredCutOff> -p <printOutUnit> -dr <allRsquaredFolder> -de <allSmallesteQTLPvalueFolder>'
#    sys.exit(2)
#for opt, arg in opts:
#    if opt == '-h':
#        print 'baggingByLD.print.arg.py -c <chromosome> -r <rsquaredCutOff> -p <printOutUnit> -dr <allRsquaredFolder> -de <allSmallesteQTLPvalueFolder>'
#        sys.exit()
#    elif opt in ("-c", "--chromosome"):
#        CHR = arg
#    elif opt in ("-r", "--rsquaredCutOff"):
#        rsquaredCutOff = arg
#    elif opt in ("-p", "--printOutUnit"):
#        printOutUnit = arg
#    elif opt in ("-dr", "--allRsquaredFolder"):
#        allRsquaredFolder = arg
#    elif opt in ("-de", "--allSmallesteQTLPvalueFolder"):
#        allSmallesteQTLPvalueFolder = arg
print 'Running Chromosome ', CHR
print 'with rsquaredCutOff ', rsquaredCutOff
print 'with printOutUnit ', printOutUnit
print 'using R2 files at', allRsquaredFolder
print 'using eQTL pvalue files at', allSmallesteQTLPvalueFolder


## 1. load eQTL, oneFile and twofile merged rsquare table data
eqtlPvalue_hash = {}
eQTLFileName = 'CHR%s.cis.gz.smallestPvalue.sorted.k1' % CHR
with open(os.path.join(allSmallesteQTLPvalueFolder, eQTLFileName)) as f:
    for line in f:
	(key, val, tmp) = line.split()
	key = re.sub(r"^.*:", "", key)
	eqtlPvalue_hash[str(key)] = float(val)

## 2. load all R2 files of the given chromosome
position1 = []
position2 = []
pfiles_1 = glob.glob('%s/*/*chr%s.*.gz' % (allRsquaredFolder, CHR))
pfiles_2 = glob.glob('%s/*/chr%s/*chr%s.*.gz' % (allRsquaredFolder, CHR, CHR))
pairwiseR2Files = pfiles_1 + pfiles_2
for file in pairwiseR2Files:
	print file
	with gzip.open(file,'r') as fin:
		for line in fin:
			(pos1, pos2, r2) = line.split()
			if r2 > rsquaredCutOff:
				position1.append(pos1)
				position2.append(pos2)

## 3. bagging
bag = [[] for i in range(0)]
pvalue = []
bagi = 0
unitBagi = 0
while len(eqtlPvalue_hash) > 1:
    if bagi >= len(bag):
            bag.append([])
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
    if len(indices_position2) > 0:
        rsquaredPositions.extend([ position1[index] for index in indices_position2 ])
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
        bagOutFile.close()
        pvalueOutFile = open('./bagPvalue.chr%s.print.txt' % CHR, 'a')
        for item in pvalue:
            pvalueOutFile.write("%s\n" % item)
        bagOutFile.close()
        pvalueOutFile.close()
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
    bagOutFile.close()
    pvalueOutFile.close()
