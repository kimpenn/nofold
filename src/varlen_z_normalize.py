#====================================================================
# Sarah Middleton
# 
# varlen_z_normalize.py
# 
#
# python varlen_z_normalize.py ../demo/demo1/demo1.bitscore ../demo/demo1/demo1.seq_info.txt ../norm/varlen.scale_means.txt ../norm/varlen.scale_sds.txt ../demo/demo1/demo1.zNorm.bitscore
#====================================================================
import utils_file_readers as ufr
import utils_analysis as ua
import sys, random, math

if len(sys.argv) == 6:
	inFile = sys.argv[1]
	lenFile = sys.argv[2]
	meanFile = sys.argv[3]
	sdFile = sys.argv[4]
	outFile = sys.argv[5]
else:
	print "Incorrect args. Exiting."
	sys.exit()


(header, scores, error) = ufr.read_bitscore(inFile)
fileHeader = "\t".join(header)
if error:
	print "Error reading bitscore.Exiting."
	sys.exit()

idOrder = []
id2Len = {}
ins = open(lenFile, 'r')
ins.readline()
for line in ins:
	things = line.split()
	id = things[0]
	length = things[1]
	id2Len[id] = int(length)
	idOrder.append(id)
ins.close()
	
# read means
means = {}
ins = open(meanFile,'r')
meanHeader = ins.readline()
splitMeanHeader = meanHeader.split()
for line in ins:
	things = line.split()
	size = int(things[0])
	meanVals = things[1:]
	means[size] = {}
	for i in range(len(splitMeanHeader)):
		cm = splitMeanHeader[i]
		if meanVals[i] == "NA":
			means[size][cm] = "NA"
		else:
			means[size][cm] = float(meanVals[i])
ins.close()

# read sds
sds = {}
ins = open(sdFile,'r')
sdHeader = ins.readline()
splitSdHeader = sdHeader.split()
for line in ins:
	things = line.split()
	size = int(things[0])
	sdVals = things[1:]
	sds[size] = {}
	for i in range(len(splitSdHeader)):
		cm = splitSdHeader[i]
		if sdVals[i] == "NA":
			sds[size][cm] = "NA"
		else:
			sds[size][cm] = float(sdVals[i])
ins.close()

maxSize = max(means.keys())

outs = open(outFile, 'w')
for cm in header:
	if cm in splitMeanHeader:
		outs.write(cm + "\t")
outs.write("\n")
for id in idOrder:
	outs.write(id)
	idLen = id2Len[id]
	if idLen > maxSize:
		print ">>> Warning: in varlen_z_normalize.py: seq for %s was too long (%snt). Using largest length parameters available (%snt)." % (id, idLen, maxSize)
		idLen = maxSize
	for i in range(len(scores[id])):
		cm = header[i]
		if cm in means[idLen]:
			newScore = (float(scores[id][i]) - float(means[idLen][cm])) / float(sds[idLen][cm])
			outStr = "\t%.4f" % newScore
			outs.write(outStr)
	outs.write("\n")
outs.close()


print ""
print "Output printed to", outFile
print ""
	
	
	
	
	
	
	
	
	
	