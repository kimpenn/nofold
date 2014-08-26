#====================================================================
# Sarah Middleton
# 
# varlen_z_normalize_pc.py
# 
#
# python varlen_z_normalize_pc.py ../demo/demo1/demo1.zNorm.pcNorm100.bitscore ../norm/varlen.zNorm.pcNorm.scale_pc.txt  ../demo/demo1/demo1.zNorm.pcNorm100.zNorm.bitscore
#====================================================================
import utils_file_readers as ufr
import utils_analysis as ua
import sys, random, math

if len(sys.argv) == 4:
	inFile = sys.argv[1]
	normFile = sys.argv[2]
	outFile = sys.argv[3]
else:
	print "Incorrect args. Exiting."
	sys.exit()


(header, scores, error) = ufr.read_bitscore(inFile)
fileHeader = "\t".join(header)
if error:
	print "Error reading bitscore.Exiting."
	sys.exit()

idOrder = []
ins = open(inFile, 'r')
ins.readline()
for line in ins:
	things = line.split()
	id = things[0]
	idOrder.append(id)
ins.close()
	
# read norm constants
means = {}
sds = {}
ins = open(normFile, 'r')
ins.readline()
for line in ins:
	things = line.split()
	id = things[0]
	means[id] = things[1]
	sds[id] = things[2]
ins.close()

outs = open(outFile, 'w')
for cm in header:
	if cm in means:
		outs.write(cm + "\t")
outs.write("\n")
for id in idOrder:
	outs.write(id)
	for i in range(len(scores[id])):
		cm = header[i]
		if cm in means:
			newScore = (float(scores[id][i]) - float(means[cm])) / float(sds[cm])
			outStr = "\t%.4f" % newScore
			outs.write(outStr)
	outs.write("\n")
outs.close()


print ""
print "Output printed to", outFile
print ""
	
	
	
	
	
	
	
	
	
	