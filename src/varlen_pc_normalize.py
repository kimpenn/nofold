#====================================================================
# Sarah Middleton
# 
# varlen_pc_normalize.py
# To use all PC's, specify 0 for last param
#
# python varlen_pc_normalize.py ../demo/demo1/demo1.zNorm.bitscore ../norm/varlen.eigs_subset.txt ../demo/demo1/demo1.zNorm.pcNorm100.bitscore 100
#====================================================================
import utils_file_readers as ufr
import utils_analysis as ua
import sys, random, math

if len(sys.argv) == 5:
	inFile = sys.argv[1]
	eigFile = sys.argv[2]
	outFile = sys.argv[3]
	numPCs = int(sys.argv[4])
else:
	print "Incorrect args. Exiting."
	sys.exit()


# read in scores (rows = samples, cols = features)
(header, scores, error) = ufr.read_bitscore(inFile)
fileHeader = "\t".join(header)
if error:
	print "Error reading bitscore. Exiting."
	sys.exit()

# get order of sequences
idOrder = []
ins = open(inFile, 'r')
ins.readline()
for line in ins:
	things = line.split()
	idOrder.append(things[0])
ins.close()
	
# read in eigenvectors (rows = eigs, cols = features)
(headerEigs, eigs, error) = ufr.read_bitscore(eigFile) #not a bitscore file, but same fmt
if error:
	print "Error reading eigFile. Exiting."
	sys.exit()

# remove uneeded PCs
removedCount = 0
if numPCs > 0:
	pcs = eigs.keys()
	for pc in pcs:
		num = int(pc[2:])
		if num > numPCs:
			del eigs[pc]
			removedCount += 1
print "Removed %s PCs" % removedCount

# check that headers are the same
for i in range(len(header)):
	if header[i] != headerEigs[i]:
		print "Error: headers do not match: %s vs %s. Exiting." % (header[i], headerEigs[i])
		sys.exit()


# get dimensions
numSamps = len(scores.keys())
numEigs = len(eigs.keys())
numFeatures = len(header)


# map (scores x t(eigs))
transf = {}
for id in scores:
	transf[id] = []
	for eig in sorted(eigs):
		sum = 0
		for i in range(numFeatures):
			sum += (eigs[eig][i] * scores[id][i])
		transf[id].append(sum)


# print results
outs = open(outFile, 'w')
for eig in sorted(eigs):
	outs.write(eig + "\t")
outs.write("\n")
for id in idOrder:
	outStr = id
	for val in transf[id]:
		outStr += "\t%.3f" % val
	outs.write(outStr + "\n")
outs.close()


print ""
print "Output printed to", outFile
print ""
	
	
	
	
	
	
	
	
	
	
