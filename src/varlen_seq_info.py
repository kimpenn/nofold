#====================================================================
# Sarah Middleton
# 
# varlen_seq_info.py
# Prints file with length and other interesting info about each sequence, for plotting with bitscore data.
# The "group" is taken to be the first part of the ID, which is separated from the rest of the ID by an underscore.
# If there is no underscore, the whole ID is used as the "group". 
#
# python varlen_seq_info.py ../demo/demo1/demo1.db ../demo/demo1/demo1.bitscore ../demo/demo1/demo1.seq_info.txt
#====================================================================
import utils_file_readers as ufr
import utils_analysis as ua
import sys, random

if len(sys.argv) == 4:
	dbFile = sys.argv[1]
	bsFile = sys.argv[2]
	outFile = sys.argv[3]
else:
	print "Incorrect args. Exiting."
	sys.exit()


(seqs, error) = ufr.read_fasta(dbFile)
if error:
	print "Error reading input fasta. Exiting."
	sys.exit()

# read bitscore IDs IN ORDER
idList = []
ins = open(bsFile, 'r')
ins.readline()
for line in ins:
	things = line.split()
	id = things[0]
	idList.append(id)
ins.close()

# get info
info = {}
for id in idList:
	if id in seqs:
		info[id] = {}
		seq = seqs[id]
		idParts = id.split("_")
		info[id]['len'] = len(seq)
		info[id]['gc'] = ua.gc_content(seq)
		info[id]['group'] = idParts[0]
	
# print results to file IN ORDER
outs = open(outFile, 'w')
outs.write("ID\tLen\tGC\tGroup\n")
for id in idList:
	outStr = "%s\t%s\t%.3f\t%s\n" % (id, info[id]['len'], info[id]['gc'], info[id]['group'])
	outs.write(outStr)
outs.close()

print ""
print "Output printed to", outFile
print ""



















