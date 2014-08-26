#=================================================================================================
# Sarah Middleton
# 
# filter_clusters2.py
# - Filters clusters based on uniqueness, avg dist, and size.
# - Does no anaylsis based on group.
# - After filtering, clusters are ranked and a non-overlapping set is chosen
# - Output is a file listing info about all clusters that passed the filters
# 
# Required args: 
#   CLUSTER_FILE = relative path to cluster results file (as output from hclust_analysis.r)
#   MATRIX_FILE = relative path to the distance matrix for these data (from hclust_analysis.r)
#   UPPER_BOUNDS = relative path to the file containing upper avg dist bounds determined by clustering
#      of random sequences using calculate_clust_thresh.py -OR- a tuple w/o spaces (e.g. ci95,ci99) to
#      use as the cutoff for all cluster sizes. 
#   OUT_FILE = relative path and file name prefix for output files:
#      <OUT_FILE>_[options]_all.txt = all unique clusters of size >= MIN_SIZE
#      <OUT_FILE>_[options]_top.txt = top non-overlapping clusters, as chosen by the filters/ranking
#      the [options] suffix is automatically generated and prevents overwriting.
#   MIN_SIZE = minimum cluster size allowed to pass.
#   RANK_METHOD = "sensitive" or "specific". Method for ranking and chosing non-overlapping final set.
#      Sensitive method: Rank firstly by size (large to small) and secondly by avg dist (small to large). 
#      pick clusters from top to bottom, throwing out any that overlap with a sequence in a previous cluster.
#      Specific: Rank by avg dist (small to large). Pick clusters as above.
#   FASTA_DB = the fasta db of sequence in these clusters or ""
#=================================================================================================
import sys, random, time, os
from utils_file_readers import *
from utils_analysis import *

# command line options (required)
if len(sys.argv) == 10:
	CLUSTER_FILE = sys.argv[1]
	MATRIX_FILE = sys.argv[2]
	UPPER_BOUNDS = sys.argv[3]
	OUT_FILE = sys.argv[4]
	MIN_SIZE = int(sys.argv[5])
	RANK_METHOD = sys.argv[6]
	FASTA_DB = sys.argv[7]
	SCI_BOUND = float(sys.argv[8])
	RATIO_BOUND = float(sys.argv[9])
	
	if RANK_METHOD == "sensitive": rm = "Sens"
	elif RANK_METHOD == "specific": rm = "Spec"
	else:
		print "Error: please choose 'sensitive' or 'specific' for RANK_METHOD."
		print "Exiting."
		sys.exit(1)
	outSuffix = "s%sr%s" % (MIN_SIZE, rm)
else:
	print "Incorrect command line args."
	print "Exiting."
	sys.exit(1)


### read in files

# read in distance matrix
print ""
print "Reading distance matrix file..."
(idList, dissim, error) = read_comp_matrix(MATRIX_FILE) #read file
totalIDs = len(idList)
if error:
	print "Could not read dist matrix. Exiting."
	sys.exit(1)

# read in bounds
bounds = {}
if os.path.isfile(UPPER_BOUNDS):
	print "Reading bounds file..."
	ins = open(UPPER_BOUNDS, 'r')
	ins.readline() #skip header
	lines = ins.readlines()
	ins.close()
	for line in lines:
		line = line.rstrip('\r\n')
		things = line.split()
		size = things[0]
		ci95 = things[-3]
		ci99 = things[-2]
		bounds[int(size)] = {'95': float(ci95), '99': float(ci99)}
	maxBound = max(bounds.keys())
	if maxBound < totalIDs:
		for i in range(maxBound, (totalIDs + 1)): #copy max bound (this is conservative, if anything, since avg dist increases with cluster size)
			bounds[i] = bounds[maxBound]
else:
	things = UPPER_BOUNDS.split(",")
	ci95 = things[-2]
	ci99 = things[-1]
	for i in range(totalIDs + 1): #set this bound for all cluster sizes
		bounds[i] = {'95': float(ci95), '99': float(ci99)}

# read in clusters 
print "Reading cluster file..."
clusters = {}
ins = open(CLUSTER_FILE, 'r')
lines = ins.readlines()
ins.close()
for line in lines:
	line = line.rstrip('\r\n')
	if line != "":
		clusters[line] = True #identical clusters are removed in this way


### filter and analyze

# calculate avg dist for each cluster, check if it passes avg dist bounds
print "Calculating avg dist..."
totalClust = 0
aboveSize = 0
below95 = 0
below99 = 0
for cluster in clusters:
	totalClust += 1
	members = cluster.split(',')
	clustSize = len(members)
	cov = 0
	sig = "no"
	passed = False
	
	if clustSize >= 2:
		cov = get_avg_dist(members, dissim) #calculate avg dist
	if clustSize >= MIN_SIZE:
		aboveSize += 1
		if cov < bounds[clustSize]['95']:
			below95 += 1
			passed = True
			sig = "*"
		if cov < bounds[clustSize]['99']:
			below99 += 1
			passed = True
			sig = "**"
	clusters[cluster] = {'pass': passed, 'size': clustSize, 'cov': cov, 'sig': sig, 'cID': totalClust}

# pick non-overlapping clusters using indicated algorithm
finalSet = []
if RANK_METHOD == "sensitive":
	finalSet = pick_nonoverlap_sensitive(clusters, None)
elif RANK_METHOD == "specific":
	finalSet = pick_nonoverlap_specific(clusters, None)


# calculate MFE of cons struct for each cluster in finalSet
mfes = {}
"""
if FOLD > 0:
	folder = os.path.dirname(CLUSTER_FILE)
	tmpOutPath = "%s/tmp_mlocarna_output" % folder
	if not os.path.exists(tmpOutPath):
		os.makedirs(tmpOutPath)
	for cluster in clusters:
		if (clusters[cluster]['size'] <= FOLD) and (clusters[cluster]['size'] > 1):
			# create a fasta file for this set
			fastaOut = "%s/tmp.fa" % tmpOutPath
			outs = open(fastaOut, 'w')
			idList = cluster.split(",")
			seqs = get_sequences(idList, FASTA_DB, cutID=25)
			if len(seqs.keys()) == 0:
				print "Problem with getting seqs from fasta_db"
			for id in seqs:
				outStr = ">%s\n%s\n" % (id, seqs[id])
				outs.write(outStr)
			outs.close()
			
			# fold seqs
			tmpOutPath2 = "%s/clust%s" % (tmpOutPath, clusters[cluster]['cID'])
			(mfe, struct, error) = align_fold(fastaOut, outPath=tmpOutPath, otherOpts="--noLP")
			if error:
				print "Error: could not fold %s, skipping." % cluster
			mfes[cluster] = mfe
"""


### print results

# print _all file
print "Printing..."
outFile = "%s_%s_all.txt" % (OUT_FILE, outSuffix)
outs = open(outFile, 'w')
outs.write("ClustID\tSigLvl\tCount\tMFE\tAvgDist\tMembers\n")
tmpCovs = {} #for sorting by avg dist
for cluster in clusters:
	tmpCovs[cluster] = clusters[cluster]['cov']
for cluster in sorted(tmpCovs, key=tmpCovs.get): #print in order of increasing avg dist
	mfe = "-"
	if cluster in mfes:
		mfe = mfes[cluster]
	if clusters[cluster]['size'] >= 2:
		outString = "%s\t%s\t%s\t%s\t%.2f\t%s\n" % (clusters[cluster]['cID'], clusters[cluster]['sig'], clusters[cluster]['size'], mfe, clusters[cluster]['cov'], cluster)
		outs.write(outString)
outs.close()
print ""
print "All clusters with %s or more members printed to %s" % (MIN_SIZE, outFile)

# print _top file
print "Printing..."
outFile = "%s_%s_top.txt" % (OUT_FILE, outSuffix)
outs = open(outFile, 'w')
outs.write("ClustID\tSigLvl\tCount\tMFE\tAvgDist\tMembers\n")
for cluster in finalSet:
	mfe = "-"
	if cluster in mfes:
		mfe = mfes[cluster]
	outString = "%s\t%s\t%s\t%s\t%.2f\t%s\n" % (clusters[cluster]['cID'], clusters[cluster]['sig'], clusters[cluster]['size'], mfe, clusters[cluster]['cov'], cluster)
	outs.write(outString)
outs.close()
print ""
print "Top clusters after filtering printed to", outFile



# print summary to screen
print ""
print "Summary:"
print "# unique clusters total:", totalClust
print "# clusters with at least",MIN_SIZE,"members:", aboveSize
print "# clusters with avg distance less than ci95:", below95
print "# clusters with avg distance less than ci99:", below99
print "# clusters in final set, after filtering:", len(finalSet)
print ""











