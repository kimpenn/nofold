#=================================================================================================
# Sarah Middleton
# 
# calculate_clust_thresh.py
# - Does multiple clusterings of random subsets, gets distribution of avg distances for each cluster size
# - This can be used to determine a cutoff avg distance to use to reduce false discovery rate
# - To get the best cutoffs for larger clusters, do many iterations and draw from a large # random seqs
# - Note: the included file of scored random seqs (bgDinucWT_50) contains 10,000 sequences. Therefore this can't
#   be used to get thresholds for datasets that are larger than 10,000 sequences (and probably won't be very accurate for 
#   datasets larger than maybe 5,000 sequences, since similar subsets will be generated each time).
#
# Required args: 
#   TARGET_DB_SIZE = size of the target database for which these bounds will be used. E.g. if you're using these bounds
#      for the clustering of a dataset of 1000 sequences, set this to 1000. It's fine to round up or down a little bit,
#      but rounding up too much will give overly stringent bounds, and rounding down too much can increase the FDR.
#
# Advanced args (must edit file.):
#   RRT_FILE = bitscore file for the "random real" seqs that will be used
#   OUT_FILE = prefix name for output file. right now, only one file will be created:
#      <OUT_FILE>_<TARGET_DB_SIZE>seq.txt - cutoffs for each size cluster, for all sizes with >100 observations
#   NUM_ITERATIONS = number of subsets to create, cluster, and analyze. more = longer runtime. Automatically calculated to yield ~10,000 clusters of size 3
#   SCALING = "scaled" or "unscaled" or "prescaled". If "scaled", the file will be scaled. Otherwise, the file is used as-is.
#   CLUST_METHOD = method to use for hclust in R, e.g. "ward" or "average"
#   SCALE_VALS = File with mean/sd constants for scaling
#   DIST_METHOD = method to use to calculate distance matrix. Options: euclidean, pearson, spearman
#
# Usage ex: 
#   python calculate_clust_thresh.py 100
#=================================================================================================
import sys, random, time, os, subprocess
from utils_file_readers import *
from utils_analysis import *

# command line options (required)
if len(sys.argv) == 2:
	TARGET_DB_SIZE = int(sys.argv[1])
else:
	print "Must indicate command line args: TARGET_DB_SIZE"
	print "Exiting."
	sys.exit(1)

# advanced args -- currently set to be consistent with the clustering 
# linkage/distance measure used for clustering real datasets
RRT_FILE = "../thresh/bgDinucWT_50.zNorm.pcNorm100.zNorm.bitscore"
OUT_FILE = "../thresh/bounds"
SCALING = "prescaled"
SCALE = "FALSE"
CLUST_METHOD = "average"
SCALE_VALS = "none"
DIST_METHOD = "spearman"

# estimate the number of iterations needed to get ~10000 clusters of size 3 across all trials
m = 0.1645	#estimated from previous tests
b = -2.9654	#estimated from previous tests
desiredNum = 10000
NUM_ITERATIONS = int(desiredNum / (m * TARGET_DB_SIZE + b))


# temp output
if not os.path.exists("../thresh/tmp"):
	os.makedirs("../thresh/tmp")
timestamp = int(time.time())
tmpSubset = "../thresh/tmp/%s_subset" % timestamp
tmpHclustOut = "../thresh/tmp/%s_hclust" % timestamp
tmpClusters = "../thresh/tmp/%s_hclust.clusters" % timestamp
tmpDist = "../thresh/tmp/%s_hclust.distMatrix" % timestamp
tmpHeights = "../thresh/tmp/%s_hclust.heights" % timestamp


# read in random seq bitscore lines
ins = open(RRT_FILE, 'r')
cmHeader = ins.readline()
rrtLines = ins.readlines()
ins.close()


# for each iteration, create a new subset, cluster, and record avg distances within each possible cluster
results = {}
heights = {}
for iter in range(NUM_ITERATIONS):
	print "Iteration %s of %s:" % ((iter + 1), NUM_ITERATIONS)
	
	# create and print random subset
	print "  Creating subset..."
	linesCopy = rrtLines[:] 
	outs = open(tmpSubset, 'w')
	outs.write(cmHeader) #print header
	for i in range(TARGET_DB_SIZE):
		randLine = random.choice(linesCopy)
		outs.write(randLine)
		linesCopy.remove(randLine)
	outs.close()
	
	# cluster subset
	print "  Clustering subset... This may take a while, depending on the size of your subset."
	start = time.time()
	command = "Rscript get_clusters.r %s %s %s %s %s %s" % (tmpSubset, tmpHclustOut, SCALE, CLUST_METHOD, SCALE_VALS, DIST_METHOD)
	print "  ", command
	p = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
	retval = p.wait()
	if retval != 0:
		print "Error: something went wrong with get_clusters.r"
		print "Exiting."
		sys.exit(1)
	print "    Finished clustering in %.2f seconds" % (time.time() - start)
	
	# read resulting dist matrix
	print "  Reading distance file..."
	(idList, dissim, error) = read_comp_matrix(tmpDist) #read file
	if error:
		print "Exiting."
		sys.exit(1)
	
	# read resulting clusters
	print "  Reading cluster file..."
	clusters = {}
	ins = open(tmpClusters, 'r')
	lines = ins.readlines()
	ins.close()
	for line in lines:
		line = line.rstrip('\r\n')
		if line != "":
			clusters[line] = True #identical clusters are removed in this way
	
	# read resulting heights file
	print "  Reading heights file..."
	heights[iter] = []
	ins = open(tmpHeights, 'r')
	for line in ins:
		ht = line.rstrip('\r\n')
		heights[iter].append(ht)
	ins.close()
	
	# get avg distance for each cluster, record
	print "  Getting avg distances..."
	for cluster in clusters:
		members = cluster.split(",")
		size = len(members)
		if size > 1:
			cov = get_avg_dist(members, dissim)
			if size not in results:
				results[size] = []
			results[size].append(cov)
		

# sort and print results
print ""
print "Printing..."


# calculate and print _bySize file
# remember: smaller avg distance is better
# won't print anything with fewer than 50 covs
outFileName2 = "%s_%sseq.txt" % (OUT_FILE, TARGET_DB_SIZE)
outs = open(outFileName2, 'w')
outs.write("ClustSize\tMean\tCI75\tCI90\tCI95\tCI99\tNumObs\n")
for size in sorted(results.keys()):
	covs = results[size]
	numCovs = len(covs)
	if numCovs < 50:
		break #stop printing	
	mean = "%.4f" % (float(sum(covs)) / numCovs)
	bot25percent = int(numCovs * 0.25)
	bot10percent =int(numCovs * 0.10)
	bot5percent = int(numCovs * 0.05) #how many results make up the smallest 5% of results (i.e. the best 5%)
	bot1percent = int(numCovs * 0.01) #how many results make up the smallest 1% of results (i.e. the best 1%)
	sortedCovs = sorted(covs)
	ci75 = "%.4f" % sortedCovs[bot25percent]
	ci90 = "%.4f" % sortedCovs[bot10percent]
	ci95 = "%.4f" % sortedCovs[bot5percent]
	ci99 = "%.4f" % sortedCovs[bot1percent]
	outs.write(str(size) + "\t" + mean + "\t" + ci75 + "\t" + ci90 + "\t" + ci95 + "\t" + ci99 + "\t" + str(numCovs) + "\n")
outs.close()



print ""
print "Output printed to", outFileName2
print ""





