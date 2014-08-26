#=================================================================================================
# Sarah Middleton
# 
# utils_analysis.py
#
# Functions for common data analyses and sequence manipulations
#=================================================================================================
from utils_file_readers import *


def run_command(command, verbose=False):
	"""
	Run the given command using the system shell
	*fix this to print output as it goes
	"""
	import subprocess
	error = False
	
	if verbose == True: 
		print command
		print ""
	
	job = subprocess.Popen(command, bufsize=0, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
	jobOutput = []
	if verbose == True:
		for line in job.stdout:
			print "   ", line,
			jobOutput.append(line)
	else:
		jobOutput = job.stdout.readlines()
	result = job.wait()
	if result != 0:
		error = True
	
	return (jobOutput, result, error)


#-----------------------------------------------------------------------------------------------
# Sequence processing
#-----------------------------------------------------------------------------------------------


def r_adjust_size(seq, size, bufferSource):
	"""
	Given a sequence, trims or buffers the sequence from the right to make it the indicated size.
	Uses sequence from bufferSource (a dictionary) for buffering. 
	Example of bufferSource:
		BUFFER_SOURCE = "../data/UCSCgenes_knownGene_mm9_WT.fasta"
		(buffers, error) = read_fasta(BUFFER_SOURCE)
		r_adjust_size(seqs, 50, buffers)
	"""
	import random
	if len(seq) > size: #trim (from right)
		seq = seq[:size]
	elif len(seq) < size: #buffer (on right)
		found = False
		needed = size - len(seq)
		while(found == False):
			rrtID = random.choice(list(bufferSource.keys()))
			rrtSeq = bufferSource[rrtID]
			if len(rrtSeq) > size:
				found = True
		rStart = int(random.randrange(0, (len(rrtSeq) - needed)))
		bufferSeq = rrtSeq[rStart:(rStart + needed)]
		seq += bufferSeq
	if len(seq) != size: # make sure the resulting sequence is the desired length
		print "Warning: in r_adjust_seq(): final seq size was incorrect."
		
	return seq



def l_adjust_size(seq, size, bufferSource):
	"""
	Given a sequence, trims or buffers the sequence from the left to make it the indicated size.
	Uses sequence from bufferSource (a dictionary) for buffering. 
	Example of bufferSource:
		BUFFER_SOURCE = "../data/UCSCgenes_knownGene_mm9_WT.fasta"
		(buffers, error) = read_fasta(BUFFER_SOURCE)
		l_adjust_size(seqs, 50, buffers) 
	"""
	import random
	if len(seq) > size: #trim (from left)
		start = len(seq) - size
		seq = seq[start:]
	elif len(seq) < size: #buffer (on left)
		found = False
		needed = size - len(seq)
		while(found == False):
			rrtID = random.choice(list(bufferSource.keys()))
			rrtSeq = bufferSource[rrtID]
			if len(rrtSeq) > size:
				found = True
		rStart = int(random.randrange(0, (len(rrtSeq) - needed)))
		bufferSeq = rrtSeq[rStart:(rStart + needed)]
		seq = bufferSeq + seq
	if len(seq) != size: # make sure the resulting sequence is the desired length
		print "Warning: in l_adjust_seq(): final seq size was incorrect."
		
	return seq



def create_sliding_window(seqs, size, slide, buffer=False, bufferSource=None, reverse=False, offset=0):
	"""
	-Takes a db of variable length sequences ("seqs") and converts it to a set of fragments of the same size
	-If size != slide, a sliding window will be used to generate fragments
	-db should be a dictionary of sequences in the fmt dict[id]->seq 
	-output is in same format as input
	-if reverse=True, sliding will start from 3' end
	-if buffer=True (requires bufferSource, a db of seqs in fmt dict[id]->seq, be supplied), add random sequence 
		to fragments that are too short (but must be at least size/2)
	"""
	fragments = {}
	
	if reverse == False: #start at 5'
		for id in seqs:
			seq = seqs[id]
			seqLen = len(seq)
			subseqNum = 0
			currSt = offset
			currEnd = currSt + size
			
			while (currEnd <= seqLen):
				subseq = seq[currSt:currEnd]
				newID = "%s_%s" % (id, subseqNum)
				fragments[newID] = subseq
				
				if len(newID) > 25:
					print "Warning: in create_sliding_window(): id %s is longer than 25 char" % newID
				if len(subseq) != size:
					print "Warning: in create_sliding_window(): id %s sequence is incorrect size" % newID
				
				subseqNum += 1
				currSt += slide
				currEnd += slide
			
			if buffer == True:
				if currSt <= seqLen:
					subseq = seq[currSt:]
					if len(subseq) >= (size / 2):
						subseq = r_adjust_size(subseq, size, bufferSource)
						newID = "%s_%s" % (id, subseqNum)
						fragments[newID] = subseq
						
	else: #start at 3'
		for id in seqs:
			seq = seqs[id]
			seqLen = len(seq)
			subseqNum = 0
			currEnd = seqLen - offset
			currSt = currEnd - size
			
			while (currSt >= 0):
				subseq = seq[currSt:currEnd]
				newID = "%s_%s" % (id, subseqNum)
				fragments[newID] = subseq
				
				if len(newID) > 25:
					print "Warning: in create_sliding_window(): id %s is longer than 25 char" % newID
				if len(subseq) != size:
					print "Warning: in create_sliding_window(): id %s sequence is incorrect size" % newID
				
				subseqNum += 1
				currSt -= slide
				currEnd -= slide
			
			if buffer == True:
				if currEnd > 0:
					subseq = seq[:currEnd]
					if len(subseq) >= (size / 2):
						subseq = l_adjust_size(subseq, size, bufferSource)
						newID = "%s_%s" % (id, subseqNum)
						fragments[newID] = subseq
	
	return fragments
	

def get_sequences(idList, fastaFile, cutID=-1):
	"""
	Given a list of sequence IDs and a fasta file of sequences, prints all the sequences that
	match any of the given IDs.
	If cutID is set to something greater than 0, the id in the fasta file will be cut to that length
	(necessary if the id was truncated by another program, such as infernal)
	"""
	ins = open(fastaFile, 'r')
	lines = ins.readlines()
	ins.close()
	
	seqs = {}
	id = ""
	for line in lines:
		line = line.rstrip('\r\n')
		if ">" in line:
			if (cutID > 0):
				id = line[1:(cutID + 1)]
			else:
				id = line[1:]
			if id in idList:
				seqs[id] = ""
		else:
			if id in idList:
				seqs[id] += line
	
	return seqs



def shuffle_nt(seq):
	"""
	Given a sequence, shuffles it while maintaining exact nt frequency
	"""
	import random
	
	strList = list(seq)
	random.shuffle(strList)
	shuffSeq = "".join(strList)
	
	return shuffSeq



def batch_shuffle_nt(inFile, outFile):
	"""
	Given an input file in fasta format, shuffles each sequence individually, maintaining
	exact nucleotide frequency, and outputs in fasta format to outFile
	"""
	from utils_file_readers import read_fasta
	
	(seqs, error) = read_fasta(inFile)
	if error:
		print "Error reading file", inFile
	
	try:
		outs = open(outFile, 'w')
	except IOError:
		print "Error writing file", outFile
		error = True
	else:
		for id in sorted(seqs.keys()):
			seq = seqs[id]
			shuffSeq = shuffle_nt(seq)
			outStr = ">%s\n%s\n" % (id, shuffSeq)
			outs.write(outStr)
		outs.close()
	
	return error


def shared_seq(seq1, seq2, pos1, pos2):
 	"""
 	Checks the seq surrounding pos1 on seq1 and pos2 on seq2 for identical sequence.
 	Returns the number of nt up and downstream that match
 	"""
 	len1 = len(seq1)
	len2 = len(seq2)
	
	# compare nt at -1, -2, ... until a mismatch, record that number "up"
	up = 0
	while ((pos1 - up) > 0) and ((pos2 - up) > 0):
		if seq1[(pos1 - up)] == seq2[(pos2 - up)]:
			up += 1
		else:
			break
	
	# compare nt at +1, +2, ... until a mismatch, record that number "down"
	down = 0
	while ((pos1 + down) < len1) and ((pos2 + down) < len2):
		if seq1[(pos1 + down)] == seq2[(pos2 + down)]:
			down += 1
		else:
			break
	
	return (up, down)



def get_cdna(seq):
	"""
	Convert an mRNA sequence to cDNA. Multi-line seqs ok. Outputs in same format as input
	(conserves line breaks, upper/lower case)
	"""
	
	
	
def get_db_size(db):
	"""
	Given a dictionary in format dict[id]->seq, calculate the total length of all seqs.
	"""
	totalLen = 0
	for id in db:
		totalLen += len(db[id])
	return totalLen
	



#-----------------------------------------------------------------------------------------------
# Cluster/Group Analysis - functions that analyze a list of genes
#-----------------------------------------------------------------------------------------------



def get_avg_dist(ids, distMatrix):
	"""
	Given a list of ids and a reference to a distance (or similarity) matrix, compute the average 
	distance between all members of the list.
	If only one id is given, returns 0
	"""
	totalDissim = 0
	totalComparisons = 0
	size = len(ids)
	if size > 1:
		for i in range(size):
			for j in range((i + 1), size): #avoid redundant comparisons
				totalDissim += distMatrix[ids[i]][ids[j]] #get dissim for this pair and add it
				totalComparisons += 1
		avgDissim = float(totalDissim) / totalComparisons
	else:
		print "Warning: in get_avg_dist(): only one id was given, so no comparison done."
		avgDissim = 0
		
	return avgDissim



def get_cluster_group(members, assignments=None):
	"""
	Given a list of the ids in a cluster and a reference to a dict containing the group
	assignments of each id, return the group with the greatest number of members in the cluster.
	* Dict must be in format assignments[seqID] -> groupID
	"""
	tally = {}
	for seqID in members:
		if assignments == None:
			things = seqID.split("_")
			group = things[0]
		else:
			group = assignments[seqID]
		if group not in tally:
			tally[group] = 0
		tally[group] += 1
	topGroup = max(tally, key=tally.get)
	topCount = tally[topGroup]
	
	return (topGroup, topCount)



def pick_nonoverlap_sensitive(clusters, enrichPvalCutoff):
	"""
	"Sensitive" mode: pick the broadest good clusters
	   * Best if used with a stringent CoV upper bound to avoid overly noisy clusters
	Given a dict with clusters ( structure should be dict[cluster] -> {'pass': Boolean, 'size': int, 'cov': float, 'enr': float} )
	return a list of the top set of non-overlapping clusters whose enrichment pval is <= enrichPvalCutoff
	Algorithm: Rank firstly by size (large to small) and secondly by CoV (small to large). Pick 
	clusters from top to bottom, throwing out any that overlap with a sequence in a previous cluster.
	"""
	finalSet = []
	usedIDs = {}
	
	tmpSizes = {}
	for cluster in clusters: #get size of each passing cluster and store in a dict for easier sorting
		if (clusters[cluster]['pass'] == True):
			if (enrichPvalCutoff == None) or (clusters[cluster]['enr'] <= enrichPvalCutoff):
				size = clusters[cluster]['size']
				if size not in tmpSizes:
					tmpSizes[size] = []
				tmpSizes[size].append(cluster) 
	
	sortedSizes = sorted(tmpSizes.keys()) #sorts by size
	sortedSizes.reverse() #re-orders from large to small
	for size in sortedSizes: 
		tmpCovs = {}
		for cluster in tmpSizes[size]: #for each cluster of this size, get CoVs
			tmpCovs[cluster] = clusters[cluster]['cov'] #store in a dict for easier sorting
		for cluster in sorted(tmpCovs, key=tmpCovs.get): #sort clusters of this size by CoV, small to large
			passed = True
			members = cluster.split(",") #get members of this cluster
			for id in members: #check if any ids overlap with those used in another chosen cluster
				if id in usedIDs:
					passed = False
			if passed == True:
				for id in members: #mark ids in this cluster as used
					usedIDs[id] = True
				finalSet.append(cluster) #add this cluster to the final set
	
	return finalSet



def pick_nonoverlap_specific(clusters, enrichPvalCutoff):
	"""
	"Specific" mode: picks the tightest clusters first. minimizes false positives.
	   * Tends to pick the smaller clusters. Sensitivity can be increased by doing cluster expansion and merging later.
	Given a dict with clusters ( of the format dict[cluster] -> {'pass': Boolean, 'size': int, 'cov': float} ),
	return a list of the top set of non-overlapping clusters whose enrichment pval is <= enrichPvalCutoff
	Algorithm: Rank by CoV (small to large). Pick clusters from top to bottom, throwing out any that overlap with a 
	sequence in a previous cluster.
	"""
	finalSet = []
	usedIDs = {}
	
	tmpCovs = {}
	for cluster in clusters:
		if (clusters[cluster]['pass'] == True):
			if (enrichPvalCutoff == None) or (clusters[cluster]['enr'] <= enrichPvalCutoff):
				tmpCovs[cluster] = clusters[cluster]['cov']
	for cluster in sorted(tmpCovs, key=tmpCovs.get): #sort by value small to large
		passed = True
		members = cluster.split(",") #get members of this cluster
		for id in members: #check if any ids overlap with those used in another chosen cluster
			if id in usedIDs:
				passed = False
		if passed == True:
			for id in members: #mark ids in this cluster as used
				usedIDs[id] = True
			finalSet.append(cluster) #add this cluster to the final set
		
	return finalSet

def pick_nonoverlap_advanced(clusters, fastaDB, minSize, sciBound=0.75, ratioBound=1, mlocarnaArgs="--free-endgaps --indel=-100 --noLP --threads=4"):
	"""
	"Advanced" mode: picks clusters by specific method, then picks any that have a SCI > sciBound
	   * Tends to pick the smaller clusters. Sensitivity can be increased by doing cluster expansion and merging later.
	Given a dict with clusters ( structure should be dict[cluster] -> {'pass': Boolean, 'size': int, 'cov': float} ),
	return a list of the top set of non-overlapping clusters whose enrichment pval is <= enrichPvalCutoff
	Algorithm: Rank by CoV (small to large). Pick clusters from top to bottom, throwing out any that overlap with a 
	sequence in a previous cluster.
	"""
	import os, time
	finalSet = []
	usedIDs = {}
	
	# get specific clusters that 'passed' bounds cutoff
	tmpCovs = {}
	for cluster in clusters:
		if (clusters[cluster]['pass'] == True) and (clusters[cluster]['size'] == minSize):
			tmpCovs[cluster] = clusters[cluster]['cov']
	for cluster in sorted(tmpCovs, key=tmpCovs.get): #sort by value small to large
		passed = True
		members = cluster.split(",") #get members of this cluster
		for id in members: #check if any ids overlap with those used in another chosen cluster
			if id in usedIDs:
				passed = False
		if passed == True:
			for id in members: #mark ids in this cluster as used
				usedIDs[id] = True
			finalSet.append(cluster) #add this cluster to the final set
	
	# get structural clusters that did not pass the bounds cutoff
	runList = []
	for cluster in clusters:
		if (clusters[cluster]['pass'] == False) and (clusters[cluster]['size'] == minSize):
			passed = True
			members = cluster.split(",") #get members of this cluster
			for id in members: #check if any ids overlap with those used in another chosen cluster
				if id in usedIDs:
					passed = False
			if passed == True:
				for id in members: #mark ids in this cluster as used
					usedIDs[id] = True
				runList.append(cluster)
	print "Finding missed clusters using mlocarna and RNAz, this may take some time."
	tmpOutPath = "tmp_mlocarna_output"
	if not os.path.exists(tmpOutPath):
		os.makedirs(tmpOutPath)
	counter = 1
	timestamp = time.time()
	for cluster in runList:
		print "  Testing cluster %s of %s" % (counter, len(runList))
		counter += 1
		# create a fasta file for this set
		fastaOut = "%s/tmp%s.fa" % (tmpOutPath, timestamp)
		outs = open(fastaOut, 'w')
		idList = cluster.split(",")
		seqs = get_sequences(idList, fastaDB, cutID=25)
		if len(seqs.keys()) == 0:
			print "Problem with getting seqs from fasta_db for %s. Skipping." % cluster
		else: #print to fasta file
			for id in seqs:
				outStr = ">%s\n%s\n" % (id, seqs[id])
				outs.write(outStr)
			outs.close()
				
			# fold seqs
			(mfe, struct, error) = align_fold(fastaOut, outPath=None, otherOpts=mlocarnaArgs)
			if error:
				print "Error: could not fold %s, skipping." % cluster
			else: #run rnaz
				alnFile = "%s/tmp%s.out/results/result.aln" % (tmpOutPath, timestamp)
				(mpi, sci, error) = run_rnaz(alnFile)
				if error:
					print "Error: could not run rnaz on %s, skipping." % cluster
				else:
					if sci >= sciBound:
						ratio = (sci * 100) / mpi
						if ratio >= ratioBound:
							finalSet.append(cluster)
	
	return finalSet



def gene_overlap(mappingFile):
	"""
	Given a file with gene ID mappings to chromosomal positions, creates and returns a dictionary
	of overlapping genes.
	Ex. file: ../data/UCSCgenes_mm9_mappings1-2.txt
	"""
	starts = {}
	ends = {}
	ins = open(mappingFile, 'r')
	ins.readline() #skip header
	lines = ins.readlines()
	ins.close()
	for line in lines:
		line = line.rstrip('\r\n')
		things = line.split()
		if (len(things) == 5) or (len(things) == 4):
			startID = things[1] + "_" + things[2]
			endID = things[1] + "_" + things[3]
			if startID not in starts:
				starts[startID] = []
			starts[startID].append(things[0])
			if endID not in ends:
				ends[endID] = []
			ends[endID].append(things[0])
		
	overlapping = {}
	for pos in starts:
		idList = starts[pos]
		if len(idList) > 1:
			for i in range(len(idList)):
				id1 = idList[i]
				if id1 not in overlapping:
					overlapping[id1] = []
				for j in range(len(idList)):
					id2 = idList[j]
					if id1 != id2:
						overlapping[id1].append(id2)
	for pos in ends:
		idList = ends[pos]
		if len(idList) > 1:
			for i in range(len(idList)):
				id1 = idList[i]
				if id1 not in overlapping:
					overlapping[id1] = []
				for j in range(len(idList)):
					id2 = idList[j]
					if (id1 != id2) and (id2 not in overlapping[id1]):
						overlapping[id1].append(id2)
	
	return overlapping
	
	
	
	

def group_cm_zAvg(idList, bitscoreFile, scale=False):
	"""
	finds the avg scores for each CM for the sequences in idList.
	assumes file is already zscored. if not, use scale=True.
	"""
	(header, scores, error) = read_bitscore(bitscoreFile)
	groupAvgs = {}
	if not error:
		normScores = scores
		if scale == True:
			normScores = scale_scores(scores)
		numSeqs = len(idList)
		for i in range(len(header)):
			cmName = header[i]
			total = 0
			for id in idList:
				total += normScores[id][i]
			avg = float(total) / numSeqs
			groupAvgs[cmName] = avg
	
	return groupAvgs
		


def group_cm_pVal(idList, bitscoreFile):
	"""
	idList can be a string (comma separated, no spaces) or a list (will be converted to str). 
	calls r script to get p-value of the difference in distribution between the subset and the
	whole for each cm. uses Welch's one-sided two-sample t-test
	"""
	import os, subprocess
	
	if isinstance(idList, list): 
		idStr = ",".join(idList) 
	else:
		idStr = idList
	
	outPath = "tmp/"
	if not os.path.exists(outPath):
		os.makedirs(outPath)
	outFile = outPath + os.path.basename(bitscoreFile) + "_subgroup_pval"
	
	command = "Rscript get_signif_features.r %s %s %s" % (bitscoreFile, outFile, idStr)
	p = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
	retval = p.wait()
	if retval != 0:
		print "Warning: in group_cm_pVal(): r script failed. May be a cluster with only one member."
	
	groupPvals = {}
	try: 
		ins = open(outFile, 'r')
	except IOError: 
		print "Error: in group_cm_pVal(): could not open r script output file", outFile
	else:
		lines = ins.readlines()
		ins.close()
		for line in lines:
			(cm, pval) = line.split()
			groupPvals[cm] = float(pval)
	
	return groupPvals
	
	

def get_idList_overlap(idList1, idList2):
	"""
	Given two lists of ids, counts the overlap between them
	Returns number of overlaps, and the fraction of each cluster that is overlapped.
	"""
	overlap = 0
	if (len(idList1) > 0) and (len(idList2) > 0):
		for id in idList1:
			if id in idList2:
				overlap += 1
		frac1 = float(overlap) / len(idList1)
		frac2 = float(overlap) / len(idList2)
	else:
		frac1 = 0
		frac2 = 0
	
	return (overlap, frac1, frac2)


def find_overlapping_clusters(clustDict1, clustDict2, thresh):
	"""
	Given two cluster dictionaries (format should be dict[clustID]['idList'] = ['id1',...'idN']) and a threshold for
	what percent overlap is required for merging, this function will compare clusters across the two 
	files and merge when necessary. 
	Returns a dictionary with info about clusters that can be merged (format: dict[clustID1][clustID2] = numOverlap)
	"""
	overlapping = {}
	for id1 in clustDict1:
		members1 = clustDict1[id1]
		for id2 in clustDict2:
			members2 = clustDict[id2]
			(numOverlap, fraction1, fraction2) = get_idList_overlap(members1, members2)
			if (fraction1 >= thresh) or (fraction2 >= thresh):
				if id1 not in overlapping:
					overlapping[id1] = {}
				overlapping[id1][id2] = numOverlap
	
	return overlapping
	
	
	

	



#-----------------------------------------------------------------------------------------------
# Misc
#-----------------------------------------------------------------------------------------------

def gc_content(seq):
	"""
	returns % gc
	"""
	seq = seq.upper()
	numG = seq.count('G')
	numC = seq.count('C')
	return (float(numG + numC) / len(seq))
	


def count_occurrences(seq, subseq):
	"""
	Count (potentially overlapping) instances of a subsequence in a string
	"""
	seq = seq.upper()
	subseq = subseq.upper()
	count = 0
	index = 0
	done = False
	while (done == False):
		index = seq.find(subseq, index)
		if (index == -1):
			done = True
		else:  
			count += 1
			index += 1 # add one so this pos won't be found again
	return count


def make_partition(numElements, numPartitions):
	"""
	Create a list indicating a random partition of numElements
	"""
	import random
	partition = []
	partSize = int(numElements / numPartitions)
	remainder = numElements % numPartitions
	
	# assign initial values
	for i in range(numPartitions):
		for j in range(partSize):
			partition.append(i)
	if remainder > 0:
		for i in range(remainder):
			partition.append(i)
	
	# shuffle values
	random.shuffle(partition)
	
	# qc
	if len(partition) != numElements:
		print "Warning in utils_analysis.make_partition(): size of partition array does not match number of elements."
	
	return partition
	
	

def euclidean_dist(list1, list2):
	"""
	Compute euclidean distance between two vectors. Must be of the same length.
	"""
	dist = -1
	if len(list1) != len(list2):
		print "In euclidean_dist(): vectors must be of same length"
	else:
		sum = 0
		for	i in range(len(list1)):
			sum += ((list1[i] - list2[i]) ** 2)
		dist = sum ** 0.5
	
	return dist


def pearson_dist(list1, list2):
	"""
	Compute pearson correlation two vectors, subtract from 1 to get pearson distance. 
	Vectors must be of the same length.
	"""
	dist = -1
	if len(list1) != len(list2):
		print "In pearson_dist(): vectors must be of same length"
	else:
		dotProd12 = 0
		dotProd11 = 0
		dotProd22 = 0
		for	i in range(len(list1)):
			dotProd12 += (list1[i] * list2[i])
			dotProd11 += (list1[i] * list1[i])
			dotProd22 += (list2[i] * list2[i])
		corr = float(dotProd12) / ((dotProd11 * dotProd22) ** 0.5)
		dist = 1 - corr #a number between 0 and 2
	
	return dist

	
def pop_corr_sim(list1, list2):
	"""
	experimental... basically pearson correlation, but assuming the total population has been standardized
	(so all means are 0, all SDs are 1)
	"""
	sim = -1
	if len(list1) != len(list2):
		print "In pop_corr_sim(): vectors must be of same length"
	else:
		dotProd12 = 0
		for	i in range(len(list1)):
			dotProd12 += (list1[i] * list2[i])
		sim = float(dotProd12) / (len(list1))
	
	return sim


def spearman_dist(list1, list2):
	"""
	Compute spearman correlation two vectors, subtract from 1 to get spearman distance. 
	Vectors must be of the same length.
	"""
	dist = -1
	if len(list1) != len(list2):
		print "In spearman_dist(): vectors must be of same length"
	else:
		order1 = sorted(range(len(list1)), key=lambda z: list1[z]) #returns order of indices of sorted list. sort ascending.
		ranks1 = sorted(range(len(order1)), key=lambda z: order1[z]) #run again to get ranks
		order2 = sorted(range(len(list2)), key=lambda z: list2[z]) 
		ranks2 = sorted(range(len(order2)), key=lambda z: order2[z]) 
		dotProd12 = 0
		dotProd11 = 0
		dotProd22 = 0
		for	i in range(len(list1)):
			dotProd12 += (ranks1[i] * ranks2[i])
			dotProd11 += (ranks1[i] * ranks1[i])
			dotProd22 += (ranks2[i] * ranks2[i])
		corr = float(dotProd12) / ((dotProd11 * dotProd22) ** 0.5)
		dist = 1 - corr #a number between 0 and 2
	
	return dist
	
	
	

def mass_dist(list1, list2, scores):
	"""
	Compute "mass" distance (as described in Yona et al. 2006; doi: 10.1093/bioinformatics/btl127)
	A measure between 0 and 1 based on the distribution of values for each feature. Disfavors similarity that is
	based too heavily on "average" values.
	"""
	dist = -1
	if len(list1) != len(list2):
		print "In mass_dist(): vectors must be of same length"
	else:
		sum = 0
		numFeatures = len(list1)
		numSeqs = len(scores.keys())
		
		# for each feature i, count number of seqs with feature value between that of list1[i] and list2[i]
		for	i in range(len(list1)):
			if list1[i] > list2[i]:
				val1 = list1[i]
				val2 = list2[i]
			else:
				val2 = list1[i]
				val1 = list2[i]
			
			count = 0
			for id in scores:
				if (scores[id][i] < val1) and (scores[id][i] > val2):
					count += 1
				
			mass = float(count) / numSeqs #convert to freq [0,1]
			sum += mass
			
		dist = float(sum) / numFeatures #normalize by total number of features to make it [0,1]
	
	return dist


def intersect_sim2(list1, list2):
	"""
	Computes 'intersect' similarity between each pair of seqs A and B as sum over i of (min(Ai,Bi))
	Adds the amount by which they vary in the same direction.
	"""
	similarity = -1
	if len(list1) != len(list2):
		print "In intersect_dist(): vectors must be of same length"
	else:
		sum = 0
		numFeats = len(list1)
		for	i in range(numFeats):
			if (list1[i] < 0) and (list2[i] < 0):
				sum += abs(max(list1[i], list2[i]))
			elif (list1[i] >= 0) and (list2[i] >= 0):
				sum += min(list1[i], list2[i])
			else:
				sum += 0
		similarity = float(sum) / numFeats
	
	return similarity
	
def intersect_sim(list1, list2):
	"""
	
	"""
	similarity = -1
	if len(list1) != len(list2):
		print "In intersect_dist(): vectors must be of same length"
	else:
		sum = 0
		numFeats = len(list1)
		for	i in range(numFeats):
			mult = list1[i] * list2[i]
			if mult > 0:
				sum += mult
		similarity = float(sum) / numFeats
	
	return similarity
	
def scale_scores(scores):
	"""
	Given a dict of scores created by read_bitscore(), z-scores the "columns"
	"""
	scaled = {}
	numCols = len(scores[scores.keys()[0]])
	numRows = len(scores.keys())
	
	# first pass: get sample means for each feature
	means = []
	for i in range(numCols):
		total = 0
		for id in scores:
			total += scores[id][i]
		means.append(float(total) / numRows)
	
	# second pass: calculate std devs for each feature
	sds = []
	for i in range(numCols):
		sumOfSq = 0
		for id in scores:
			sumOfSq += ((scores[id][i] - means[i]) ** 2)
		sds.append((sumOfSq / (numRows - 1)) ** 0.5)
	
	# third pass: calculate zscores
	for id in scores:
		scaled[id] = []
		for i in range(numCols):
			zscore = (scores[id][i] - means[i]) / sds[i]
			scaled[id].append(zscore)
	
	# final pass: error check
	for id in scaled:
		if len(scaled[id]) != len(scores[id]):
			print "Error: in scale_scores(): something went wrong with scaling."
	
	return scaled
	



def compute_dist_matrix(scores, method="euclidean"):
	"""
	Creates a dist matrix from a given scores (output from read_bitscore)
	Since this will be a symmetrical matrix, compute one triangle and then copy it
	method can be euclidean, pearson, or spearman
	"""
	distMatrix = []
	idList = scores.keys()
	
	# create one triangle of dist scores
	for i in range(0,len(idList)):
		if (i % 100) == 0:
			print "  finished dist for %s examples" % i
		distMatrix.append([])
		for j in range(0, len(idList)):
			if j > i:
				if method == "euclidean":
					dist = euclidean_dist(scores[idList[i]], scores[idList[j]])
				elif method == "pearson":
					dist = pearson_dist(scores[idList[i]], scores[idList[j]])
				elif method == "spearman":
					dist = pearson_dist(scores[idList[i]], scores[idList[j]])
				distMatrix[i].append(dist)
			else:
				distMatrix[i].append(0)
	
	# copy over triangle
	for i in range(0,len(idList)):
		for j in range(0, len(idList)):
			if j < i:
				distMatrix[i][j] = distMatrix[j][i]
	
	return distMatrix, idList
	
	


def knn_predict(distMatrix, realClasses, k):
	"""
	Find the k nearest neighbors for each example in distMatrix, predict class based on majority
	vote. Return list of predictions.
	* Note: current version gets rid of all exact matches
	"""
	predictions = []
	
	# for each row, sort by score, take top K hits, find max class
	for i in range(0,len(distMatrix)):
		current = list(distMatrix[i])
		if len(current) != len(distMatrix):
			print "j dimension not equal to i dimension"
		currentSorted = sorted(range(len(current)), key=lambda z: current[z]) #returns order of indices of sorted list. sort ascending.
		votes = {}
		neighborOrder = []
		for j in range(0,k):
			valid = False
			while valid == False: 
				if current[currentSorted[j]] == 0:
					popped = currentSorted.pop(j)
				else:
					valid = True
			
			classVote = realClasses[currentSorted[j]]
			neighborOrder.append(classVote)
			if classVote not in votes:
				votes[classVote] = 1
			else:
				votes[classVote] += 1
		
		# pick the best class based on majority vote
		# if tie, pick the class that was the nearest neighbor
		votesSorted = sorted(votes, key=votes.get, reverse=True) # sort descending.
		maxVote = votes[votesSorted[0]]
		candidates = []
		for cl in votes:
			if votes[cl] == maxVote:
				candidates.append(cl)
		if len(candidates) == 1:
			bestClass = candidates[0]
		elif len(candidates) < 1:
			print "In knn_predict(): Something wrong with picking best class."
			bestClass = "unknown"
		else:
			bestClass = ""
			bestRank = k + 1 #ensures this will be lower rank than any candidate
			for cl in candidates:
				rank = neighborOrder.index(cl)
				if rank < bestRank:
					bestRank = rank
					bestClass = cl
		
		predictions.append(bestClass) 
		

	return predictions




def rna_fold(seq, psfile = True):
	"""
	Folds the given sequence using RNAfold. (Assumes ViennaRNA is installed.)
	"seq" can not include a header line.
	"""
	import subprocess
	
	if psfile:
		command = "echo %s | RNAfold -p -noLP" % seq
	else:
		command = "echo %s | RNAfold -p -noLP -noPS" % seq
	
	p = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
	lines = p.stdout.readlines()
	strucLine = lines[1]
	retval = p.wait()
	
	mfe = structLine[(structLine.find(' (')+2):structLine.find(')\n')]
	struct = structLine[structLine.find('\n')+1:structLine.find(' ')]
	
	return struct, mfe



def bulk_fold(seqDict, outPath):
	"""
	Folds each sequence in the dictionary, renaming each output .ps file and moving
	it to the outPath so it is not overwritten
	"""
	return False
	
	
	

def align_fold(inputFile, outPath="tmp_mlocarna_output", otherOpts=""):
	"""
	Requires locARNA be installed.
	Aligns and folds seqs in given fasta file. Returns cons struct and mfe.
	"""
	import os, subprocess
	error = False
	
	if outPath == None:
		command = "mlocarna %s %s" % (inputFile, otherOpts)
	else:
		if not os.path.exists(outPath):
			os.makedirs(outPath)
		command = "mlocarna %s --tgtdir=%s %s" % (inputFile, outPath, otherOpts)
	p = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
	resultLines = p.stdout.readlines()
	for line in resultLines:
		if "alifold" in line:
			finalStruct = line
	retval = p.wait()
	if retval != 0:
		print "Error: problem while folding. Output:"
		for line in resultLines:
			print "    ", line
		error = True
	
	try:
		things = finalStruct.split()
	except UnboundLocalError:
		print "UnboundLocalError."
		things = ["","",""]
	struct = things[1]
	mfe = things[2][1:]
	
	return (mfe, struct, error)


def run_rnaz(alnFile, progPath="RNAz", progArgs=""):
	"""
	Run RNAz on provided alignment file and return MPI and SCI
	"""
	import subprocess
	error = False
	mpi = 0
	sci = 0
	
	command = "%s %s %s" % (progPath, alnFile, progArgs)
	p = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
	for line in p.stdout.readlines():
		if "Mean pairwise identity" in line:
			things = line.split()
			mpi = float(things[-1])
		if "Structure conservation index" in line:
			things = line.split()
			if things[-1] == "range.":
				print "  RNAz: 'range' found:", line
				sci = 0
			else:
				sci = float(things[-1])
	retval = p.wait()
	if retval != 0:
		print "Error: problem while running RNAz on", alnFile
		error = True
	
	return (mpi, sci, error)


def cm_calibrate_forecast(cmFile, cmcalibratePath="cmcalibrate"):
	"""
	Runs cmcalibrate --forecast on the given cm to estimate how long it will take to calibrate.
	Returns expected time to run
	"""
	import subprocess
	
	error = False
	timeToRun = 0
	command = "%s --forecast 1 %s" % (cmcalibratePath, cmFile)
	p = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
	for line in p.stdout.readlines():
		if "# all" == line[:5]:
			things = line.split()
			timeToRun = things[-1]
	retval = p.wait()
	if retval != 0:
		print "Error: problem while forecasting cm", cmFile
		error = True
	
	return (timeToRun, error)



def cm_calibrate(params):
	"""
	Runs cmcalibrate on the given cm. Unpacks a list of parameters (to faciliate multiprocessing).
	Returns elapsed time
	"""
	import subprocess, time
	
	cmFile = params[0]
	cmcalibratePath = params[1]
	otherOpts = params[2]
	error = False
	timeStarted = time.time()
	command = "%s %s %s" % (cmcalibratePath, cmFile, otherOpts)
	p = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
	retval = p.wait()
	if retval != 0:
		print "Error: problem while calibrating cm", cmFile
		error = True
	timeTaken = time.time() - timeStarted
	print "Finished calibrating %s in %s" % (cmFile, timeTaken)
	
	return (timeTaken, error)



def cm_search(params):
	"""
	Runs cmsearch on the given cm and db. Unpacks a list of parameters (to faciliate multiprocessing).
	"""
	import subprocess
	
	cmFile = params[0]
	dbFile = params[1]
	cmsearchPath = params[2]
	options = params[3]
	error = False
	print "Searching", cmFile
	command = "%s %s %s %s" % (cmsearchPath, options, cmFile, dbFile)
	p = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
	lines = p.stdout.readlines()
	retval = p.wait()
	if retval != 0:
		print "Error: in cm_search(): problem searching", cmFile
		print "Output:"
		for line in lines:
			print "  ", line
		error = True
	return error
	














