#=================================================================================================
# Sarah Middleton
# 
# expand_clusters.py
# Expands clusters using original database.
# - given a list of top clusters, the original fasta db, a "bg"  fasta db (optional), and a cluster overlap threshold for merging:
# - searches each CM from the top clusters on the db
# - parses output of search
# - checks if original seqs were found by the search
# - (if bg db provided) tests for enrichment of the non-bg seqs
# - checks for overlapping clusters, merges if found
# - prints results in cluster format
# - use -h for help
# - note: assumes CLUSTER_FILE is has been annotated and turned into CMs, and that CLUSTER_FILE is
#   in the directory where the rest of the annotation results are stored. If you are using --ethresh,
#   assumes you have calibrated these CMs and they are names *.c.cm.
# 
# Required args: 
#   CLUSTER_FILE = any output file from filter_clusters2.py
#   FASTA_DB = original fasta database for the clustered sequences
#=================================================================================================
if __name__ == '__main__':
	import sys, random, time, os, subprocess, glob
	from optparse import OptionParser
	from multiprocessing import Pool, cpu_count
	from utils_file_readers import *
	from utils_analysis import *
	
	usageMsg = "Usage: %prog CLUSTER_FILE FASTA_DB [options]"
	parser = OptionParser(usage=usageMsg)
	parser.add_option("--dist-matrix", action="store", dest="DIST_MATRIX", help="If you want to calculate the avg. dissim. in the new clusters, supply a dissim matrix for all seqs.")
	parser.add_option("--bg", action="store", dest="BG_FASTA", help="requires --bg-id as well. a fasta db containing background sequences to search in addition to original fasta db. this is required for doing enrichment analysis.")
	parser.add_option("--bg-id", action="store", dest="BG_ID", help="group id of the bg seqs. must be first part of the name for each bg seq.")
	parser.add_option("--merge", action="store", dest="MERGE_THRESH", help="if this is specified, clusters overlapping by the specified fraction (e.g. '0.5') will be merged into one cluster in output.")
	parser.add_option("--ethresh", action="store", default=1, dest="E_THRESH", help="E-value threshold to use with cmsearch. Default is 1")
	parser.add_option("--bit-thresh", action="store", default=None, dest="BIT_THRESH", help="Bitscore threshold to use with cmsearch. This can be used if CMs were not calibrated. If specified, this overrides --ethresh.")
	parser.add_option("--glocal", action="store_true", default=False, dest="DO_GLOCAL", help="Use 'glocal' alignment for search (global w/ respect to CM, local w/ respect to seq)")
	parser.add_option("--qc", action="store_true", default=False, dest="DO_QC", help="When all sequences are of known groups, sens/spec/acc can be computed and output to the cluster file. GroupID must be first part of ID (e.g.: GROUP_gene123_3")
	parser.add_option("--skip-search", action="store_true", default=False, dest="SKIP_SEARCH", help="For when you've already done a search and just want to run different options.")
	parser.add_option("--cpu", action="store", default=1, dest="NUM_CPU", help="Number of cpu to use. Default is 1.")
	parser.add_option("--infernal-path", action="store", default="", dest="INFERNAL_PATH", help="Absolute path to infernal executables, if not in default path.")
	
	
	# read/process args
	(opts, args) = parser.parse_args()
	if len(args) == 2:
		CLUSTER_FILE = args[0]
		FASTA_DB = args[1]
	else:
		print "Required args: CLUSTER_FILE FASTA_DB"
		print "Use -h for help. Exiting."
		sys.exit(1)
	
	CMSEARCH_PATH = "%scmsearch" % opts.INFERNAL_PATH
	TIMESTAMP = int(time.time())
	tmpFisher = "tmp/tmp_fisher_%s.txt" % TIMESTAMP
	bgID = "None"
	if opts.BG_ID != None:
		bgID = opts.BG_ID
	gloc = "Loc"
	if opts.DO_GLOCAL == True:
		gloc = "Gloc"
	if opts.BIT_THRESH != None:
		thresh = opts.BIT_THRESH
		marker = "bs"
		FN_END = -3
	else:
		thresh = opts.E_THRESH
		marker = "e"
		FN_END = -5
	if not os.path.exists("tmp"):
		os.makedirs("tmp")
	RUNSTAMP = "%s%sbg%s%s" % (marker, thresh, bgID, gloc)
	
	# read dist matrix
	if opts.DIST_MATRIX != None:
		(idList, dissimMatrix, error) = read_comp_matrix(opts.DIST_MATRIX)
	
	# read cluster file
	clusters = {}
	ins = open(CLUSTER_FILE, 'r')
	ins.readline() #skip header
	lines = ins.readlines()
	ins.close()
	for line in lines:
		things = line.split()
		cID = things[0]
		members = things[-1]
		clusters[cID] = members.split(",")
	
	# get list of cms
	fn = os.path.basename(CLUSTER_FILE)
	fnParts = fn.split(".")
	dbName = fnParts[0]
	workingDir = os.path.dirname(CLUSTER_FILE)
	cmDir = "%s/%s_clust_cms" % (workingDir, dbName)
	cmList = []
	if opts.BIT_THRESH != None: #use non-calibrated CMs
		cmExt = cmDir + "/*.cm"
		fileList = glob.glob(cmExt)
		for name in fileList:
			base = os.path.basename(name)
			if ".c." not in base:
				cmID = base[5:FN_END] # [5,-3]
				if cmID in clusters:
					cmList.append(name)
	else: #use only calibrated CMs
		cmExt = cmDir + "/*.c.cm"
		fileList = glob.glob(cmExt)
		for name in fileList:
			base = os.path.basename(name)
			cmID = base[5:FN_END] #[5,-5]
			if cmID in clusters:
				cmList.append(name)
	
	
	
	# if bg db provided, merge with FASTA_DB
	if opts.BG_FASTA != None:
		mergedFilename = "%s_withBG.db" % FASTA_DB
		dbToUse = mergedFilename
		error = merge_files(FASTA_DB, opts.BG_FASTA, mergedFilename, header=False)
		if error:
			print "Could not merge target database with bg database. Exiting."
			sys.exit(1)
	else:
		dbToUse = FASTA_DB
	
	# for each cm, search the database
	execList = []
	print "Predicted time to search:"
	for i in range(len(cmList)):
		command = "%s --forecast 1 %s %s" % (CMSEARCH_PATH, cmList[i], dbToUse)
		p = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
		for line in p.stdout.readlines():
			if "  all" == line[:5]:
				things = line.split()
				timeToRun = things[-1]
				print "%s: %s" % (cmList[i], timeToRun)
		retval = p.wait()
		
		newFolder = "%s/cmsearch_%s" % (os.path.dirname(cmList[i]), RUNSTAMP)
		if not os.path.exists(newFolder):
			os.makedirs(newFolder)
		cmName = os.path.basename(cmList[i][:FN_END])
		alnOut = "%s/%s_aln.txt" % (newFolder, cmName)
		tabOut = "%s/%s_tab.txt" % (newFolder, cmName)
		gcOut = "%s/%s_gc.txt" % (newFolder, cmName)
		if opts.BIT_THRESH != None:
			threshType = "-T"
		else:
			threshType = "-E"
		if opts.DO_GLOCAL == True:
			gloc = "-g"
		else:
			gloc = ""
		
		options = "--toponly --informat FASTA %s %s -o %s --tabfile %s --gcfile %s %s" % (threshType, thresh, alnOut, tabOut, gcOut, gloc)
		params = [cmList[i], dbToUse, CMSEARCH_PATH, options]
		execList.append(params)
	print ""
	
	if not opts.SKIP_SEARCH:
		# create a pool of processes to run CM searching in parallel
		print "Creating %s processes to split %s jobs." % (int(opts.NUM_CPU), len(cmList))
		print ""
		pool = Pool(processes = int(opts.NUM_CPU))
		searchStart = time.time()
		result = pool.map_async(cm_search, execList)
		ret = result.wait()
		#print "Results: (false = no error)"
		#print result.get()
		elapsedTime = time.time() - searchStart
		
	
	# parse _tab output for each cm and create new cluster file
	print "Processing output and printing new cluster file"
	newClustFile = "%s_expanded_%s.txt" % (CLUSTER_FILE, RUNSTAMP)
	newClusters = {}
	outs = open(newClustFile, 'w')
	outs.write("ClustID\tSigLvl\tCount\tAvgDissim\tMembers\n")
	for cmfile in cmList:
		justFilename = os.path.basename(cmfile)
		cID = justFilename[5:FN_END]
		if cID not in clusters:
			print "cID not in orig cluster file:", cID
		else:
			tabFileFolder = "%s/cmsearch_%s" % (os.path.dirname(cmfile), RUNSTAMP)
			tabFile = "%s/%s_tab.txt" % (tabFileFolder, justFilename[:FN_END])
			(results, error) = read_cmsearch_tabfile(tabFile)
			if error:
				print "Error: problem reading tabfile:", tabFile
				print "Exiting."
				sys.exit(1)
			
			# check if all original seqs are in new cluster
			notFound = 0
			for id in clusters[cID]:
				if id not in results:
					#results[id] = {}
					print "[clust%s] Warning: orig id %s not found in expanded cluster." % (cID, id)
					notFound += 1
			if notFound > 0:
				print "[clust%s] %s / %s orig ids were not recovered in search.\n" % (cID, notFound, len(clusters[cID]))
			
			# get cluster info and print
			expIdList = results.keys()
			expIdStr = ",".join(expIdList)
			if opts.DIST_MATRIX == None:
				outStr = "%s\t%s\t%s\t%s\t%s\n" % (cID, "-", len(expIdList), "-", expIdStr)
			else:
				avgDist = get_avg_dist(expIdList, dissimMatrix)
				outStr = "%s\t%s\t%s\t%.2f\t%s\n" % (cID, "-", len(expIdList), avgDist, expIdStr)
			outs.write(outStr)	
			newClusters[cID] = expIdList[:]
	outs.close()
	print ""
	print "New cluster file printed to", newClustFile
	
	
	# also output a "merged" cluster file, if --merge option is used
	# "chain merging" is allowed to occur (where two non overlapping clusters are merged due to common overlapping cluster)
	if opts.MERGE_THRESH != None:
		print "Merging highly overlapping clusters"
		mergedOutfile = "%s_expanded_merged_%s.txt" % (CLUSTER_FILE, RUNSTAMP)
		cids = newClusters.keys()
		numClust = len(newClusters.keys())
		mergeList = {}
		
		# create 'graph' of connections between clusters based on overlap
		# when comparing two clusters, if either overlaps the other by more than the thresh, create an edge between them.
		for i in range(numClust):
			mergeList[cids[i]] = []
		for i in range(numClust):
			for j in range(i, numClust): #note: this includes comparing the cluster to itself
				cid1 = cids[i]
				cid2 = cids[j]
				(num, frac1, frac2) = get_idList_overlap(newClusters[cid1], newClusters[cid2])
				if (frac1 >= float(opts.MERGE_THRESH)) or (frac2 >= float(opts.MERGE_THRESH)): 
					if cid2 not in mergeList[cid1]:
						mergeList[cid1].append(cid2)
					if cid1 not in mergeList[cid2]:
						mergeList[cid2].append(cid1)
		
		# do a BFS of the graph to find connected components. these will be the merged clusters.
		mergedClusters = {}
		usedList = {} #so we don't repeat ourselves
		for cid in cids:
			usedList[cid] = False
		newCid = 1
		for id in mergeList: #do a BFS of the network
			if usedList[id] == False:
				usedList[id] = True
				queue = mergeList[id][:] #start a queue containing the clusters connected to this cluster
				front = 0
				end = len(queue)
				while front < end:
					currentId = queue[front]
					for subid in mergeList[currentId]:
						if subid not in queue:
							queue.append(subid)
							end += 1
					usedList[currentId] = True
					front += 1
				mergedClusters[("M" + str(newCid))] = queue[:]
				newCid += 1
		
		outs = open(mergedOutfile, 'w')
		outs.write("ClustID\tMergedIDs\tCount\tAvgDissim\tMembers\n")
		mergedClustersList = mergedClusters.keys()
		for newcid in mergedClustersList:
			print "Merging the following clusters into [clust%s]:" % newcid,
			finalIdList = []
			for oldcid in mergedClusters[newcid]:
				print "%s," % oldcid,
				for seqID in newClusters[oldcid]:
					if seqID not in finalIdList:
						finalIdList.append(seqID)
			print ""
			
			if len(finalIdList) == 0:
				print "  > Removing this cluster, as it has no sequences."
				del mergedClusters[newcid]
			else:
				finalIdStr = ",".join(finalIdList)
				mergedIdStr = ",".join(mergedClusters[newcid])
				if opts.BG_FASTA != None:
					outStr = "%s\t%s\t%s\t%s\t%s\n" % (cID, "-", len(expIdList), "-", expIdStr)
				else:
					avgDist = get_avg_dist(finalIdList, dissimMatrix)
					outStr = "%s\t%s\t%s\t%.2f\t%s\n" % (newcid, mergedIdStr, len(finalIdList), avgDist, finalIdStr)
				outs.write(outStr)	
				mergedClusters[newcid] = finalIdList
		
		outs.close()
			
				
	
	# also output a "QC" cluster file with sens/spec/acc if the --qc option is used (including a merged one if specified)
	if opts.DO_QC == True:
		if opts.BG_FASTA != None:
			print "Can not do QC with bg fasta file."
		else:
			if opts.MERGE_THRESH != None:
				# get group totals
				allGroupsTotal = len(idList)
				groupTotals = {}
				for id in idList:
					idParts = id.split("_")
					groupID = idParts[0]
					if groupID not in groupTotals:
						groupTotals[groupID] = 0
					groupTotals[groupID] += 1
				
				mergedQCOut = "%s_QC_expanded_merged_%s.txt" % (CLUSTER_FILE, RUNSTAMP)
				outs = open(mergedQCOut, 'w')
				outs.write("ClustID\tCount\tTopGroup\tSens\tSpec\tF1\tFDR\tMembers\n")
				for newcid in mergedClusters:
					seqIdList = mergedClusters[newcid]
					seqIdStr = ",".join(seqIdList)
					clusterSize = len(seqIdList)
					(topGroup, tp) = get_cluster_group(seqIdList)
					fp = clusterSize - tp
					tn = (allGroupsTotal - groupTotals[topGroup]) - fp
					fn = groupTotals[topGroup] - tp
					sens = float(tp) / (tp + fn)
					spec = float(tn) / (fp + tn)
					f1 = float(2 * tp) / (2 * tp + fp + fn)
					fdr = float(fp) / (fp + tp)
					
					outStr = "%s\t%s\t%s\t%.3f\t%.3f\t%.3f\t%.3f\t%s\n" % (newcid, clusterSize, topGroup, sens, spec, f1, fdr, seqIdStr)
					outs.write(outStr)
				outs.close()
	
	
	
	# checks for enrichment of each cluster for target seqs if --bg is provided
	# DOES consider original cluster seqs in calculation
	if opts.BG_ID != None:
		# get fg/bg totals
		fgTotal = 0
		bgTotal = 0
		(seqs, error) = read_fasta(dbToUse)
		if error:
			print "Error: problem reading fasta db", dbToUse
			print "Exiting."
			sys.exit(1) 
		for id in seqs:
			idParts = id.split("_")
			if idParts[0] == opts.BG_ID:
				bgTotal += 1
			else:
				fgTotal += 1
		if bgTotal == 0:
			print "Warning: no bg seqs found."
		if fgTotal == 0:
			print "Warning: no fg seqs found."
		
		enrOut = "%s_expanded_enrich_%s.txt" % (CLUSTER_FILE, RUNSTAMP)
		outs = open(enrOut, 'w')
		outs.write("ClustID\tCount\tTP\tFP\tSig\tEnrichP\tMembers\n")
		for cid in newClusters:
			newIdList = []
			seqIdList = newClusters[cid]
			for id in seqIdList:
				#if id not in clusters[cid]:
				newIdList.append(id)
			
			clusterSize = len(newIdList)
			#oldClusterSize = len(clusters[cid])
			sortedNewIdList = sorted(newIdList)
			seqIdStr = ",".join(sortedNewIdList)
			tp = 0
			fp = 0
			for id in newIdList:
				idParts = id.split("_")
				if idParts[0] == opts.BG_ID:
					fp += 1
				else:
					tp += 1
			#fn = (fgTotal - oldClusterSize) - tp
			fn = fgTotal - tp
			tn = bgTotal - fp
			command = "Rscript get_enrichment.r %s %s %s %s %s %s" % (tp, fp, fn, tn, "greater", tmpFisher)
			p = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
			retval = p.wait()
			ins = open(tmpFisher, 'r')
			pval = ins.readline()
			pval = float(pval.rstrip('\r\n'))
			ins.close()
			
			sig = "-"
			if pval <= 0.05:
				sig = "*"
			if pval <= 0.01:
				sig = "**"
			if pval <= 0.001:
				sig = "***"
			
			tpp, fpp = 0, 0
			if clusterSize > 0:
				tpp = "%.2f" % (float(tp) / clusterSize)
				fpp = "%.2f" % (float(fp) / clusterSize)
			outStr = "%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (cid, clusterSize, tpp, fpp, sig, pval, seqIdStr)
			outs.write(outStr)
		outs.close()
		
		if opts.MERGE_THRESH != None:
			mergedEnrOut = "%s_expanded_merged_enrich_%s.txt" % (CLUSTER_FILE, RUNSTAMP)
			outs = open(mergedEnrOut, 'w')
			outs.write("ClustID\tCount\tTP\tFP\tSig\tEnrichP\tMembers\n")
			for newcid in mergedClusters:
				seqIdList = mergedClusters[newcid]
				clusterSize = len(seqIdList)
				seqIdStr = ",".join(seqIdList)
				tp = 0
				fp = 0
				for id in seqIdList:
					idParts = id.split("_")
					if idParts[0] == opts.BG_ID:
						fp += 1
					else:
						tp += 1
				fn = fgTotal - tp
				tn = bgTotal - fp
				command = "Rscript get_enrichment.r %s %s %s %s %s %s" % (tp, fp, fn, tn, "greater", tmpFisher)
				p = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
				retval = p.wait()
				ins = open(tmpFisher, 'r')
				pval = ins.readline()
				pval = float(pval.rstrip('\r\n'))
				ins.close()
				
				sig = "-"
				if pval <= 0.05:
					sig = "*"
				if pval <= 0.01:
					sig = "**"
				if pval <= 0.001:
					sig = "***"
				
				tpp, fpp = 0, 0
				if clusterSize > 0:
					tpp = "%.2f" % (float(tp) / clusterSize)
					fpp = "%.2f" % (float(fp) / clusterSize)
				outStr = "%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (newcid, clusterSize, tpp, fpp, sig, pval, seqIdStr)
				outs.write(outStr)
			outs.close()
	
	
	
	
