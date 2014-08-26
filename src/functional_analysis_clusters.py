#=================================================================================================
# Sarah Middleton
# 
# functional_analysis_clusters.py
# - Main output is <CLUSTER_FILE>.summary and <CLUSTER_FILE>.details
# - use -h for help
# 
# 
#=================================================================================================
import sys, random, time, os, subprocess, shutil
from optparse import OptionParser
from utils_file_readers import *
from utils_analysis import *

usageMsg = "Usage: %prog CLUSTER_FILE FASTA_DB BITSCORE_FILE [options]"
parser = OptionParser(usage=usageMsg)
parser.add_option("--maxSize", action="store", default=-1, dest="MAX_SIZE", help="skip any clusters larger than this size. by default does everything.")
parser.add_option("--cmA", action="store_true", default=False, dest="DO_CMA", help="get top cms for cluster using avg zscore")
parser.add_option("--cmP", action="store_true", default=False, dest="DO_CMP", help="get top cms for cluster using p-values (not as useful, in general)")
parser.add_option("--seqsOut", action="store", dest="FASTA_OUTPATH", help="print fasta seqs for the cluster to OUT_PATH")
parser.add_option("--struct", action="store", dest="STRUCT", help="get structure using locARNA in global mode using specified extra arguments. requires locARNA installation and --seqsOut")
parser.add_option("--cmbuild", action="store", dest="CM_OUT_PATH", help="build a cm for each cluster from the locARNA output. requires Infernal installation and --struct. This could instead be done using train_cms.py")
parser.add_option("--rnaz", action="store_true", default=False, dest="DO_RNAZ", help="use RNAz to get MPI, SCI, and z-score of cluster alignment. can only be used with --struct")
parser.add_option("--rnaz-path", action="store", default="", dest="RNAZ_LOC", help="location of rnaz executable, if not in default path")
parser.add_option("--locarna-path", action="store", default="", dest="LOCARNA_LOC", help="location of locarna executable, if not in default path")
parser.add_option("--aln", action="store_true", default=False, dest="PRINT_ALN", help="print cluster mult. alignment to the summary file?")
parser.add_option("--scale", action="store_true", default=False, dest="SCALE", help="Should the bitscore file be scaled when analyzing top CMs? (note that most bitscore files will have already been normalized, so this is not needed)")
#parser.add_option("--gene", action="store", dest="CHR_LOC", help="get overlapping genes using chrom locations")
#parser.add_option("--glistOut", action="store", dest="GLIST_OUTFILE", help="print genes in each cluster to indicated file for use with DAVID.")



# read/process args
(opts, args) = parser.parse_args()
#print "Options:", opts
#print "Positional args:", args
if len(args) < 2:
	print "Error: Must specify CLUSTER_FILE and FASTA_DB"
	print "Exiting"
	sys.exit(1)
if (opts.STRUCT != None) and (opts.FASTA_OUTPATH == None):
	print "Error: --struct and --seqsOut must be used together"
	print "Exiting."
	sys.exit(1)
if (opts.CM_OUT_PATH != None) and (opts.STRUCT == None):
	print "Error: --cmbuild and --struct must be used together"
	print "Exiting."
	sys.exit(1)
if (opts.DO_RNAZ != None) and (opts.STRUCT == None):
	print "Error: --rnaz and --struct must be used together"
	print "Exiting."
	sys.exit(1)
CLUSTER_FILE = args[0]
FASTA_DB = args[1]
BITSCORE_FILE = args[2]
RFAM_ANNOTS = "rfid_info_table.txt" # pre-made file; should be in src directory

# read in rfam annotations
rfamAnnots = {}
if ((opts.DO_CMA == True) or (opts.DO_CMP == True)) and os.path.exists(RFAM_ANNOTS):
	ins = open(RFAM_ANNOTS, 'r')
	ins.readline() # skip header
	for line in ins:
		lineParts = line.split()
		id = lineParts[0]
		shortName = lineParts[1]
		longName = lineParts[2]
		rfamAnnots[shortName] = longName
	ins.close()

# read in clusters
ins = open(CLUSTER_FILE, 'r')
header = ins.readline() 
lines = ins.readlines()
ins.close()

# open details output file
outFile1 = CLUSTER_FILE + ".details"
outs = open(outFile1, 'w')
numClusts = len(lines)
currClust = 0

# dict for summary
summary = {}

# analyze each cluster
for line in lines:
	currClust += 1
	print "Analyzing cluster %s of %s:" % (currClust, numClusts)
	
	line = line.rstrip('\r\n')
	things = line.split()
	cluster = things[-1]
	dissim = things[-2]
	seqCount = things[2]
	sigLvl = things[1]
	clustID = things[0]
	idList = cluster.split(",")
	
	if len(idList) < 2:
		print "  Cluster had less than 2 sequences. Skipping."
		continue
	
	summary[clustID] = {}
	summary[clustID]['sig'] = sigLvl
	summary[clustID]['size'] = seqCount
	summary[clustID]['members'] = cluster
	summary[clustID]['dissim'] = dissim
	summary[clustID]['mfe'] = "-"
	summary[clustID]['mpi'] = "-"
	summary[clustID]['sci'] = "-"
	summary[clustID]['ratio'] = "-" 
	summary[clustID]['covar'] = "-"
	summary[clustID]['zscore'] = "-"
	summary[clustID]['topCMs'] = []
	
	outs.write(header)
	outs.write(line)
	outs.write("\n")
	
	
	# Get highest scoring CMs by avg zscore
	if opts.DO_CMA == True:
		print "  Getting best CMs by z-score avg (--cmA)..."
		outs.write("Best CMs (by avg z-score):\n")
		cmAvgs = group_cm_zAvg(idList, BITSCORE_FILE, opts.SCALE)
		cmRanks = sorted(cmAvgs, key=cmAvgs.get)
		cmRanks.reverse()
		for i in range(10):
			cmName = cmRanks[i]
			cmNameLong = ""
			if cmName in rfamAnnots:
				cmNameLong = rfamAnnots[cmName]
			outStr = "  %2s. (Z = %.2f) %s - %s\n" % ((i + 1), cmAvgs[cmName], cmName, cmNameLong)
			outStr2 = "%s,%.3f" % (cmName, cmAvgs[cmName])
			outs.write(outStr)
			summary[clustID]['topCMs'].append(outStr2)
		outs.write("\n")
	else:
		summary[clustID]['topCMs'].append("-")
	
	# Get highest scoring CMs by pval
	if opts.DO_CMP == True:
		print "  Getting best CMs by pval (--cmP)..."
		outs.write("Best CMs (by 1-sided pvalue):\n")
		cmPvals = group_cm_pVal(idList, BITSCORE_FILE)
		cmRanks = sorted(cmPvals, key=cmPvals.get)
		for i in range(10):
			cmName = cmRanks[i]
			outStr = "  %s. %s (%s)\n" % ((i + 1), cmName, cmPvals[cmName])
			outs.write(outStr)
		outs.write("\n")		
	
	# print original fasta format seqs for each cluster (prints to NEW file. this is needed for folding.)
	if opts.FASTA_OUTPATH != None:
		if not os.path.exists(opts.FASTA_OUTPATH):
			os.makedirs(opts.FASTA_OUTPATH)
		print "  Printing sequences to %s (--seqsOut)..." % opts.FASTA_OUTPATH
		fastaOut = "%sclust%s.fa" % (opts.FASTA_OUTPATH, clustID)
		outs2 = open(fastaOut, 'w')
		seqs = get_sequences(idList, FASTA_DB, cutID=25)
		for id in seqs:
			outStr = ">%s\n%s\n" % (id, seqs[id])
			outs2.write(outStr)
		outs2.write("\n")
		outs2.close()
	
	# Get multiple alignment and structure
	if (opts.STRUCT != None):
		print "  Aligning and folding using mlocarna (--struct)..."
		outs.write("Structure and mfe from mlocarna:\n")
		structOutputFolder = "%sclust%s_global.out/" % (opts.FASTA_OUTPATH, clustID)
		inputFile = fastaOut #from above
		otherOpts = opts.STRUCT
		command = "%smlocarna %s --tgtdir=%s %s" % (opts.LOCARNA_LOC, inputFile, structOutputFolder, otherOpts)
		p = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
		resultLines = p.stdout.readlines()
		passed = True
		finalStruct = ""
		totalMFE = "NA"
		covContribution = "NA"
		structPart = ""
		for rline in resultLines:
			if "alifold" in rline:
				finalStruct = rline
				things = finalStruct.split()
				structPart = things[1]
				totalMFE = things[-5].lstrip("(")
				covContribution = things[-1].rstrip(")")
				if float(totalMFE) > -1:
					passed = False
		retval = p.wait()
		if retval != 0:
			print "Error: problem while folding. Command tried:", command
			print "Exiting."
			sys.exit(1)
		if (passed == False) and ("--noLP" not in otherOpts):
			# try folding with --LP
			print "   > Poor folding, trying with --LP..."
			command = "%smlocarna %s --tgtdir=%s %s --LP" % (opts.LOCARNA_LOC, inputFile, structOutputFolder, otherOpts)
			p = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
			resultLines = p.stdout.readlines()
			for rline in resultLines:
				if "alifold" in rline:
					finalStruct = rline
					things = finalStruct.split()
					structPart = things[1]
					totalMFE = things[-5].lstrip("(")
					covContribution = things[-1].rstrip(")")
			retval = p.wait()
			if retval != 0:
				print "Error: problem while folding. Command tried:", command
				print "Exiting."
				sys.exit(1)
		# print this struct to a file in the structOutputFolder
		structOutputFile = "%scons_struct.txt" % structOutputFolder
		outs2 = open(structOutputFile, 'w')
		outs2.write(finalStruct)
		outs2.close()
		outs.write("  Command: " + command + "\n")
		outs.write("  " + structPart + "\n")
		outStr = "  MFE = %s (Covariance contribution: %s)\n" % (totalMFE, covContribution)
		outs.write(outStr)
		
		# make a copy of the resulting structure .ps file with a proper name and in the parent folder for easier access
		pathToPic = "%sresults/alirna.ps" % structOutputFolder
		newPathToPic = "%scons_struc_pics/" % opts.FASTA_OUTPATH
		newPicFilename = "%sclust%s.ps" % (newPathToPic, clustID)
		if not os.path.exists(newPathToPic):
			os.makedirs(newPathToPic)
		shutil.copyfile(pathToPic, newPicFilename)
		outs.write("\n")
		summary[clustID]['mfe'] = "%.2f" % float(totalMFE)
		summary[clustID]['covar'] = "%.2f" % float(covContribution)
	
	# Get RNAz stats for multiple alignment
	if opts.DO_RNAZ == True:
		if (opts.MAX_SIZE == -1) or (int(seqCount) <= int(opts.MAX_SIZE)):
			print "  Getting RNAz stats (--rnaz)..."
			outs.write("RNAz stats:\n")
			structOutputFolder = "%sclust%s_global.out/" % (opts.FASTA_OUTPATH, clustID)
			alnFile = "%sresults/result.aln" % structOutputFolder
			command = "%sRNAz %s" % (opts.RNAZ_LOC, alnFile)
			p = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
			resultLines = p.stdout.readlines()
			(mpi, sci, zscore, ratio) = ("NA","NA","NA","NA")
			if len(resultLines) > 21:
				if "Mean pairwise identity" in resultLines[6]:
					things = resultLines[6].split()
					mpi = things[-1]
					try:
						mpi = "%.2f" % float(mpi)
					except ValueError:
						mpi = "NF"
				if "Structure conservation index" in resultLines[15]:
					things = resultLines[15].split()
					sci = things[-1]
					try:
						sci = "%.2f" % float(sci)
					except ValueError:
						sci = "NF"
				if "Mean z-score" in resultLines[14]:
					things = resultLines[14].split()
					zscore = things[-1]
					try:
						zscore = float(zscore)
					except ValueError:
						zscore = "NF"
				outStr = "  MPI: %s, SCI: %s, ZSCORE: %s\n" % (mpi, sci, zscore)
				outs.write(outStr)
				outs.write("  Full summary:\n")
				for i in range(3, 21):
					outs.write("  " + resultLines[i])
			else:
				for rline in resultLines:
					print rline
			retval = p.wait()
			if retval != 0:
				print "Error: problem while running RNAz. Command tried:", command
				print "Exiting."
				sys.exit(1)
			outs.write("\n")
			summary[clustID]['mpi'] = mpi
			summary[clustID]['sci'] = sci
			summary[clustID]['zscore'] = zscore
			
			try:
				ratio = "%.2f" % ((float(sci) * 100) / float(mpi))
			except ValueError:
				ratio = "NA"
			summary[clustID]['ratio'] = ratio
			
	
	# Print multiple alignment in the summary file
	if opts.PRINT_ALN == True:
		print "  Printing multiple alignment (--aln)..."
		outs.write("Multiple alignment from mlocarna:\n")
		structOutputFolder = "%sclust%s_global.out/" % (opts.FASTA_OUTPATH, clustID)
		alnFile = "%sresults/result.aln" % structOutputFolder
		(seqs, error) = read_clustalw(alnFile)
		for id in seqs:
			outStr = "  %s\t%s\n" % (id, seqs[id])
			outs.write(outStr)
		outs.write("\n")
	
	# Train new CM for this cluster
	if (opts.CM_OUT_PATH != None):
		print "  Building cm from structural alignment (--cmbuild)"
		if not os.path.exists(opts.CM_OUT_PATH):
			os.makedirs(opts.CM_OUT_PATH)
		cmStoFile = "%sclust%s.sto" % (opts.CM_OUT_PATH, clustID)
		cmOutFile = "%sclust%s.cm" % (opts.CM_OUT_PATH, clustID)
		mlocarna_to_stockholm(structOutputFolder, cmStoFile) #creates stockholm format file for this cluster from mlocarna results
		command = "cmbuild -F %s %s" % (cmOutFile, cmStoFile)
		p = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
		resultLines = p.stdout.readlines()
		retval = p.wait()
		if retval != 0:
			print "Error: problem while building cm."
			print "cmbuild output:"
			for rline in resultLines:
				print "  ", rline
			print ""
			outs.write("Note: this caused a CMbuild error.\n")
		
		
	
	outs.write("//")
	outs.write("\n\n\n\n\n\n\n")
	outs.flush()


outs.close()

# open summary output file
outFile2 = CLUSTER_FILE + ".summary"
outs = open(outFile2, 'w')
outs.write("ClustID\tInfo\tCount\tMFE\tCovContr\tMPI\tSCI\tSCI/MPI\tZscore\tAvgDissim\tTopCMs\tMembers\n")
for id in summary:
	outStr = "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (id, summary[id]['sig'], summary[id]['size'], summary[id]['mfe'], summary[id]['covar'], summary[id]['mpi'], summary[id]['sci'], summary[id]['ratio'] , summary[id]['zscore'], summary[id]['dissim'], ";".join(summary[id]['topCMs']), summary[id]['members'])
	outs.write(outStr)
outs.close()


'''
# print DAVID input file
if opts.GLIST_OUTFILE != None:
	output = {}
	for line in lines:
		things = line.split()
		cid = "clust" + things[0]
		gidList = things[-1].split(",")
		output[cid] = []
		for fullID in gidList:
			(dbname, gid) = fullID.split("_", 2)
			if gid not in output[cid]:
				output[cid].append(gid)
	
	outs = open(opts.GLIST_OUTFILE, 'w')
	header = "\t".join(sorted(output.keys()))
	outs.write(header + "\n")
	maxLen = 0
	for cid in output:
		if len(output[cid]) > maxLen:
			maxLen = len(output[cid])
	for i in range(maxLen):
		for cid in sorted(output.keys()):
			try: gid = output[cid][i] + "\t"
			except IndexError: gid = "\t"
			outs.write(gid)
		outs.write("\n")
	outs.close()

'''

print "Details printed to", outFile1
print "Summary printed to", outFile2


