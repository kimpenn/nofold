#====================================================================
# Sarah Middleton
# 
# score_and_normalize.py
# Score a database in FASTA format for structures using Infernal.
# Splits the job into the indicated number of processes. 
# Raw score files are output to a folder called cmscore_results_rfam that 
# is created within the same folder as the fasta file.
# Score files are consolidated into a single .bitscore file and normalized.
# ==> Use *.zNorm.pcNorm100.zNorm.bitscore for further analysis/clustering.
#
# Required args: a fasta file of sequences
# Usage examples: 
#    python score_and_normalize.py ../demo/demo1/demo1.db --cpus=4
#    python score_and_normalize.py ../demo/demo1/demo1.db --cpus=4 --infernal-path=../../../TOOLS/Infernal-1.0.2/infernal-1.0.2/src/
#====================================================================
import subprocess, sys, os, time, glob
import utils_file_readers as ufr
from multiprocessing import Pool, cpu_count
from optparse import OptionParser

#------------------------------------------------------------------
# Runs the db through Infernal and scores with the given CM.
# Takes a list of params and unpacks them.
#------------------------------------------------------------------
def score_cm(params):
	seqDB = params[0]
	cmFile = params[1]
	outFile = params[2]
	progPath = params[3]
	command = "%s --search -a --infile %s %s > %s" % (progPath, seqDB, cmFile, outFile)
	
	start = time.time() # time the scoring for this CM
	job = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
	jobOutput = job.stdout.read()
	if "Error" in jobOutput:
		print "  ", jobOutput
	result = job.wait()
	totalTime = time.time() - start
	
	if result == 0:
		res = "SUCCESS"
	else:
		res = "*****ERROR*****"
	print "%s -- result: %s -- time: %.3f" % (cmFile, res, totalTime)
		
	return totalTime
	
#------------------------------------------------------------------
# Main method
# Spawns the indicated number of processes and splits the CM
# scoring jobs among them.
#------------------------------------------------------------------
if __name__ == '__main__':
	
	usageMsg = "Usage: %prog FASTA.db [options]"
	parser = OptionParser(usage=usageMsg)
	parser.add_option("--cpus", action="store", type='int', default=1, dest="MAX_CPU", help="Maximum number of CPUs to use. Default is [%default].")
	parser.add_option("--infernal-path", action="store", default=None, dest="INFERNAL_PATH", help="Path to folder containing Infernal executables. Default is to assume they have been added to your PATH and can be called by name.")
	parser.add_option("--cm-folder", action="store", default="../models/rfam_cms/", dest="CM_FOLDER", help="Folder containing the covariance models you would like to use. Default is to use the 1,973 Rfam CMs included in [%default].")
	
	
	# read/process args
	(opts, args) = parser.parse_args()
	scriptDir = os.path.dirname(sys.argv[0])
	if scriptDir != "":
		scriptDir += "/"
		
	if len(args) == 1:
		DB_FILE = os.path.abspath(os.path.expanduser(args[0]))
		NUM_PROCESSES = int(opts.MAX_CPU)
		if opts.INFERNAL_PATH == None:
			PROG_PATH = "cmscore"
		else:
			PROG_PATH = os.path.expanduser(opts.INFERNAL_PATH).rstrip('/') + "/cmscore"
		if opts.CM_FOLDER == "../models/rfam_cms/":
			CM_FOLDER = os.path.abspath("%s%s" % (scriptDir, opts.CM_FOLDER))
		else:
			CM_FOLDER = os.path.abspath(os.path.expanduser(opts.CM_FOLDER))
	else:
		print "Incorrect args."
		print "Use -h for help. Exiting."
		sys.exit()
	
		
	# set up additional file names
	workingDir = os.path.dirname(DB_FILE)
	filename = os.path.basename(DB_FILE)
	fnParts = filename.split(".")
	dbName = fnParts[0]
	scoresOut = "%s/cmscore_results_rfam/" % (workingDir)
	bitOut = "%s/%s.bitscore" % (workingDir, dbName)
	
	# check that everything exists
	if not os.path.exists(DB_FILE):
		print ">> Error: could not find indicated fasta file."
		print "   (tried: %s)" % DB_FILE
		print ">> Exiting."
		sys.exit()
	if not os.path.exists(CM_FOLDER):
		print ">> Error: could not find CM folder."
		print "   (tried: %s)" % CM_FOLDER
		print ">> Exiting."
		sys.exit()
	if (opts.INFERNAL_PATH != None) and (not os.path.exists(PROG_PATH)):
		print ">> Error: could not find cmscore."
		print "   (tried: %s)" % PROG_PATH 
		print "   Ensure that you have the correct version of Infernal (1.0.2) and that either cmscore is added to your PATH or you have supplied the correct path to the folder containing the executable."
		print ">> Exiting."
		sys.exit()
		
	# check number of sequences
	(seqs, error) = ufr.read_fasta(DB_FILE)
	if error:
		print ">> Error: could not read fasta file. Check that it is in the correct format (see README)."
		print ">> Exiting."
		sys.exit()
	numSeqs = len(seqs.keys())
		
	print ""
	print "Files/paths to be used:\n"
	print "    output directory:", workingDir
	print "       sequence file:", DB_FILE
	print "           cm folder:", CM_FOLDER
	print "     cmscore command:", PROG_PATH
	print ""
		
	
	# create output directory if necessary
	if os.path.exists(scoresOut):
		print ""
		print "Output directory", scoresOut, "already exists."
		response = raw_input("Ok to overwrite existing CM output files? (y/n) ")
		if response != "y":
			print "Exiting."
			sys.exit()
	else:
		print "Output directory", scoresOut, "does not exist, creating."
		os.makedirs(scoresOut)
	
	# get list of CMs, create list of parameters for splitting jobs
	execList = []
	count = 0
	cmList = glob.glob(CM_FOLDER + "/*.cm")
	for cmFile in cmList:
		count += 1
		cmName = os.path.basename(cmFile)
		outFile = scoresOut + "scores_" + cmName + ".txt"		
		params = [DB_FILE, cmFile, outFile, PROG_PATH]
		execList.append(params)

	print ""
	print "======================================================"
	print "Beginning scoring process."
	print NUM_PROCESSES, "processes will be created."
	print count, "jobs will be assigned to the process pool."
	print ""
	
	# create a pool of processes to run CM scoring in parallel
	pool = Pool(processes = NUM_PROCESSES)
	timeStart = time.time()
	result = pool.map_async(score_cm, execList)
	result.wait()
	elapsedTime = time.time() - timeStart
	
	
	# calculate avg time
	timings = result.get()
	totalTime = 0
	for t in timings:
		totalTime += t
	avgTime = float(totalTime) / count
	cpuHrs = (float(totalTime) / 60) / 60
	elapsedHrs = (float(elapsedTime) / 60) / 60
	
	
	print ""
	print "Finished scoring %s sequences on %s CMs in %.2fs / %.4fhr on %s core(s) (%.4fcpu-s / %.4fcpu-hrs, or ~ %.2fs/cm)" % (numSeqs, count, elapsedTime, elapsedHrs, NUM_PROCESSES, totalTime, cpuHrs, avgTime)
	print "Output printed in", scoresOut
	
	# print bitscore file (table of scores, rows=seqs, cols=CMs)
	print ""
	print "======================================================"
	print "Printing score table (bitscore file)..."
	command = "python parse_cmscore_output.py %s %s" % (scoresOut, bitOut)
	result = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
	for line in result.stdout.readlines():
		print "   ", line
	retval = result.wait()
	if retval != 0:
		sys.exit()
	
	
	# normalize bitscore file
	print ""
	print "======================================================"
	print "Normalizing bitscore file..."
	command = "python pca_normalize_pipeline.py %s %s --verbose" % (bitOut, DB_FILE)
	result = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
	for line in result.stdout.readlines():
		print line,
	retval = result.wait()
	if retval != 0:
		sys.exit()
	
	print ""
	print "======================================================"
	print "Finished."
	

