#=================================================================================================
# Sarah Middleton
# 
# train_cms.py
# - creates CMs for each cluster in the given file.
# - optionally calibrates the CMs. this takes a while.
# - outputs to a folder called <DBNAME>_clust_cms
# 
#=================================================================================================
import sys, subprocess, time, shutil, os
from optparse import OptionParser
from multiprocessing import Pool, cpu_count
from utils_analysis import *

if __name__ == '__main__': 
	
	usageMsg = "Usage: %prog CLUSTER_FILE [options]"
	parser = OptionParser(usage=usageMsg)
	parser.add_option("--cpu", action="store", type='int', default=1, dest="NUM_CPU", help="number of cores to use.")
	parser.add_option("--forecast", action="store_true", default=False, dest="DO_FORECAST", help="do forecasting for runtime of calibration.")
	parser.add_option("--build", action="store_true", default=False, dest="DO_BUILD", help="do building step")
	parser.add_option("--calibrate", action="store_true", default=False, dest="DO_CALIB", help="do calibration step")
	parser.add_option("--infernal-path", action="store", default="", dest="PROG_PATH", help="path to infernal executables, if not in default path.")
	(opts, args) = parser.parse_args()
	if len(args) == 1:
		CLUSTER_FILE = args[0]
	else:
		print "Must indicate required command line args: CLUSTER_FILE"
		print "Exiting."
		sys.exit(1)
	
	# create output folder
	fn = os.path.basename(CLUSTER_FILE)
	fnParts = fn.split(".")
	DB_NAME = fnParts[0]
	outDir = "%s/%s_clust_cms" % (os.path.dirname(CLUSTER_FILE), DB_NAME)
	if not os.path.exists(outDir):
		os.makedirs(outDir)
	
	# read in list of cluster ids
	clusterList = []
	ins = open(CLUSTER_FILE, 'r')
	ins.readline()
	lines = ins.readlines()
	ins.close()
	for line in lines:
		(cid, rest) = line.split(None, 1)
		clusterList.append(cid)
	
	# create CMs
	cmList = []
	alnDir = "%s/%s_structs" % (os.path.dirname(CLUSTER_FILE), DB_NAME)
	cmbuildPath = "%scmbuild" % opts.PROG_PATH
	if opts.DO_BUILD == True:
		print ""
		print "Creating cms in %s" % outDir
		for cid in clusterList:
			structDir = "%s/%s_structs/clust%s_global.out/" % (os.path.dirname(CLUSTER_FILE), DB_NAME, cid)
			stockholmFile = "%s/clust%s.sto" % (outDir, cid)
			error = mlocarna_to_stockholm(structDir, stockholmFile) #creates stockholm format file for this cluster from mlocarna results
			if error: 
				print "Error while creating stockholm file for %s" % cid
				print "Exiting."
				sys.exit(1)
			else:
				cmFile = "%s/clust%s.cm" % (outDir, cid)
				command = "%s -F %s %s" % (cmbuildPath, cmFile, stockholmFile)
				p = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
				resultLines = p.stdout.readlines()
				retval = p.wait()
				if retval != 0:
					print "Error: problem while building cm %s" % cid
					print "cmbuild output:"
					for rline in resultLines:
						print "  ", rline
					print "Skipping."
					print ""
					if os.path.exists(cmFile):
						os.remove(cmFile)
					
					#print "Exiting."
					#sys.exit(1)
				else:
					cmList.append(cmFile)
	else: #assume the cms have already been built
		for cid in clusterList:
			cmFile = "%s/clust%s.cm" % (outDir, cid)
			cmList.append(cmFile)
				
	# check times
	total = 0
	cmcalibratePath = "%scmcalibrate" % opts.PROG_PATH
	if opts.DO_FORECAST == True:
		print ""
		print "Calibration forecast:"
		for cmfile in cmList:
			(timeToRun, error) = cm_calibrate_forecast(cmfile, cmcalibratePath)
			if error:
				print "Could not forecast cm %s" % cmfile
				print "Exiting."
				sys.exit(1)
			print "%s: %s" % (cmfile, timeToRun)
			(hours, minutes, seconds) = timeToRun.split(":")
			convertedTime = int(seconds) + (60 * int(minutes)) + (3600 * int(hours))
			total += convertedTime
		print ""
		totalHours = float(total) / 3600
		print "Total expected time: %.2f hours" % totalHours
		print "(so on %s core(s), about %.2f hours)" % (opts.NUM_CPU, (totalHours / opts.NUM_CPU))
		print ""
	
	if opts.DO_CALIB == True:
		print "Calibrating CMs"
		# prepare input for process pool
		execList = []
		for cmfile in cmList:
			newCmfile = "%s.c.cm" % cmfile[:-3] #full path
			shutil.copyfile(cmfile, newCmfile) # copy cm file so it's not overwritten
			params = [newCmfile, cmcalibratePath, ""]
			execList.append(params)
		
		# create a pool of processes to run CM scoring in parallel
		print "Creating %s processes to split %s jobs." % (opts.NUM_CPU, len(cmList))
		pool = Pool(processes = opts.NUM_CPU)
		timeStart = time.time()
		result = pool.map_async(cm_calibrate, execList)
		ret = result.wait()
		elapsedTime = time.time() - timeStart
		
		print ""
		print "Finished."
		print "Total time elapsed: %.2f sec" %  elapsedTime
		print ""
	
