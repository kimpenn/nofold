#=================================================================================================
# Sarah Middleton
# 
# pca_normalize_pipeline.py
# Pipeline automates steps from raw bitscore file -> pca-normalized bitscore file
# Everything will be output in the same folder as the bitscore file.
# Requires pre-existing scaling constants.
# 
# Required args: 
#    BITSCORE_FILE = path to raw .bitscore file 
#    DB_FILE = path to .db file
# 
# Notes:
#   
#    
# Example usage:  
#    python pca_normalize_pipeline.py ../demo/demo1/demo1.bitscore ../demo/demo1/demo1.db --verbose
#=================================================================================================
import sys, os
from optparse import OptionParser
from utils_file_readers import *
from utils_analysis import *

print ""

usageMsg = "Usage: %prog BITSCORE_FILE DB_FILE [options]"
parser = OptionParser(usage=usageMsg)

parser.add_option("--num-pc", action="store", dest="NUM_PC", default="100", help="Number of PC axes to retain. Default is [%default].")
parser.add_option("--out-prefix", action="store", dest="OUT_PREFIX", default=None, help="Prefix to use for output files. Default is to use the name of the input bitscore file up to the first period.")
parser.add_option("--mean-scale", action="store", dest="SCALE_MEANS", default="../norm/varlen2.scale_means.txt", help="File containing mean values of each original feature. Default is [%default].")
parser.add_option("--sd-scale", action="store", dest="SCALE_SDS", default="../norm/varlen2.scale_sds.txt", help="File containing SD values of each original feature. Default is [%default].")
parser.add_option("--pc-scale", action="store", dest="SCALE_PCS", default="../norm/varlen2.zNorm.pcNorm.scale_pc.txt", help="File containing scaling constants for PC axes. Default is [%default].")
parser.add_option("--eigs", action="store", dest="EIGS", default="../norm/varlen2.eigs_subset.txt", help="File containing mean values of each original feature. Default is [%default].")
parser.add_option("--verbose", action="store_true", default=False, dest="VERBOSE", help="Print output from called scripts. May give better indications of progress.")


# read/process args
(opts, args) = parser.parse_args()
if len(args) == 2:
	BITSCORE_FILE = os.path.abspath(args[0])
	DB_FILE = os.path.abspath(args[1])
else:
	print "Incorrect args."
	print "Use -h for help. Exiting."
	sys.exit()

# create file names
outFolder = os.path.dirname(BITSCORE_FILE)
fn = os.path.basename(BITSCORE_FILE)
fnParts = fn.split(".")
dbName = fnParts[0]
if opts.OUT_PREFIX != None:
	dbName = opts.OUT_PREFIX
SEQINFO_FILE = "%s/%s.seq_info.txt" % (outFolder, dbName)
NORM1_FILE = "%s/%s.zNorm.bitscore" % (outFolder, dbName)
NORM2_FILE = "%s/%s.zNorm.pcNorm%s.bitscore" % (outFolder, dbName, opts.NUM_PC)
NORM3_FILE = "%s/%s.zNorm.pcNorm%s.zNorm.bitscore" % (outFolder, dbName, opts.NUM_PC)


# get seq info

print ""
print "---------------------"
print "Step 1: Get seq info."
command = "python varlen_seq_info.py %s %s %s" % (DB_FILE, BITSCORE_FILE, SEQINFO_FILE)
(output, res, error) = run_command(command, verbose=opts.VERBOSE)
if error:
	print ">>Error detected. Exiting."
	sys.exit()


# length normalize rfam

print ""
print "---------------------"
print "Step 2: Standardize raw bitscore file, taking into account sequence length effects."
command = "python varlen_z_normalize.py %s %s %s %s %s" % (BITSCORE_FILE, SEQINFO_FILE, opts.SCALE_MEANS, opts.SCALE_SDS, NORM1_FILE)
(output, res, error) = run_command(command, verbose=opts.VERBOSE)
if error:
	print ">>Error detected. Exiting."
	sys.exit()

	
# pcNormalize

print ""
print "---------------------"
print "Step 3: Convert to PCA space, retain top axes."
command = "python varlen_pc_normalize.py %s %s %s %s" % (NORM1_FILE, opts.EIGS, NORM2_FILE, opts.NUM_PC)
(output, res, error) = run_command(command, verbose=opts.VERBOSE)
if error:
	print ">>Error detected. Exiting."
	sys.exit()


# re-normalize

print ""
print "---------------------"
print "Step 4: Standardize PC axes."
command = "python varlen_z_normalize_pc.py %s %s %s" % (NORM2_FILE, opts.SCALE_PCS, NORM3_FILE)
(output, res, error) = run_command(command, verbose=opts.VERBOSE)
if error:
	print ">>Error detected. Exiting."
	sys.exit()

print ""
print "Normalized bitscore printed to", NORM3_FILE
print "Finished."



















