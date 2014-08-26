#=================================================================================================
# Sarah Middleton
# 
# utils_file_readers.py
# Functions for reading in various file types 
#=================================================================================================


def read_comp_matrix(fileName):
	"""
	Given a pairwise comparison matrix (e.g. a dissimilarity matrix from R), returns a dict 
	containing the data and a list containing the unique ids. Assumes rows and columns are labeled.
	"""
	dissim = {}
	idList = []
	error = False
	
	try:
		ins = open(fileName, 'r')
	except IOError:
		error = True
		print "Error: in read_comp_matrix(): could not open", fileName
	else:
		header = ins.readline().split()
		lines = ins.readlines()
		ins.close()
		
		# process header
		for id in header:
			id = id.strip('\"') #get rid of quotes that R puts around the id
			idList.append(id)
		
		# process scores
		for line in lines:
			line = line.rstrip('\r\n')
			(id, rest) = line.split(None, 1)
			id = id.strip('\"') #get rid of quotes that R puts around the id
			scores = rest.split()
			if id in dissim: #error check
				print "Error: in read_comp_matrix(): repeat id (%s)" % id
				error = True
				break
			dissim[id] = {}
			if len(scores) != len(idList): #file format check
				print "Error: in read_comp_matrix(): number of scores (%s) does not match number of ids in header (%s)" % (len(scores), len(idList))
				error = True
				break
			for i in range(len(scores)):
				otherID = idList[i]
				score = float(scores[i])
				if (id == otherID) and (score != 0): #file format check
					print "Warning: in read_comp_matrix(): non-zero score for distance between an id and itself:", id
				dissim[id][otherID] = score
	
	return (idList, dissim, error)
	


def read_fasta(fileName):
	"""
	Given a fasta file, returns a dict with the sequence ID as the key and the concatonated sequence as the value.
	Header lines must start with ">"
	"""
	seqs = {}
	error = False
	
	try:
		ins = open(fileName, 'r')
	except IOError:
		error = True
		print "Error: in read_fasta(): could not open", fileName
	else:
		lines = ins.readlines()
		ins.close()
		
		id = ""
		for line in lines:
			line = line.rstrip('\r\n')
			if (len(line) > 0) and (">" == line[0]):
				id = line[1:]
				if id in seqs:
					print "Warning: in read_fasta(): repeat id (%s) in file. Overwriting." % id
				seqs[id] = ""
			else:
				seqs[id] += line.upper()
	
	return (seqs, error)
	

def read_fasta_UCSC(fileName):
	"""
	Same as above, but does extra header processing for UCSC-style fasta files. Extracts just the id and seq.
	Header should something look like this: 
		>mm9_knownGene_uc007aet.1 range=chr1:3195985-3205713 5'pad=0 3'pad=0 strand=- repeatMasking=none
	"""
	seqs = {}
	error = False
	
	try:
		ins = open(fileName, 'r')
	except IOError:
		error = True
		print "Error: in read_fasta_UCSC(): could not open", fileName
	else:
		lines = ins.readlines()
		ins.close()
		
		id = ""
		for line in lines:
			line = line.rstrip('\r\n')
			if ">" == line[0]:
				things = line.split()
				(assemb, track, id) = things[0].split("_")
				if id in seqs:
					print "Warning: in read_fasta(): repeat id (%s) in file. Overwriting." % (id, fileName)
				seqs[id] = ""
			else:
				seqs[id] += line.upper()
		
	return (seqs, error)	




def read_bitscore(fileName):
	"""
	Reads a "bitscore" format file (output from parse_cmscore_output.py). Returns list of column names, 
	dict with format scores[id] -> [score1, score2, ..., scoreN]
	"""
	scores = {}
	header = ""
	error = False
	
	try:
		ins = open(fileName, 'r')
	except IOError:
		error = True
		print "Error: in read_bitscore(): could not open", fileName
	else:
		header = ins.readline().split()
		lines = ins.readlines()
		ins.close()
		for line in lines:
			line = line.rstrip('\r\n')
			(id, scoreString) = line.split(None, 1)
			scoreList = scoreString.split()
			scores[id] = []
			for score in scoreList:
				scores[id].append(float(score))
	
	return (header, scores, error)


def read_clustalw(fileName):
	"""
	Reads a ClustalW format alignment file. Concatonates any split lines, and returns a dict with seqIDs for keys and 
	concatonated sequences as values.
	"""
	seqs = {}
	error = False
	
	# open alignment file
	try:
		ins = open(fileName, 'r')
	except IOError:
		error = True
		print "Error: in read_clustalw(): could not open", fileName
	else:
		ins.readline() # skip header
		lines = ins.readlines()
		ins.close()
		
		for line in lines:
			line = line.rstrip('\r\n')
			if line != "":
				(id, alnSeq) = line.split()
				if id not in seqs:
					seqs[id] = ""
				seqs[id] += alnSeq
		
	return (seqs, error)


def read_cmsearch_tabfile(fileName):
	"""
	Given an output file from cmscore in "tab" output format, reads the relevant info into a dict.
	Dict format: dict[seqID]['bit'/'eval'/'gc'] --> value
	"""
	info = {}
	error = False
	#HEADER_END = 23 #line where read info starts
	
	try:
		ins = open(fileName, 'r')
	except IOError:
		error = True
		print "Error: in read_cmsearch_tabfile(): could not open", fileName
	else:
		ins.readline() # skip header
		lines = ins.readlines()
		ins.close()
		for line in lines:
			if (line != "") and ("#" != line[0]):
				things = line.split()
				id = things[1]
				if (id not in info) or ((id in info) and (info[id]['bit'] < things[6])):
					info[id] = {}
					info[id]['bit'] = things[6]
					info[id]['eval'] = things[7]
					info[id]['gc'] = things[8]
					info[id]['start'] = things[2]
					info[id]['stop'] = things[3]

	return (info, error)


def mlocarna_to_stockholm(mlocFolder, stockOutFile):
	"""
	Given a folder containing the standard output of mlocarna, creates a stockholm-format file
	from the alifold.out and result.aln files.
	"""
	error = False

	# get cons struct from alifold.out
	fileName = "%scons_struct.txt" % mlocFolder
	try:
		ins = open(fileName, 'r')
	except IOError:
		error = True
		print "Error: in mlocarna_to_stockholm(): could not open", fileName
	else:
		lines = ins.readlines()
		consStruc = lines[-1]
		things = consStruc.split()
		consStruc = things[1]
		ins.close()
		
	# get alignment from result.aln
	fileName = "%sresults/result.aln" % mlocFolder
	(alnLines, alnError) = read_clustalw(fileName) #reads alignment into one line per seq
	if alnError:
		error = True
	
	# print in stockholm format
	if not error:
		try:
			outAln = open(stockOutFile, 'w')
		except IOError:
			error = True
			print "Error: in mlocarna_to_stockholm(): could not open outFile:", stockOutFile
		else:
			outAln.write("# STOCKHOLM 1.0\n\n")
			for id in alnLines:
				fixedStr = alnLines[id].replace("~", "-")
				outStr = "%s\t%s\n" % (id, fixedStr)
				outAln.write(outStr)
			outAln.write("#=GC SS_cons\t" + consStruc + "\n")
			outAln.write("//\n")
			outAln.close()
			
	return error





def merge_files(file1, file2, outFile, header=False):
	"""
	Appends file2 to file1 and writes result to outFile. If header, maintains the header from
	file1 and discards the first line of file2.
	"""
	error = False
	try:
		ins1 = open(file1, 'r')
	except IOError:
		error = True
		print "Error: in merge_files(): could not open", file1
	else:
		try:
			ins2 = open(file2, 'r')
		except IOError:
			error = True
			print "Error: in merge_files(): could not open", file2
		else:
			if header == True:
				lines2.readline() #skip header
			lines1 = ins1.readlines()
			lines2 = ins2.readlines()
			ins1.close()
			ins2.close()
			try:
				outs = open(outFile, 'w')
			except IOError:
				error = True
				print "Error: in merge_files(): could not open", outFile
			else:
				for line in lines1:
					outs.write(line)
				for line in lines2:
					outs.write(line)
				outs.close()
	
	return error
	



	