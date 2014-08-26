#====================================================================
# Sarah Middleton
# 
# get_signif_features.r
# given a bitscore file and a subset of ids (e.g. from a cluster),
# determines if any features are significantly higher for the
# subset than the whole set. does a one-sided Welch's t-test
#====================================================================

#get command line args
args<-commandArgs(TRUE) 
bitscoreFile <- args[1]
outputFile <- args[2]
idStr <- args[3] # a string, comma separated, no spaces


# read bitscore file
wholeset <- read.table(bitscoreFile)

# scale scores
wholeset <- scale(wholeset)

# convert scoreList to vector
idList <- unlist(strsplit(idStr, ","))

# get a separate frame that contains only scores from ids in idList
subset <- wholeset[rownames(wholeset) %in% idList,]

# compare corresponding columns in the two sets
cms <- colnames(wholeset)
numCols <- length(cms)
results <- rep(NA, numCols)
for (i in 1:numCols) 
{
	tt <- t.test(subset[,i], wholeset[,i], alternative="g")
	pval <- tt["p.value"]
	outStr <- paste(cms[i], pval, sep="\t")
	results[i] <- outStr
}
cat(results, file=outputFile, sep="\n") #overwrite if exists
