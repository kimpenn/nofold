#====================================================================
# Sarah Middleton
# 
# get_clusters.r
# given a bitscore file, clusters using hclust and outputs the 
# distance matrix and clusters to file.
# 
# Usage:
# Rscript get_clusters.r filename.bitscore cluster_output.txt True "%s" "%s" "%s"' % (BITSCORE_FILE, CLUSTER_FILE, opts.SCALE, opts.HCLUST_METHOD, opts.SCALE_FILE, opts.HCLUST_DIST
#====================================================================


# scaling by custom-defined mean/var
# 'data' is a dataframe of values to be scaled (e.g. bitscore format)
# 'means' must correspond to the columns in 'data'
# 'vars' must correspond to the columns in 'data'
customScale <- function(data, means, vars)
{
	datMat <- as.matrix(data)
	if ((length(colnames(data)) != length(means)) | (length(colnames(data)) != length(vars))) {
		print("In customScale(): number of data features did not match size of mean/var vectors")
	} else {
		datMat <- apply(datMat,1,function(x) (x-means)/vars)
	}
	return(as.data.frame(t(datMat)))
}

#get command line args
args<-commandArgs(TRUE) 
bitscoreFile <- args[1]
outputFile <- args[2]
scaleData <- args[3] #should the data be scaled?
clustMethod <- args[4]
customScaleVals <- args[5] #ignored if scaleData==False, required o/w
distMethod <- args[6] #added 6/4/13. valid args are 'euclidean', 'pearson', 'spearman'

library("fastcluster")

# read file and scale if necessary
scores <- read.table(bitscoreFile)

if ((scaleData == TRUE) | (scaleData == "TRUE") | (scaleData == "True")) { 
	scaleVals <- read.table(customScaleVals)
	scores <- customScale(scores,scaleVals$Mean,scaleVals$Var) #z-score normalize bitscores 
} 

# dist()
distMatrix <- 0
if (distMethod == 'euclidean') {
	distMatrix <- dist(scores, method="euclidean")
} else if (distMethod == 'pearson') {
	distMatrix <- as.dist(1-cor(t(as.matrix(scores)), method='pearson'))
} else if (distMethod == 'spearman') {
	distMatrix <- as.dist(1-cor(t(as.matrix(scores)), method='spearman'))
} else {
	print("Invalid distance parameter.")
}

# print distance matrix
matrixOut <- paste(outputFile, ".distMatrix", sep="") #create filename
printMatrix <- as.matrix(distMatrix)
write.table(printMatrix, file=matrixOut) #overwrite if exists

# hclust()
hc <- hclust(distMatrix, method=clustMethod)

# cut tree into all possible clusters and print results
maxClusts <- length(row.names(scores))
clustersOut <- paste(outputFile, ".clusters", sep="")
totalExpected <- (maxClusts * (maxClusts + 1)) / 2 #total number of clusters, including redundant ones
results <- rep(NA, totalExpected)
clustCount <- 0
for (i in 1:maxClusts) 
{
	branches <- cutree(hc, k=i) #cut tree into i clusters
	numClusts <- length(unique(branches)) #check how many clusters were created
	for (j in 1:numClusts)
	{
		clustCount <- clustCount + 1
		clustNames <- names(which(branches == j)) #get list of seqs in this cluster
		results[clustCount] <- paste(clustNames, collapse = ",") #combine into one string and add to results[]
	}
}
cat(results, file=clustersOut, sep="\n") #overwrite if exists


# print tree heights
heightsOut <- paste(outputFile, ".heights", sep="")
cat(hc$height, file=heightsOut, sep="\n") #overwrite if exists


