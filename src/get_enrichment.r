#====================================================================
# Sarah Middleton
# 
# get_enrichment.r
# finds enrichment using fisher's exact test
#====================================================================

#get command line args
args<-commandArgs(TRUE) 
tp <- as.integer(args[1])
fp <- as.integer(args[2])
fn <- as.integer(args[3])
tn <- as.integer(args[4])
testType <- args[5]
outputFile <- args[6]

# create matrix
m1 <- rbind(c(tp, fp), c(fn, tn))

# run test
fe <- fisher.test(m1, alternative=testType)
pval <- as.character(fe["p.value"])

cat(pval, file=outputFile, sep="\n") #overwrite if exists
