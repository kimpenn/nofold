# NoFold: A motif finder for RNA secondary structures


NoFold is an approach for characterizing and clustering RNA secondary structures without computational folding or alignment. It works by mapping each RNA sequence of interest to a structural feature space, where each coordinate within the space corresponds to the probabilistic similarity of the sequence to an empirically defined structure model (e.g. Rfam family covariance models). NoFold provides scripts for mapping sequences to this structure space, extracting any robust clusters that are formed, and annotating those clusters with structural and functional information.

This repo hosts the development version of NoFold, which may differ from the version described in the paper. You can download the paper version [here](http://kim.bio.upenn.edu/software/nofold.shtml).

Questions: contact sarahmid@mail.med.upenn.edu.


<br>
# Download

**This repo contains the following:**
- All code needed to run NoFold (scoring, clustering, annotation)
- 1,973 Rfam covariance models (will be used automatically)
- Pre-made threshold files appropriate for datasets of up to ~4,000 sequences
- A script for generating thresholds specific to your dataset size, if necessary
- A demo dataset for testing your installation

**In addition, you will need the following external programs:**
- [Python 2.X](https://www.python.org/downloads/)
- [Infernal (v.1.0.2)](http://eddylab.org/software/infernal/infernal-1.0.2.tar.gz) (NOT the newest version, which doesn't contain cmscore)
- [R](https://www.r-project.org/) with [fastcluster](http://danifold.net/fastcluster.html) package 
- [LocARNA](http://www.bioinf.uni-freiburg.de/Software/LocARNA/) 

We recommend also installing [RNAz](http://www.tbi.univie.ac.at/~wash/RNAz/) for additional annotation of clusters. Use --rnaz option when clustering if available.


<br>
# Installation


1. Install required external software (see above)
2. Unzip the tar file: tar -zxvf nofold.tar.gz
3. NoFold is now ready to go. Try running the demo dataset to verify everything 
   is working properly (see below).



<br>
# Quick-start example with demo datatset


Run the following from within the `/src/` directory. You may need to change the number of CPUs used, or provide the paths to your Infernal/LocARNA installations. See "Script Details" for more information.

1. Get normalized structural coordinates of sequences:
```
python score_and_normalize.py ../demo/demo1/demo1.db --cpus=4
```

2. Extract any clusters that form within the structure space:
```
python nofold_pipeline.py ../demo/demo1/demo1.zNorm.pcNorm100.zNorm.bitscore ../demo/demo1/demo1.db --cpus=4 --bounds-file=../thresh/bounds_30seq.txt --verbose
```

The second script outputs a detailed annotation file of the identified clusters:
```
../demo/demo1/demo1.clusters_s3rSpec_top.txt_expanded_merged_bs11.50bgNoneGloc.txt.details
```


<br>
# General usage guide



## 1\. (Prep) Create FASTA-format input file.

### Sequence IDs
The starting input to NoFold is simply a file of RNA sequences. First, ensure that your sequences in are in FASTA format with the following header format:    
```
>seqID
```
or, for background files:
```
>groupID_seqID
```

- The header line should be a single unique identifier (no spaces) for each sequence. 
- This identifier must be <= 25 characters, otherwise it will be clipped by Infernal in the score output file.
- This identifier can optionally contain a "groupID" (e.g. "test" or "bg") which should be separated  from the rest of the ID by an underscore. This information is used in situations where a background database is being used to test enrichment, in which case all sequences belonging to the background set should have the same group ID. 


### FASTA file name
The file name of the fasta file will be used to prefix all output files, and should therefore be chosen so that this will not cause any other files in the same directory to be overwritten. To determine this prefix (the "DBNAME"), NoFold uses only the part of the file name up to the first period. 

For example:
```
mydataset.db   --> DBNAME will be "mydataset"
my_dataset.db  --> DBNAME will be "my_dataset"
my.dataset.db  --> DBNAME will be "my"
```
- The .db extension is not required, this is just the convention we use here to indicate the sequence file.
- We recommend choosing a short DBNAME, since this will be used to prefix all output files/folders.
- We recommend placing this file in a new directory, since several output files and directories will be generated within the directory the sequence file is located in.



<br>
## 2\. (Prep) Check if appropriate thresholds are available. 

The thresholds used for filtering clusters are dependent on the size of your input database. We have provided threshold files for datasets of up to about 4,000 sequences, which can be found in the /thresh/ directory. Choose the file `"bounds_<DATASET_SIZE>seq.txt"` such that `<DATASET_SIZE>` is closest to the number of sequences in your dataset. In general, rounding up will result in slightly more stringent thresholds, while rounding down may result in less stringency. We usually round to the nearest 100.

### (Optional) Generating your own threshold files
If there are no appropriate threshold files, or you would like to generate your own, use the script `calculate_clust_thresh.py` as follows (from within the `/src/` directory):
```    
python calculate_clust_thresh.py <SIZE>
```
For example, to create a threshold set for a 550-sequence dataset, you would run:
```    
python calculate_clust_thresh.py 550
``` 
This will output a new bounds file in `/thresh/`, which you can use for clustering (see directions for clustering below). Note, creating thresholds usually takes a few minutes because it's basically creating and clustering many random datasets. Getting thresholds for large datasets (>3,000) may take over an hour to complete. Please note that since the random dataset used as a source for this process only contains 10,000 sequences, you will not be able to get thresholds for datasets larger than that unless you generate your own larger dataset.

Note: we are looking into implementing a model-based estimation of these thresholds, which will make it easier to get thresholds for larger datasets. This will be introduced in a future update.



<br>
## 3\. Score and normalize sequences.

In this step, your sequences will be scored against the 1,973 Rfam CMs and the scores will be normalized as described in the [NoFold paper](http://rnajournal.cshlp.org/content/20/11/1671.long). To do this, run the following command from within the `/src/` directory:
```    
python score_and_normalize.py <DB_FILE> --cpus=<NUM_CPU>
```
For example, to run the demo1 dataset on 4 cpu cores:
```    
python score_and_normalize.py ../demo/demo1/demo1.db --cpus=4
```
The final output of this script will be a file called `*.zNorm.pcNorm100.zNorm.bitscore` which contains all of the normalized scores in a matrix format. Make sure to use this file, and NOT any of the other `.bitscore` files, as the input to the clustering script (below). The other `.bitscore` files are intermediate output from the normalization process.

### Note on time
The scoring step takes ~0.012 seconds per CM per 50nt sequence on one core. This scales linearly with increasing numbers of sequences, but scales quadratically with increasing sequence length. Therefore, we recommend that you convert long sequences into a sliding window dataset (preferably with a "slide" of ~70% of the sequence length). NoFold can theoretically parallelize scoring up to 1,973 CPUs (i.e. the number of CMs). We are looking into adding MPI support in a future update.



<br>
## 4\. Cluster scored sequences.

In this step, your scored sequences will be clustered, filtered, and annotated. To do this, run the following command from within the `/src/` directory:
```    
python nofold_pipeline.py <NORMALIZED_BITSCORE_FILE> <DB_FILE> --cpus=<NUM_CPU> --bounds-file=<THRESH_FILE> --verbose
```
For example, to run the demo1 dataset on 4 cpu cores:
```    
python nofold_pipeline.py ../demo/demo1/demo1.zNorm.pcNorm100.zNorm.bitscore ../demo/demo1/demo1.db --cpus=4 --bounds-file=../thresh/bounds_30seq.txt --verbose
```
If you have RNAz installed, you can also add the --rnaz option to output useful RNAz annotations in the `*.details` files. Type -h for additional options and usage.

There are several output files from this (see usage details below) but the main file of interest is:
```
*.clusters_s3rSpec_top.txt_expanded_merged_bs##bgNoneGloc.txt.details
```
which contains annotation of all expanded clusters that passed the filters. Consensus structures for these clusters (predicted by LocARNA) can be found in  `/*_structs/cons_struc_pics/` folder (they are labeled by cluster ID).

### Enrichment testing
If you would like to test your clusters for enrichment, you need to prepare a separate fasta file of "background" sequences. Ensure that each sequence in the file has a common "groupID" in the header, e.g. `bg_seq1`, `bg_seq2`, etc. (see fasta format notes above). This is used to distinguish between query and background sequences during the test. Enrichment testing is integrated into the NoFold pipeline script, so just add the following options onto your clustering command:
```
--bg-db=<BG_DB> --bg-id=<GROUP_ID>
```
where `<BG_DB>` is the path to your background fasta file (`*.db`) and `<GROUP_ID>` is the groupID portion of the sequence IDs. For example, here is how we can cluster `demo1.db` with enrichment testing:
```
python nofold_pipeline.py ../demo/demo1/demo1.zNorm.pcNorm100.zNorm.bitscore ../demo/demo1/demo1.db --cpus=4 --bounds-file=../thresh/bounds_30seq.txt --verbose --bg-db=../demo/demo1/demoBG.db --bg-id=DEMOBG
```
This causes an additional file to be printed that contains the enrichment p-value of each cluster. If any sequences from the background were found to have the same structure, those sequence's IDs will be listed under "Members":
```
/demo/demo1/demo1.clusters_s3rSpec_top.txt_expanded_merged_bs11.50bgNoneGloc.txt.summary_expanded_enrich_bs11.50bgDEMOBGGloc.txt
```


<br>
# Script details



> ### `score_and_normalize.py`

    External software requirements:
        - Python 2.X
        - Infernal v. 1.0.2 
    
    Input files: 
        - A fasta file of query sequences to score 
        - A set of CMs to use for scoring
    
    Output files:
        - A folder called /cmscore_results_rfam/ containing score files for each CM, which list the bitscores of each query sequence (can be deleted)
        - *.bitscore - A file containing raw scores in a matrix format. Rows = sequences, Columns = CMs.
        - *.seq_info.txt - A file containing some stats about each sequence
        - *.zNorm.bitscore - intermediate normalization file, used for cluster annotation
        - *.zNorm.pcNorm100.bitscore - intermediate normalization file, not currently used
        => *.zNorm.pcNorm100.zNorm.bitscore - final normalized file, main file used for clustering

    Usage:
        python score_and_normalize.py FASTA_FILE [options]


    Options:
        -h, --help          show this help message and exit
        
        --cpus=MAX_CPU      Maximum number of CPUs to use. Default is [1].
        
        --infernal-path=INFERNAL_PATH
                            Path to Infernal executables. Default is to assume
                            they have been added to your PATH and can be called by
                            name.
        
        --cm-folder=CM_FOLDER
                            Folder containing the covariance models you would like
                            to use. Default is to use the 1,973 Rfam CMs included
                            in [../models/rfam_cms/].

                        
<br>
> ### `nofold_pipeline.py`
    
    External software requirements:
        - Python 2.X
        - R with fastcluster package
        - LocARNA (v.1.7.2 tested, others may work)
        - Infernal v.1.0.2
        - (Optional) RNAz
    
    Input files:
        - the normalized bitscore file output from score_and_normalize.py (*.zNorm.pcNorm100.zNorm.bitscore)
        - the original fasta file used for scoring
        - a file with appropriate thresholds (i.e. bounds_#seq.txt, where # is roughly the number of query sequences). Note: you must specify either this file (recommended), or a sigle pair of thresholds as --bounds=THRESH1,THRESH2
    
    Output files:
        -Intermediate output:
            - *.hierarchical.distMatrix - distance matrix between all query sequences
            - *.hierarchical.heights - join heights of dendrogram
            - *.hierarchical.clusters - list of all possible clusters from dendrogram
            - *.clusters_s3rSpec_top.txt - list of clusters that passed initial thresholds
            - *.clusters_s3rSpec_all.txt - list of all clusters
            - *.clusters_s3rSpec_top.txt.details - annotation of clusters that passed thresh
            - *.clusters_s3rSpec_top.txt.summary - summary stats of clusters that passed thresholds
            - *.clusters_s3rSpec_top.txt_expanded_bs##bgNoneGloc.txt - summary of clusters after expansion 
                               (## = bitscore threshold used for expansion)
            - *.clusters_s3rSpec_top.txt_expanded_merged_bs##bgNoneGloc.txt - summary of clusters after merging 
                               (## = bitscore threshold used for expansion)
            
        -Main output:
            => *.clusters_s3rSpec_top.txt_expanded_merged_bs##bgNoneGloc.txt.details - MAIN CLUSTER ANNOTATION FILE
                               detailed annotation of each final cluster (## = bitscore threshold used for expansion)
            => *.clusters_s3rSpec_top.txt_expanded_merged_bs##bgNoneGloc.txt.summary - MAIN CLUSTER ANNOTATION FILE
                               summary stats for each final cluster (## = bitscore threshold used for expansion)
            => folder /*_clust_cms/ - cluster-CMs, a CM generated for each threshold-passing cluster
            => folder /*_structs/ - predicted structures of each threshold-passing cluster

    Usage:
        python nofold_pipeline.py BITSCORE_FILE FASTA_FILE --bounds-file=THRESHOLD_FILE [options]
    
    
    Options:
          -h, --help            show this help message and exit

          Common options:
            --cpus=MAX_CPUS     Maximum number of threads to use. Default is [1].
            --bounds-file=BOUNDS_FILE
                                Bounds/threshold file, as generated by 
                                calculate_cluster_thresh.py. Overrides --bounds.
            --rnaz              Use RNAz to annotate clusters. Default is [False].
            --verbose           Print output from called scripts. May give 
                                better indications of progress.
            --infernal-path=INFERNAL_PATH
                                Path to Infernal executables. Default is to 
                                assume they have been added to your PATH and can 
                                be called by name.
            --locarna-path=LOCARNA_PATH
                                Path to LocARNA executables. Default is to assume 
                                they have been added to your PATH and can be 
                                called by name.
            --rnaz-path=RNAZ_PATH
                                Path to RNAz executable. Default is to assume it 
                                has been added to your PATH and can be called by 
                                name. This is ignored unless --rnaz is also 
                                specified.

          Enrichment testing:
            --bg-id=BG_ID       Common id of bg seqs. Should be first part of ID, 
                                separated from the rest of the ID by an 
                                underscore. Needed for use of --bg-db
            --bg-db=BG_DB       Fasta file of bg seqs. Only needed if doing 
                                enrichment test.

          Clustering tweaks:
            --min-clust-size=MIN_CLUST_SIZE
                                Minimum size cluster to be considered. All 
                                smaller clusters will be discarded. Default 
                                is [3].
            --bit-thresh=BIT_THRESH
                                Bitscore threshold to use while expanding 
                                clusters using Infernal cmsearch. Default is 
                                log2(db size).
            --merge=MERGE_FRAC  Fraction of two clusters that must overlap for 
                                the two clusters to be merged. Default is [0.5].

          Advanced options:
            --bounds=BOUNDS     Two upper bounds to use for filtering by similarity.
                                Must be specified if --bounds-file is not. Fmt:
                                int,int. First number should be the less stringent
                                (higher) thresh, second is more stringent.
            --orig-bs=ORIG_BS   Original scaled bitscore file, e.g.
                                db_name.zNorm.bitscore. Used for annotating which CMs
                                scored highly in a cluster. Default is to assume this
                                file is in the same folder as the supplied bitscore
                                file.
            --clust-method=CLUST_METHOD
                                Method for clustering sequences. Currently, the only
                                choice is 'hierarchical'
            --hclust-method=HCLUST_METHOD
                                Method for hierarchical clustering. Choices: anything
                                that works with hclust in R. Note that if you change
                                this, you must generate an appropriate threshold file
                                using the same settings.
            --hclust-dist=HCLUST_DIST
                                Distance measure for hierarchical clustering. Choices:
                                euclidean, pearson, spearman (default). For the corr
                                measures, (1-corr) is used. Note that if you change
                                this, you must generate an appropriate threshold file
                                using the same settings.
            --filter=FILTER     Type of filtering for clusters. Choices: specific,
                                sensitive. Default is [specific].
            --calibrate         Calibrate cluster-CMs before expanding, allowing the
                                used of an E-value threshold. Warning: this can take a
                                VERY long time. Default is [False]
            --e-val=E_THRESH    E-value threshold to use while expanding clusters
                                using Infernal cmsearch. Default is [1]. Only used if
                                --calibrate is specified.
            --scale             Standardize the input bitscore file? Default is False,
                                since the input bitscore files are usually already
                                scaled by score_and_normalize.py.
            --scale-file=SCALE_FILE
                                A file containing pre-determined means and variances
                                for each feature, as generated in
                                calculate_universal_scale.py. Default is none. This
                                should only be used with the --scale option.
                                

<br>
> ### `calculate_clust_thresh.py`

    External software requirements:
        - Python 2.X
        - R with fastcluster package
    
    Input files:
        - None (automatically uses the included bitscore file with random sequences, 
          bgDinucWT_50.zNorm.pcNorm100.zNorm.bitscore)
    
    Output files:
        => bounds_#seq.txt - calculated thresholds for each size cluster, used for clustering
    
    Usage:
        python caclulate_clust_thresh.py <SIZE>



<br>
# Reference


Middleton, SA and Kim, J. 2014. NoFold: RNA structure clustering without folding or alignment. RNA 20:1671-1683. 
[online](http://rnajournal.cshlp.org/content/20/11/1671.long)