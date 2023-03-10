# tronko
A rapid phylogeny-based method for accurate community profiling of large-scale metabarcoding datasets

In the tronko package there are two modules: `tronko-build` and `tronko-assign`. `tronko-build` is for building custom reference databases that tronko-assign uses as input. We have two reference databases currently available for download with `tronko-assign`. Cytochrome oxidase I (COI) which was custom built with <a href="https://github.com/limey-bean/CRUX_Creating-Reference-libraries-Using-eXisting-tools">CRUX</a> using forward primer `GGWACWGGWTGAACWGTWTAYCCYCC` and reverse primer `TANACYTCnGGRTGNCCRAARAAYCA`. 16S which was custom built with <a href="https://github.com/limey-bean/CRUX_Creating-Reference-libraries-Using-eXisting-tools">CRUX</a> using forward primer `GTGCCAGCMGCCGCGGTAA` and reverse primer `GACTACHVGGGTATCTAATCC`. 
<a href="https://doi.org/10.5281/zenodo.7407318"><img src="https://zenodo.org/badge/DOI/10.5281/zenodo.7407318.svg" alt="DOI"></a>

Alignment-based and composition-based assignment methods calculate the lowest common ancestor (LCA) using data only in the leaf nodes of a phylogeny (A). The advantage of Tronko is that it stores fractional likelihoods in all nodes of a phylogeny and calculates the LCA based on all nodes in the tree (B).
<img src="https://github.com/lpipes/tronko/blob/main/Overview_Figure.jpg?raw=true">

# tronko-build
`tronko-build` is for building custom reference databases to be used with `tronko-assign`.

	tronko-build [OPTIONS]
	
		-h,				usage:
		-y, use a partition directory (you want to partition or you have multiple clusters)
		-l, use only single tree (do not partition)
		-t [FILE], rooted phylogenetic tree [FILE: Newick]
		-m [FILE], multiple sequence alignment [FILE: FASTA]
		-d [DIRECTORY], REQUIRED, output directory
		-x [FILE], taxonomy file [FILE: FASTA_header	domain;phylum;class;order;family;genus;species]
		-e [DIRECTORY], directory for reading multiple clusters
		-n [INT], number of partitions in read directory
		-b [INT], restart partitions with partition number [default: 0]
		-s, partition using sum-of-pairs score [can't use with -f]
		-u [FLOAT], minimum threshold for sum of pairs score [default: 0.5]
		-v, partition using minimum number of leaf nodes [can't use with -s, use with -f]
		-f [INT], don't partition less than the minimum number of leaf nodes [can't use with -s, use with -v]
		-g, don't change missing data

# tronko-assign
`tronko-assign` is for species assignment of queries. It requires a `tronko-build` database.

	tronko-assign [OPTIONS]
	
		-h, usage:
		-p, use paired reads
		-s, use single reads
		-r, use a reference
		-v, when using single reads, reverse-complement it
		-z, when using paired-end reads,  reverse-complement the second read
		-f [FILE], path to reference database file
		-o [FILE], path to output file
		-g [FILE], path to single-end reads file
		-1 [FILE], path to paired-end forward read file
		-2 [FILE], path to paired-end reverse read file
		-a [FILE], path to fasta file (for bwa database)
		-c [INT], LCA cut-off to use [default:5]
		-C [INT], number of cores
		-L [INT], number of lines to read for assignment
		-P, print alignments to stdout
		-w, use Needleman-Wunsch Alignment (default: WFA)
		-q, Query is FASTQ [default is FASTA]
		-e, Use only a portion of the reference sequences
		-n [INT], Padding (Number of bases) to use in the portion of the reference sequences
		-5 [FILE], Print tree number and leaf number
		-6, Skip the bwa build if already exists
		-u, Score constant [default: 0.01]

Tronko uses the <a href="https://github.com/smarco/WFA2-lib">Wavefront Alignment Algorithm (version 2)</a> or <a href="https://github.com/noporpoise/seq-align">Needleman-Wunsch Algorithm</a> for semi-global alignments. It uses <a href="https://github.com/lh3/bwa">bwa</a> for alignment to leaf nodes, and uses <a href="https://github.com/DavidLeeds/hashmap">David Leeds' hashmap</a> for hashmap implementation in C. `tronko-assign` does not reverse complement your reads automatically. You must use options `-v` or `-z` to reverse complement your read for better alignment to the reference database. For more information on the direction of your reads based on your library prep, please refer to this helpful blog here: <a href="http://onetipperday.blogspot.com/2012/07/how-to-tell-which-library-type-to-use.html">http://onetipperday.blogspot.com/2012/07/how-to-tell-which-library-type-to-use.html</a>.

## Example output

The output file is a tab-delimited text file where only the forward readname is retained (if using paired-end reads). The output displays the taxonomic path for assignment, the score, the number of forward read mismatches with the `bwa` hit, the number of reverse read mismatches with the `bwa` hit, the tree number for the best assignment (0 if using 1 tree), and the node number the read (or reads in the case of paired-end reads) was assigned to. For single-end reads, the `Reverse_Mismatch` will always be 0 and the `Forward_Mismatch` is the number of read mismatches with the `bwa` hit.

```
Readname	Taxonomic_Path	Score	Forward_Mismatch	Reverse_Mismatch	Tree_Number	Node_Number
GU572157.1_8_1	Eukaryota;Chordata;Aves;Charadriiformes;Alcidae;Uria;Uria aalge	-54.258690	5.000000	4.00000	0	1095
GU572157.1_7_1	Eukaryota;Chordata;Aves;Charadriiformes;Alcidae;Uria;Uria aalge	-42.871226	1.000000	6.00000	0	1095
GU572157.1_6_1	Eukaryota;Chordata;Aves;Charadriiformes;Alcidae;Uria;Uria aalge	-59.952407	7.000000	3.00000	0	1098
GU572157.1_5_1	Eukaryota;Chordata;Aves;Charadriiformes;Alcidae;Uria;Uria aalge	-31.483761	2.000000	4.00000	0	1095
GU572157.1_4_1	Eukaryota;Chordata;Aves;Charadriiformes;Alcidae;Uria;Uria aalge	-31.483761	2.000000	3.00000	0	1095
GU572157.1_3_1	Eukaryota;Chordata;Aves;Charadriiformes;Alcidae;Uria;Uria aalge	-54.258690	2.000000	7.00000	0	1095
```

# INSTALLATION

`tronko-assign` has no dependencies. `tronko-build` only has dependencies if using the partition procedure for the initial trees. These are <a href="https://github.com/stamatak/standard-RAxML">`raxmlHPC-PTHREADS`</a>, <a href="https://github.com/refresh-bio/FAMSA">`famsa`</a>, `nw_reroot` from <a href="https://anaconda.org/bioconda/newick_utils/files">Newick utilties</a>, <a href="https://raw.githubusercontent.com/lpipes/tronko/main/scripts/fasta2phyml.pl">`fasta2phyml.pl`</a>, and <a href="https://ftp.gnu.org/gnu/sed/">`sed`</a>, which must be installed in your path.

	cd tronko/tronko-build
	make
	../tronko-assign
	make

# SINGULARITY CONTAINER

To use the container download:

	singularity pull library://lpipes/tronko/tronko:1.0

To run `tronko-assign` with the container:

	singularity exec --bind <root-dir-to-bind-to-container> tronko.sif tronko-assign

To run `tronko-build` with the container:

	singularity exec --bind <root-dir-to-bind-to-container> tronko.sif tronko-build

# `tronko-assign` Usage
Tronko does not detect the correct orientation of the reads. If your reverse read needs to be reverse complemented use the option `-z`. The default options of Tronko assume that your reads are in FASTA format. If you want to assign reads in FASTQ format, use the option `-q`. You will also need a FASTA file (not gzipped) of all of your reference sequences in the reference database (use the option `-a`). To skip the `bwa index` build use `-6`. The reads (and reference database file) can be gzipped or not gzipped. Assigning paired-end reads in FASTA format:
```
tronko-assign -r -f [tronko-build REFERENCE DB FILE] -p -1 [FORWARD READS FASTA] -2 [REVERSE READS FASTA] -a [REFERENCE SEQUENCES FASTA] -o [OUTPUT FILE]
```
Assigning single-end reads in FASTA format:
```
tronko-assign -r -f [tronko-build REFERENCE DB FILE] -s -g [READS FASTA] -a [REFERENCE SEQUENCES FASTA] -o [OUTPUT FILE]
```

# `tronko-build` Usage

## `tronko-build` Simple Usage (using 1 phylogenetic tree)
```
tronko-build -l -t [Rooted Newick Tree] -m [Multiple Sequence Alignment FASTA] -x [TAXONOMY FILE] -d [OUTPUT DIRECTORY] 
```
The taxonomy file is a `.txt` file that has the following format:
```
FASTA_header\tdomain;phylum;class;order;family;genus;species
```
The tree file, MSA file, and the taxonomy file must all contain identical corresponding names. The MSA file should not contain any line breaks. To remove line breaks we recommend
```
sed -i ':a; $!N; /^>/!s/\n\([^>]\)/\1/; ta; P; D' test.fasta
```
## `tronko-build` example datasets

An example dataset for a single tree is provided in `tronko-build/example_datasets/single_tree`. The dataset includes 1466 COI sequences from the order Charadriiformes. The MSA is `tronko-build/example_datasets/single_tree/Charadriiformes_MSA.fasta`, the taxonomy file is `tronko-build/example_datasets/single_tree/Charadriiformes_taxonomy.txt`, and the tree file is `tronko-build/example_datasets/single_tree/RAxML_bestTree.Charadriiformes.reroot`. To build the reference database for `tronko-assign` for this dataset, run `tronko-build` with the following command (`-l` specifies using a single tree):
```
tronko-build -l -m tronko-build/example_datasets/single_tree/Charadriiformes_MSA.fasta -x tronko-build/example_datasets/single_tree/Charadriiformes_taxonomy.txt -t tronko-build/example_datasets/single_tree/RAxML_bestTree.Charadriiformes.reroot -d tronko-build/example_datasets/single_tree
```
An example dataset for multiple trees is provided in `tronko-build/example_datasets/multiple_trees/one_MSA`. The dataset includes 99 COI sequences from different species. The MSA is `tronko-build/example_datasets/multiple_trees/one_MSA/1_MSA.fasta`, the taxonomy file is `tronko-build/example_datasets/multiple_trees/one_MSA/1_taxonomy.txt`, and the tree file is `tronko-build/example_datasets/multiple_trees/one_MSA/RAxML_bestTree.1.reroot`. To build the reference database and partition it using the Sum-of-Pairs score (see manuscript), run `tronko-build` with the following command [-d is REQUIRED and must contain the full path NOT relative path):
```
tronko-build -y -e tronko-build/example_datasets/multiple_trees/one_MSA -n 1 -d [FULL PATH TO OUTPUT DIRECTORY] -s
```

An example dataset to not partition multiple tree is provided in `tronko-build/example_datasets/multiple_trees/multiple_MSA`. The dataset include 99 COI sequences from different species. To build the reference database and not partition the database, set the `-f` parameter higher than the number of sequences in the dataset. Run `tronko-build` with the following command:
```
tronko-build -y -e tronko-build/example_datasets/multiple_trees/multiple_MSA -n 5 -d outdir_multiple_MSA -v -f 500
```

A successful run will produce a `tronko-assign` reference database named `reference_tree.txt` in the output directory. The `reference_tree.txt` database file is what is used to run `tronko-assign`.

## `tronko-assign` example datasets

An example dataset for `tronko-assign` is provided in `example_datasets/single_tree`. This dataset contains single-end reads (`example_datasets/single_tree/missingreads_singleend_150bp_2error.fasta`) and paired-end reads (forward read: `example_datasets/single_tree/missingreads_pairedend_150bp_2error_read1.fasta` and reverse read: `example_datasets/single_tree/missingreads_pairedend_150bp_2error_read2.fasta`). For the single-end reads, there are 164 150bp reads with 2% simulated error/polymorphisms from the following sequence (taxonomically classified on NCBI as Uria aalge):

```
>GU572157.1
CCTGGCTGGTAATCTAGCCCATGCCGGAGCTTCAGTGGATTTAGCAATCTTCTCCCTTCACTTAGCAGGTGTATCATCTATTCTAGGCGCTATCAACTTTATCACAACAGCCATCAACATAAAGCCTCCAGCCCTCTCACAATACCAAACCCCCCTATTCGTATGATCAGTACTTATCACTGCTGTCCTACTACTACTCTCACTCCCAGTACTTGCTGCTGGTATCACTATATTACTAACAGATCGAAACTTAAACACAACATTCTTTGATCCAGCTGGAGGTGGTGACCCAGTACTTTACCAACACCTCTTC
``` 

This particular sequence, `GU572157.1`, has been removed from the single tree example database. To obtain assignments, run `tronko-assign` with the single tree example database from `tronko-build` with the following command for the single-end reads (using the Needleman-Wunsch alignment and a default LCA cut-off of 5):

```
tronko-assign -r -f tronko-build/example_datasets/single_tree/reference_tree.txt -a tronko-build/example_datasets/single_tree/Charadriiformes.fasta -s -g example_datasets/single_tree/missingreads_singleend_150bp_2error.fasta -o example_datasets/single_tree/missingreads_singleend_150bp_2error_results.txt -w
```

To obtain assignments, run `tronko-assign` with the single tree example database from `tronko-build` with the following command for the paired-end reads (using the Needleman-Wunsch alignment and a default LCA cut-off of 5):

```
tronko-assign -r -f tronko-build/example_datasets/single_tree/reference_tree.txt -a tronko-build/example_datasets/single_tree/Charadriiformes.fasta -p -1 example_datasets/single_tree/missingreads_pairedend_150bp_2error_read1.fasta -2 example_datasets/single_tree/missingreads_pairedend_150bp_2error_read2.fasta -o example_datasets/single_tree/missingreads_pairedend_150bp_2error_results.txt -w
```

An example dataset for `tronko-assign` with multiple trees is provided in `example_datasets/multiple_trees`. This dataset contains single-end reads (`example_datasets/multiple_trees/missingreads_singleend_150bp_2error.fasta`) and paired-end reads (forward read: `example_datasets/multiple_trees/missingreads_pairedend_150bp_2error_read1.fasta` and reverse read: `example_datasets/multiple_trees/missingreads_pairedend_150bp_2error_read2.fasta`). For the single-end reads, there are 84 150bp reads with 2% simulated error/polymorphisms from the following sequence (taxonomically classified on NCBI as Sagitta elegans):

```
>KP857443.1
TTTGAGCACTGTGGGACATAGGGGTGGAGCAGTGGATTTGGGTATCTTCTCTTTGCACCTGGCTGGCGTTAGAAGAATCTTGGGGAGAGCTAATTTTATTACCACTATCACCAATATAAAAGGGGAAGGTATGACTATAGAACTCATGCCTTTATTCGTGTGGGCGGTGCTCCTCACGGCTGTCTTACTTTTACTCTCTCTACCTGTATTAGCTGGGGCTATCACAATGTTAC
```

This particular sequence, `KP857443.1`, has been removed from the multiple trees example databases. To obtain assignments, run `tronko-assign` with the multiple trees example database with partitions from `tronko-build` with the following command for the single-end reads (using the Needleman-Wunsch alignment `-w` and a default LCA cut-off of 5 `-c 5`):

```
tronko-assign -r -f out_oneMSA/reference_tree.txt -s -g example_datasets/multiple_trees/missingreads_singleend_150bp_2error.fasta -o example_datasets/multiple_trees/missingreads_singleend_150bp_2error_partition_results.txt -a tronko-build/example_datasets/multiple_trees/one_MSA/1.fasta
```

To assign paired-end reads on the same reference database run `tronko-assign`:
```
tronko-assign -r -f out_oneMSA/reference_tree.txt -p -1 example_datasets/multiple_trees/missingreads_pairedend_150bp_2error_read1.fasta -2 example_datasets/multiple_trees/missingreads_pairedend_150bp_2error_read2.fasta -o example_datasets/multiple_trees/missingreads_pairedend_150bp_2error_partition_results.txt -a -a tronko-build/example_datasets/multiple_trees/one_MSA/1.fasta
```

To obtain assignments, run `tronko-assign` with the multiple trees example database with multiple MSAs from `tronko-build` with the following command for the single-end reads (using the Needleman-Wunsch alignment `-w` and a default LCA cut-off of 5 `-c 5`):

```
tronko-assign -r -f outdir_multiple_MSA/reference_tree.txt -s -g example_datasets/multiple_trees/missingreads_singleend_150bp_2error.fasta -o example_datasets/multiple_trees/missingreads_singleend_150bp_2error_multiple_MSA_results.txt -a tronko-build/example_datasets/multiple_trees/one_MSA/1.fasta
```

To assign paired-end reads on the same reference database run `tronko-assign`:
```
tronko-assign -r -f outdir_multiple_MSA/reference_tree.txt -p -1 example_datasets/multiple_trees/missingreads_pairedend_150bp_2error_read1.fasta -2 example_datasets/multiple_trees/missingreads_pairedend_150bp_2error_read2.fasta -o example_datasets/multiple_trees/missingreads_pairedend_150bp_2error_multiple_MSA_results.txt -a -a tronko-build/example_datasets/multiple_trees/one_MSA/1.fasta
```

## More on `tronko-build` Usage with multiple trees

`tronko-build` requires a multiple sequence alignment (FASTA format), rooted phylogenetic tree (Newick format), and a corresponding taxonomy file for each cluster build. All of the files should be in one directory and specify the directory with `-e` with each cluster being designated by a number. MSA files should be named `[Number]_MSA.fasta`, taxonomy files should be named `[Number]_taxonomy.txt`, and tree files should be named `RAxML_bestTree.[Number].reroot`. Example of the contents of a directory containing 3 clusters:
```
1_MSA.fasta
2_MSA.fasta
3_MSA.fasta
1_taxonomy.txt
2_taxonomy.txt
3_taxonomy.txt
RAxML_bestTree.1.reroot
RAxML_bestTree.2.reroot
RAxML_bestTree.3.reroot
```
Once you have the cluster files prepared and you have created an output directory (the program assumes the output directory already exists), an example command is
```
tronko-build -y -e [DIRECTORY CONTAINING MSA, TAX, and TREE FILES] -n [NUMBER OF PARTITIONS] -d [OUTPUT DIRECTORY] -s
```
For the example with 3 clusters, an example command partitioning by sum-of-pairs score would be:
```
tronko-build -y -e [DIRECTORY CONTAINING MSA, TAX, and TREE FILES] -n 3 -d output -s
```
The reference database file will be output to `[OUTPUT DIRECTORY]/reference_tree.txt`. The `reference_tree.txt` file is the reference database file that `tronko-assign` requires for assignment.

To partition the reference database further (see manuscript for details) the following dependencies are needed to be installed in your path (no dependencies needed for `tronko-assign` and for building a reference database with a single tree or without partitioning the database): <a href="https://github.com/stamatak/standard-RAxML">`raxmlHPC-PTHREADS`</a>, <a href="https://github.com/refresh-bio/FAMSA">`famsa`</a>, `nw_reroot` from <a href="https://anaconda.org/bioconda/newick_utils/files">Newick utilties</a>, <a href="https://raw.githubusercontent.com/lpipes/tronko/main/scripts/fasta2phyml.pl">`fasta2phyml.pl`</a>, and <a href="https://ftp.gnu.org/gnu/sed/">`sed`</a>. Partitioning the database further is only needed when the underlying MSA is unreliable. An example command to create the reference database and partition a database that contains 100 initial clusters using the sum-of-squares approach:

```
tronko-build -y -e initial_clusters_directory -d outdir -n 100 -s
```

The `reference_tree.txt` file will be output to the `outdir` directory. An example command to create the reference database and partition a database that contains 100 initial clusters using a threshold for the number of leaf nodes (i.e., 500 leaf nodes):

```
tronko-build -y -e initial_clusters_directory -d outdir -n 100 -v -f 500
```

The `reference_tree.txt` file will be output to the `outdir` directory.

# `tronko-assign` using pre-built 16S and COI databases

First download the databases from Zenodo (<a href="https://doi.org/10.5281/zenodo.7407318"><img src="https://zenodo.org/badge/DOI/10.5281/zenodo.7407318.svg" alt="DOI"></a>). To assign FASTQ (`-q`) paired-end reads (`-p`) using the 16S database (and reverse-complement the reverse read `-z`) with Needleman-Wunsch (`-w`), 16 cores (`-C 16`), and an LCA cut-off of 5 (`-c 5`):

```
tronko-assign -r -q -p -z -w -C 16 -c 5 -f 16S_tronko_build.txt.gz -a 16S.fasta -1 16S_TW-DR-1-S88_F_filt.fastq.gz -2 16S_TW-DR-1-S88_R_filt.fastq.gz -o 16S_TW-DR-1-S88_results.txt 
```  

With the CO1 database and same parameters:
```
tronko-assign -r -q -p -z -w -C 16 -c 5 -f CO1_tronko_build.txt.gz -a CO1.fasta -1 CO1_TW-DR-1-S88_F_filt.fastq.gz -2 CO1_TW-DR-1-S88_R_filt.fastq.gz -o CO1_TW-DR-1-S88_results.txt
```

# Performance

We performed a leave-one-species-out test comparing Tronko (with LCA cut-offs for the score of 0, 5, 10, 15, and 20 with Needleman-Wunsch alignment) to kraken2, metaphlan2, and MEGAN for 1,467 COI sequences from 253 species from the order Charadriiformes using 150bp x 2 paired-end sequences and 150bp and 300bp single-end sequences using 0, 1, and 2% error/polymorphism.
<img src="https://github.com/lpipes/tronko/blob/main/LSO.png?raw=true">
Using leave-one-species-out and simulating reads (both paired-end and single-end) with a 0-2% error (or polymorphism), Tronko detected the correct genus more accurately than the other methods even when using an aggressive cut-off (i.e., when cut-off=0) (D and G).

# Citation

Pipes L, and Nielsen R (2022) A rapid phylogeny-based method for accurate community profiling of large-scale metabarcoding datasets. bioRXiv.
https://www.biorxiv.org/content/10.1101/2022.12.06.519402v1 
