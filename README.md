# tronko
A rapid phylogeny-based method for accurate community profiling of large-scale metabarcoding datasets

In the tronko package there are two modules: `tronko-build` and `tronko-assign`. `tronko-build` is for building custom reference databases that tronko-assign uses as input. We have two reference databases currently available for download with `tronko-assign`. Cytochrome oxidase I (COI) which was custom built with <a href="https://github.com/limey-bean/CRUX_Creating-Reference-libraries-Using-eXisting-tools">CRUX</a> using forward primer `GGWACWGGWTGAACWGTWTAYCCYCC` and reverse primer `TANACYTCnGGRTGNCCRAARAAYCA`. 16S which was custom built with <a href="https://github.com/limey-bean/CRUX_Creating-Reference-libraries-Using-eXisting-tools">CRUX</a> using forward primer `GTGCCAGCMGCCGCGGTAA` and reverse primer `GACTACHVGGGTATCTAATCC`. 
<a href="https://doi.org/10.5281/zenodo.7407318"><img src="https://zenodo.org/badge/DOI/10.5281/zenodo.7407318.svg" alt="DOI"></a>

Alignment-based and composition-based assignment methods calculate the lowest common ancestor (LCA) using data only in the leaf nodes of a phylogeny (A). The advantage of Tronko is that it stores fractional likelihoods in all nodes of a phylogeny and calculates the LCA based on all nodes in the tree (B).
<img src="https://github.com/lpipes/tronko/blob/main/Overview_Figure.jpg?raw=true">

# tronko-build
`tronko-build` is for building custom reference databases to be used with `tronko-assign`.

	tronko-build [OPTIONS]
	
		-h, --help				usage:
		-y, --partition-directory		use a partition directory (you have multiple clusters)
		-l, --single-tree			use only single tree (do not partition)
		-t, --tree-file				rooted phylogenetic tree [FILE: Newick]
		-m, --msa-file				multiple sequence alignment [FILE: FASTA]
		-d, --partitions-directory		output directory for partitions
		-x, --tax-file				taxonomy file [FILE: FASTA_header\tdomain;phylum;class;order;family;genus;species]
		-e, --read-directory			directory for multiple cluster
		-n, --number-of-partitions		number of partitions in read directory
		-b, --where-to-restart-partitions	restart partitions with partition number
		-s, --use-spscore			partition using sum-of-pairs score [can't use with -f]
		-u, --sum-of-pairs			minimum threshold for sum of pairs score [default: 0.1]
		-v, --use-minleaves			partition using minimum number of leaf nodes [can't use with -s, use with -f]
		-f, --minimum-leaf-nodes to retain	don't partition less than the minimum number of leaf nodes [can't use with -s, use with -v]
		-g, --no-change-missingdata		don't flag missing data

# tronko-assign
`tronko-assign` is for species assignment of queries. It requires a `tronko-build` database.

	tronko-assign [OPTIONS]
	
		-h, --help			usage: [-paired] [-single] [-reference] [-ntree]
		-p, --paired			use paired reads
		-s, --single			use single reads
		-r, --reference			use a reference
		-v, --reverse-single-read	when using single read reverse it
		-z, --reverse-paired-read	when using pairs reverse second read
		-f, --reference-file		path to reference file
		-o, --results			path to output file
		-g, --single-read-file		path to single read file
		-1, --paired-read-file1		path to paired read 1 file
		-2, --paired-read-file2		path to paired read 2 file
		-a, --fasta-file		path to fasta file (for bwa database)
		-c, --cinterval			score cut-off to use [default:5]
		-C, --number-of-cores		number of cores
		-L, --number-of-lines-to-read	number of lines to read for assignment
		-P, --print-alignments		print alignments to stdout
		-w, --use-nw			use Needleman-Wunsch
		-q, --fastq			Query is FASTQ [default is FASTA]
		-e, --use-leaf-portion		Use only a portion of leaf
		-n, --padding [INT]		Padding to use in leaf portion
		-5, --print-node-info		[FILE] Print tree number and leaf number
		-6, --skip-bowtie2-build	Skip the bowtie2 build
		-u, --score-constant		Score constant [default: 0.01]

# INSTALLATION

	cd tronko-build
	make
	../tronko-assign
	make

# `tronko-assign` Usage
Assigning paired-end reads in FASTA format
```
tronko-assign -r -f [tronko-build REFERENCE DB FILE] -p -1 [FORWARD READS FASTA] -2 [REVERSE READS FASTA] -a [REFERENCE SEQUENCES FASTA] -o [OUTPUT FILE]
```
Assigning single-end reads in FASTA format
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

## `tronko-build` Usage with multiple trees

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
Once you have the cluster files prepared, an example command is
```
tronko-build -y -e [DIRECTORY CONTAINING MSA, TAX, and TREE FILES] -n [NUMBER OF PARTITIONS] -d [OUTPUT DIRECTORY] -s
```
For the example with 3 clusters, an example command partitioning by sum-of-pairs score would be:
```
tronko-build -y -e [DIRECTORY CONTAINING MSA, TAX, and TREE FILES] -n 3 -d output -s
```

# Performance
<img src="https://github.com/lpipes/tronko/blob/main/LSO.png?raw=true">

# Citation
