# tronko
A rapid phylogeny-based method for accurate community profiling of large-scale metabarcoding datasets

In the tronko package there are two modules: `tronko-build` and `tronko-assign`. `tronko-build` is for building custom reference databases that tronko-assign uses as input. We have two reference databases currently available for download with `tronko-assign`. Cytochrome oxidase I (COI) which was custom built with <a href="https://github.com/limey-bean/CRUX_Creating-Reference-libraries-Using-eXisting-tools">CRUX</a> using forward primer `GGWACWGGWTGAACWGTWTAYCCYCC` and reverse primer `TANACYTCnGGRTGNCCRAARAAYCA`. 16S which was custom built with <a href="https://github.com/limey-bean/CRUX_Creating-Reference-libraries-Using-eXisting-tools">CRUX</a> using forward primer `GTGCCAGCMGCCGCGGTAA` and reverse primer `GACTACHVGGGTATCTAATCC`. 

Alignment-based and composition-based assignment methods calculate the lowest common ancestor (LCA) using data only in the leaf nodes of a phylogeny (A). The advantage of Tronko is that it stores fractional likelihoods in all nodes of a phylogeny and calculates the LCA based on all nodes in the tree (B).
<img src="https://github.com/lpipes/tronko/blob/main/Overview_Figure.jpg?raw=true">

# tronko-build
`tronko-build` is for building custom reference databases to be used with `tronko-assign`.

	tronko_build [OPTIONS]
	
		-h, --help				usage: [-paired] [-single] [-reference] [-ntree]
		-y, --partition-directory		use a partition directory
		-u, --sum-of-pairs			minimum sum of pairs score
		-f, --minimum-leaf-nodes to retain	don't partition less than the minimum
		-l, --single-tree			use only single tree (do not partition)
		-t, --tree-file				path to treefile
		-m, --msa-file				path to msa file
		-d, --partitions-directory		path to output partitions files
		-x, --tax-file				path to tax file
		-e, --read-directory			path to directory to read
		-n, --number-of-partitions		number of partitions in read directory
		-b, --where-to-restart-partitions
		-s, --use-spscore
		-v, --use-minleaves
		-g, --no-change-missingdata


# tronko-assign
`tronko-assign` is for species assignment of queries. It requires a `tronko-build` database.

	tronko-assign [OPTIONS]
	
		-h, --help			usage: [-paired] [-single] [-reference] [-ntree]
		-p, --paired			use paired reads
		-s, --single			use single reads
		-r, --reference			use a reference
		-q, --fastq			reads to assign are fastq
		-y, --partition-directory	use a partition directory
		-v, --reverse-single-read	when using single read reverse it
		-z, --reverse-paired-read	when using pairs reverse second read
		-f, --reference-file		path to reference file
		-t, --tree-file			path to treefile
		-m, --msa-file			path to msa file
		-d, --partitions-directory	path to output partitions files
		-o, --results			path to output file
		-x, --tax-file			path to tax file
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
		-U, --print-unassigned		Print unassigned reads
		-e, --use-leaf-portion		Use only a portion of leaf
		-n, --padding [INT]		Padding to use in leaf portion
		-5, --print-node-info		[FILE] Print tree number and leaf number
		-6, --skip-bowtie2-build	Skip the bowtie2 build
		-u, --score-constant		Score constant [default: 0.01]

# Performance

# Citation
