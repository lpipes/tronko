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
	
		-h, --help			usage:
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
		-6, --skip-bwa-build		Skip the bwa build
		-u, --score-constant		Score constant [default: 0.01]

Tronko uses the <a href="https://github.com/smarco/WFA2-lib">Wavefront Alignment Algorithm (version 2)</a> or <a href="https://github.com/noporpoise/seq-align">Needleman-Wunsch Algorithm</a> for semi-global alignments. It uses <a href="https://github.com/lh3/bwa">bwa</a> for alignment to leaf nodes, and uses <a href="https://github.com/DavidLeeds/hashmap">David Leeds' hashmap</a> for hashmap implementation in C.

## Example output

The output file is a tab-delimited text file where only the forward readname is retained (if using paired-end reads). The output displays the taxonomic path for assignment, the score, the number of forward read mismatches with the `bwa` hit, the number of reverse read mismatches with the `bwa` hit, the tree number for the best assignment (0 if using 1 tree), and the node number the read (or reads in the case of paired-end reads) was assigned to.

```
Readname	Taxonomic_Path	Score	Forward_Mismatch	Reverse_Mismatch	Tree_Number	Node_Number
M00160:15:000000000-JHG8V:1:1101:9131:1243_1:N:0:GTGCAGA+TACCATC	Eukaryota;Arthropoda;Insecta;Ephemeroptera;Baetidae;Fallceon;Fallceon sp. BOLD:AAL8084	0.000000	0.000000	1.000000	5718	44
M00160:15:000000000-JHG8V:1:1101:21631:1259_1:N:0:GTGCAGA+TACCATC	Eukaryota;Arthropoda;Insecta;Ephemeroptera;Baetidae;Fallceon;Fallceon sp. BOLD:AAL8084	0.000000	0.000000	1.000000	5718	44
M00160:15:000000000-JHG8V:1:1101:15032:1369_1:N:0:GTGCAGA+TACCATC	Eukaryota;Arthropoda;Insecta;Ephemeroptera;Baetidae;Baetis;Baetis adonis	-4.605170	1.000000	1.000000	8901	66
M00160:15:000000000-JHG8V:1:1101:15129:1391_1:N:0:GTGCAGA+TACCATC	Eukaryota;Arthropoda;Insecta;Ephemeroptera;Baetidae;Baetis;Baetis adonis	0.000000	0.000000	1.000000	8901	66
```

# INSTALLATION

	cd tronko-build
	make
	../tronko-assign
	make

# `tronko-assign` Usage
Tronko does not detect the correct orientation of the reads. If your reverse read needs to be reverse complemented use the option `-z`. The default options of Tronko assume that your reads are in FASTA format. If you want to assign reads in FASTQ format, use the option `-q`. The reads (and reference database file) can be gzipped or unzipped. Assigning paired-end reads in FASTA format:
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
Once you have the cluster files prepared and you have created an output directory (the program assumes the output directory already exists), an example command is
```
tronko-build -y -e [DIRECTORY CONTAINING MSA, TAX, and TREE FILES] -n [NUMBER OF PARTITIONS] -d [OUTPUT DIRECTORY] -s
```
For the example with 3 clusters, an example command partitioning by sum-of-pairs score would be:
```
tronko-build -y -e [DIRECTORY CONTAINING MSA, TAX, and TREE FILES] -n 3 -d output -s
```
The reference database file will be output to `[OUTPUT DIRECTORY]/reference_tree.txt`. The `reference_tree.txt` file is the reference database file that `tronko-assign` requires for assignment.

# Performance

We performed a leave-one-species-out test comparing Tronko (with LCA cut-offs for the score of 0, 5, 10, 15, and 20 with Needleman-Wunsch alignment) to kraken2, metaphlan2, and MEGAN for 1,467 COI sequences from 253 species from the order Charadriiformes using 150bp x 2 paired-end sequences and 150bp and 300bp single-end sequences using 0, 1, and 2% error/polymorphism.
<img src="https://github.com/lpipes/tronko/blob/main/LSO.png?raw=true">
Using leave-one-species-out and simulating reads (both paired-end and single-end) with a 0-2% error (or polymorphism), Tronko detected the correct genus more accurately than the other methods even when using an aggressive cut-off (i.e., when cut-off=0) (D and G).

# Citation

Pipes L, and Nielsen R (2022) A rapid phylogeny-based method for accurate community profiling of large-scale metabarcoding datasets. bioRXiv.
https://www.biorxiv.org/content/10.1101/2022.12.06.519402v1 
