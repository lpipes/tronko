# tronko
A rapid phylogeny-based method for accurate community profiling of large-scale metabarcoding datasets

In the tronko package there are two modules: `tronko-build` and `tronko-assign`. `tronko-build` is for building custom reference databases that tronko-assign uses as input. We have two reference databases currently available for download with `tronko-assign`. Cytochrome oxidase I (COI) which was custom built with <a href="https://github.com/limey-bean/CRUX_Creating-Reference-libraries-Using-eXisting-tools">CRUX</a> using forward primer `GGWACWGGWTGAACWGTWTAYCCYCC` and reverse primer `TANACYTCnGGRTGNCCRAARAAYCA`. 16S which was custom built with <a href="https://github.com/limey-bean/CRUX_Creating-Reference-libraries-Using-eXisting-tools">CRUX</a> using forward primer `GTGCCAGCMGCCGCGGTAA` and reverse primer `GACTACHVGGGTATCTAATCC`. 

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
