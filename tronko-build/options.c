#include "options.h"

static struct Options long_options[]=
{
	{"help", no_argument, 0, 'h'},
	{"partition-directory", no_argument, 0, 'y'},
	{"sum-of-pairs", required_argument, 0, 'u'},
	{"single-tree", no_argument, 0, 'l'},
	{"tree-file", required_argument, 0, 't'},
	{"msa-file", required_argument, 0, 'm'},
	{"tax-file", required_argument, 0, 'x'},
	{"partitions-directory", required_argument, 0, 'd'},
	{"read-directory", required_argument, 0, 'e'},
	{"number-of-partitions", required_argument, 0, 'n'},
	{"where-to-restart-partitions", required_argument, 0, 'b'},
	{"minimum-leaf-nodes", required_argument, 0, 'f'},
	{"use-spscore", no_argument, 0, 's'},
	{"use-minleaves", no_argument, 0, 'v'},
	{"no-change-missingdata", no_argument, 0, 'g'},
};

char usage[] = "\ntronko-build [OPTIONS]\n\
	\n\
	-h, --help				usage:\n\
	-y, --partition-directory		use a partition directory (you have multiple clusters)\n\
	-l, --single-tree			use only single tree (do not partition)\n\
	-t, --tree-file				rooted phylogenetic tree [FILE: Newick]\n\
	-m, --msa-file				multiple sequence alignment [FILE: FASTA]\n\
	-d, --partitions-directory		output directory for partitions\n\
	-x, --tax-file				taxonomy file [FILE: FASTA_header\tdomain;phylum;class;order;family;genus;species]\n\
	-e, --read-directory			directory for multiple cluster\n\
	-n, --number-of-partitions		number of partitions in read directory\n\
	-b, --where-to-restart-partitions	restart partitions with partition number [default: 0]\n\
	-s, --use-spscore			partition using sum-of-pairs score [can't use with -f]\n\
	-u, --sum-of-pairs			minimum threshold for sum of pairs score [default: 0.5]\n\
	-v, --use-minleaves			partition using minimum number of leaf nodes [can't use with -s, use with -f]\n\
	-f, --minimum-leaf-nodes to retain	don't partition less than the minimum number of leaf nodes [can't use with -s, use with -v]\n\
	-g, --no-change-missingdata		don't flag missing data\n\
	\n";

void print_help_statement(){
	printf("%s", &usage[0]);
	return;
}

void parse_options(int argc, char **argv, Options *opt){
	int option_index, success;
	char c;
	if (argc==1){
		print_help_statement();
		exit(0);
	}
	while(1){
		c=getopt_long(argc,argv,"hspvlgryu:t:m:d:o:x:1:2:a:c:e:n:b:D:f:",long_options, &option_index);
		if (c==-1) break;
		switch(c){
			case 'h':
				print_help_statement();
				exit(0);
				break;
			case 'r': //--reference
				opt->reference_mode = 1;
				break;
			case 'l': //--single-tree
				opt->number_of_trees = 1;
				break;
			case 's':
				opt->use_spscore = 1;
				break;
			case 'v':
				opt->use_min_leaves = 1;
				break;
			case 'f':
				success = sscanf(optarg, "%d", &(opt->min_leaves));
				if (!success)
					fprintf(stderr, "Could not read min leaves\n");
				break;
			case 'y':
				opt->use_partitions = 1;
				break;
			case 'u':
				success = sscanf(optarg, "%lf", &(opt->sp_score));
				if (!success)
					fprintf(stderr, "Could not read sum of pairs score\n");
				break;
			case 't':
				success = sscanf(optarg, "%s", opt->tree_file);
				if (!success)
					fprintf(stderr, "Invalid tree file.\n");
				break;
			case 'm':
				success = sscanf(optarg, "%s", opt->msa_file);
				if (!success)
					fprintf(stderr, "Invalid MSA file\n");
				break;
			case 'd':
				success = sscanf(optarg, "%s", opt->partitions_directory);
				if (!success)
					fprintf(stderr, "Invalid partitions output directory\n");
				break;
			case 'o':
				success = sscanf(optarg, "%s", opt->results_file);
				if (!success)
					fprintf(stderr, "Invalid output results file\n");
				break;
			case 'x':
				success = sscanf(optarg, "%s", opt->taxonomy_file);
				if (!success)
					fprintf(stderr, "Invalid taxonomy file\n");
				break;
			case 'g':
				opt->missing_data = 0;
				break;
			case '1':
				success = sscanf(optarg, "%s", opt->read1_file);
				if (!success)
					fprintf(stderr, "Invalid read 1 file.\n");
				break;
			case '2':
				success = sscanf(optarg, "%s", opt->read2_file);
				if (!success)
					fprintf(stderr, "Invalid read 2 file.\n");
				break;
			case 'a':
				success = sscanf(optarg, "%s", opt->fasta_file);
				if (!success)
					fprintf(stderr, "Invalid fasta file.\n");
				break;
			case 'c':
				success = sscanf(optarg, "%lf", &(opt->cinterval));
				if (!success)
					fprintf(stderr, "Invalid interval\n");
				break;
			case 'e':
				success = sscanf(optarg, "%s", opt->readdir);
				if (!success)
					fprintf(stderr, "Invalid directory");
				break;
			case 'D':
				success = sscanf(optarg, "%s", opt->rmrefdir);
				if (!success)
					fprintf(stderr, "Invalid directory");
				break;
			case 'n':
				success = sscanf(optarg, "%d", &(opt->number_of_partitions));
				if (!success)
					fprintf(stderr, "Invalid directory");
				break;
			case 'b':
				success = sscanf(optarg, "%d", &(opt->restart));
				if (!success)
					fprintf(stderr, "Invalid number\n");
				break;
		}
	}
}
