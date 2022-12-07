#include "options.h"

static struct Options long_options[]=
{
	{"help", no_argument, 0, 'h'},
	{"paired", no_argument, 0, 'p'},
	{"single", no_argument, 0, 's'},
	{"reference", no_argument, 0, 'r'},
	{"partition-directory", no_argument, 0, 'y'},
	{"reference-file", required_argument, 0, 'f'},
	{"tree-file", required_argument, 0, 't'},
	{"msa-file", required_argument, 0, 'm'},
	{"tax-file", required_argument, 0, 'x'},
	{"partitions-directory", required_argument, 0, 'd'},
	{"results", required_argument, 0, 'o'},
	{"single-read-file", required_argument, 0, 'g'},
	{"paired-read-file1", required_argument, 0, '1'},
	{"paired-read-file2", required_argument, 0, '2'},
	{"fasta-file", required_argument, 0, 'a'},
	{"Cinterval", required_argument, 0, 'c'},
	{"reverse-single-read",no_argument, 0, 'v'},
	{"reverse-paired-read",no_argument, 0, 'z'},
	{"number-of-cores",required_argument, 0,'C'},
	{"number-of-lines-to-read",required_argument, 0, 'L'},
	{"print-alignments",required_argument, 0, 'P'},
	{"use-nw",no_argument, 0, 'w'},
	{"print-unassigned",no_argument, 0, 'U'},
	{"use-leaf-portion",no_argument, 0, 'e'},
	{"padding",required_argument, 0, 'n'},
	{"fastq",no_argument, 0, 'q'},
	{"print-node-info",required_argument, 0, '5'},
	{"skip-bowtie2-build",no_argument,0, '6'},
	{"score-constant",required_argument,0, 'u'}
};

char usage[] = "\ntronko-assign [OPTIONS]\n\
	\n\
	-h, --help			usage: [-paired] [-single] [-reference] [-ntree]\n\
	-p, --paired			use paired reads\n\
	-s, --single			use single reads\n\
	-r, --reference			use a reference\n\
	-v, --reverse-single-read	when using single read reverse it\n\
	-z, --reverse-paired-read	when using pairs reverse second read\n\
	-f, --reference-file		path to reference file\n\
	-o, --results			path to output file\n\
	-g, --single-read-file		path to single read file\n\
	-1, --paired-read-file1		path to paired read 1 file\n\
	-2, --paired-read-file2		path to paired read 2 file\n\
	-a, --fasta-file		path to fasta file (for bwa database)\n\
	-c, --cinterval			score cut-off to use [default:5]\n\
	-C, --number-of-cores		number of cores\n\
	-L, --number-of-lines-to-read	number of lines to read for assignment\n\
	-P, --print-alignments		print alignments to stdout\n\
	-w, --use-nw			use Needleman-Wunsch\n\
	-q, --fastq			Query is FASTQ [default is FASTA]\n\
	-e, --use-leaf-portion		Use only a portion of leaf\n\
	-n, --padding [INT]		Padding to use in leaf portion\n\
	-5, --print-node-info		[FILE] Print tree number and leaf number\n\
	-6, --skip-bowtie2-build	Skip the bowtie2 build\n\
	-u, --score-constant		Score constant [default: 0.01]\n\
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
		c=getopt_long(argc,argv,"hpsrqw6yevUzP5:f:u:t:m:d:o:x:g:1:2:a:c:n:3:4:C:L:",long_options, &option_index);
		if (c==-1) break;
		switch(c){
			case 'h':
				print_help_statement();
				exit(0);
				break;
			case 'w':
				opt->use_nw=1;
				break;
			case 'p': //--paired
				strcpy(opt->paired_or_single,"paired");
				break;
			case 's': //--single
				strcpy(opt->paired_or_single,"single");
				break;
			case 'q': //--fastq
				opt->fastq=1;
				break;
			case '6':
				opt->skip_build=1;
				break;
			case '5': //print node info
				success = sscanf(optarg, "%s", opt->print_node_info);
				if (!success)
					fprintf(stderr, "Invalid file\n");
				break;
			case 'U':
				opt->unassigned=1;
				break;
			case 'r': //--reference
				opt->reference_mode = 1;
				break;
			case 'P': //--print-alignments
				opt->print_alignments = 1;
				break;
			case 'y':
				opt->use_partitions = 1;
				break;
			case 'v':
				opt->reverse_single_read=1;
				break;
			case 'z':
				opt->reverse_second_of_paired_read=1;
				break;	
			case 'f':
				success = sscanf(optarg, "%s", opt->reference_file);
				if (!success)
					fprintf(stderr, "Invalid reference file.\n");
				break;
			case 'u':
				success = sscanf(optarg, "%lf", opt->score_constant);
				if (!success)
					fprintf(stderr, "Could not read score constant\n");
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
				success = sscanf(optarg, "%s", opt->read1_file);
				if (!success)
					fprintf(stderr, "Invalid single read file\n");
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
			case '3':
				success = sscanf(optarg, "%s", opt->print_alignments_dir);
				opt->print_alignments_to_file = 1;
				if (!success)
					fprintf(stderr, "Invalid directory.\n");
				break;
			case '4':
				success = sscanf(optarg, "%s", opt->treedir);
				if (!success)
					fprintf(stderr, "Invalid directory.\n");
				break;
			case 'c':
				success = sscanf(optarg, "%lf", &(opt->cinterval));
				if (!success)
					fprintf(stderr, "Invalid interval\n");
				break;
			case 'e':
				opt->use_leaf_portion = 1;
				break;
			case 'n':
				success = sscanf(optarg, "%d", &(opt->padding));
				if (!success)
					fprintf(stderr, "Invalid directory");
				break;
			case 'C':
				success = sscanf(optarg, "%d", &(opt->number_of_cores));
				if (!success)
					fprintf(stderr, "Invalid number");
				break;
			case 'L':
				success = sscanf(optarg, "%d", &(opt->number_of_lines_to_read));
				if (!success)
					fprintf(stderr, "Invalid number");
				break;
		}
	}
}
