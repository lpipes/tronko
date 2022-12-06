#ifndef _OPTIONS_
#define _OPTIONS_
#include <stdlib.h>
#include <stdint.h>
#include <stdio.h>
#include <getopt.h>
#include "global.h"

void print_help_statement();
void parse_options(int argc, char **argv, Options *opt);

#endif /* _OPTIONS_ */
