
#include "granet.h"
#include <getopt.h>
#include <string.h> 
#include <errno.h>
#include <sys/stat.h>

/* MACRO used to change the color of the standard output */
#define TXT_RED                     "\033[0;31m"
#define TXT_NO_COLOR                "\033[0m"
#define TXT_GREEN                   "\033[0;32m"
#define TXT_BLUE                    "\033[0;34m"
#define TXT_YELLOW                  "\033[0;33m"
#define TXT_CYAN                    "\033[0;36m"
#define TXT_PURPLE                  "\033[0;35m"



#define USAGE "\nUsage: tutorial_graphs --graph <FILE_NAME> --format <TYPE> --graph-type <GRAPH_TYPE> --output-dir <DIR_NAME>\n\n"\
			  "\t-f, --graph <FILE_NAME>        Read the graph from file <FILE_NAME>.\n"\
			  "\t-o, --output-dir <DIR_NAME>    Write output to dir <DIR_NAME>.\n"\
			  "\t-t, --format <TYPE>            The format of the file containing the graph, admitted types are:\n"\
			  "\t                               0 - (Edges List) The first line of the file must contain the \n"\
			  "\t                                   number of nodes and edges of the graph in the following format: \n"\
			  "\t                                   \"# Nodes: <NUMBER> Edges: <NUMBER>\" \n"\
			  "\t-g, --graph-type <GRAPH_TYPE>  Admitted types are:\n"\
			  "\t                               D - Directed graph\n"\
			  "\t                               U - Undirected graph\n\n"
				  


int main(int argc, char **argv){

	int file_type = -1;					/* Format of the file containing the graph */
	char ch;                    		/* Option read by getopt_long */ 
	char *file_name = NULL;				/* Name of the file containing the graph */   
	char type = 0;
	char *output_dir = NULL;

	static struct option long_options[] ={
	    {"graph", required_argument, NULL, 'f'},
	    {"format", required_argument, NULL, 't'},
	    {"graph-type", required_argument, NULL, 'g'},
	    {"output-dir", required_argument, NULL, 'o'},
		{"help", no_argument, NULL, 'h'},
	    {NULL, 0, NULL, 0}
	};

	/* getopt_long() returns the option character when a short option is recognized. 
	 * For a long option, it returns <val> if <flag> is NULL, and 0 otherwise. */
	while ((ch = getopt_long(argc, argv, "f:t:g:o:h", long_options, NULL)) != -1){
	    switch (ch){
		    case 'f':
	            file_name = strdup(optarg);
	            break;		    
			case 'o':
	            output_dir = strdup(optarg);
				if ( mkdir(output_dir, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH ) != 0 ){
					if ( errno != EEXIST){
						exit(EXIT_FAILURE);
					}
				}
	            break;
			case 't':
	            file_type = atoi(optarg); /* Check that the file format is one of the allowed types */
				if ( file_type != EDGE_LIST ){
	                printf(USAGE);
	                exit(EXIT_FAILURE);
				}
	            break;
			case 'g':
			    type = optarg[0];
				if ( type != 'D' && type != 'U'){
	                printf(USAGE);
	                exit(EXIT_FAILURE);
				}
	            break;
			case 'h':
	           default:
	               printf(USAGE);
	               exit(EXIT_FAILURE);
	               break;
	    }
	}
	/* Check that all needed information have been provided through the command line */
	if ( file_type == -1 || file_name == NULL || type == 0 || output_dir == NULL) {
	      printf(USAGE);
	      exit(EXIT_FAILURE);
	}
	if (type == 'D'){
	    printf("This program does not support Undirected Graphs\n");
		exit(EXIT_FAILURE);
	}else if(type != 'U') { 
	      printf(USAGE);
	      exit(EXIT_FAILURE);
	}
	
	
	graph_type *g, *ng;
	
	/* Read the graph */
	fprintf(stderr, "(%s) Read Graph from file: %s\n", get_time_string(), file_name);
	g = get_graph_from_file( file_name, file_type);
	if (g == NULL){ 
		exit(EXIT_FAILURE); 
	}
	print_graph(g, 0);
	printf("\n");
	
	/* Write digraph */
	char out_file[256] = "graph.out.edgelist";
	const char *path;
	path = get_path(output_dir, out_file);
	write_graph_edgelist(g, path);
	fprintf(stderr, "(%s) Write Graph to file: %s\n", get_time_string(), path);
	
	/* Remove Nodes */
	fprintf(stderr, "(%s) Delete nodes from graph\n", get_time_string());
	vertex_type vts[2] = {1,2};
	ng = delete_nodes_from_graph ( g, vts, 2);
	print_graph(ng, 0);
	print_graph(ng, 1);
	free_graph(ng);
	printf("\n");
	
	/* Compute Degree*/
	fprintf(stderr, "(%s) Compute nodes degree\n", get_time_string());
	degree_type *deg;
	deg = get_degree(g);
	print_degree(deg, g->nodes, g->nodes);
	printf("\n");
		
	sort_degree(deg, g->nodes, ASC);
	print_degree(deg, g->nodes, 3);
	printf("\n");
	
	sort_degree(deg, g->nodes, DESC);
	print_degree(deg, g->nodes, 3);
	printf("\n");


	/* Compute Connected Components*/
	fprintf(stderr, "(%s) Compute Connected Components\n", get_time_string());
	cc_type *cc;
	cc = get_graph_cc(g);
	print_cc(cc);
	free_cc(cc);
	printf("\n");


	free(file_name);
	free(output_dir);
	return 0;

}














