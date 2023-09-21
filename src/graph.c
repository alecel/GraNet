#include "graph.h"


#define BUFFSIZE     1024
/*
 * Function: get_edge_list_from_file
 * ------------------------------------------------
 * The function read a graph file in an edge list format. It returns an array containing
 * the list of edges of the graph read from the file. The function also fill the two parameters
 * <n_nodes> and <n_edges> with the number of nodes and edges of the graph respectivally.
 *
 * @param file_name		The name of the file
 * @param n_nodes     	The parameter will contains the number of nodes in the graph
 * @param n_edges     	The parameter will contains the number of edges in the graph
 *
 * @return          	The function allocates and returns a pointer to an array containing the 
 *						list of edges of the graph. Returns a NULL pointer if an error occurs.
 */
static vertex_type *get_edge_list_from_file(const char *file_name, vertex_type *n_nodes, vertex_type *n_edges){
	FILE *fp=NULL;
	char buffer[BUFFSIZE];
	vertex_type *edge_list=NULL;
	vertex_type edge_idx=0;
	vertex_type self_loops=0;
	int err=0, n_line=0, n_comments=0;
	vertex_type nodes=0; 
	vertex_type edges=0;
	
	fp = Fopen(file_name, "r", __FILE__, __LINE__, "Opening an edge list file.");
	
	while ( fgets(buffer, BUFFSIZE, fp)){ /* READ a LINE */
		int buflen = 0;
		buflen = strlen(buffer);
        if ( buffer[buflen-1] != '\n' && !feof(fp) ) { // Check that the buffer is big enough to read a line
            fprintf(stderr, "(%s)[ERROR] File %s. The line is too long, increase the BUFFSIZE! Exit\n",get_time_string(), file_name);
			fclose(fp);
			return NULL;
        }
		n_line++;
		
		if (strchr(buffer, '#') != NULL) { // The line it's a comment
			n_comments++;
            if (strstr(buffer, "Nodes:") ) { // The line contains the number of Nodes and Edges
                err = sscanf(buffer, "# Nodes: %"PRIvertex" Edges: %"PRIvertex"\n", &nodes, &edges);
				if (err != 2 || err == EOF){
					fprintf(stderr, "(%s)[ERROR] Error reading file %s.\n",get_time_string(), file_name);
					fclose(fp);
					return NULL;
				} /* Allocate the memory for the list of edges */
				edge_list = (vertex_type*) Malloc(sizeof(vertex_type)*edges*2,__FILE__, __LINE__, "Allocating the edge list.");
            }
		}else{ // Read a new edge
			vertex_type i, j;
			if( nodes == 0 || edges == 0){
				fprintf(stderr, "(%s)[ERROR] File %s. You must specify the number of Nodes and Edges in the graph before listing the edges.\n",get_time_string(), file_name);
				fclose(fp);
				free(edge_list);
				return NULL;
			}
			sscanf(buffer, "%"PRIvertex" %"PRIvertex"\n", &i, &j);
			if (err != 2 || err == EOF){
				fprintf(stderr, "(%s)[ERROR] Error reading file %s.\n", get_time_string(), file_name);
				fclose(fp);
				free(edge_list);
				return NULL;
			}
			if (i > nodes || j > nodes){
				fprintf(stderr, "(%s)[ERROR] File %s. Vertex id bigger than the number of vertices in the graph: edge %"PRIvertex" %"PRIvertex", number of nodes %"PRIvertex".\n",get_time_string(), file_name, i, j, nodes);
				fclose(fp);
				free(edge_list);
				return NULL;
			}
			if (i == 0 || j == 0){
				fprintf(stderr, "(%s)[ERROR] File %s. Vertex id 0 not allowed, vertex id must start from 1.\n",get_time_string(), file_name);
				fclose(fp);
				free(edge_list);
				return NULL;
			}
			
			if ( i != j) { // Ignore self-loops
				if( edge_idx == (edges*2) ){
					fprintf(stderr, "(%s)[ERROR] File %s. Too many edges listed in the file (%"PRIvertex"), only %"PRIvertex" were expected\n",get_time_string(), file_name, edge_idx/2, edges );
					fclose(fp);
					free(edge_list);
					return NULL;
				}
				edge_list[edge_idx]=i;
				edge_list[edge_idx+1]=j;
				edge_idx +=2;
			} else { self_loops++; }
		}
	}	
 
    if (ferror(fp)){ /* Check reading errors */
		perror("Reading error");
		fclose(fp);
		return NULL;
	}
		
	if( ((edge_idx/2)+self_loops) > edges ){
		fprintf(stderr, "(%s)[ERROR] File %s. Too many edges listed in the file, only %"PRIvertex" were expected\n",get_time_string(), file_name, edges );
		fclose(fp);
		free(edge_list);
		return NULL;
	}
	fclose(fp);
	
	if ( ((edge_idx/2)+self_loops) < edges ){
		fprintf(stderr, "(%s)[ERROR] File %s. Not enough edges listed in the file, %"PRIvertex" were expected\n",get_time_string(), file_name, edges );
		free(edge_list);
		return NULL;
	} 
	
	*n_edges=edges-self_loops;
	*n_nodes=nodes;
	
	return edge_list;
}




/*
 * Function: get_digraph_from_edge_list
 * ------------------------------------------------
 * The function create a DIGRAPH <digraph_type> from an edge list. The edge list is an array
 * containing the edges of the graph, position i and i+1 are the src and dst of each edge.
 *
 * @param edge_list     	The array containing the list of edges.
 * @param n_nodes			Number of nodes in the graph.
 * @param n_edges			Number of edges in the graph.
 *
 * @return					The function returns a pointer to a digraph of type digraph_type containing
 * 							the graph. Use the function free_digraph to free the memory.       
 */
static digraph_type *get_digraph_from_edge_list(vertex_type *edge_list, vertex_type n_nodes, vertex_type n_edges){
	vertex_type	*in_idx;			
	vertex_type	*in_neighbours;		
	vertex_type *out_idx;			
	vertex_type *out_neighbours;	
	vertex_type *labels;
	digraph_type *g;
	vertex_type *next_in;
	vertex_type *next_out;
	
	// Allocate the memory for the graph
	g = (digraph_type *) Malloc(sizeof(digraph_type),__FILE__, __LINE__, "Graph allocation.");
	g->nodes = n_nodes;
	g->edges = n_edges;
	labels = (vertex_type *) Malloc(sizeof(vertex_type)*(n_nodes+1),__FILE__, __LINE__, "Graph labels.");
	in_idx = (vertex_type *) Calloc((n_nodes+2),sizeof(vertex_type),__FILE__, __LINE__, "Graph in_idx.");
	out_idx = (vertex_type *) Calloc((n_nodes+2),sizeof(vertex_type),__FILE__, __LINE__, "Graph out_idx.");
	in_neighbours = (vertex_type *) Malloc(sizeof(vertex_type)*(n_edges+1),__FILE__, __LINE__, "Graph in_neighbours.");
	out_neighbours = (vertex_type *) Malloc(sizeof(vertex_type)*(n_edges+1),__FILE__, __LINE__, "Graph out_neighbours.");
	// Allocate the memory for temporary structures
	next_in = (vertex_type *) Calloc(n_nodes+1,sizeof(vertex_type),__FILE__, __LINE__, "next_in");
	next_out = (vertex_type *) Calloc(n_nodes+1,sizeof(vertex_type),__FILE__, __LINE__, "next_out");
	
	in_idx[1] = out_idx[1] = 1;
	in_neighbours[0] = out_neighbours[0] = 0;
	
	// We write only indexes between 2 and n+1 of arrays <in_idx> and <out_idx>.
	for (vertex_type i = 0; i < n_edges*2; i=i+2) { // Create in_idx and out_idx arrays 
		in_idx[ edge_list[i+1] + 1 ]++;
		out_idx[ edge_list[i] + 1 ]++;
	}
	for (vertex_type i = 0; i < n_nodes+1; i++) { // Create in_idx and out_idx arrays 
		in_idx[i+1] += in_idx[i];
		out_idx[i+1] += out_idx[i];
		next_in[i] = in_idx[i];
		next_out[i] = out_idx[i];
		labels[i] = i;
	}

	for (vertex_type i = 0; i < n_edges*2; i=i+2) { // Create in_neighbours and out_neighbours arrays
		in_neighbours[ next_in[edge_list[i+1]]++ ] = edge_list[i];
		out_neighbours[ next_out[ edge_list[i]]++ ] = edge_list[i+1];
	}
	free(next_in);
	free(next_out);
	
	g->labels = labels;
	g->in_idx = in_idx;
	g->out_idx = out_idx;
	g->in_neighbours = in_neighbours;
	g->out_neighbours = out_neighbours;
    g->del_nodes_list = NULL;
    g->del_nodes_num = 0;
    g->org_nodes_num = g->nodes;
	
	return g;
}



/*
 * Function: compare_edge
 * ------------------------------------------------     
 * Compare two edges based on the ID of their verticies. 
 *
 * @param edge1   First edge to compare
 * @param edge2   Second edge to compare
 *
 * @return		0 if the edges are equal, 1 if edge1 is greater than edge2 
 *              and -1 if edge2 is greater than edge1
 */
static int compare_edge(const void *edge1, const void *edge2) {
    vertex_type *ed1 = (vertex_type *) edge1;
    vertex_type *ed2 = (vertex_type *) edge2;

    if (ed1[0] < ed2[0]) return -1;
    if (ed1[0] > ed2[0]) return  1;
    if (ed1[1] < ed2[1]) return -1;
    if (ed1[1] > ed2[1]) return  1;
    return 0;
}
/*
 * Function: remove_multiple_edges
 * ------------------------------------------------     
 * The function removes multiple edges, keeping only one in case of duplicates
 * between nodes.
 *
 * @param edge_list   An array containing the list of edges, each edge is stored 
 *                    in two adjacent position of the array
 * @param n_edges     Number of edges
 *
 * @return		Number of edges kept after removing duplicates.
 */
static vertex_type remove_multiple_edges(vertex_type *edge_list, vertex_type n_edges){
	vertex_type kept, i;
		
	kept=0;	
	qsort(edge_list, n_edges, sizeof(vertex_type[2]), compare_edge);
	
    for( i=0; i < (n_edges-1)*2; i=i+2){
		if ( edge_list[i]==edge_list[i+2] && edge_list[i+1]==edge_list[i+3]){ continue; }		
        //fprintf(stderr,"%"PRIvertex" - %"PRIvertex"\n", i, kept);
		edge_list[kept] = edge_list[i];
		edge_list[kept+1] = edge_list[i+1];
		kept = kept+2;
	}
    //fprintf(stderr,"%"PRIvertex" - %"PRIvertex"\n", i, kept);
	edge_list[kept] = edge_list[i];
	edge_list[kept+1] = edge_list[i+1];
	kept = kept+2;
	
	//fprintf(stderr, "(%s)[INFO] Number of multiple edges removed %"PRIvertex"\n", get_time_string(), n_edges-(kept/2));
	return kept/2;
}

/*
 * Function: get_digraph_from_file
 * ------------------------------------------------
 * The function create a DIGRAPH <digraph_type> from a file. The parameter <file_format>
 * specifies the format of the file. If the file format is not supported a NULL pointer is returned
 *
 * @param file_name     	The name of the file.
 * @param file_format		The format of the file.
 *
 * @return					The function returns a pointer to a digraph of type digraph_type containing
 * 							the graph. Use the function free_digraph to free the memory. If the file format 
 *                          is not supported or an error occurs a NULL pointer is returned      
 */
digraph_type *get_digraph_from_file(const char *file_name, unsigned char file_format){
	vertex_type n_nodes; 
	vertex_type n_edges;
	vertex_type *edge_list;
	digraph_type *g;
	
	//fprintf(stderr, "(%s)[INFO] Opening file %s\n",get_time_string(), file_name);
	switch (file_format){
		case EDGE_LIST:
			edge_list = get_edge_list_from_file(file_name, &n_nodes, &n_edges);
			break;
		default:
		    fprintf(stderr, "(%s)[ERROR] Unknown input file format %u\n",get_time_string(), file_format);
			return NULL;
	}
	if (edge_list == NULL){ /* An error occurred while reading the file */
		return NULL;
	}
	//fprintf(stderr, "(%s)[INFO] File %s. The graph has %"PRIvertex" Nodes and %"PRIvertex" Edges\n",get_time_string(), file_name, n_nodes, n_edges);
	n_edges = remove_multiple_edges(edge_list, n_edges);
	//fprintf(stderr, "(%s)[INFO] File %s. The graph has %"PRIvertex" Nodes and %"PRIvertex" Edges\n",get_time_string(), file_name, n_nodes, n_edges);
	g = get_digraph_from_edge_list(edge_list, n_nodes, n_edges);
	free(edge_list);
	//fprintf(stderr, "(%s)[DEBUG] Graph building completed. \n",get_time_string());
	
	return g;
}



/*
 * Function: get_graph_from_edge_list
 * ------------------------------------------------
 * The function create a GRAPH <graph_type> from an edge list. The edge list is an array
 * containing the edges of the graph, position i and i+1 are the src and dst of each edge.
 *
 * @param edge_list     	The array containing the list of edges.
 * @param n_nodes			Number of nodes in the graph.
 * @param n_edges			Number of edges in the graph.
 *
 * @return					The function returns a pointer to a graph of type graph_type containing
 * 							the graph. Use the function free_graph to free the memory.       
 */
static graph_type *get_graph_from_edge_list(vertex_type *edge_list, vertex_type n_nodes, vertex_type n_edges){
	vertex_type *idx;			
	vertex_type *neighbours;	
	vertex_type *labels;
	graph_type *g;
	vertex_type *next;
	
	// Allocate the memory for the graph
	g = (graph_type *) Malloc(sizeof(graph_type),__FILE__, __LINE__, "Graph allocation.");
	g->nodes = n_nodes;
	g->edges = n_edges;
	labels = (vertex_type *) Malloc(sizeof(vertex_type)*(n_nodes+1),__FILE__, __LINE__, "Graph labels.");
	idx = (vertex_type *) Calloc((n_nodes+2),sizeof(vertex_type),__FILE__, __LINE__, "Graph idx.");
	neighbours = (vertex_type *) Malloc(sizeof(vertex_type)*(2*n_edges+1),__FILE__, __LINE__, "Graph neighbours.");
	// Allocate the memory for temporary structures
	next = (vertex_type *) Calloc(n_nodes+1,sizeof(vertex_type),__FILE__, __LINE__, "next_in");
	
	idx[1] = 1;
	neighbours[0] = 0;
	
	// We write only indexes between 2 and n+1 of arrays <in_idx> and <out_idx>.
	for (vertex_type i = 0; i < n_edges*2; i=i+2) { // Create idx array 
		idx[ edge_list[i+1] + 1 ]++;
		idx[ edge_list[i] + 1 ]++;
	}
	for (vertex_type i = 0; i < n_nodes+1; i++) { // Create idx array
		idx[i+1] += idx[i];
		next[i] = idx[i];
		labels[i] = i;
	}
	for (vertex_type i = 0; i < n_edges*2; i=i+2) { // Create neighbours arrays
		neighbours[next[edge_list[i+1]]++] = edge_list[i];
		neighbours[next[edge_list[i]]++] = edge_list[i+1];
	}
	free(next);
	
	g->labels = labels;
	g->idx = idx;
	g->neighbours = neighbours;
    g->del_nodes_list = NULL;
    g->del_nodes_num = 0;
    g->org_nodes_num = g->nodes;
	
	return g;
}



/*
 * Function: get_graph_from_file
 * ------------------------------------------------
 * The function create a GRAPH <graph_type> from a file. The parameter <file_format>
 * specifies the format of the file. If the file format is not supported a NULL pointer is returned
 *
 * @param file_name     	The name of the file.
 * @param file_format		The format of the file.
 *
 * @return					The function returns a pointer to a graph of type graph_type containing
 * 							the graph. Use the function free_graph to free the memory. If the file format 
 *                          is not supported or an error occurs a NULL pointer is returned      
 */
graph_type *get_graph_from_file(const char *file_name, unsigned char file_format){
	vertex_type n_nodes; 
	vertex_type n_edges;
	vertex_type *edge_list;
	graph_type *g;
	
	//fprintf(stderr, "(%s)[INFO] Opening file %s\n", get_time_string(),file_name);
	switch (file_format){
		case EDGE_LIST:
			edge_list = get_edge_list_from_file(file_name, &n_nodes, &n_edges);
			break;
		default:
			fprintf(stderr, "(%s)[ERROR] Unknown input file format %u\n", get_time_string(),file_format);
			return NULL;
	}
	if (edge_list == NULL){ /* An error occurred while reading the file */
		return NULL;
	}
	//fprintf(stderr, "(%s)[INFO] File %s. The graph has %"PRIvertex" Nodes and %"PRIvertex" Edges\n", get_time_string(), file_name, n_nodes, n_edges);
	n_edges = remove_multiple_edges(edge_list, n_edges);
	//fprintf(stderr, "(%s)[INFO] File %s. The graph has %"PRIvertex" Nodes and %"PRIvertex" Edges\n", get_time_string(), file_name, n_nodes, n_edges);
	g = get_graph_from_edge_list(edge_list, n_nodes, n_edges);
	free(edge_list);
	
	return g;
}


/*
 * Function: print_graph
 * ------------------------------------------------
 * The function prints the information of the graph: number of nodes, number of edges
 * list of edges and list of labels. 
 *
 * @param g     	    The pointer to the graph to print. 
 * @param lable			If 1 print the node original id, if 0 print the current node id.   
 */
void print_graph(graph_type *g, char label){
	if (g == NULL){ return; }
	fprintf(stdout, TXT_GREEN"GRAPH:"TXT_YELLOW" Nodes: "TXT_NO_COLOR"%"PRIvertex" "TXT_YELLOW" Edges: "TXT_NO_COLOR"%"PRIvertex"\n", g->nodes, g->edges);
	fprintf(stdout, TXT_GREEN"\tEdges: \n\t"TXT_NO_COLOR);
	for (vertex_type i = 1; i < g->nodes+1; i++) {
		for(vertex_type j = g->idx[i]; j < g->idx[i+1]; j++){
			if (g->neighbours[j] > i){
				if ( label == 0 ){
					fprintf(stdout, "(%"PRIvertex",%"PRIvertex")", i, g->neighbours[j]);
				}else{
					fprintf(stdout, "(%"PRIvertex",%"PRIvertex")", g->labels[i], g->labels[g->neighbours[j]]);
				}
				fprintf(stdout, " ");
			}
		}
	}
	fprintf(stdout, "\n");
	
	fprintf(stdout, TXT_GREEN"\tConnections: \n\t"TXT_NO_COLOR);
	for (vertex_type i = 1; i < g->nodes+1; i++) {
		if ( label == 0 ){
			fprintf(stdout, "  %"PRIvertex"->{", i);
		}else{
			fprintf(stdout, "  %"PRIvertex"->{", g->labels[i]);
		}
		for(vertex_type j = g->idx[i]; j < g->idx[i+1]; j++){
			if ( label == 0 ){
				fprintf(stdout, "%"PRIvertex, g->neighbours[j]);
			}else{
				fprintf(stdout, "%"PRIvertex, g->labels[g->neighbours[j]]);
			}
			if ( (j+1) < g->idx[i+1] ) {
				fprintf(stdout, ",");
			}
		}
		fprintf(stdout, "}");
	}
	fprintf(stdout, "\n");
	
	fprintf(stdout, TXT_GREEN"\tLabels: \n\t"TXT_NO_COLOR);
	for (vertex_type i = 1; i < g->nodes+1; i++) {
		fprintf(stdout, "  %"PRIvertex"=>%"PRIvertex, i, g->labels[i]);
	}
	fprintf(stdout, "\n");
}


/*
 * Function: print_digraph
 * ------------------------------------------------
 * The function prints the information of the digraph: number of nodes, number of edges
 * list of edges and list of labels. 
 *
 * @param g     	    The pointer to the digraph to print. 
 * @param label			If 1 print the node original id, if 0 print the current node id.   
 */
void print_digraph(digraph_type *g, char label){
	if (g == NULL){ return; }
	fprintf(stdout, TXT_GREEN"DIGRAPH:"TXT_YELLOW" Nodes: "TXT_NO_COLOR"%"PRIvertex" "TXT_YELLOW" Edges: "TXT_NO_COLOR"%"PRIvertex"\n", g->nodes, g->edges);
	fprintf(stdout, TXT_GREEN"\tEdges: \n\t"TXT_NO_COLOR);
	for (vertex_type i = 1; i < g->nodes+1; i++) {
		for(vertex_type j = g->out_idx[i]; j < g->out_idx[i+1]; j++){
			if ( label == 0 ){
				fprintf(stdout, "(%"PRIvertex", %"PRIvertex") ", i, g->out_neighbours[j]);
			}else{
				fprintf(stdout, "(%"PRIvertex", %"PRIvertex") ", g->labels[i], g->labels[g->out_neighbours[j]]);
			}
		}
	}
	fprintf(stdout, "\n");
	fprintf(stdout, TXT_GREEN"\tConnections: \n\t"TXT_NO_COLOR);
	for (vertex_type i = 1; i < g->nodes+1; i++) {
		if ( label == 0 ){
			fprintf(stdout, "  %"PRIvertex"->{", i);
		}else{
			fprintf(stdout, "  %"PRIvertex"->{", g->labels[i]);
		}
		for(vertex_type j = g->out_idx[i]; j < g->out_idx[i+1]; j++){
			if ( label == 0 ){
				fprintf(stdout, "%"PRIvertex, g->out_neighbours[j]);
			}else{
				fprintf(stdout, "%"PRIvertex, g->labels[g->out_neighbours[j]]);
			}
			if ( (j+1) < g->out_idx[i+1] ) {
				fprintf(stdout, ",");
			}
		}
		fprintf(stdout, "}");
	}
	fprintf(stdout, "\n");
	
	fprintf(stdout, TXT_GREEN"\tLabels: \n\t"TXT_NO_COLOR);
	for (vertex_type i = 1; i < g->nodes+1; i++) {
		fprintf(stdout, "  %"PRIvertex"=>%"PRIvertex, i, g->labels[i]);
	}
	fprintf(stdout, "\n");
}


/*
 * Function: write_graph_edgelist
 * ------------------------------------------------
 * The function writes the graph <g> on file <file_name> in edge list format.
 *
 * @param g     	    The pointer to the graph to print.   
 * @param file_name		File used to write the graph in edge list format.
 */
void write_graph_edgelist(graph_type *g, const char* file_name){
	FILE *fp = Fopen(file_name, "w", __FILE__, __LINE__, "Opening file to write.");
	
	fprintf(fp, "# Nodes: %"PRIvertex" Edges: %"PRIvertex, g->nodes, g->edges);
	for (vertex_type i = 1; i < g->nodes+1; i++) {
		for(vertex_type j = g->idx[i]; j < g->idx[i+1]; j++){
			if (g->neighbours[j] > i){
				fprintf(fp, "\n%"PRIvertex"  %"PRIvertex, i, g->neighbours[j]);
			}
		}
	}
	fclose(fp);
}


/*
 * Function: write_digraph_edgelist
 * ------------------------------------------------
 * The function writes the digraph <g> on file <file_name> in edge list format.
 *
 * @param g     	    The pointer to the digraph to print.   
 * @param file_name		File used to write the graph in edge list format.
 */
void write_digraph_edgelist(digraph_type *g, const char* file_name){
	FILE *fp = Fopen(file_name, "w", __FILE__, __LINE__, "Opening file to write.");
	
	fprintf(fp, "# Nodes: %"PRIvertex" Edges: %"PRIvertex, g->nodes, g->edges);
	for (vertex_type i = 1; i < g->nodes+1; i++) {
		for(vertex_type j = g->out_idx[i]; j < g->out_idx[i+1]; j++){
			fprintf(fp, "\n%"PRIvertex"  %"PRIvertex, i, g->out_neighbours[j]);
		}
	}
	fclose(fp);
}

/*
 * Function: write_digraph_dot
 * ------------------------------------------------
 * The function writes the digraph <g> on file <file_name> in dot format.
 *
 * @param g     	    The pointer to the digraph to print.   
 * @param file_name		File used to write the graph in dot format.
 */
void write_digraph_dot(digraph_type *g, const char* file_name){
	FILE *fp = Fopen(file_name, "w", __FILE__, __LINE__, "Opening file to write.");
	
	fprintf(fp, "digraph G {\n");
	for (vertex_type i = 1; i < g->nodes+1; i++) {
		for(vertex_type j = g->out_idx[i]; j < g->out_idx[i+1]; j++){
			fprintf(fp, "\t%"PRIvertex" -> %"PRIvertex" ;\n", i, g->out_neighbours[j]);
		}
	}
	fprintf(fp, "}\n");
	fclose(fp);
}


/*
 * Function: write_graph_dot
 * ------------------------------------------------
 * The function writes the graph <g> on file <file_name> in dot format.
 *
 * @param g     	    The pointer to the graph to print.   
 * @param file_name		File used to write the graph in dot format.
 */
void write_graph_dot(graph_type *g, const char* file_name){
	FILE *fp = Fopen(file_name, "w", __FILE__, __LINE__, "Opening file to write.");
	
	fprintf(fp, "graph G {\n");
	for (vertex_type i = 1; i < g->nodes+1; i++) {
		for(vertex_type j = g->idx[i]; j < g->idx[i+1]; j++){
			if (g->neighbours[j] > i){
				fprintf(fp, "\t%"PRIvertex" -- %"PRIvertex" ;\n", i, g->neighbours[j]);
			}
		}
	}
	fprintf(fp, "}\n");
	fclose(fp);
}


/*
 * Function: free_digraph
 * ------------------------------------------------
 * The function frees the memory allocated for a digraph. 
 *
 * @param g     	The pointer to the digraph to free.       
 */
void free_digraph(digraph_type *g){
	if (g != NULL){
		free(g->labels);
		free(g->in_idx);
		free(g->out_idx);
		free(g->in_neighbours);
		free(g->out_neighbours);
		free(g->del_nodes_list);
		free(g);
	}
	return;
}

/*
 * Function: free_graph
 * ------------------------------------------------
 * The function frees the memory allocated for a graph. 
 *
 * @param g     	The pointer to the graph to free.       
 */
void free_graph(graph_type *g){
	if (g != NULL){
		free(g->labels);
		free(g->idx);
		free(g->neighbours);
		free(g->del_nodes_list);
		free(g);
	}
	return;
}


/*
 * Function: free_forest
 * ------------------------------------------------
 * The function frees the memory allocated for the Forest f.
 *
 * @param f		The Forest to free.       
 */
void free_forest(forest_type *f){
	if (f != NULL){
		free(f->vertex);
		free(f->parent);
		free(f->order);
		free(f);
	}
	return;
}


/*
 * Function: print_forest
 * ------------------------------------------------
 * The function prints on the standard output a Forest <f>
 *
 * @param f		The forest to print   
 */
void print_forest(forest_type *f){
	if (f == NULL){ return; }
	fprintf(stdout, TXT_GREEN"TREE:"TXT_YELLOW"\n  Nodes: "TXT_NO_COLOR"%"PRIvertex" "TXT_YELLOW" Edges: "TXT_NO_COLOR"%"PRIvertex"\n", f->nodes, f->edges);
	for (vertex_type i = 1, j = 0; i < f->all_nodes+1; i++) {
		vertex_type w;
		w = f->parent[i];
		if( w != 0 && w != i){ /* We visited the node and i is not a root of a tree */
			if ( j != 0 ){
				fprintf(stdout, ",");
			}
			j++;
			fprintf(stdout, " %"PRIvertex"->%"PRIvertex, w, i);
		}
	}
	fprintf(stdout, "\n");
}

/*
 * Function: write_forest_dot
 * ------------------------------------------------
 * The function writes on a file the Forest <f> in dot format.
 *
 * @param f				The Forest to write   
 * @param file_name		The name of the file.
 * @param type			The type of graph to write, possible values are: DIGRAPH or GRAPH.
 */
void write_forest_dot(forest_type *f, const char* file_name, char type){
	FILE *fp = Fopen(file_name, "w", __FILE__, __LINE__, "Opening file to write.");
	
	if ( type == DIGRAPH ){
		fprintf(fp, "digraph G {\n");
		for (vertex_type i = 1; i < f->all_nodes+1; i++) {
			vertex_type v;
			v = f->parent[i];
			if( v != 0 && v != i){ /* We visited the node and i is not a root of a tree */
				fprintf(fp, "\t%"PRIvertex" -> %"PRIvertex" ;\n", v, i);
			}
		}
		fprintf(fp, "}\n");
	}else if ( type == GRAPH ){
		fprintf(fp, "graph G {\n");
		for (vertex_type i = 1; i < f->all_nodes+1; i++) {
			vertex_type v;
			v = f->parent[i];
			if( v != 0 && v != i){ /* We visited the node and i is not a root of a tree */
				fprintf(fp, "\t%"PRIvertex" -- %"PRIvertex" ;\n", v, i);
			}
		}
		fprintf(fp, "}\n");
	}
	fclose(fp);
}

/*
 * Function: write_forest_edge_list
 * ------------------------------------------------
 * The function writes on a file the Forest <f> in edges list format.
 *
 * @param f				The Forest to write.
 * @param file_name		The name of the file.
 */
void write_forest_edge_list(forest_type *f, const char* file_name){
	FILE *fp = Fopen(file_name, "w", __FILE__, __LINE__, "Opening file to write.");
	
	fprintf(fp, "# Nodes: %"PRIvertex" Edges: %"PRIvertex, f->nodes, f->edges);
	for (vertex_type i = 1; i < f->all_nodes+1; i++) {
		vertex_type v;
		v = f->parent[i];
		if( v != 0 && v != i){ /* We visited the node and i is not a root of a tree */
			fprintf(fp, "\n%"PRIvertex" %"PRIvertex, v, i);
		}
	}
	fclose(fp);
}


/*
 * Function: get_graph_dfs
 * ------------------------------------------------
 * The function does a Depth First Search (DFS) of a graph <g>. If <unconn> is 0
 * the search stop when no other nodes can be reached from the <root>. If <unconn>
 * is 1 the search keep visiting nodes until all nodes have been visited.
 *
 * @param g     	 The pointer to the graph.
 * @param root		 The vertex to use as the root for the search.
 * @param unconn	 To decide how to proceed if the graph is disconnected.
 *    
 * @return           The function returns a forest representing the DFS visit of the graph.  
 */
forest_type *get_graph_dfs(graph_type *g, vertex_type root, char unconn){
	char *visited;						/* List of visited nodes */
	vertex_type *stack;					/* Stack structure used for the visit of the graph */
	vertex_type stack_idx, dfs_idx;
	vertex_type n_nodes, n_edges;
	vertex_type last_unvisited;			/* The node with the lowest index we haven't yet visited */
	char proceed;
	vertex_type v, w;
	forest_type *dfs_forest;
	
	n_nodes = g->nodes;
	dfs_forest = (forest_type *) Malloc(sizeof(forest_type),__FILE__, __LINE__, "dfs_tree.");
	dfs_forest->vertex = (vertex_type *) Malloc(sizeof(vertex_type)*(n_nodes+1),__FILE__, __LINE__, "dfs_tree->vertex.");
	dfs_forest->parent = (vertex_type *) Calloc(n_nodes+1,sizeof(vertex_type),__FILE__, __LINE__, "dfs_tree->parent.");
	dfs_forest->order = (vertex_type *) Calloc(n_nodes+1,sizeof(vertex_type),__FILE__, __LINE__, "dfs_tree->order.");
	dfs_forest->all_nodes = n_nodes;
	visited = (char *) Calloc(n_nodes+1, sizeof(vertex_type),__FILE__, __LINE__, "Visited array.");
	stack = (vertex_type *) Malloc(sizeof(vertex_type)*(g->edges*2),__FILE__, __LINE__, "Stack array.");
	stack_idx = 1;
	dfs_idx = 1;
	n_edges = 0;
	stack[stack_idx] = root;
	proceed = 1;
	last_unvisited = root;
	dfs_forest->parent[root] = root;
	
	while (proceed == 1){ 
		v = stack[stack_idx];
		stack_idx--;
		if ( visited[v] == 0){ /* New node, not yet visited */
			visited[v] = 1;
			dfs_forest->order[v] = dfs_idx;
			dfs_forest->vertex[dfs_idx] = v;
			dfs_idx++;
			/* Put the neighbours of the node in the stack */
			for(vertex_type j = g->idx[v+1]-1; j >= g->idx[v]; j--){
				w = g->neighbours[j];
				if( visited[w] == 0 ){
					stack_idx++;
					stack[stack_idx] = w;
					if ( dfs_forest->parent[w] == 0) { n_edges++;}
					dfs_forest->parent[w] = v;
				}
			}
		}
		if( stack_idx == 0){ /* This was the last node in the stack */
			if (unconn == 0){ break;}
			if (dfs_idx == n_nodes+1){ /* check if we visited all nodes */
				proceed = 0;
			}else{
				for(vertex_type i = last_unvisited; i < n_nodes+1; i++ ){ 
					if(visited[i] == 0){
						stack_idx++;
						stack[stack_idx] = i;
						last_unvisited = i+1;
						dfs_forest->parent[i] = i;
						break;
					}
				}
			}
		}
	}
	free(visited);
	free(stack);
	dfs_forest->nodes = dfs_idx-1;
	dfs_forest->edges = n_edges;
	return dfs_forest;
}

/*
 * Function: get_digraph_dfs
 * ------------------------------------------------
 * The function does a Depth First Search (DFS) of a digraph <g>. If <unconn> is 0
 * the search stop when no other nodes can be reached from the <root>. If <unconn>
 * is 1 the search keep visiting nodes until all nodes have been visited.
 *
 * @param g     	 The pointer to the graph.
 * @param root		 The vertex to use as the root for the search.
 * @param unconn	 To decide how to proceed if the graph is disconnected.
 *    
 * @return           The function returns a forest representing the DFS visit of the digraph.  
 */
forest_type *get_digraph_dfs(digraph_type *g, vertex_type root, char unconn){
	char *visited;				/* List of visited nodes */
	vertex_type *stack;					/* Stack structure used for the visit of the graph */
	vertex_type stack_idx, dfs_idx;
	vertex_type n_nodes, n_edges;
	vertex_type last_unvisited;			/* The node with the lowest index we haven't yet visited */
	char proceed;
	vertex_type v, w;
	forest_type *dfs_forest;
	
	n_nodes = g->nodes;
	dfs_forest = (forest_type *) Malloc(sizeof(forest_type),__FILE__, __LINE__, "dfs_tree.");
	dfs_forest->vertex = (vertex_type *) Malloc(sizeof(vertex_type)*(n_nodes+1),__FILE__, __LINE__, "dfs_tree->vertex.");
	dfs_forest->parent = (vertex_type *) Calloc(n_nodes+1,sizeof(vertex_type),__FILE__, __LINE__, "dfs_tree->parent.");
	dfs_forest->order = (vertex_type *) Calloc(n_nodes+1,sizeof(vertex_type),__FILE__, __LINE__, "dfs_tree->order.");
	dfs_forest->all_nodes = n_nodes;
	visited = (char *) Calloc(n_nodes+1, sizeof(vertex_type),__FILE__, __LINE__, "Visited array.");
	stack = (vertex_type *) Malloc(sizeof(vertex_type)*(g->edges*2),__FILE__, __LINE__, "Stack array.");
	stack_idx = 1;
	dfs_idx = 1;
	n_edges = 0;
	stack[stack_idx] = root;
	proceed = 1;
	last_unvisited = root;
	dfs_forest->parent[root] = root;
	
	while (proceed == 1){ 
		v = stack[stack_idx];
		stack_idx--;
		if ( visited[v] == 0){ /* New node, not yet visited */
			visited[v] = 1;
			dfs_forest->order[v] = dfs_idx;
			dfs_forest->vertex[dfs_idx] = v;
			dfs_idx++;
			/* Put the neighbours of the node in the stack*/
			for(vertex_type j = g->out_idx[v+1]-1; j >= g->out_idx[v]; j--){
				w = g->out_neighbours[j];
				if( visited[w] == 0 ){
					stack_idx++;
					stack[stack_idx] = w;
					if ( dfs_forest->parent[w] == 0) { n_edges++;}
					dfs_forest->parent[w] = v;
				}
			}
		}
		if( stack_idx == 0){ /* This was the last node in the stack */
			if (unconn == 0){ break;}
			if (dfs_idx == n_nodes+1){ /* check if we visited all nodes */
				proceed = 0;
			}else{
				for(vertex_type i = last_unvisited; i < n_nodes+1; i++ ){ 
					if(visited[i] == 0){
						stack_idx++;
						stack[stack_idx] = i;
						last_unvisited = i+1;
						dfs_forest->parent[i] = i;
						break;
					}
				}
			}
		}
	}
	free(visited);
	free(stack);
	dfs_forest->nodes = dfs_idx-1;
	dfs_forest->edges = n_edges;
	return dfs_forest;
}

/*
 * Function: get_digraph_bfs
 * ------------------------------------------------
 * The function does a Breadth First Search (BFS) of a digraph <g>. If <unconn> is 0
 * the search stop when no other nodes can be reached from the <root>. If <unconn>
 * is 1 the search keep visiting nodes until all nodes have been visited.
 *
 * @param g     	 The pointer to the graph.
 * @param root		 The vertex to use as root for the search.
 * @param unconn	 To decide how to proceed if the graph is disconnected.
 *    
 * @return           The function returns a forest representing the BFS visit of the digraph. .   
 */
forest_type *get_digraph_bfs(digraph_type *g, vertex_type root, char unconn){
	vertex_type *queue;					/* Queue structure used to visit the graph */
	vertex_type head, tail, bfs_idx;
	vertex_type n_nodes, n_edges;
	vertex_type last_unvisited;			/* The node with the lowest index we haven't yet visited */
	char proceed, *visited_in_queue;
	vertex_type v, w;
	forest_type *bfs_forest;
	
	n_nodes = g->nodes;
	bfs_forest = (forest_type *) Malloc(sizeof(forest_type),__FILE__, __LINE__, "bfs_tree.");
	bfs_forest->vertex = (vertex_type *) Malloc(sizeof(vertex_type)*(n_nodes+1),__FILE__, __LINE__, "bfs_tree->vertex.");
	bfs_forest->parent = (vertex_type *) Calloc(n_nodes+1,sizeof(vertex_type),__FILE__, __LINE__, "bfs_tree->parent.");
	bfs_forest->order = (vertex_type *) Calloc(n_nodes+1,sizeof(vertex_type),__FILE__, __LINE__, "bfs_tree->order.");
	bfs_forest->all_nodes = n_nodes;
	visited_in_queue = (char *) Calloc(n_nodes+1, sizeof(vertex_type),__FILE__, __LINE__, "Visited array.");
	queue = (vertex_type *) Malloc(sizeof(vertex_type)*n_nodes,__FILE__, __LINE__, "Queue array.");
	head = 0;
	tail = 1;
	bfs_idx = 1;
	n_edges = 0;
	queue[head] = root;
	last_unvisited = root;
	proceed = 1;
	visited_in_queue[root]=1;
	bfs_forest->parent[root] = root;
	
	while( proceed == 1){
		v = queue[head];
		head++;
		if (visited_in_queue[v] == 1 ){ /* New node to visit */
			visited_in_queue[v] = 2;
			bfs_forest->order[v]=bfs_idx;
			bfs_forest->vertex[bfs_idx]=v;
			bfs_idx++;
			for(vertex_type j = g->out_idx[v]; j <= g->out_idx[v+1]-1; j++){
				w = g->out_neighbours[j];
				if( visited_in_queue[w] == 0 ){
					queue[tail] = w;
					tail++;
					n_edges++;
					bfs_forest->parent[w] = v;
					visited_in_queue[w] = 1;
				}
			}
		}
		
		if (head == tail){
			if (unconn == 0) {break;}
			if (bfs_idx == n_nodes+1){ /* check if we visited all nodes */
				proceed = 0;
			}else{
				for(vertex_type i = last_unvisited; i < n_nodes+1; i++ ){ 
					if(visited_in_queue[i] == 0){
						queue[tail] = i;
						tail++;
						last_unvisited = i+1;
						visited_in_queue[i] = 1;
						bfs_forest->parent[i] = i;
						break;
					}
				}
			}	
		}
	}
	free(visited_in_queue);
	free(queue);
	bfs_forest->nodes = bfs_idx-1;
	bfs_forest->edges = n_edges;
	return bfs_forest;
}

/*
 * Function: get_graph_bfs
 * ------------------------------------------------
 * The function does a Breadth First Search (BFS) of a graph <g>. If <unconn> is 0
 * the search stop when no other nodes can be reached from the <root>. If <unconn>
 * is 1 the search keep visiting nodes until all nodes have been visited.
 *
 * @param g     	 The pointer to the graph.
 * @param root		 The vertex to use as root for the search.
 * @param unconn	 To decide how to proceed if the graph is disconnected.
 *    
 * @return           The function returns a forest representing the BFS visit of the graph. .   
 */
forest_type *get_graph_bfs(graph_type *g, vertex_type root, char unconn){
	vertex_type *queue;					/* Queue structure used to visit the graph */
	vertex_type head, tail, bfs_idx;
	vertex_type n_nodes, n_edges;
	vertex_type last_unvisited;			/* The node with the lowest index we haven't yet visited */
	char proceed, *visited_in_queue;
	vertex_type v, w;
	forest_type *bfs_forest;
	
	n_nodes = g->nodes;
	bfs_forest = (forest_type *) Malloc(sizeof(forest_type),__FILE__, __LINE__, "bfs_tree.");
	bfs_forest->vertex = (vertex_type *) Malloc(sizeof(vertex_type)*(n_nodes+1),__FILE__, __LINE__, "bfs_tree->vertex.");
	bfs_forest->parent = (vertex_type *) Calloc(n_nodes+1,sizeof(vertex_type),__FILE__, __LINE__, "bfs_tree->parent.");
	bfs_forest->order = (vertex_type *) Calloc(n_nodes+1,sizeof(vertex_type),__FILE__, __LINE__, "bfs_tree->order.");
	bfs_forest->all_nodes = n_nodes;
	visited_in_queue = (char *) Calloc(n_nodes+1, sizeof(vertex_type),__FILE__, __LINE__, "Visited array.");
	queue = (vertex_type *) Malloc(sizeof(vertex_type)*n_nodes,__FILE__, __LINE__, "Queue array.");
	head = 0;
	tail = 1;
	bfs_idx = 1;
	n_edges = 0;
	queue[head] = root;
	last_unvisited = root;
	proceed = 1;
	visited_in_queue[root]=1;
	bfs_forest->parent[root] = root;
	
	while( proceed == 1){
		v = queue[head];
		head++;
		if (visited_in_queue[v] == 1 ){ /* New node to visit */
			visited_in_queue[v] = 2;
			bfs_forest->order[v]=bfs_idx;
			bfs_forest->vertex[bfs_idx]=v;
			bfs_idx++;
			for(vertex_type j = g->idx[v]; j <= g->idx[v+1]-1; j++){
				w = g->neighbours[j];
				if( visited_in_queue[w] == 0 ){
					queue[tail] = w;
					tail++;
					n_edges++;
					bfs_forest->parent[w] = v;
					visited_in_queue[w] = 1;
				}
			}
		}
		
		if (head == tail){
			if (unconn == 0) {break;}
			if (bfs_idx == n_nodes+1){ /* check if we visited all nodes */
				proceed = 0;
			}else{
				for(vertex_type i = last_unvisited; i < n_nodes+1; i++ ){ 
					if(visited_in_queue[i] == 0){
						queue[tail] = i;
						tail++;
						last_unvisited = i+1;
						visited_in_queue[i]=1;
						bfs_forest->parent[i] = i;
						break;
					}
				}
			}	
		}
	}
	free(visited_in_queue);
	free(queue);
	bfs_forest->nodes = bfs_idx-1;
	bfs_forest->edges = n_edges;
	return bfs_forest;
}

/*
 * Function: get_digraph_scc
 * ------------------------------------------------
 * The function computes the Strongly Connected Components (SCC) of a digraph g. 
 * The function implements the iterative version of the Tarjan algorithm, presented in
 * Tarjan, R. E. (1972), "Depth-first search and linear graph algorithms". 
 *
 * @param g     	The pointer to the digraph. 
 *
 * @return			The function returns the pointer to a cc data stracture containing the SCC
 *					composition of the digraph g. 
 */
cc_type *get_digraph_scc(digraph_type *g){
	vertex_type *lowlink, *number;	
	vertex_type *scc_stack, *stack;	
	char *node_in_scc_stack;
	cc_type *scc;	
	vertex_type stack_idx, last_unvisited, dfs_idx, scc_stack_idx, n_nodes;
	char proceed;
	
	n_nodes = g->nodes; 
	lowlink = (vertex_type *) Malloc(sizeof(vertex_type)*(n_nodes+1), __FILE__, __LINE__, "Lowlink array.");
	number = (vertex_type *) Calloc( n_nodes+1, sizeof(vertex_type), __FILE__, __LINE__, "Number array.");
	stack = (vertex_type *) Malloc(sizeof(vertex_type)*(g->edges*2),__FILE__, __LINE__, "Stack array.");
	scc_stack = (vertex_type *) Malloc(sizeof(vertex_type)*(n_nodes+1),__FILE__, __LINE__, "SCC Stack array.");
	node_in_scc_stack = (char *) Calloc((n_nodes+1), sizeof(char), __FILE__, __LINE__, "Check SCC Stack array.");
	scc = (cc_type*) Malloc(sizeof(cc_type), __FILE__, __LINE__, "SCC." );
	scc->vertex = (vertex_type *) Malloc(sizeof(vertex_type)*(n_nodes+1), __FILE__, __LINE__, "scc->vertex array.");
	scc->idx = (vertex_type *) Malloc(sizeof(vertex_type)*(n_nodes+1), __FILE__, __LINE__, "scc->idx array.");
	
	proceed = 1;
	last_unvisited = 1;
	stack_idx = 1;
	dfs_idx = 1;
	scc_stack_idx = 0;
	stack[stack_idx] = 1;
	scc->number = 0;
	scc->idx[0] = 0;
	
	while(proceed == 1){
		vertex_type v;
		vertex_type w;
		v = stack[stack_idx];
		if(number[v] == 0){ /* _ New node to visit, we keep it in the stack _ */
			scc_stack_idx++;	/* put the node in the stack of points */
			scc_stack[scc_stack_idx]=v;
			node_in_scc_stack[v]=1;
			lowlink[v]=dfs_idx;	/* set lowlink and number of the node */
			number[v]=dfs_idx;
			dfs_idx++;
			for(vertex_type j = g->out_idx[v+1]-1; j >= g->out_idx[v]; j--){ /* For each neighbour of the node */
				w = g->out_neighbours[j];
				if( number[w] == 0 ){ /* Put in the stack the nodes we didn't visit yet */
					stack_idx++;
					stack[stack_idx] = w;	
				}else if ( number[w] < number[v] && node_in_scc_stack[w] == 1){  /* Update the lowlink of v if the neighbour has been visited and it is in the stack of points  */
					lowlink[v] = MIN( lowlink[v], number[w]);
				}
			}
		}else if(node_in_scc_stack[v] == 1) { /* _ A node we visited and it is still in the stack of points _ */
			stack_idx--;	/* Remove the node from the stack */
			for(vertex_type j = g->out_idx[v+1]-1; j >= g->out_idx[v]; j--){ /* For each neighbour of the node */
				w = g->out_neighbours[j];
				if ( node_in_scc_stack[w] == 1 ){ /* Update the lowlink of v if the neighbour w has been visited and w is in the stack of points  */
					lowlink[v] = MIN( lowlink[v], lowlink[w]);
				}
			}
			if(lowlink[v] == number[v] ){ /* v is the root of a tree, v is the root of a SCC */
				vertex_type scc_size = 0;
				vertex_type scc_idx = scc->idx[scc->number];
				scc->number++;
				w = scc_stack[scc_stack_idx];
				while( number[w] >= number[v] ){
					scc_size++;
					node_in_scc_stack[w] = 0;
					scc->vertex[scc_idx] = w; 	/* update scc data structure */
					scc_idx++;
					if ( scc_stack_idx == 1 ) { break; }
					scc_stack_idx--; 	/* read next vertex that is possibly in the scc */
					w = scc_stack[scc_stack_idx];
				}	
				scc->idx[scc->number] = scc_idx;
			}
		}else{ /* ___ A node we visited and it is NO MORE in the stack of points (scc_stack) ___ */
			stack_idx--;	/* Remove the node from the stack */
		}		
		if( stack_idx == 0 ){ /* This was the last node in the stack */
			if ( dfs_idx == n_nodes+1 ){ /* We visited all nodes */
				proceed = 0; 
			}else{
				for(vertex_type i = last_unvisited; i < n_nodes+1; i++ ){  /* There are other nodes to visit */
					if(number[i] == 0){
						stack_idx++;
						stack[stack_idx] = i;
						last_unvisited = i+1;
						break;
					}
				}
			}
		}	
	}
	scc->idx = (vertex_type*) Realloc(scc->idx, sizeof(vertex_type)*(scc->number+1), __FILE__, __LINE__, "scc->idx array.");
    free(stack);
    free(scc_stack);
    free(lowlink);
    free(number);
    free(node_in_scc_stack);
	return scc;	
}


/*
 * Function: print_cc
 * ------------------------------------------------
 * The function prints the CC decomposition of a graph. 
 *
 * @param cc     	 The pointer to the CC to print.  
 */
void print_cc (cc_type *cc){
	for(vertex_type i = 0; i < cc->number; i++){
		fprintf(stdout, TXT_GREEN"CC %"PRIvertex":"TXT_NO_COLOR, i);
		for(vertex_type j = cc->idx[i]; j < cc->idx[i+1]; j++){
			fprintf(stdout, " %"PRIvertex, cc->vertex[j]);
		}
		fprintf(stdout,  " - size: "TXT_YELLOW"%"PRIvertex""TXT_NO_COLOR"\n", cc->idx[i+1]-cc->idx[i]);
	}
}

/*
 * Function: write_cc
 * ------------------------------------------------
 * The function writes the CC decomposition of a graph. 
 *
 * @param cc     	 The pointer to the CC to write.
 * @param file_name  The name of file to write the output.  
 */
void write_cc (cc_type *cc, const char *file_name){
	FILE *fp = Fopen(file_name, "w", __FILE__, __LINE__, "Opening file to write.");

	fprintf(fp, "{ \n");
	for(vertex_type i = 0; i < cc->number; i++){
		fprintf(fp, "\t\"%"PRIvertex"\": { \n",i);
		fprintf(fp, "\t\t\"nodes\": [");
		for(vertex_type j = cc->idx[i]; j < cc->idx[i+1]-1; j++){
			fprintf(fp, " %"PRIvertex",", cc->vertex[j]);
		}
		fprintf(fp, " %"PRIvertex"],\n", cc->vertex[ cc->idx[i+1]-1  ]);
		fprintf(fp, "\t\t\"size\" : %"PRIvertex"\n", cc->idx[i+1]-cc->idx[i]);
		if (i != cc->number-1){
			fprintf(fp, "\t},\n");
		}else{
			fprintf(fp, "\t}\n");
		}
	}
	fprintf(fp, "}");
	fclose(fp);
	return;
}


/*
 * Function: free_cc
 * ------------------------------------------------
 * The function frees the memory allocated for a CC. 
 *
 * @param cc     	The pointer to the CC to free.       
 */
void free_cc (cc_type *cc){
	if(cc != NULL){
		free(cc->vertex);
		free(cc->idx);
		free(cc);
	}
}


/*
 * Function: get_inoutdegree
 * ------------------------------------------------
 * The function computes the in-degree, out-degree of each node in the digraph <g>.
 *
 * @param g     	The digraph.   
 *
 * @return			The function returns a pointer to a didegree_type array containing the degrees of 
 *					the nodes. A NULL pointer is returned if an error occurs.
 */
didegree_type *get_inoutdegree(digraph_type *g){
	didegree_type *nodes_degree;
	
	if (g == NULL){ return NULL; }
	nodes_degree = (didegree_type *) Malloc( g->nodes*sizeof(didegree_type), __FILE__, __LINE__, "malloc nodes_degree.");

	for (vertex_type i = 1; i < g->nodes+1; i++) {
		nodes_degree[i-1].v = i;
		nodes_degree[i-1].in_deg = g->in_idx[i+1]-g->in_idx[i];
		nodes_degree[i-1].out_deg = g->out_idx[i+1]-g->out_idx[i];
		nodes_degree[i-1].all_deg = nodes_degree[i-1].in_deg+nodes_degree[i-1].out_deg;
	}
	return nodes_degree;
}


/*
 * Function: write_inoutdegree
 * ------------------------------------------------
 * The function prints the degree of the nodes in <d>.
 *
 * @param d     		The degree of the nodes. 
 * @param n_nodes  		Number of graph nodes.
 * @param top_n  		Number of node to print.
 * @param file			The name of the file used to write in/out degree information
 */
void write_inoutdegree(didegree_type *d, vertex_type n_nodes, vertex_type top_n, const char *file_name){
	FILE *file = Fopen(file_name, "w", __FILE__, __LINE__, "Opening file to write.");
	
	if (d == NULL){ return; }
	vertex_type max_nodes = MIN(n_nodes, top_n);
	fprintf(file, "Node, InDeg, OutDeg, AllDeg\n");
	for(vertex_type i=0; i < max_nodes; i++){
		fprintf(file, "%"PRIvertex", %"PRIvertex", %"PRIvertex", %"PRIvertex "\n", d[i].v, d[i].in_deg, d[i].out_deg,d[i].all_deg);
	}
	fclose(file);
}

/*
 * Function: print_inoutdegree
 * ------------------------------------------------
 * The function prints the degree of the nodes in <d>.
 *
 * @param d     		The degree of the nodes. 
 * @param n_nodes  		Number of graph nodes.
 * @param top_n  		Number of node to print.
 */
void print_inoutdegree(didegree_type *d, vertex_type n_nodes, vertex_type top_n){
	if (d == NULL){ return; }
	vertex_type max_nodes = MIN(n_nodes, top_n);
	printf("Node, InDeg, OutDeg, AllDeg\n");
	for(vertex_type i=0; i < max_nodes; i++){
		printf("%"PRIvertex", %"PRIvertex", %"PRIvertex", %"PRIvertex"\n", d[i].v, d[i].in_deg, d[i].out_deg,d[i].all_deg);
	}
}


static int outdeg_asc (const void * a, const void * b) {
   return ( (*(didegree_type*)a).out_deg - (*(didegree_type*)b).out_deg );
}
static int outdeg_desc (const void * a, const void * b) {
   return ( (*(didegree_type*)b).out_deg - (*(didegree_type*)a).out_deg );
}
static int indeg_asc (const void * a, const void * b) {
   return ( (*(didegree_type*)a).in_deg - (*(didegree_type*)b).in_deg );
}
static int indeg_desc (const void * a, const void * b) {
   return ( (*(didegree_type*)b).in_deg - (*(didegree_type*)a).in_deg );
}
static int alldeg_asc (const void * a, const void * b) {
   return ( (*(didegree_type*)a).all_deg - (*(didegree_type*)b).all_deg );
}
static int alldeg_desc (const void * a, const void * b) {
   return ( (*(didegree_type*)b).all_deg - (*(didegree_type*)a).all_deg );
}
/*
 * Function: sort_inoutdegree
 * ------------------------------------------------
 * The function sorts the degree values in <d>.
 *
 * @param d			The array of didgree_type to sort.  
 * @param n_nodes   Number of nodes in the graph.
 * @param field		The field of the array to sort, possible values are: INDEG or OUTDEG.
 * @param order		The order of the sort, possible values are: ASC or DESC.
 */
void sort_inoutdegree(didegree_type *d, vertex_type n_nodes, char field, char order){
	switch (field){
		case INDEG:
			if (order == ASC){
				qsort(d, n_nodes, sizeof(didegree_type), indeg_asc);
			}else if (order == DESC){
				qsort(d, n_nodes, sizeof(didegree_type), indeg_desc);
			}
			break;
		case OUTDEG:
			if (order == ASC){
				qsort(d, n_nodes, sizeof(didegree_type), outdeg_asc);
			}else if (order == DESC){
				qsort(d, n_nodes, sizeof(didegree_type), outdeg_desc);
			}
			break;
		case ALLDEG:
			if (order == ASC){
				qsort(d, n_nodes, sizeof(didegree_type), alldeg_asc);
			}else if (order == DESC){
				qsort(d, n_nodes, sizeof(didegree_type), alldeg_desc);
			}
			break;
		default:
			break;
	}	
}



/*
 * Function: get_pairs_connectivity
 * -----------------------------------
 * Computes the connectivity value of a graph G.
 * Let G a directed or undirected graph, let C_1,C_2,...,C_l be its 
 * (weakly or strongly) connected components.
 * We define the connectivity value of G as:
 *
 *      f(G) = \sum_{i=1}^{l} \binom{C_i}{2}
 *
 * f(G) equals the number of vertex pairs in G that are (weakly or strongly) 
 * connected (i.e., pairwise connectivity value)
 *
 * @param scc       the scc decomposition of the digraph
 *
 * @return          the connectivity value of the digraph
 */
uint64_t get_pairs_connectivity(cc_type *cc) {
    uint64_t conn_val;		  /* the connectivity value */

    conn_val = 0;
    for ( vertex_type i = 1; i <= cc->number; i++) {
		uint64_t c_size = 0;
		c_size = cc->idx[i] - cc->idx[i - 1];
        conn_val += c_size*(c_size-1)/2 ;
    }
    return conn_val;
}




/*
 * Function: delete_nodes_from_graph
 * --------------------------------------
 * The function removes a set of nodes from a graph.
 *
 * @param g			the graph 
 * @param vts		the set of nodes to remove
 * @param n_vts		the number of nodes to remove
 *
 * @return          a copy of the original graph without the nodes to be removed    	
 */
graph_type *delete_nodes_from_graph(graph_type *g, vertex_type *vts, vertex_type n_vts) {
    graph_type *g_minus;				/* The new graph G\A */
	vertex_type n_nodes, n_edges, n_kept_nodes;
	char *to_remove;
	vertex_type *new_id;

	n_nodes = g->nodes;
	n_edges = g->edges;
	n_kept_nodes = n_nodes - n_vts;
    g_minus = (graph_type *) Calloc( 1, sizeof(graph_type),__FILE__, __LINE__, "Calloc of g_minus");
    g_minus->nodes = n_kept_nodes;
	g_minus->edges = n_edges;
    g_minus->labels = (vertex_type *) Malloc( (n_kept_nodes + 1)*sizeof(vertex_type),__FILE__, __LINE__, "Malloc of g_minus->labels"); 
	g_minus->idx = (vertex_type *) Malloc( (n_kept_nodes + 2)*sizeof(vertex_type), __FILE__, __LINE__,"Malloc of g_minus->idx");  
	g_minus->neighbours = (vertex_type *) Malloc( (n_edges*2 + 1)*sizeof(vertex_type),__FILE__, __LINE__, "Malloc of g_minus->neighbours");
    g_minus->del_nodes_num = g->del_nodes_num + n_vts;
    g_minus->del_nodes_list = (vertex_type *) Malloc( (g_minus->del_nodes_num)*sizeof(vertex_type),__FILE__, __LINE__, "Malloc of g_minus->del_nodes_list");
    g_minus->org_nodes_num = g->org_nodes_num;
	to_remove = (char *) Calloc(n_nodes+1, sizeof(char), __FILE__, __LINE__,"Calloc to_remove array");
	
	if (g->del_nodes_num > 0 ){
		memcpy(g_minus->del_nodes_list, g->del_nodes_list, (g->del_nodes_num)*sizeof(vertex_type));
    	for(vertex_type x = g->del_nodes_num, y=0; x < g_minus->del_nodes_num; x++, y++){
			g_minus->del_nodes_list[x] = g->labels[vts[y]];
    	}
	}else{
		memcpy(g_minus->del_nodes_list, vts, n_vts*sizeof(vertex_type));
	}
	
    for(vertex_type x = 0; x < n_vts; x++){
		if(vts[x] > n_nodes){
			return NULL;
		}
		to_remove[vts[x]] = 1;
    }
	new_id = (vertex_type *) Calloc(n_nodes+1, sizeof(vertex_type),__FILE__, __LINE__, "Calloc to_remove array");
	
	/* Remove edges */
	vertex_type j, h, del_edges, kept_edges;
    for (j = 1, kept_edges = 1, h = 1, del_edges = 0; j < n_nodes+1; j++) { 
		/* Checks all vertices, scans idx and neighbours */
        if ( to_remove[j] == 1) { /* We remove <j> and all edges starting from <j> */
			del_edges += g->idx[j+1] - g->idx[j];
        }else{   /* We do not remove <j>, <h> is the ID of <j> in the new graph. */
            g_minus->labels[h] = g->labels[j];
            g_minus->idx[h] = g->idx[j] - del_edges;
			new_id[j] = h;
            h++;
			/* We do not remove <j>, we checks if we can copy its neighbors in the new graph */
            for (vertex_type i = g->idx[j]; i < g->idx[j + 1]; i++) { /* Scan all neighbors <u> of <j> */
				vertex_type u;
                u = g->neighbours[i];
				if (to_remove[u] == 0){ /* we keep the neighbour */
					g_minus->neighbours[kept_edges] = u;
					kept_edges++;
				}else{ del_edges++; }
            }
        }
    }
    g_minus->idx[h] = g->idx[j] - del_edges; /* Update idx last index [n+1]*/
    for(vertex_type x = 1; x < kept_edges; x++){ /* update nodes id */
		g_minus->neighbours[x] = new_id[g_minus->neighbours[x]];
    }

    g_minus->edges -= del_edges/2;
    g_minus->neighbours = Realloc( g_minus->neighbours, (kept_edges) * sizeof(vertex_type), __FILE__, __LINE__, "Realloc of g_minus->neighbours" );
	free(to_remove);
    free(new_id);

    return g_minus;
}




/*
 * Function: delete_nodes_from_digraph
 * --------------------------------------
 * The function removes a set of nodes from a digraph.
 *
 * @param g			the digraph 
 * @param vts		the set of nodes to remove
 * @param n_vts		the number of nodes to remove
 *
 * @return          a copy of the original digraph without the nodes to be removed    	
 */
digraph_type *delete_nodes_from_digraph(digraph_type *g, vertex_type *vts, vertex_type n_vts) {
    digraph_type *g_minus;				/* The new graph G\A */
	vertex_type n_nodes, n_edges, n_kept_nodes;
	char *to_remove;
	vertex_type *new_id;

	n_nodes = g->nodes;
	n_edges = g->edges;
	n_kept_nodes = n_nodes - n_vts;
    g_minus = (digraph_type *) Calloc( 1, sizeof(digraph_type),__FILE__, __LINE__, "Calloc of g_minus");
    g_minus->nodes = n_kept_nodes;
	g_minus->edges = n_edges;
    g_minus->labels = (vertex_type *) Malloc( (n_kept_nodes + 1)*sizeof(vertex_type),__FILE__, __LINE__, "Malloc of g_minus->labels");
    g_minus->in_idx = (vertex_type *) Malloc( (n_kept_nodes + 2)*sizeof(vertex_type), __FILE__, __LINE__,"Malloc of g_minus->in_idx");
    g_minus->out_idx = (vertex_type *) Malloc( (n_kept_nodes + 2)*sizeof(vertex_type),__FILE__, __LINE__, "Malloc of g_minus->out_idx");
    g_minus->in_neighbours = (vertex_type *) Malloc( (n_edges + 1)*sizeof(vertex_type),__FILE__, __LINE__, "Malloc of g_minus->in_neighbours");
    g_minus->out_neighbours = (vertex_type *) Malloc( (n_edges + 1)*sizeof(vertex_type), __FILE__, __LINE__,"Malloc of g_minus->out_neighbours");
    g_minus->del_nodes_num = g->del_nodes_num + n_vts;
    g_minus->del_nodes_list = (vertex_type *) Malloc( (g_minus->del_nodes_num)*sizeof(vertex_type),__FILE__, __LINE__, "Malloc of g_minus->del_nodes_list");
    g_minus->org_nodes_num = g->org_nodes_num;
	to_remove = (char *) Calloc(n_nodes+1, sizeof(char), __FILE__, __LINE__,"Calloc to_remove array");
	
	if (g->del_nodes_num > 0 ){
		memcpy(g_minus->del_nodes_list, g->del_nodes_list, (g->del_nodes_num)*sizeof(vertex_type));
    	for(vertex_type x = g->del_nodes_num, y=0; x < g_minus->del_nodes_num; x++, y++){
			g_minus->del_nodes_list[x] = g->labels[vts[y]];
    	}
	}else{
		memcpy(g_minus->del_nodes_list, vts, n_vts*sizeof(vertex_type));
	}
	
    for(vertex_type x = 0; x < n_vts; x++){
		if(vts[x] > n_nodes){
			return NULL;
		}
		to_remove[vts[x]] = 1;
    }
	new_id = (vertex_type *) Calloc(n_nodes+1, sizeof(vertex_type),__FILE__, __LINE__, "Calloc to_remove array");
	
	/* Remove out edges */
	vertex_type j, h, del_edges, kept_edges;
    for (j = 1, kept_edges = 1, h = 1, del_edges = 0; j < n_nodes+1; j++) { 
		/* Checks all vertices, scans out_idx and out_neighbours */
        if ( to_remove[j] == 1) { /* We remove <j> and all edges starting from <j> */
			del_edges += g->out_idx[j+1] - g->out_idx[j];
        }else{   /* We do not remove <j>, <h> is the ID of <j> in the new graph. */
            g_minus->labels[h] = g->labels[j];
            g_minus->out_idx[h] = g->out_idx[j] - del_edges;
			new_id[j] = h;
            h++;
			/* We do not remove <j>, we checks if we can copy its neighbors in the new graph */
            for (vertex_type i = g->out_idx[j]; i < g->out_idx[j + 1]; i++) { /* Scan all neighbors <u> of <j> */
				vertex_type u;
                u = g->out_neighbours[i];
				if (to_remove[u] == 0){ /* we keep the neighbour */
					g_minus->out_neighbours[kept_edges] = u;
					kept_edges++;
				}else{ del_edges++; }
            }
        }
    }
    g_minus->out_idx[h] = g->out_idx[j] - del_edges; /* Update out_idx last index [n+1]*/
    for(vertex_type x = 1; x < kept_edges; x++){ /* update nodes id */
		g_minus->out_neighbours[x] = new_id[g_minus->out_neighbours[x]];
    }
	
	/* Remove in edges */
    for (j = 1, kept_edges = 1, h = 1, del_edges = 0; j < n_nodes+1; j++) { /* Checks all vertices */
        if ( to_remove[j] == 1) { /* We remove <j> and all edges starting from <j> */
            del_edges += g->in_idx[j + 1] - g->in_idx[j];
        } else {  /* We do not remove <j>, <h> is the ID of <j> in the new graph. */
            g_minus->in_idx[h] = g->in_idx[j] - del_edges;
            h++;
            for (vertex_type i = g->in_idx[j]; i < g->in_idx[j + 1]; i++) { /* Scan all neighbors <u> of <j> */
				vertex_type u;
                u = g->in_neighbours[i];
				if (to_remove[u] == 0){ /* we keep the neighbour */
					g_minus->in_neighbours[kept_edges] = new_id[u];
					kept_edges++;
				}else{ del_edges++; }
            }
        }
    }
    g_minus->in_idx[h] = g->in_idx[j] - del_edges; /* Update in_idx last index [n+1]*/

    g_minus->edges -= del_edges;
    g_minus->in_neighbours = Realloc( g_minus->in_neighbours, (kept_edges) * sizeof(vertex_type), __FILE__, __LINE__, "Realloc of g_minus->in_neighbours" );
    g_minus->out_neighbours = Realloc( g_minus->out_neighbours,  (kept_edges) * sizeof(vertex_type), __FILE__, __LINE__, "Realloc of g_minus->out_neighbours" );
	free(to_remove);
    free(new_id);
	
    return g_minus;
}



/*
 * Function: get_digraph_critical_nodes_bruteforce
 * --------------------------------------
 * The function checks which set of nodes of size <k> is critical in the graph <g> and returns
 * a list of them. A set of nodes is critical if it minimaizes the connectivity of <g>.
 * Let G a digraph, let C_1, C_2, ..., C_l be its strongly connected components.
 * We define the connectivity value of G as:
 *
 *      f(G) = \sum_{i=1}^{l} \binom{C_i}{2}
 *
 * f(G) equals the number of vertex pairs in G that are strongly connected
 * (i.e., pairwise strong connectivity value).
 * For each set of vertecies A of size <k> in V the function creates G\A and computes f(G\A).
 *
 * @param g             	the digraph
 * @param k					the number of verticies to remove 
 * @param max_set			maximum number of critical set to return
 * @param end_conn			residual connectivity after removing any critical set
 * @param tot_sets			total number of critical sets found
 *
 * @return              	An array of critical sets or NULL if <max_set> is set to 0. If an error occurs
 *                          both <end_conn> and <tot_sets> are set to 0 and a NULL pointer is returned.
 */
vertex_type ** get_digraph_critical_nodes_bruteforce ( digraph_type *g, vertex_type k, uint64_t max_set, uint64_t *end_conn, uint64_t *tot_sets ) {
	vertex_type n;					/* Number of verticies in G 								*/
	vertex_type **opt_set;          /* List of optimal critical sets found 						*/
    uint64_t n_crtl_set;            /* Number of critical sets found 							*/
    uint64_t crtl_set_idx;          /* Index of the next critical set to insert in <opt_set>	*/    

    n_crtl_set = 0;
    crtl_set_idx = 0;
	n = g->nodes;
	if( k > n ){ *tot_sets = 0; *end_conn=0; return NULL; }
	if (max_set != 0){ // number of optimal set to store
    	opt_set = (vertex_type **) Malloc( max_set*sizeof(vertex_type *),__FILE__, __LINE__, "Malloc of opt_set" );
    	for (vertex_type i = 0; i < max_set; i++){
        	opt_set[i] = (vertex_type *) Malloc( k*sizeof(vertex_type),__FILE__, __LINE__, "Malloc of opt_set[i]" );
    	}
	}else{ opt_set = NULL; }

    /* Compute initial values SCC and Connectivity */
    cc_type *scc;                          /* Number of scc found in G */
	uint64_t min_conn;                      /* Min connectivity value for set of size <k> */
    vertex_type *set_gen, l, j, lst_e;
    digraph_type *diff_g;                   /* G after verticies removal */
    uint64_t conn, idx_s;
	
    /* Start Generating the subsets */
	set_gen = Malloc( k*sizeof(vertex_type), __FILE__, __LINE__,"Malloc of set_gen" );
	for ( j = 0; j < k; j++){
		set_gen[j] = j+1;
	}
	l = k-1;
	lst_e = k-1;
	idx_s = 0;		// Try to remove the first subset
    diff_g = delete_nodes_from_digraph( g, set_gen, k );
    scc = get_digraph_scc(diff_g);
    conn = get_pairs_connectivity(scc);
    min_conn = conn;
	if (max_set != 0){
    	memcpy(opt_set[crtl_set_idx], set_gen, k*sizeof(vertex_type));
	}
    n_crtl_set++;
    crtl_set_idx++;
    free_digraph(diff_g);
    free_cc(scc);
	idx_s++;
	
	if ( k != n ){ 
        while( (set_gen[l]+(k-l) <= n || l != 0)){ /* while there are valid elments to increment */
            set_gen[lst_e]++;
            if( set_gen[lst_e] > n ){ /* increment the next incremental element */
                for ( j = lst_e; j > l; j--) { 
                    /* Check if there is an incrementable element between <l> and the last element
                     * of <set_gen> starting from the end of <set_gen>, otherwise we use <l> */
                    if ( set_gen[j]+(k-j) <= n ){
                        l = j;
                        break;
                    }
                }
                if ( set_gen[l]+(k-l) > n ){ /* Check if element <l> can be incremented, otherwise we pass to the next <l--> */
                    l--;
                }
                set_gen[l]++;
                for ( j = l+1; j < k; j++){ /* set the values of the new starting indexes from <l> */
                    set_gen[j] = set_gen[j-1]+1;
                }	
            }
            diff_g = delete_nodes_from_digraph( g, set_gen, k );
            scc = get_digraph_scc(diff_g);
            conn = get_pairs_connectivity(scc);
            if(conn < min_conn){ /* Save the set as the new optimal solution */
                min_conn = conn;
                n_crtl_set = 1;
                crtl_set_idx = 0;
                if (max_set != 0){
                    memcpy(opt_set[crtl_set_idx], set_gen, k*sizeof(vertex_type));
                    crtl_set_idx++;
                }
            }else if ( min_conn == conn ){ /* Add the set to the optimal solution */
                n_crtl_set++;
                if (crtl_set_idx < max_set){
                    memcpy(opt_set[crtl_set_idx], set_gen, k*sizeof(vertex_type));
                    crtl_set_idx++;
                }
            }
            free_digraph(diff_g);
            free_cc(scc);
            idx_s++;
        }
        free(set_gen);
    }
	*end_conn = min_conn;
	*tot_sets = n_crtl_set;

	return opt_set;
}


/*
 * Function: get_reverse_digraph
 * ------------------------------------------------
 * Return a pointer to the reverse digraph of <g>.
 *
 * @param g		The digraph to revert.
 *
 * @return		A pointer to the reverse digraph.
 */
digraph_type *get_reverse_digraph(digraph_type *g) {
    digraph_type *g_rev;
	g_rev = (digraph_type *) Malloc(sizeof(digraph_type),__FILE__, __LINE__, "Reverse digraph");
    g_rev->out_idx  = g->in_idx;
    g_rev->out_neighbours = g->in_neighbours;
    g_rev->in_idx   = g->out_idx;
    g_rev->in_neighbours  = g->out_neighbours;
    g_rev->labels = g->labels;
    g_rev->nodes = g->nodes;
    g_rev->edges = g->edges;
    g_rev->del_nodes_list = g->del_nodes_list;
    g_rev->del_nodes_num = g->del_nodes_num;
    g_rev->org_nodes_num = g->org_nodes_num;
    return g_rev;
}


static int id_asc (const void * a, const void * b) {
   return ( *(vertex_type*)a - *(vertex_type*)b );
}
/*
 * Function: get_subgraph_graph
 * ------------------------------------------------
 * Create a new graph from a graph <g> containing only the nodes in the list <nodes>.
 *
 * @params g		 The graph to use as source
 * @params nodes	 The list of nodes to keep
 * @params n_nodes   The number of nodes to keep, i.e. the size of <nodes>
 *
 * @return	The subgraph of <g> containing only the nodes in <nodes>  
 */
graph_type *get_subgraph_graph(graph_type *g, vertex_type *nodes, vertex_type n_nodes){
	
	graph_type *sub_g;
	vertex_type *r_nodes, to_remove;
	
	to_remove = g->nodes-n_nodes;
	r_nodes = (vertex_type *) Malloc(sizeof(vertex_type)*(to_remove),__FILE__, __LINE__, "Nodes to remove");
	
	//fprintf(stderr, "(%s)[DEBUG] Sort nodes to keep.\n",get_time_string());
	qsort(nodes, n_nodes, sizeof(vertex_type), id_asc);
	//fprintf(stderr, "(%s)[DEBUG] Identifying nodes to remove.\n",get_time_string());
	for (vertex_type i=1, j=0, k=0; i<=g->nodes; i++ ){
		if (i == nodes[j] && j<n_nodes){
			j++;
		}else{
			r_nodes[k]=i;
			k++;
		}
	}
	//fprintf(stderr, "(%s)[DEBUG] Deleting nodes from graph.\n",get_time_string());
	sub_g = delete_nodes_from_graph(g, r_nodes , to_remove);
	
	free(r_nodes);
	return sub_g;
}


/*
 * Function: get_subgraph_digraph
 * ------------------------------------------------
 * Create a new graph from a graph <g> containing only the nodes in the list <nodes>.
 *
 * @params g		 The graph to use as source
 * @params nodes	 The list of nodes to keep
 * @params n_nodes   The number of nodes to keep, i.e. the size of <nodes>
 *
 * @return		A new graph, the subgraph <sub_g> of <g> containing only the nodes in <nodes>  
 */
digraph_type *get_subgraph_digraph(digraph_type *g, vertex_type *nodes, vertex_type n_nodes){
	
	digraph_type *sub_g;
	vertex_type *r_nodes, to_remove;
	
	to_remove = g->nodes-n_nodes;
	r_nodes = (vertex_type *) Malloc(sizeof(vertex_type)*(to_remove),__FILE__, __LINE__, "Nodes to remove");
	
	//fprintf(stderr, "(%s)[DEBUG] Sort nodes to keep.\n",get_time_string());
	qsort(nodes, n_nodes, sizeof(vertex_type), id_asc);
	//fprintf(stderr, "(%s)[DEBUG] Identifying nodes to remove.\n",get_time_string());
	for (vertex_type i=1, j=0, k=0; i<=g->nodes; i++ ){
		if (j < n_nodes && i == nodes[j] ){
			j++;
		}else{
			r_nodes[k]=i;
			k++;
		}
	}
	//fprintf(stderr, "(%s)[DEBUG] Deleting nodes from digraph.\n",get_time_string());
	sub_g = delete_nodes_from_digraph(g, r_nodes , to_remove);
	
	free(r_nodes);
	return sub_g;
}


/*
 * Function: get_dominators_naive
 * ------------------------------------------------
 * Compute the dominators Tree using a sequenze of DFS for each node v in V wiht
 * w != r, where r is the root of the digraph g. The funciont assumes that
 * it is possible to reach any vertex in the graph from the root.
 *
 * Determine, by means of a search from r, the set S of vertices reachable from r 
 * by paths which avoid v. The verticesin V - {v} - S are exactly those which v dominates.
 * Knowing the set of vertices dominated by each vertex, construct the dominator tree.
 *
 * @param g			The flow graph for which we compute the dominator tree.
 * @param root		The node to use as root for the flow graph.
 *
 * @return			The dominator tree <dom> represented as a parents array, i.e.
 *					dom[v]=w means that v dominates w, and there is an edge v->w
 *					in the dominator tree.
 */
vertex_type * get_dominators_pm(digraph_type *g, vertex_type root ){
	vertex_type *dom;
	vertex_type *visited;				/* List of visited nodes */
	vertex_type *stack;					/* Stack structure used for the visit of the graph */
	vertex_type stack_idx, n_nodes, n_visited;
	char proceed;
	vertex_type *tmp_dom_s;
	vertex_type **tmp_dom;

	n_nodes = g->nodes;
	visited = (vertex_type *) Malloc(sizeof(vertex_type)*(n_nodes+1),__FILE__, __LINE__, "Visited array.");
	stack = (vertex_type *) Malloc(sizeof(vertex_type)*(g->edges*2),__FILE__, __LINE__, "Stack array.");
	dom = (vertex_type *) Malloc(sizeof(vertex_type)*(n_nodes+1),__FILE__, __LINE__, "Dominators dom");
	// list of nodes dominated by v tmp_dom[v]
	tmp_dom = (vertex_type **) Malloc(sizeof(vertex_type*)*(n_nodes+1),__FILE__, __LINE__, "Dominators tmp_dom");
	// number of nodes dominated by v tmp_dom_s[v]
	tmp_dom_s = (vertex_type *) Malloc(sizeof(vertex_type)*(n_nodes+1),__FILE__, __LINE__, "Dominators tmp_dom_s");

	// Init Dominators tree
	for ( vertex_type i=1; i<=n_nodes; i++ ){  dom[i] = root;  }
	dom[root] = 0;
	tmp_dom_s[root] = 0;

	// Compute dominance node by node
	for ( vertex_type i=1; i<=n_nodes; i++ ){
		digraph_type *gd;
		if( i == root ){ continue; }
		// Remove node <i> to search the set of nodes dominated by <i>
		vertex_type r_nodes[] = {i};
		gd = delete_nodes_from_digraph(g, r_nodes , 1);
		
		memset(visited, 0, (n_nodes+1)*sizeof(vertex_type));
		// find the label of the root node in the new graph <gd>
		for ( vertex_type t=1; t < n_nodes; t++ ){  
			  if (gd->labels[t] == root){
				  stack_idx = 1;
				  stack[stack_idx] = t;
				  break;
			  }
		}
		proceed = 1;
		n_visited = 0;
		while (proceed == 1){ 
			vertex_type v;
			v = stack[stack_idx];
			stack_idx--;
			if ( visited[v] == 0){ /* New node, not yet visited */
				visited[v] = 1;
				n_visited++;
				/* Put the neighbours of the node in the stack*/
				for(vertex_type j = gd->out_idx[v+1]-1; j >= gd->out_idx[v]; j--){
					vertex_type w;
					w = gd->out_neighbours[j];
					if( visited[w] == 0 ){
						stack_idx++;
						stack[stack_idx] = w;
					}
				}
			}
			if( stack_idx == 0){ /* This was the last node in the stack */
				proceed = 0;
			}
		}
		
		tmp_dom_s[i] = n_nodes-n_visited-1;
		if (tmp_dom_s[i] > 0 ){ // Only if i dominates at least one node
			tmp_dom[i] = (vertex_type *) Malloc(sizeof(vertex_type)*tmp_dom_s[i],__FILE__, __LINE__, "Dominators dom");
			for ( vertex_type j=1, idx = 0; j<=n_nodes-1; j++ ){ 
				vertex_type j_lab;
				j_lab = gd->labels[j];
				if( i != j_lab && visited[j] == 0){
					tmp_dom[i][idx] = j_lab;
					idx++;
				} 
			}
		}
		free_digraph(gd);
	}
	free(stack);
	
	for( vertex_type j=1 ; j<=n_nodes; j++ ){ // per ogni nodo j del grafo
		if( j == root ){ continue; }
		memset(visited, 0, (n_nodes+1)*sizeof(vertex_type));
		for( vertex_type i=0; i < tmp_dom_s[j]; i++ ){ // per ogni nodo v dominato da j
			vertex_type v = tmp_dom[j][i];
			visited[v]++;
			for( vertex_type h=0; h < tmp_dom_s[v]; h++ ){ // per ogni nodo w dominato da v
				vertex_type w = tmp_dom[v][h];
				visited[w]++;
			}
		}
		for( vertex_type v=1; v <= n_nodes; v++ ){
			if (visited[v] == 1) {
				dom[v] = j;
			}
		}
	}
	
	for( vertex_type j=1 ; j<=n_nodes; j++ ){
		if ( tmp_dom_s[j] > 0 ){
			free(tmp_dom[j]);
		}
	}
	free(tmp_dom_s);
	free(tmp_dom);
	free(visited);
	return dom;
}



/*
 * Function: eval
 * ------------------------------------------------
 * The function evaluate a vertex v, if v is the root of a tree in the forest, 
 * the function return v. Otherwise, let r be the root of the tree in the forest 
 * which contains v. Return any vertex u != r of minimum semi(u) on the path r->v.
 *
 * @param v         The vertex to evaluate
 * @param ancestor  The Forest used by the function.
 * @param label     The vertex u to return while evaluating v.
 * @param dfs		The DFS tree.	
 * @param stack		The stack used to visit nodes.
 * @param sdom		The semidominator tree. 
 *
 *
 * @return		If v is the root of a tree in the forest, return v. Otherwise, 
 * 				let r be the root of the tree in the forest which contains v.
 *				Return any vertex u != r of minimum semi(u) on the path r->v.
 */
static vertex_type eval(vertex_type v, vertex_type *ancestor, vertex_type *label, forest_type *dfs, vertex_type *stack, vertex_type *sdom){ 
	vertex_type z, root;
	vertex_type stack_idx;
	
	stack_idx = 1;	
	stack[stack_idx] = v;
	z = v;
	// determine the sequence v = v_k, v_h, ..., v_0 = r
	while( ancestor[z] != 0 ){
		z = ancestor[z];
		stack_idx++;
		stack[stack_idx] = z;
	}
	root = stack[stack_idx];
	while ( stack_idx != 0 ) {
		vertex_type q;
		q = stack[stack_idx];
		stack_idx--;
		if (ancestor[q] != 0 ){ // there is an edge z->q
			z = ancestor[q];
			ancestor[q] = root;
			if( dfs->order[sdom[label[z]]] < dfs->order[sdom[label[q]]] ){
				label[q] = label[z];
			}
		}
	}
	return label[v];	
} 

/*
 * Function: get_dominators
 * ------------------------------------------------
 * Compute the dominators Tree using the lengauer-tarjan algorithm. The algorithm
 * assumes that the dfs is a tree and not a forest, thus it is possible to reach
 * any vertex in the graph from the root.
 *
 * @param g			The graph to use for the computation.
 * @param dfs       The DFS visit of the graph <g>
 *
 * @return			The dominator tree <dom> represented as a parents array, i.e.
 *					dom[v]=w means that v dominates w, and there is an edge v->w
 *					in the dominator tree.
 */
vertex_type *get_dominators_lt(digraph_type *g, forest_type *dfs){
	vertex_type *sdom, *dom;
	vertex_type n_nodes, u;
	vertex_type *label, *ancestor;
	vertex_type *stack, *bucket;
	
	n_nodes = g->nodes;
	sdom = (vertex_type *) Malloc(sizeof(vertex_type)*(n_nodes+1),__FILE__, __LINE__, "sdom");
	dom = (vertex_type *) Malloc(sizeof(vertex_type)*(n_nodes+1),__FILE__, __LINE__, "sdom");
	stack = (vertex_type *) Malloc(sizeof(vertex_type)*(n_nodes+1),__FILE__, __LINE__, "stack");
	// Initially ancestor(v) = 0 and label(v) = v for each vertex v -- ancestor(v) = 0 only if v is a root
    ancestor = (vertex_type *) Calloc( n_nodes+1, sizeof(vertex_type),__FILE__, __LINE__, "ancestor");
	label = (vertex_type *) Malloc(sizeof(vertex_type)*(n_nodes+1),__FILE__, __LINE__, "label");
	// bucket(w): A set of vertices whose semidominator is w. Each vertex has only one semidonimator.
	bucket = (vertex_type *) Calloc( n_nodes+1, sizeof(vertex_type),__FILE__, __LINE__, "ancestor");

	for ( vertex_type i=0; i<=n_nodes; i++ ){ // initialize sdom structure
		sdom[i] = i;
		label[i] = i;
	}
	// Visit nodes in reverse preorder, skip dfs->vertex[v==1] that is the root
	for ( vertex_type i=n_nodes; i > 1; i-- ){ 
		vertex_type w;
		w = dfs->vertex[i];
		// Evaluate each predecessor of w 
		for(vertex_type j = g->in_idx[w+1]-1; j >= g->in_idx[w]; j--){
			vertex_type v;
			v = g->in_neighbours[j];
			u = eval(v, ancestor, label, dfs, stack, sdom); 
			if ( dfs->order[sdom[u]] < dfs->order[sdom[w]] ) {
				sdom[w] = sdom[u];
			}
		}
		// Add the link to the forest used by the <eval> procedure
		ancestor[w] = dfs->parent[w];
		//add w to bucket(vertex(semi(w))
		bucket[w] = sdom[w];
		
		//for each v in bucket(parent(w))
		for( vertex_type v = n_nodes; v > 0; v-- ){
			if( bucket[v] == dfs->parent[w]){ // v is in the bucket(parent(w))
				bucket[v] = 0;
				u = eval(v, ancestor, label, dfs, stack, sdom); 
				if ( dfs->order[sdom[u]] < dfs->order[sdom[v]] ) {
					dom[v] = u;
				}else{
					dom[v] = dfs->parent[w];
				}
				
			}
		}
	}
	
	free(label);
	free(bucket);
	free(stack);
	free(ancestor);
	
	for ( vertex_type i=2; i <= n_nodes; i++ ){ 
		vertex_type w;
		w = dfs->vertex[i];
		if( dom[w] != sdom[w]){
			dom[w] = dom[dom[w]];
		}
	}
	dom[dfs->vertex[1]] = 0;
	
	free(sdom);
	return dom;
}


/*
 * Function: write_dominators_tree_dot
 * ------------------------------------------------
 * The function writes on a file a Dominators Tree in dot format.
 *
 * @param parent		The Dominators Tree to write   
 * @param nodes			The number of node in the Tree
 * @param file_name		The name of the file.
 * @param type			The type of graph to write, possible values are: DIGRAPH or GRAPH.
 */
void write_dominators_tree_dot(vertex_type *parent, vertex_type nodes, const char* file_name, char type){
	FILE *fp = Fopen(file_name, "w", __FILE__, __LINE__, "Opening file to write.");
	if ( type == DIGRAPH ){
		fprintf(fp, "digraph G {\n");
		for (vertex_type i = 1; i < nodes+1; i++) {
			vertex_type v;
			v = parent[i];
			if( v != 0 ){
				fprintf(fp, "\t%"PRIvertex" -> %"PRIvertex" ;\n", v, i);
			}
		}
		fprintf(fp, "}\n");
	}else if ( type == GRAPH ){
		fprintf(fp, "graph G {\n");
		for (vertex_type i = 1; i < nodes+1; i++) {
			vertex_type v;
			v = parent[i];
			if( v != 0 ){ /* We visited the node and i is not a root of a tree */
				fprintf(fp, "\t%"PRIvertex" -- %"PRIvertex" ;\n", v, i);
			}
		}
		fprintf(fp, "}\n");
	}
	fclose(fp);
}

/*
 * Function: write_dominators_tree_edge_list
 * ------------------------------------------------
 * The function writes on a file a Dominators Tree in edges list format.
 *
 * @param parent		The parent array representing the Tree.
 * @param nodes			The number of nodes in the Tree.
 * @param file_name		The name of the file.
 */
void write_dominators_tree_edge_list(vertex_type *parent, vertex_type nodes, const char* file_name){
	FILE *fp = Fopen(file_name, "w", __FILE__, __LINE__, "Opening file to write.");
	fprintf(fp, "# Nodes: %"PRIvertex" Edges: %"PRIvertex, nodes, nodes-1);
	for (vertex_type i = 1; i < nodes+1; i++) {
		vertex_type v;
		v = parent[i];
		if( v != 0 ){ 
			fprintf(fp, "\n%"PRIvertex" %"PRIvertex, v, i);
		}
	}
	fclose(fp);
}


/*
 * Function: get_graph_from_digraph
 * ------------------------------------------------
 * The function return an undirected graph from a digraph <g>.
 * 
 * @param g		The digraph to transform into an undirected graph.
 *
 * @return		A new graph from the directed graph, NULL if an error occurred.
 *
 */
graph_type * get_graph_from_digraph(digraph_type *g){
	vertex_type *edge_list;
	vertex_type edge_idx;
	
	edge_list = (vertex_type*) Malloc(sizeof(vertex_type)*g->edges*2, __FILE__, __LINE__, "Allocating edge_list");	
	edge_idx = 0;
	
	for(vertex_type v = 1; v <= g->nodes ; v++){
		for(vertex_type j = g->out_idx[v]; j <= g->out_idx[v+1]-1; j++){
			char check;
			vertex_type w;
			check = 0;
			w = g->out_neighbours[j];
			// CHECK if w has an id smaller than v, check if v is an out neighbour of w, 
			// because in that case the edge has already been added
			if( w < v ){
				for(vertex_type k = g->out_idx[w]; k <= g->out_idx[w+1]-1; k++){
					if ( g->out_neighbours[k] == v){
						check = 1;
						//printf("DUPLICATE %"PRIvertex" -> %"PRIvertex"\n", v, w);
						break;
					}
				}
			}
			if (check == 0){
				edge_list[edge_idx] = v;
				edge_list[edge_idx+1] = w;
				edge_idx += 2;
				//printf("%"PRIvertex" -> %"PRIvertex"\n", v, w);
			}
		}
		
	}
		
	graph_type *ug = get_graph_from_edge_list( edge_list,  g->nodes, edge_idx/2);
	free(edge_list);
	
	return ug;
}



/*
 * Function: get_degree
 * ------------------------------------------------
 * The function computes degree of each node in the graph <g>.
 *
 * @param g     	The graph.   
 *
 * @return			The function returns a pointer to a degree_type array containing the degrees of 
 *					the nodes. A NULL pointer is returned if an error occurs.
 */
degree_type *get_degree(graph_type *g){
	degree_type *nodes_degree;
	
	if (g == NULL){ return NULL; }
	nodes_degree = (degree_type *) Malloc( g->nodes*sizeof(degree_type), __FILE__, __LINE__, "malloc nodes_degree.");

	for (vertex_type i = 1; i < g->nodes+1; i++) {
		nodes_degree[i-1].v = i;
		nodes_degree[i-1].deg = g->idx[i+1]-g->idx[i];
	}
	return nodes_degree;
}


/*
 * Function: write_degree
 * ------------------------------------------------
 * The function prints the degree of the nodes in <d>.
 *
 * @param d     		The degree of the nodes. 
 * @param n_nodes  		Number of graph nodes.
 * @param top_n  		Number of node to print.
 * @param file			The name of the file used to write degree information
 */
void write_degree(degree_type *d, vertex_type n_nodes, vertex_type top_n, const char *file_name){
	FILE *file = Fopen(file_name, "w", __FILE__, __LINE__, "Opening file to write.");
	
	if (d == NULL){ return; }
	vertex_type max_nodes = MIN(n_nodes, top_n);
	fprintf(file, "Node,Deg\n");
	for(vertex_type i=0; i < max_nodes; i++){
		fprintf(file, "%"PRIvertex",%"PRIvertex"\n", d[i].v, d[i].deg);
	}
	fclose(file);
}

/*
 * Function: print_degree
 * ------------------------------------------------
 * The function prints the degree of the nodes in <d>.
 *
 * @param d     		The degree of the nodes. 
 * @param n_nodes  		Number of graph nodes.
 * @param top_n  		Number of node to print.
 */
void print_degree(degree_type *d, vertex_type n_nodes, vertex_type top_n){	
	if (d == NULL){ return; }
	vertex_type max_nodes = MIN(n_nodes, top_n);
	printf( "Node,Deg\n");
	for(vertex_type i=0; i < max_nodes; i++){
		printf( "%"PRIvertex",%"PRIvertex"\n", d[i].v, d[i].deg);
	}
}

static int deg_asc (const void * a, const void * b) {
   return ( (*(degree_type*)a).deg - (*(degree_type*)b).deg );
}
static int deg_desc (const void * a, const void * b) {
   return ( (*(degree_type*)b).deg - (*(degree_type*)a).deg );
}
/*
 * Function: sort_degree
 * ------------------------------------------------
 * The function sorts the degree values in <d>.
 *
 * @param d			The array of dgree_type to sort.  
 * @param n_nodes   Number of nodes in the graph.
 * @param order		The order of the sort, possible values are: ASC or DESC.
 */
void sort_degree(degree_type *d, vertex_type n_nodes, char order){
	if (order == ASC){
		qsort(d, n_nodes, sizeof(degree_type), deg_asc);
	}else if (order == DESC){
		qsort(d, n_nodes, sizeof(degree_type), deg_desc);
	}
}



/*
 * Function: get_digraph_wcc
 * ------------------------------------------------
 * The function computes the Weakly Connected Components (WCC) of a digraph g. 
 *
 * @param g     	The pointer to the digraph. 
 *
 * @return			The funciton returns the pointer to a cc data stracture containing the WCC
 *					decomposition of the digraph g. 
 */
cc_type *get_digraph_wcc(digraph_type *g){
	vertex_type *queue;					      /* Queue structure used to visit the graph */
	vertex_type head, tail;       			  /* Head and Tail index of the queue */
	vertex_type n_nodes, n_visited;           /* Number of nodes in the graph */
	vertex_type last_unvisited;			      /* The node with the lowest index we haven't yet visited */
	char proceed;
	vertex_type v, w;
	cc_type *wcc;
	vertex_type n_wcc;
	/* visited_in_queue[v] has value 1 if v is in the queue, has value 2 if v has been visited, i.e.
	 * v has already been extracted from the queue */
	char *visited_in_queue;					  
	
	
	n_nodes = g->nodes;
	wcc = (cc_type*) Malloc(sizeof(cc_type), __FILE__, __LINE__, "CC." );
	wcc->vertex = (vertex_type *) Malloc(sizeof(vertex_type)*(n_nodes+1), __FILE__, __LINE__, "cc->vertex array.");
	wcc->idx = (vertex_type *) Malloc(sizeof(vertex_type)*(n_nodes+1), __FILE__, __LINE__, "cc->idx array.");
	visited_in_queue = (char *) Calloc(n_nodes+1, sizeof(char),__FILE__, __LINE__, "Visited array.");
	queue = (vertex_type *) Malloc(sizeof(vertex_type)*n_nodes,__FILE__, __LINE__, "Queue array.");
	head = 0;
	n_visited = 0;
	tail = 1;
	queue[head] = 1;
	last_unvisited = 1;
	proceed = 1;
	wcc->idx[0] = 0;
	n_wcc = 1;
	visited_in_queue[1] = 1;
	
	while( proceed == 1){
		v = queue[head];
		head++;
		if (visited_in_queue[v] == 1 ){ /* New node to visit */
			visited_in_queue[v] = 2;
			wcc->vertex[n_visited] = v;
			n_visited++;
			for(vertex_type j = g->out_idx[v]; j <= g->out_idx[v+1]-1; j++){
				w = g->out_neighbours[j];
				if( visited_in_queue[w] == 0 ){
					queue[tail] = w;
					tail++;
					visited_in_queue[w] = 1;
				}
			}
			for(vertex_type j = g->in_idx[v]; j <= g->in_idx[v+1]-1; j++){
				w = g->in_neighbours[j];
				if( visited_in_queue[w] == 0 ){
					queue[tail] = w;
					tail++;
					visited_in_queue[w] = 1;
				}
			}
		}
		
		if (head == tail){
			if (n_visited == n_nodes){ /* check if we visited all nodes */
				proceed = 0;
				wcc->idx[n_wcc] = n_visited;
			}else{
				for(vertex_type i = last_unvisited; i < n_nodes+1; i++ ){ 
					if(visited_in_queue[i] == 0){
						wcc->idx[n_wcc] = n_visited;
						n_wcc++; 
						queue[tail] = i;
						visited_in_queue[i] = 1;
						tail++;
						last_unvisited = i+1;
						break;
					}
				}
			}	
		}
	}
	wcc->number = n_wcc;
	free(visited_in_queue);
	free(queue);
	return wcc;
	
}


/*
 * Function: get_graph_cc
 * ------------------------------------------------
 * The function computes the Connected Components (CC) of a graph g.
 *
 * @param g     	The pointer to the digraph. 
 *
 * @return			The funciton returns the pointer to a cc data structure containing the CC
 *					decomposition of the digraph g. 
 */
cc_type *get_graph_cc(graph_type *g){
	vertex_type *queue;					      /* Queue structure used to visit the graph */
	vertex_type head, tail;       			  /* Head and Tail index of the queue */
	vertex_type n_nodes, n_visited;           /* Number of nodes in the graph */
	vertex_type last_unvisited;			      /* The node with the lowest index we haven't yet visited */
	char proceed, *visited_in_queue;
	vertex_type v, w;
	cc_type *cc;
	vertex_type n_cc;
	
	n_nodes = g->nodes;
	cc = (cc_type*) Malloc(sizeof(cc_type), __FILE__, __LINE__, "CC." );
	cc->vertex = (vertex_type *) Malloc(sizeof(vertex_type)*(n_nodes+1), __FILE__, __LINE__, "cc->vertex array.");
	cc->idx = (vertex_type *) Malloc(sizeof(vertex_type)*(n_nodes+1), __FILE__, __LINE__, "cc->idx array.");
	visited_in_queue = (char *) Calloc(n_nodes+1, sizeof(vertex_type),__FILE__, __LINE__, "Visited array.");
	queue = (vertex_type *) Malloc(sizeof(vertex_type)*n_nodes,__FILE__, __LINE__, "Queue array.");
	head = 0;
	n_visited = 0;
	tail = 1;
	queue[head] = 1;
	last_unvisited = 1;
	proceed = 1;
	cc->idx[0] = 0;
	n_cc = 1;
	visited_in_queue[1]=1;
	
	while( proceed == 1){
		v = queue[head];
		head++;
		if (visited_in_queue[v] == 1 ){ /* New node to visit */
			visited_in_queue[v] = 2;
			cc->vertex[n_visited] = v;
			n_visited++;
			for(vertex_type j = g->idx[v]; j <= g->idx[v+1]-1; j++){
				w = g->neighbours[j];
				if( visited_in_queue[w] == 0 ){
					queue[tail] = w;
					tail++;
					visited_in_queue[w] = 1;
				}
			}
		}
		
		if (head == tail){
			if (n_visited == n_nodes){ /* check if we visited all nodes */
				proceed = 0;
				cc->idx[n_cc] = n_visited;
			}else{
				for(vertex_type i = last_unvisited; i < n_nodes+1; i++ ){ 
					if(visited_in_queue[i] == 0){
						cc->idx[n_cc] = n_visited;
						n_cc++; 
						queue[tail] = i;
						tail++;
						last_unvisited = i+1;
						visited_in_queue[i] = 1;
						break;
					}
				}
			}	
		}
	}
	cc->number = n_cc;
	free(visited_in_queue);
	free(queue);
	return cc;
}



/*
 * Function: get_cc_distribution
 * ------------------------------------------------
 * The function returns the distribution of the sizes of the CC. 
 *
 * @param cc     	The CC decomposition of the graph.
 * @param size		The size of the array returned by the function.  
 *
 * @return			An array containing the sizes of the CC.
 */
uint64_t * get_cc_distribution (cc_type *cc, uint64_t *size){
	uint64_t *cc_dist;
	
	*size = cc->number;
	cc_dist = (uint64_t *) Malloc( sizeof(uint64_t)*(cc->number), __FILE__, __LINE__, "cc_dist array");
	for(vertex_type i = 0; i < cc->number; i++){
		cc_dist[i] = cc->idx[i+1]-cc->idx[i];
	}
	return cc_dist;
}



/*
 * Function: write_cc_distribution
 * ------------------------------------------------
 * The function writes the distribution of the sizes of the CC. 
 *
 * @param cc_dist       Array containing the sizes of the CC.
 * @param size 			Lenght of <cc_dis>.
 * @param file_name		The name of the file.
 */
void write_cc_distribution (uint64_t *cc_dist, uint64_t size, const char *file_name){
	FILE *file = Fopen(file_name, "w", __FILE__, __LINE__, "Opening file to write.");	
	for (uint64_t i = 0; i<size-1; i++){
		fprintf(file, "%"PRIu64", ", cc_dist[i]);
	}
	fprintf(file, "%"PRIu64"\n", cc_dist[size-1]);
    fclose(file);
	return;
}



/*
 * Function: print_cc_distribution
 * ------------------------------------------------
 * The function prints on the standard output the distribution of the sizes of the CC. 
 *
 * @param cc_dist       Array containing the sizes of the CC.
 * @param size 			Lenght of <cc_dis>.
 */
void print_cc_distribution (uint64_t *cc_dist, uint64_t size){
	for (uint64_t i = 0; i<size-1; i++){
		printf("%"PRIu64", ", cc_dist[i]);
	}
	printf( "%"PRIu64"\n", cc_dist[size-1]);	
	return;
}



/*
 * Function: get_giant_cc_vertices
 * ------------------------------------------------
 * The function returns the list of vertices in the giant CC of the graph.
 *
 * @param cc        The CC decomposition of the graph.
 * @param size		The size of the returned array.
 *
 * @return			An array containing the list of vertices in the giant CC.      	 
 */
vertex_type *get_giant_cc_vertices (cc_type *cc, uint64_t *size){
	vertex_type id_cc;
	vertex_type *vertices, i, j;
	
	*size = 0;
	id_cc = 0;
	vertices = NULL;	
	for( i = 0; i < cc->number; i++){
		if ( (*size) <  (cc->idx[i+1]-cc->idx[i]) ) {
			*size = cc->idx[i+1]-cc->idx[i];
			id_cc = i;
		}
	}
	vertices = (vertex_type *) Malloc(sizeof(vertex_type)*(*size),__FILE__, __LINE__, "Giant CC vertices array");
	for(i=0, j = cc->idx[id_cc]; j < cc->idx[id_cc+1]; j++,i++){
		vertices[i] = cc->vertex[j];
	}
	
	return vertices;
}


/*
 * Function: write_digraph_labels
 * ------------------------------------------------
 * The function writes to a file the original labels of the nodes in a 
 * csv format: NODE_ID, ORIGINAL_ID
 *
 * @param g        		The digraph
 * @param file_name		The name of the file.
 *   	 
 */
void write_digraph_labels (digraph_type *g, const char *file_name){
	FILE *file = Fopen(file_name, "w", __FILE__, __LINE__, "Opening file to write.");	
	fprintf(file, "NodeID, OriginalID\n");
	for(vertex_type i=1; i < g->nodes; i++ ){
		fprintf(file, "%"PRIvertex", %"PRIvertex" \n", i, g->labels[i]);
	}
	fclose(file);
	return;
}



/*
 * Function: write_graph_labels
 * ------------------------------------------------
 * The function writes to a file the original labels of the nodes in a 
 * csv format: NODE_ID, ORIGINAL_ID
 *
 * @param g        		The graph
 * @param file_name		The name of the file.
 *   	 
 */
void write_graph_labels (graph_type *g, const char *file_name){
	FILE *file = Fopen(file_name, "w", __FILE__, __LINE__, "Opening file to write.");	
	fprintf(file, "NodeID, OriginalID\n");
	for(vertex_type i=1; i < g->nodes; i++ ){
		fprintf(file, "%"PRIvertex", %"PRIvertex" \n", i, g->labels[i]);
	}
	fclose(file);
	return;
}





