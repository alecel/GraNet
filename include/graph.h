#ifndef CRANIC_GRAPH_H_
#define CRANIC_GRAPH_H_


#include "utility.h"
#include <inttypes.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>


#ifdef LARGE
	typedef uint64_t vertex_type;
	#define PRIvertex PRIu64
#else
	typedef uint32_t vertex_type;
	#define PRIvertex PRIu32
#endif
	
/* Supported graph file formats */	
#define EDGE_LIST 0
	
/* Supported graph types */
#define DIGRAPH 	0
#define FLOWGRAPH	1
#define GRAPH		2



/*
 * Type: cc_type
 * ------------------------------------------------
 * The CC decomposition of a graph/digraph is implemented using two arrays:
 * 		- <vertex> contains the nodes of the graph;
 *		- <idx> contains the pointers to the nodes in <vertex> (in the same CC).
 *
 * The vertices of CC <i> are in the vector <vertex> from index idx[i] included to idx[i+1] excluded. 
 *	
 * Example:
 * 		vertex[idx[i]], ..., vertex[idx[i+1]-1] are in CC <i>.
 *
 *		idx[i+1]-idx[i] is the size of CC <i>.
 *	
 */
typedef struct cc {
    vertex_type *vertex;      /* List of vertices ordered by CC, size |V|=n */
    vertex_type *idx;         /* Indexes of the vertices in the same CC, size <number>+1 */
    vertex_type number;       /* Number of CC */
} cc_type;



/*
 * Type: digraph_type
 * ------------------------------------------------
 * A Digraph is stored using two adjacent lists representing IN and OUT edges.
 * Each list is implemented using 2 arrays:
 * 		- <in_idx> and <in_neighbours> for IN edges;
 *		- <out_idx> and <out_neighbours> for OUT edges.
 * Nodes' ids start from 1.
 *
 * Given a graph G=(V,E) with |V|=n and |E|=m. 
 * The arrays <in_idx> and <out_idx> have size n+2, while the arrays 
 * <in_neighbours> and <out_neighbours> have size m+1. 
 * The positions with index 0 of these 4 arrays are never used for data but only 
 * for implementation purposes.
 * in_idx[0] and out_idx[0] are always set to 0.
 * in_idx[1] and out_idx[1] are always set to 1.
 * The arrays <in_idx> and <out_idx> store the start and end index of the
 * IN and OUT neighbours of each node v of G.
 * Given a node v of G, in_idx[v] is the start index of its IN neighbours
 * in <in_neighbours> and in_idx[v+1] is the end index of its IN neighbours.
 *
 * Example:
 *
 *		out_idx[v+1] - out_idx[v] is the number of OUT neighbours of node v.
 *
 *		out_neighbours[out_idx[v]], out_neighbours[out_idx[v]+1], ..., 
 * 		out_neighbours[out_idx[v+1]-1] are the OUT neighbours of node v.
 *
 * The array <labels> has size n+1 and contains the original ID of the nodes, usually v = labels[v].
 * In case of node removal it could happen that v != labels[v] for some v, because nodes
 * are relabelled after removal.
 */
typedef struct digraph {
	vertex_type	*in_idx;			/* indices of neighbours in <in_neighbours>, size n+2 */
	vertex_type	*in_neighbours;		/* list of IN neighbours (edges), size m+1  */
	vertex_type *out_idx;			/* indices of neighbours in <out_neighbours>, size n+2 */
	vertex_type *out_neighbours;	/* list of OUT neighbours (edges), size m+1  */
    vertex_type *labels;			/* nodes labels, size n+1, labels[i] is the original ID of nodes i */
    vertex_type nodes;				/* number of nodes, |V|=n */
    vertex_type edges;				/* number of edges, |E|=m */
	vertex_type org_nodes_num;		/* number of nodes in the original graph */	
	vertex_type	del_nodes_num;		/* number of deleted nodes */
	vertex_type	*del_nodes_list;    /* list of deleted nodes with their original label */
} digraph_type;


/*
 * Type: graph_type
 * ------------------------------------------------
 * A Graph (undirected graph) is stored using an adjacent list.
 * The list is implemented using 2 arrays: <idx> and <neighbours> edges;
 * Nodes' ids start from 1.
 *
 * Given a graph G=(V,E) with |V|=n and |E|=m. The array <idx> has size n+2, while 
 * the array <neighbours> has size 2m+1. The positions with index 0 of these 2 
 * arrays are never used for data but only for implementation purposes.
 * idx[0] is always set to 0, in_idx[1] is always set to 1.
 * The array <idx> stores the start and end index of the neighbours of each node v of G.
 * Given a node v of G, idx[v] is the start index of its neighbours
 * in <neighbours> and idx[v+1] is the end index of its neighbours.
 *
 * Example:
 *
 *		idx[v+1] - idx[v] is the number of neighbours of node v.
 *
 *		neighbours[idx[v]], neighbours[idx[v]+1], ..., 
 * 		neighbours[idx[v+1]-1] are the neighbours of node v.
 *
 * The array <labels> has size n+1 and contains the original ID of the nodes, usually v = labels[v].
 * In case of node removal it could happen that v != labels[v] for some v, because nodes
 * are relabelled after removal.
 */
typedef struct graph {
	vertex_type	*idx;				/* indices of neighbours size n+2 */
	vertex_type	*neighbours;		/* list of neighbours (edges), size 2m+1  */
    vertex_type *labels;			/* nodes labels, size n+1, labels[i] is the original ID of nodes i */
    vertex_type nodes;				/* number of nodes, |V|=n */
    vertex_type edges;				/* number of edges, |E|=m */
	vertex_type org_nodes_num;		/* number of nodes in the original graph */	
	vertex_type	del_nodes_num;		/* number of deleted nodes */
	vertex_type	*del_nodes_list;    /* list of deleted nodes with their original label */
} graph_type;

/*
 * Type: didegree_type
 * ------------------------------------------------
 * The didegree_type structure contains the IN and OUT degree of a vertex v.
 */
typedef struct dideg{
	vertex_type v;					/* vertex id */
	vertex_type in_deg;				/* in degree of vertex <v> */
	vertex_type out_deg;			/* out degree of vertex <v> */
	vertex_type all_deg;			/* out degree of vertex <v> */
} didegree_type;


/*
 * Type: degree_type
 * ------------------------------------------------
 * The degree_type structure contains degree of a vertex v.
 */
typedef struct deg{
	vertex_type v;					/* vertex id */
	vertex_type deg;				/* degree of vertex <v> */
} degree_type;


/*
 * Type: forest_type
 * ------------------------------------------------
 * The forest_type structure represents a DFS or BFS visit of a graph.
 * Given a forest F=(V,E) with |V|=n and |E|=m, three arrays are used to store the forest:
 *
 *     1. <vertex> contains the nodes in visit order. vertex[i] = v means
 *       that v was the i-th node to be visited.
 *   	  
 *     2. <order> contains the order of the nodes. order[v] = i means 
 *	     that v was the i-th to be visited.
 *
 *     3. <parent> containing the parent id of a node, parent[v] = w means that
 *		 w is the parent of v in T.
 *
 * The variables nodes and edges contains the number of nodes and edges in the forest. While
 * the variable all_nodes contains the orginal number of nodes in the graph before the visit.
 *
 * For each v in V if parent[v] == v than v is a root of a tree T in forest F.
 * For each v in V if parent[v] == 0 than v has not been visited.
 */
typedef struct forest{
	vertex_type *order;
	vertex_type *vertex;
	vertex_type *parent;
	vertex_type all_nodes;
	vertex_type nodes;
	vertex_type edges;
} forest_type;


/*
 * Function: get_digraph_from_file
 * ------------------------------------------------
 * The function create a DIGRAPH <digraph_type> from a file. The parameter <file_format>
 * specifies the format of the file. If the file format is not supported a NULL pointer is returned.
 *
 * @param file_name     	The name of the file.
 * @param file_format		The format of the file.
 *
 * @return					The function returns a pointer to a digraph of type digraph_type containing
 * 							the graph. Use the function free_digraph to free the memory. If the file format 
 *                          is not supported or an error occurs a NULL pointer is return      
 */
digraph_type *get_digraph_from_file(const char *file_name, unsigned char file_format);

/*
 * Function: print_digraph
 * ------------------------------------------------
 * The function prints the information of the digraph: number of nodes, number of edges
 * list of edges and list of labels. 
 *
 * @param g     	    The pointer to the digraph to print. 
 * @param label			If 1 print the node original id, if 0 print the current node id.   
 */
void print_digraph(digraph_type *g, char label);

/*
 * Function: write_digraph_dotformat
 * ------------------------------------------------
 * The function writes the digraph <g> on file <file_name> in dot format.
 *
 * @param g     	    The pointer to the digraph to print.   
 * @param file_name		File used to write the graph in dot format.
 */
void write_digraph_dot(digraph_type *g, const char* file_name);

/*
 * Function: write_digraph_edgelist
 * ------------------------------------------------
 * The function writes the digraph <g> on file <file_name> in edge list format.
 *
 * @param g     	    The pointer to the digraph to print.   
 * @param file_name		File used to write the graph in dot format.
 */
void write_digraph_edgelist(digraph_type *g, const char* file_name);

/*
 * Function: free_digraph
 * ------------------------------------------------
 * The function frees the memory allocated for a digraph. 
 *
 * @param g     	The pointer to the digraph to free.       
 */
void free_digraph(digraph_type *g);

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
forest_type *get_digraph_dfs(digraph_type *g, vertex_type root, char unconn);

/*
 * Function: get_digraph_bfs
 * ------------------------------------------------
 * The function does a Breadth First Search (BFS) of a digraph <g>. If <unconn> is 0
 * the search stop when no other nodes can be reached from the <root>. If <unconn>
 * is 1 the search keep visiting nodes until all nodes have been visited.
 *
 * @param g     	 The pointer to the graph.
 * @param root		 The vertex to use as the root for the search.
 * @param unconn	 To decide how to proceed if the graph is disconnected.
 *    
 * @return           The function returns a forest representing the BFS visit of the digraph.   
 */
forest_type *get_digraph_bfs(digraph_type *g, vertex_type root, char unconn);

/*
 * Function: write_forest_dot
 * ------------------------------------------------
 * The function writes on a file the Forest <f> in dot format.
 *
 * @param f				The Forest to write   
 * @param file_name		The name of the file.
 * @param type			The type of graph to write, possible values are: DIGRAPH or GRAPH.
 */
void write_forest_dot(forest_type *f, const char* file_name, char type);

/*
 * Function: write_forest_edge_list
 * ------------------------------------------------
 * The function writes on a file the Forest <f> in edges list format.
 *
 * @param f				The Forest to write.
 * @param file_name		The name of the file.
 */
void write_forest_edge_list(forest_type *f, const char* file_name);

/*
 * Function: print_forest
 * ------------------------------------------------
 * The function prints on the standard output a Forest <f>
 *
 * @param f		The forest to print   
 */
void print_forest(forest_type *f);

/*
 * Function: free_forest
 * ------------------------------------------------
 * The function frees the memory allocated for the Forest <f>.
 *
 * @param f		The Forest to free.       
 */
void free_forest(forest_type *f);


/*
 * Function: get_digraph_scc
 * ------------------------------------------------
 * The function computes the Strongly Connected Components (SCC) of a digraph g. 
 * The function implements the iterative version of the Tarjan algorithm, presented in
 * Tarjan, R. E. (1972), "Depth-first search and linear graph algorithms". 
 *
 * @param g     	The pointer to the digraph. 
 *
 * @return			The funciton returns the pointer to a cc data stracture containing the SCC
 *					composition of the digraph g. 
 */
cc_type *get_digraph_scc(digraph_type *g);

/*
 * Function: print_scc
 * ------------------------------------------------
 * The function prints the CC decomposition of a graph. 
 *
 * @param cc     	 The pointer to the CC to print.  
 */
void print_cc(cc_type *cc);

/*
 * Function: free_cc
 * ------------------------------------------------
 * The function frees the memory allocated for a cc. 
 *
 * @param cc     	The pointer to the cc to free.       
 */
void free_cc(cc_type *cc);

/*
 * Function: get_inoutdegree
 * ------------------------------------------------
 * The function computes the in-degree, out-degree and total degree of each node in the digraph <g>.
 *
 * @param g     	The digraph.   
 *
 * @return			The function returns a pointer to a degree_type array containing the degrees of 
 *					the nodes. A NULL pointer is returned if an error occurs.
 */
didegree_type *get_inoutdegree(digraph_type *g);

/*
 * Function: write_inoutdegree
 * ------------------------------------------------
 * The function prints the degree of the nodes in <d>.
 *
 * @param d     		The degree of the nodes. 
 * @param n_nodes  		Number of graph nodes.
 * @param top_n  		Number of node to print.
 * @param file_name		The name of the file used to write in/out degree information
 */
void write_inoutdegree(didegree_type *d, vertex_type n_nodes, vertex_type top_n, const char *file_name);


/* Sort Order */
#define ASC		0
#define DESC	1
/* Sort Field */
#define INDEG	0
#define OUTDEG	1
#define DEG		2
#define ALLDEG	3
/*
 * Function: sort_inoutdegree
 * ------------------------------------------------
 * The function sorts the degree values in <d>.
 *
 * @param d			The array of didgree_type to sort   
 * @param field		The field of the array to sort, possible values are: INDEG or OUTDEG
 * @param order		The order of the sort, possible values are: ASC or DESC 
 */
void sort_inoutdegree(didegree_type *d, vertex_type n_nodes, char field, char order);

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
uint64_t get_pairs_connectivity(cc_type *cc);


/*
 * Function: get_digraph_critical_nodes_bruteforce
 * --------------------------------------
 * The function checks which set of nodes of size <k> is critical in the graph <g> and returns
 * a list of them. A set of nodes is critical if it minimaizes the connectivity of <g>.
 * Let G be a digraph, let C_1, C_2, ..., C_l be its strongly connected components.
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
vertex_type ** get_digraph_critical_nodes_bruteforce ( digraph_type *g, vertex_type k, uint64_t max_set, uint64_t *end_conn, uint64_t *tot_sets ); 


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
digraph_type *delete_nodes_from_digraph ( digraph_type *g, vertex_type *vts, vertex_type n_vts);


/*
 * Function: get_dominators_naive
 * ------------------------------------------------
 * Compute the dominators Tree using a sequenze of DFS for each node v in V with
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
vertex_type *get_dominators_pm ( digraph_type *g, vertex_type root );


/*
 * Function: get_reverse_digraph
 * ------------------------------------------------
 * Return a pointer to the reverse digraph of <g>.
 *
 * @param g		The digraph to revert.
 *
 * @return		A pointer to the reverse digraph.
 */
digraph_type *get_reverse_digraph ( digraph_type *g );


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
vertex_type *get_dominators_lt(digraph_type *g, forest_type *dfs );



/*
 * Function: get_subgraph_digraph
 * ------------------------------------------------
 * Create a new graph from a graph <g> containing only the nodes in the list <nodes>.
 *
 * @params g		 The graph to use as source
 * @params nodes	 The listo of nodes to keep
 * @params n_nodes   The number of nodes to keep, i.e. the size of <nodes>
 *
 * @return	The subgraph of <g> containing only the nodes in <nodes>  
 */
digraph_type *get_subgraph_digraph(digraph_type *g, vertex_type *nodes, vertex_type n_nodes);


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
graph_type *get_graph_from_file(const char *file_name, unsigned char file_format);

/*
 * Function: print_graph
 * ------------------------------------------------
 * The function prints the information of the graph: number of nodes, number of edges
 * list of edges and list of labels. 
 *
 * @param g     	    The pointer to the graph to print. 
 * @param lable			If 1 print the node original id, if 0 print the current node id.   
 */
void print_graph(graph_type *g, char lable);

/*
 * Function: free_graph
 * ------------------------------------------------
 * The function frees the memory allocated for a graph. 
 *
 * @param g     	The pointer to the graph to free.       
 */
void free_graph(graph_type *g);


/*
 * Function: write_graph_edgelist
 * ------------------------------------------------
 * The function writes the graph <g> on file <file_name> in edge list format.
 *
 * @param g     	    The pointer to the graph to print.   
 * @param file_name		File used to write the graph in edge list format.
 */
void write_graph_edgelist(graph_type *g, const char* file_name);


/*
 * Function: write_graph_dot
 * ------------------------------------------------
 * The function writes the graph <g> on file <file_name> in dot format.
 *
 * @param g     	    The pointer to the graph to print.   
 * @param file_name		File used to write the graph in dot format.
 */
void write_graph_dot(graph_type *g, const char* file_name);


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
 * @return           The function returns a forest representing the BFS visit of the graph.   
 */
forest_type *get_graph_bfs(graph_type *g, vertex_type root, char unconn);


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
forest_type *get_graph_dfs(graph_type *g, vertex_type root, char unconn);


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
graph_type *get_graph_from_digraph(digraph_type *g);


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
degree_type *get_degree(graph_type *g);


/*
 * Function: sort_degree
 * ------------------------------------------------
 * The function sorts the degree values in <d>.
 *
 * @param d			The array of dgree_type to sort.  
 * @param n_nodes   Number of nodes in the graph.
 * @param order		The order of the sort, possible values are: ASC or DESC.
 */
void sort_degree(degree_type *d, vertex_type n_nodes, char order);


/*
 * Function: write_degree
 * ------------------------------------------------
 * The function prints the degree of the nodes in <d>.
 *
 * @param d     		The degree of the nodes. 
 * @param n_nodes  		Number of graph nodes.
 * @param top_n  		Number of node to print.
 * @param file_name		The name of the file used to write degree information
 */
void write_degree(degree_type *d, vertex_type n_nodes, vertex_type top_n, const char *file_name);


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
cc_type *get_digraph_wcc(digraph_type *g);

/*
 * Function: get_graph_cc
 * ------------------------------------------------
 * The function computes the Connected Components (CC) of a graph g.
 *
 * @param g     	The pointer to the digraph. 
 *
 * @return			The funciton returns the pointer to a cc data stracture containing the CC
 *					decomposition of the digraph g. 
 */
cc_type *get_graph_cc(graph_type *g);

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
graph_type *delete_nodes_from_graph(graph_type *g, vertex_type *vts, vertex_type n_vts);


/*
 * Function: print_inoutdegree
 * ------------------------------------------------
 * The function prints the degree of the nodes in <d>.
 *
 * @param d     		The degree of the nodes. 
 * @param n_nodes  		Number of graph nodes.
 * @param top_n  		Number of node to print.
 */
void print_inoutdegree(didegree_type *d, vertex_type n_nodes, vertex_type top_n);
	
/*
 * Function: print_degree
 * ------------------------------------------------
 * The function prints the degree of the nodes in <d>.
 *
 * @param d     		The degree of the nodes. 
 * @param n_nodes  		Number of graph nodes.
 * @param top_n  		Number of node to print.
 */
void print_degree(degree_type *d, vertex_type n_nodes, vertex_type top_n);

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
graph_type *get_subgraph_graph(graph_type *g, vertex_type *nodes, vertex_type n_nodes);


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
uint64_t * get_cc_distribution (cc_type *cc, uint64_t *size);

/*
 * Function: write_cc_distribution
 * ------------------------------------------------
 * The function writes the distribution of the sizes of the CC. 
 *
 * @param cc_dist       Array containing the sizes of the CC.
 * @param size 			Lenght of <cc_dis>.
 * @param file_name		The name of the file.
 */
void write_cc_distribution (uint64_t *cc_dist, uint64_t size, const char *file_name);

/*
 * Function: print_cc_distribution
 * ------------------------------------------------
 * The function prints on the standard output the distribution of the sizes of the CC. 
 *
 * @param cc_dist       Array containing the sizes of the CC.
 * @param size 			Lenght of <cc_dis>.
 */
void print_cc_distribution (uint64_t *cc_dist, uint64_t size);


/*
 * Function: write_cc
 * ------------------------------------------------
 * The function writes the CC decomposition of a graph. 
 *
 * @param cc     	 The pointer to the CC to write.
 * @param file_name  The name of file to write the output.  
 */
void write_cc (cc_type *cc, const char *file_name);


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
vertex_type *get_giant_cc_vertices (cc_type *cc, uint64_t *size);



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
void write_digraph_labels (digraph_type *g, const char *file_name);



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
void write_graph_labels (graph_type *g, const char *file_name);


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
void write_dominators_tree_dot(vertex_type *parent, vertex_type nodes, const char* file_name, char type);


/*
 * Function: write_dominators_tree_edge_list
 * ------------------------------------------------
 * The function writes on a file a Dominators Tree in edges list format.
 *
 * @param parent		The parent array representing the Tree.
 * @param nodes			The number of nodes in the Tree.
 * @param file_name		The name of the file.
 */
void write_dominators_tree_edge_list(vertex_type *parent, vertex_type nodes, const char* file_name);



#endif



