# <a name="top"></a> GraNet
* [GraNet](#granet)
* [Tutorial](#tutorial)
* [Documentation](#documentation)


---

## <a name="granet"></a> GraNet ([top](#top))

GraNet is a C library providing supports for directed and undirected graphs manipulation.

##### Import the library header

```C
#include "granet.h"
``` 
### <a name="deps"></a> Dependencies ([top](#top))

You probably need to install `libbsd-dev` or `libbsd-devel` on your system to get access to `strlcat` and `strlcpy` functions.


---


# <a name="tutorial"></a> Tutorial ([top](#top))
* [Digraphs](#tdigraph)
* [Graphs](#graph)

## <a name="tdigraph"></a> Digraph ([tutorial](#tutorial))

### Create a Digraph

You can create a new digraph reading its structure from a file. Currently the only file format supported is the `EDGE_LIST` format.

```sh
digraph_type *dg;
dg = get_digraph_from_file( file_name, EDGE_LIST);
```

The function `get_digraph_from_file` returns a pointer to a digraph of type `digraph_type`. If the file format is not supported or an error occurs a NULL pointer is returned.  
> 
> Currently the only file format supported to import a digraph is the edges list format (*EDGE_LIST*). A file in this format contains the following information:
> 
>  - A line starting with '#' specifying the number of Nodes and Edges of the graph *(any other line starting with '#' is interpreted as a comment and skipped*).
> -  A line for each edge of the graph containing a pair of numbers, that are the nodes connected by the edge.
> 
> The numbers used to identify nodes must be greater then 0, i.e. the first node id is 1. Below there is an example of a graph in edges list format.  
> 
> ```sh
> # Nodes: 4 Edges: 3
> 1 2
> 2 3
> 3 4
>``` 

### Print the Digraph on the standard output

You can print a digraph to the standard output as follows:

```sh
print_digraph(dg, 0);
```

### Write the Digraph to a File

You can print a digraph to a File in the edgelist format as follows:

```sh
/* Write digraph */
char out_file[256] = "digraph.out.edgelist";
write_digraph_edgelist(dg, get_path(output_dir,out_file));
```

### Remove nodes from the digraph

We can remove any subset of nodes from the digraph. 

```sh
vertex_type vts[2] = {1,2};
ndg = delete_nodes_from_digraph ( dg, vts, 2);
```
The ID of the nodes in the new digraph will start from 1, to retrive the previous ID of each node you can use the lables. To print the digraph with the original ID labels specify 1 as the second parametr of the print\_digraph function.

```sh
print_digraph(dg, 1);
```

### Free

Remember to release the memory allocated once you don't need the digraph anymore, for example after removing some nodes.

```sh
free_digraph(ndg);
```
 
### Compute IN/OUT degree

You can compute the IN and OUT degree of the nodes with the function `get_inoutdegree`.

```sh
didegree_type *deg;
deg = get_inoutdegree(g);
```
The function returns a pointer to a `didegree_type` array containing the degrees of the nodes. A NULL pointer is returned if an error occurs.

You can sort the `didegree_type` array in ascending (ASC) or descending (DESC) order by in or out degree (INDEG or OUTDEG).

```sh
sort_inoutdegree(deg, dg->nodes, OUTDEG, DESC);
sort_inoutdegree(deg, dg->nodes, INDEG, ASC);
```

You can write the `didegree_type` array to a file (write\_inoutdegree) or print them to the standard output (print\_inoutdegree). The output is in [CSV](https://en.wikipedia.org/wiki/Comma-separated_values) format.

```sh
print_inoutdegree( deg, dg->nodes, dg->nodes);
write_inoutdegree( deg, dg->nodes, dg->nodes, file_name);
```

### Compute Strongly Connected Components

You can compute the Strongly Connected Components of the digraph as follows:

```sh
cc_type *scc;
scc = get_digraph_scc(dg);
```

### Compute Weakly Connected Components

You can compute the Weakly Connected Components of the digraph as follows:

```sh
cc_type *wcc;
wcc = get_digraph_wcc(dg);
```

## <a name="graph"></a> Graph ([tutorial](#tutorial))

### Create a Graph

You can create a new graph reading its structure from a file. Currently the only file format supported is the `EDGE_LIST` format.

```sh
graph_type *g;
g = get_graph_from_file( file_name, EDGE_LIST);
```

The function `get_graph_from_file` returns a pointer to a graph of type `graph_type`. If the file format is not supported or an error occurs a NULL pointer is returned. 

### Print the Graph on the standard output

You can print a digraph to the standard output as follows:

```sh
print_graph(g, 0);
```

### Write the Graph to a File

You can print a graph to a File in the edgelist format as follows:

```sh
/* Write graph */
char out_file[256] = "graph.out.edgelist";
write_graph_edgelist(g, get_path(output_dir,out_file));
```

### Remove nodes from the digraph

We can remove any subset of nodes from the graph. 

```sh
vertex_type vts[2] = {1,2};
ndg = delete_nodes_from_graph ( g, vts, 2);
```
The ID of the nodes in the new graph will start from 1, to retrive the previous ID of each node you can use the lables. To print the graph with the original ID labels specify 1 as the second parametr of the print\_graph function.

```sh
print_graph(g, 1);
```

### Free

Remember to release the memory allocated once you don't need the graph anymore, for example after removing some nodes.

```sh
free_graph(ndg);
```

### Compute degree

You can compute the degree of the nodes with the function `get_degree`.

```sh
degree_type *deg;
deg = get_degree(g);
```
The function returns a pointer to a `degree_type` array containing the degrees of the nodes. A NULL pointer is returned if an error occurs.

You can sort the `degree_type` array in ascending (ASC) or descending (DESC) order.

```sh
sort_degree(deg, g->nodes, ASC);
sort_degree(deg, g->nodes, DESC);
```

You can write the `degree_type` array to a file (write\_degree) or print them to the standard output (print\_degree). The output is in [CSV](https://en.wikipedia.org/wiki/Comma-separated_values) format.

```sh
print_degree(deg, g->nodes, 3);
write_degree( deg, dg->nodes, dg->nodes, file_name);
```

### Compute Connected Components

You can compute the Connected Components of the graph as follows:

```sh
cc_type *cc;
cc = get_graph_cc(g);
```



---

# <a name="documentation"></a> Documentation ([top](#top))

* [Digraph](#digraph)
    * [Data Structure](#digraphDS)
    * [Read Digraph](#digraphGet) 
    * [Write Digraph](#digraphWrite)
    * [Digraph Operations](#digraphOp)
* [Graph](#graph)
    * [Data Structure](#graphDS)
    * [Read Graph](#graphGet) 
    * [Write Graph](#graphWrite)
    * [Graph Operations](#graphOp)
* [File Format](#file)
    * [Edge List](#edgelist)
* [Graph Traversal](#traversal)
    * [Depth-First Search](#dfs)
    * [Breadth-First Search](#bfs)
* [Forest](#forest)
* [Connected Components](#cc)
    * [CC](#ccc)	
    * [WCC](#wcc)
    * [SCC](#scc)
    * [Giant Component](#giantcc)
    * [Distribution](#ccdist)
* [Degree](#degree)
    * [Get Degree](#degreeGet)
    * [Sort Degree](#degreeSort)
    * [Write Degree](#degreeWrite)
* [Shortest Paths](#sp)
* [Connectivity](#connectivity)
* [Dominators](#dominators)


## <a name="digraph"></a> Digraph ([documentation](#documentation))

### <a name="digraphDS"></a> Data Structure ([documentation](#documentation))

A digraph is stored in a `digraph_type` data structure. A `digraph_type` contains: the number of nodes (`nodes`), the number of edges (`edges`), the list of IN edges (`in_neighbours`), the list of OUT edges (`out_neighbours`), the number of deleted nodes (`del_nodes_num`), the list of deleted nodes (`del_nodes_list`), the number of nodes in the original digraph (`org_nodes_num`), the original nodes' ID (`labels`). A valid nodes ID is in the range 1 to n. 


```C
typedef struct digraph {
	vertex_type	*in_idx;	        /* indices of neighbours in <in_neighbours>, size n+2 */
	vertex_type	*in_neighbours;	    /* list of IN neighbours (edges), size m+1  */
	vertex_type *out_idx;           /* indices of neighbours in <out_neighbours>, size n+2 */
	vertex_type *out_neighbours;    /* list of OUT neighbours (edges), size m+1  */
    vertex_type *labels;            /* nodes labels, size n+1, labels[i] is the original ID of nodes i */
    vertex_type nodes;              /* number of nodes, |V|=n */
    vertex_type edges;              /* number of edges, |E|=m */
	vertex_type org_nodes_num;	    /* number of nodes in the original graph */	
	vertex_type	del_nodes_num;	    /* number of deleted nodes */
	vertex_type	*del_nodes_list;    /* list of deleted nodes with their original label */
} digraph_type;

```

 Given a digraph G=(V,E) with |V|=n and |E|=m. Arrays `in_idx` and `out_idx` have size n+2, while the arrays `in_neighbours` and `out_neighbours` have size m+1. The positions with index 0 of these 4 arrays are never used to store data but only for implementation purposes. In particular `in_idx[0]` and `out_idx[0]` are always set to 0, while `in_idx[1]` and `out_idx[1]` are always set to 1.
Arrays `in_idx` and `out_idx` store the start and end index of the IN and OUT neighbours of each node v in V. Given a node v in V, `in_idx[v]` is the start valid index of its IN neighbours (stored in `in_neighbours`) and `in_idx[v+1]-1` is the end valid index of its IN neighbours.

> The number of OUT neighbours of node v is: out\_idx[v+1] - out\_idx[v]
> 
> The list of OUT neighbours of node v is: out\_neighbours[out\_idx[v]], out\_neighbours[out\_idx[v]+1], ..., out\_neighbours[out\_idx[v+1]-1]
> 

The array `labels` has size n+1 and contains the original ID of the nodes, usually v = labels[v]. In case of nodes removal it could happen that v != labels[v] for some v, because nodes are relabelled after removal, starting from 1 to n where n is the current number of nodes in the graph.


##### Release Memory
The function `free_digraph` is used to free the memory allocated for a digraph.

```sh
void free_digraph(digraph_type *g);
``` 

### <a name="digraphGet"></a> Read Digraph ([documentation](#documentation))

The function `get_digraph_from_file` is used to read a digraph from a file. A digraph is represented by a `digraph_type` variable.  

```sh
digraph_type *get_digraph_from_file(const char *file_name, unsigned char file_format, char check)
```

The first parameter of the funciton is the name of the file containing the graph, the second parameter is the format of the file *(currently the only supported value is EDGE_LIST, see [Edge List](#edgelist))*.

```sh
digraph_type *g;
g = get_digraph_from_file( file_name, EDGE_LIST);
``` 



### <a name="digraphWrite"></a> Write Digraph ([documentation](#documentation))

##### Standard Output
 
The function `print_digraph` is used to write a digraph to the standard output. The parameter `label` can be used to show the original node id or the current node id: if label is set to `1` the function prints the node original id, if label is set to `0` the funciton prints the current node id.

```sh
void print_digraph(digraph_type *g, char label)
```
>Expected output:
>
>```sh
>DIGRAPH:
>	Nodes: 4 Edges: 3
>GRAPH EDGES:
>	1->{2}  2->{3}  3->{4}  
>LABELS:
>	1=>1  2=>2  3=>3  4=>4  
>```

##### File
 
The function `write_digraph_edgelist` is used to write a digraph to a file in edge list format.

```sh
void write_digraph_edgelist(digraph_type *g, const char* file_name)
```

The function `write_digraph_dot` is used to write a digraph to a file in [dot](https://en.wikipedia.org/wiki/DOT_%28graph_description_language%29) format.

```sh
void write_digraph_dot(digraph_type *g, const char* file_name)
```

The function `write_digraph_labels` is used to write the labels of the graph in a CSV format, the columns are: NodeID, OriginalID.

```sh
write_digraph_labels (digraph_type *g, const char *file_name)
```

### <a name="digraphOp"></a> Digraph Operations ([documentation](#documentation))

#### Nodes Removal

The function `delete_nodes_from_digraph` is used to remove a set of nodes from a digraph.  It returns a copy of the original digraph without the nodes. The `vts` array contains the set of nodes to remove, `n_vts` is the number of nodes to remove, i.e. the size of `vts`

```sh
digraph_type *delete_nodes_from_digraph(digraph_type *g, vertex_type *vts, vertex_type n_vts)
```

#### Reverse Digraph

The function `get_reverse_digraph` returns a pointer to the reverse digraph of `g`. The reverse digraph of a digraph `g` is the digraph `g` with the direction of all its edges reverted.

Let G = (V,E) be a strongly connected graph. The reverse digraph of G, denoted by G<sup>R</sup> = (V, E<sup>R</sup>), is the digraph that results from G by reversing the direction of all edges.

```
digraph_type *get_reverse_digraph ( digraph_type *g )
```
#### Subgraph

The function `get_subgraph_digraph` returns a copy of a digraph `g` containing only the nodes listed in `nodes`.

```C
digraph_type *get_subgraph_digraph(digraph_type *g, vertex_type *nodes, vertex_type n_nodes)
```

#### Get Graph from Digraph

The function `get_graph_from_digraph` returns an undirected copy of a digraph `g`.

```C
graph_type * get_graph_from_digraph(digraph_type *g)
```

## <a name="graph"></a> Graph ([documentation](#documentation))

### <a name="graphDS"></a> Data Structure ([documentation](#documentation))

A (undirected) graph is stored in a `graph_type` data structure. A `graph_type` contains: the number of nodes (`nodes`), the number of edges (`edges`), the list of edges (`neighbours`), the number of deleted nodes (`del_nodes_num`), the list of deleted nodes (`del_nodes_list`), the number of nodes in the original digraph (`org_nodes_num`), the original nodes' ID (`labels`). A valid nodes ID is in the range 1 to n. 


```C
typedef struct graph {
	vertex_type	*idx;				/* indices of neighbours size n+2 */
	vertex_type	*neighbours;		/* list of neighbours (edges), size m+1  */
    vertex_type *labels;			/* nodes labels, size n+1, labels[i] is the original ID of nodes i */
    vertex_type nodes;				/* number of nodes, |V|=n */
    vertex_type edges;				/* number of edges, |E|=m */
	vertex_type org_nodes_num;		/* number of nodes in the original graph */	
	vertex_type	del_nodes_num;		/* number of deleted nodes */
	vertex_type	*del_nodes_list;    /* list of deleted nodes with their original label */
} graph_type
```

Given a graph G=(V,E) with |V|=n and |E|=m. Array `idx` has size n+2, while the array `neighbours` has size 2m+1. The positions with index 0 of these 2 arrays are never used to store data but only for implementation purposes. In particular `idx[0]` is always set to 0, while `idx[1]` is always set to 1.
Array `idx` stores the start and end index of the neighbours of each node v in V. Given a node v in V, `idx[v]` is the start valid index of its neighbours (stored in `neighbours`) and `idx[v+1]-1` is the end valid index of its neighbours.

> The number of neighbours of node v is: idx[v+1] - idx[v]
> 
> The list of neighbours of node v is: neighbours[idx[v]], neighbours[idx[v]+1], ..., neighbours[idx[v+1]-1]
> 

The array `labels` has size n+1 and contains the original ID of the nodes, usually v = labels[v]. In case of nodes removal it could happen that v != labels[v] for some v, because nodes are relabelled after removal, starting from 1 to n where n is the current number of nodes in the graph.

##### Release Memory
The function `free_graph` is used to free the memory allocated for a graph.

```sh
void free_graph(graph_type *g)
``` 

###  <a name="graphGet"></a> Read Graph ([documentation](#documentation))

The function `get_graph_from_file` is used to read a graph from a file. A graph is represented by a `graph_type` variable.  

```sh
graph_type *get_graph_from_file(const char *file_name, unsigned char file_format, char check)
```

The first parameter of the funciton is the name of the file containing the graph, the second parameter is the format of the file *(currently the only supported value is EDGE_LIST, see [Edge List](#edgelist))*.

```sh
graph_type *g;
g = get_graph_from_file( file_name, EDGE_LIST, 0 )
``` 

###  <a name="graphWrite"></a> Write Graph ([documentation](#documentation))

##### Standard Output
 
The function `print_graph` is used to write a graph to the standard output. The parameter `label` can be used to show the original node id or the current node id: if label is set to `1` the function prints the node original id, if label is set to `0` the funciton prints the current node id.

```sh
void print_graph(graph_type *g, char label)
```
>Expected output:
>
>```sh
>GRAPH: Nodes: 3  Edges: 3
>	Edges: 
>	(1,2) (1,3) (2,3) 
>	Connections: 
>	  1->{2,3}  2->{1,3}  3->{1,2}
>	Labels: 
>
>	  1=>1  2=>2  3=>3
>```

##### File
 
The function `write_graph_edgelist` is used to write a graph to a file in edge list format.

```sh
void write_graph_edgelist(graph_type *g, const char* file_name)
```

The function `write_graph_dot` is used to write a digraph to a file in [dot](https://en.wikipedia.org/wiki/DOT_%28graph_description_language%29) format.

```sh
void write_graph_dot(graph_type *g, const char* file_name)
```

The function `write_graph_labels` is used to write the labels of the graph in a CSV format, the columns are: NodeID, OriginalID.

```sh
write_graph_labels (graph_type *g, const char *file_name)
```

###  <a name="graphOp"></a> Graph Operations ([documentation](#documentation))

#### Nodes Removal

The function `delete_nodes_from_graph` is used to remove a set of nodes from a graph.  It returns a copy of the original graph without the nodes. The `vts` array contains the set of nodes to remove, `n_vts` is the number of nodes to remove, i.e. the size of `vts`

```C
graph_type *delete_nodes_from_graph(graph_type *g, vertex_type *vts, vertex_type n_vts)
```
  
#### Subgraph

The function `get_subgraph_graph` returns a copy of a graph `g` containing only the nodes listed in `nodes`.

```C
graph_type *get_subgraph_graph(graph_type *g, vertex_type *nodes, vertex_type n_nodes)
```

## <a name="file"></a> File Format ([documentation](#documentation))

### <a name="edgelist"></a> Edge List ([documentation](#documentation))

A file in edge list format contains the following information:

 - A line starting with '#' specifying the number of Nodes and Edges of the graph *(any other line starting with '#' is interpreted as a comment and skipped*).
-  A line for each edge of the graph containing a pair of numbers, that are the nodes connected by the edge.

The numbers used to identify nodes must be greater then 0, i.e. the first node id is 1. Below there is an example of graph in the edges list format.

```sh
# Nodes: 4 Edges: 3
1 2
2 3
3 4
```

## <a name="traversal"></a> Graph Traversal ([documentation](#documentation))

### <a name="dfs"></a> Depth-First Search (DFS) ([documentation](#documentation))

The function `get_digraph_dfs` is used to traverse a digraph in DFS order. The function return a [forest](https://en.wikipedia.org/wiki/Tree_%28graph_theory%29#Forest) represented by a `forest_type` variable. The `root` parameter specify the node used to start the digraph traversal.

```sh
forest_type *get_digraph_dfs(digraph_type *g, vertex_type root)
```

The function `get_graph_dfs` is used to traverse a graph in DFS order. The function return a [forest](https://en.wikipedia.org/wiki/Tree_%28graph_theory%29#Forest) represented by a `forest_type` variable. The `root` parameter specify the node used to start the graph traversal.

```sh
forest_type *get_graph_dfs(graph_type *g, vertex_type root)
```

### <a name="bfs"></a> Breadth-First Search (BFS) ([documentation](#documentation))

The function `get_digraph_bfs` is used to traverse a digraph in BFS order. The function return a [forest](https://en.wikipedia.org/wiki/Tree_%28graph_theory%29#Forest) represented by a `forest_type` variable. The `root` parameter specify the node used to start the digraph traversal.

```sh
forest_type *get_digraph_bfs(digraph_type *g, vertex_type root)
```
The function `get_graph_bfs` is used to traverse a graph in BFS order. The function return a [forest](https://en.wikipedia.org/wiki/Tree_%28graph_theory%29#Forest) represented by a `forest_type` variable. The `root` parameter specify the node used to start the graph traversal.

```sh
forest_type *get_graph_bfs(graph_type *g, vertex_type root)
```
## <a name="forest"></a> Forest ([documentation](#documentation))

### <a name="forestDS"></a>  Data Structure ([documentation](#documentation))

A forest is stored in the `forest_type` data structure. It contains the number of nodes (`nodes`) and edges (`edges`) of the forest, the visit order of the nodes (`order` and `vertex`) and the trees (`parent`) of the forest. The `vertex` array contains the nodes in visit order: `vertex[i] = v` means that v was the i-th node to be visited, where v is the noded ID. The `order` array contains the order in which nodes have been visited: `order[v] = i` means that v was the i-th to be visited. The `parent` array contains the parent node of each node in the forest, `parent[v] = w` means that w is the parent of v in T and `parent[v] = 0` only if v is the root of T.

```sh
typedef struct forest{
	vertex_type *order;
	vertex_type *vertex;
	vertex_type *parent;
	vertex_type nodes;
	vertex_type edges;
} forest_type
```

##### Release Memory
The function `free_forest` is used to free the memory allocated for a forest.

```sh
void free_forest(forest_type *f)
``` 

### <a name="forestWrite"></a> Write Forest ([documentation](#documentation))

##### Standard Output

The function `print_forest` is used to write a forest to the standard output.

```sh
void print_forest(vertex_type *parent, vertex_type n_nodes, vertex_type n_edges)
```

##### File

The function `write_forest_edge_list` is used to write a forest to a file in edge list format.

```sh
void write_forest_edge_list(vertex_type *parent, vertex_type n_nodes, vertex_type n_edges, const char* file_name)
```

The function `write_forest_dot` is used to write a forest to a file in [dot](https://en.wikipedia.org/wiki/DOT_%28graph_description_language%29) format.

```sh
void write_forest_dot(vertex_type *parent, vertex_type n_nodes, const char* file_name, char type)
```


## <a name="cc"></a> Connected Components ([documentation](#documentation))


##### <a name="ccDS"></a>  Data Structure 

A set of CCs is stored in a `cc_type` data structure. A `cc_type` contains: the number of CC `number`, the list of vertices ordered by CC `vertex` and the indexes of the vertices in the same CC `idx`.

```sh
typedef struct cc {
    vertex_type *vertex;      /* List of vertices ordered by CC, size |V|=n */
    vertex_type *idx;         /* Indexes of the vertices in the same CC, size <number>+1 */
    vertex_type number;       /* Number of CC */
} cc_type
```

The vertices of CC i are in the vector `vertex` starting from index `idx[i]` to `idx[i+1]-1`. 

> The list of vertices of CC i is: vertex[idx[i]], ..., vertex[idx[i+1]-1]
>
> The size of CC i is: idx[i+1]-idx[i]

##### Release Memory 
The function `free_cc` is used to free the memory allocated for a CC.

```sh
void free_cc (cc_type *cc)
``` 

### <a name="ccc"></a>  Connected Components (CCs) ([documentation](#documentation))
The function `get_graph_cc` returns the CCs of a graph `g`. The decomposition in CC is represented by a `cc_type` variable.

```sh
cc_type *get_graph_cc(graph_type *g)
``` 

### <a name="wcc"></a> Weakly Connected Components (WCCs) ([documentation](#documentation))

#####  Data Structure 
A set of WCCs is stored in a `cc_type` data structure. A `cc_type` contains: the number of WCC `number`, the list of vertices ordered by WCC `vertex` and the indexes of the vertices in the same WCC `idx`.

#####  Get WCC 
The function `get_digraph_wcc` returns the WCCs of a digraph `g`. 

```C
cc_type *get_digraph_wcc(digraph_type *g)
``` 

### <a name="scc"></a> Strongly Connected Components (SCCs) ([documentation](#documentation))

#####  Data Structure 
A set of [SCCs](https://en.wikipedia.org/wiki/Strongly_connected_component) is stored in a `cc_type` data structure. A `cc_type` contains: the number of SCC `number`, the list of vertices ordered by SCC `vertex` and the indexes of the vertices in the same SCC `idx`.


#####  Get SCC 
The function `get_digraph_scc` returns the SCCs of a digraph `g`.  The funciton uses the Tarjan algorithm [3].

```sh
cc_type *get_digraph_scc(digraph_type *g)
``` 

### <a name="giantcc"></a> Giant Component ([documentation](#documentation))

The function `get_giant_cc_vertices` is used to get the list of vertices in the giant component (CC/SCC/WCC). Functions `get_subgraph_graph` and `get_subgraph_digraph` can be used to get the giant component subgraph from the original graph. 

```sh
vertex_type *get_giant_cc_vertices (cc_type *cc, uint64_t *size)
``` 

### Write CC/WCC/SCC

#### Standard Output

The function `print_cc` is used to write the CC decomposition of a graph to the standard output.

```sh
void print_cc(cc_type *cc)
```
>Expected output:
>
>```
>CC 0: NODE_1 ... NODE_N - size: X
>CC 1: NODE_1 ... NODE_N - size: Y
>...
>CC N: NODE_1 ... NODE_N - size: Z
>```

#### File

The function `write_cc` is used to write the CC decomposition of a graph to a file in [json](https://en.wikipedia.org/wiki/JSON) format. 

```sh
void write_cc (cc_type *cc, const char *file_name)
```

> Expected output:
>
> ```
>{ 
>	CC_0: { 
>		"nodes": [ NODE_1, ... , NODE_N],
>		"size" : NUMBER
>	},
>	...
>	CC_N: { 
>		"nodes": [ NODE_1, ... , NODE_N],
>		"size" : NUMBER
>	}
>}
> ```


### <a name="ccdist"></a> Connected Components Distribution (CC/WCC/SCC) ([documentation](#documentation))

The function `get_cc_distribution` is used to compute the distribution of the sizes of the CCs/WCCs/SCC of a graph/digraph.

```sh
uint64_t * get_cc_distribution (cc_type *cc, uint64_t *size)
```

### Write Distribution

#### Standard Output

The function `print_cc_distribution` is used to print on the standard output the distribution of the sizes of the CCs of a graph.

```sh
void print_cc_distribution (uint64_t *cc_dist, uint64_t size)
```
	
> Expected output for print and write:
>
> ```
> CC_1_SIZE, CC_2_SIZE, ..., CC_N_SIZE
> ```

#### File

The function `write_cc_distribution` is used to write the distribution of the sizes of the CCs of a graph.

```sh
void write_cc_distribution (uint64_t *cc_dist, uint64_t size, const char *file_name)
```




## <a name="degree"></a> Degree ([documentation](#documentation))

### <a name="degreeGet"></a>  Get Degree ([documentation](#documentation))
##### IN/OUT Degree 
The function `get_inoutdegree` is used to compute the IN and OUT degree of the nodes of a digraph. The values are stored in an array of  `didegree_type` variables.

```sh
didegree_type *get_inoutdegree(digraph_type *g)
``` 

##### Degree 
The function `get_degree` is used to compute the degree of the nodes of a graph. The values are stored in an array of  `degree_type` variables.

```sh
degree_type *get_degree(graph_type *g)
``` 
### <a name="degreeSort"></a> Sort Degree ([documentation](#documentation))
##### IN/OUT Degree
The function `sort_inoutdegree` is used to sort the IN and OUT degree values stored in an array  of  `didegree_type` variables.

```sh  
void sort_inoutdegree(didegree_type *d, vertex_type n_nodes, char field, char order)
``` 

- Values can be sorted for *IN* or *OUT* degree using the third parameter of the function `char field` specifying: ***INDEG*** or ***OUTDEG***.
- Values can be sorted in *ascendent* or *descendent* order using the fourth parameter of the function `char order` specifying: ***ASC*** or ***DESC***.

##### Degree
The function `sort_degree` is used to sort the degree values stored in an array  of  `degree_type` variables.

```sh  
void sort_degree(degree_type *d, vertex_type n_nodes, char order)
``` 

- Values can be sorted in *ascendent* or *descendent* order using the fourth parameter of the function `char order` specifying: ***ASC*** or ***DESC***.


### <a name="degreeWrite"></a> Write Degree ([documentation](#documentation))

##### IN/OUT Degree
The function `write_inoutdegree` is used to write IN and OUT degree values to a file. The third parameter of the function is used to specify how many values (rows) to write.

```sh
void write_inoutdegree(didegree_type *d, vertex_type n_nodes, vertex_type top_n, const char *file_name)
```

The function `write_inoutdegree` is used to print on standard output IN and OUT degree values.

```sh
void print_inoutdegree(didegree_type *d, vertex_type n_nodes, vertex_type top_n)
```

> Values are written/printed in CSV format as follows:
> 
>```sh
>Node,  InDeg,  OutDeg
>2,  3,  3
>5,  3,  3
>...
>```

##### Degree
The function `write_degree` is used to write degree values to a file. The third parameter of the function is used to specify how many values (rows) to write.

```sh
void write_degree(degree_type *d, vertex_type n_nodes, vertex_type top_n, const char *file_name)
```
The function `print_degree` is used to print on standard output degree values.

```sh
void print_degree(degree_type *d, vertex_type n_nodes, vertex_type top_n)
```
> Values are written/printed in CSV format as follows:
> 
>```sh
>Node,  Deg
>2,  4
>5,  3
>...
>```


## <a name="connectivity"></a> Connectivity ([documentation](#documentation))

### Connectivity
The function `get_pairs_connectivity` is used to compute the nodes pairs connectivity. 

```sh
uint64_t get_pairs_connectivity(scc_type *scc)
```
Nodes pairs connectivity is defined as follows:


>  Let G be a directed graph, let C<sub>1</sub>,C<sub>2</sub>,...,C<sub>l</sub> be its strongly connected components.  We define the connectivity value of G as:
>  
> f(G) = &sum;<sup>l</sup><sub>i=1</sub></i>  \binom{C<sub>i</sub>}{2}
> 
> f(G) equals the number of vertex pairs in G that are strongly
> connected  (i.e., pairwise strong connectivity value)

### Critical Nodes

The function `get_critical_nodes_bruteforce` is used to compute the critical nodes sets of a digraph.   The function checks which set of nodes of size `k` is critical in the digraph `g` and returns a list of them whose maximum size is `max_set`. A set of nodes is critical if it minimaizes the connectivity of `g`. 

The function returns: 1) an array of critical sets or NULL if `max_set` is set to 0; 2) the initial connectivity of the graph `init_conn`; 3) the final connectivity `end_conn` after removing a critical set of size `k`; 4) the total number of critical sets `tot_sets` of size `k`. 

```sh
vertex_type ** get_critical_nodes_bruteforce ( digraph_type *g, vertex_type k, uint64_t max_set, uint64_t *init_conn, uint64_t *end_conn, uint64_t *tot_sets )
```

## <a name="dominators"></a> Dominators ([documentation](#documentation))

### Flow Graphs and Dominators Trees
A *flow graph* is a directed graph with a distinguished start vertex s such that every vertex is reachable from s. Let G<sub>s</sub> be a *flow graph* with start vertex s. A vertex u is a dominator of a vertex v (u dominates v) if every path from s to v in G<sub>s</sub> contains u. The dominator relation is reflexive and transitive. Its transitive reduction is a rooted tree, the dominator tree D: u dominates v if and only if u is an ancestor of v in D.

The functions both return the dominator tree represented as a parents array (`vertex_type *dom`), *i.e.* `dom[v]=w` means that `v` dominates `w`, and there is an edge `v->w` in the dominator tree.

### Computing Dominators
The function `get_dominators_pm` returns the dominator tree of digraph `g` using the Purdom-Moore [1] algorithm.

**Precondition:** the digraph `g` whose start vertex is `root` is a flow graph.

```sh
vertex_type * get_dominators_pm ( digraph_type *g, vertex_type root )
```
The function `get_dominators_lt` returns the dominator tree of digraph `g` using the Lengauer-Tarjan [2] simple algorithm whose time complexity is *O( m log(n) )*

**Precondition:** the digraph `g` whose start vertex is `root` is a flow graph.

```sh
vertex_type * get_dominators_lt(digraph_type *g, forest_type *dfs )
```

## References

[1] P. W. Purdom, Jr. and E. F. Moore. "Immediate predominators in a directed graph". Communications of the ACM, 15(8):777–778, 1972.

[2] T. Lengauer and R. E. Tarjan. "A fast algorithm for finding dominators in a flowgraph". ACM Transactions on Programming Languages and Systems, 1(1):121–41, 1979.

[3] R. Tarjan. "Depth-first search and linear graph algorithms". 12th Annual Symposium on Switching and Automata Theory (swat 1971), East Lansing, MI, USA, 1971, pp. 114-121, doi: 10.1109/SWAT.1971.10.