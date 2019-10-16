#ifndef SPARSE_GRAPH_H
#define SPARSE_GRAPH_H

/* Prototypes for dynamically allocated sparse graph.  For use in
 * parallel, the assembler only pays attention to edges which involve
 * nodes in some [lo, hi] range.  All true node numbers should be
 * non-negative; negative node indices will be ignored.
 */

typedef struct sparse_graph_t sparse_graph_t;


/* Create a new dynamic sparse graph object
 *   lo - Index of first node for this processor
 *   hi - Index of the last node for this processor
 */
sparse_graph_t* make_sparse_graph(int lo, int hi);


/* Free an existing sparse graph object 
 *   graph - Object to be destroyed
 */
void destroy_sparse_graph(sparse_graph_t* graph);


/* Add an edge to a sparse graph (no self-loops, no non-positive indices) 
 *   graph - Graph object to be modified
 *   src   - source node for edge
 *   dest  - destination node for edge
 */
void sparse_graph_add(sparse_graph_t* graph, int src, int dest);


/* Add an element contribution to a sparse graph 
 *   graph          - Graph object to be modified
 *   ix(nen,numelt) - element connectivity array (column-major)
 */
void sparse_graph_assemble(sparse_graph_t* graph, int* ix, 
                           int numelt, int nen);

/* Convert the linked-list adjacency representation to a CSR index structure. 
 *   graph - Graph object to be converted
 *   pia   - Newly-allocated row pointer array for CSR format
 *   pja   - Newly-allocated column index array for CSR format
 */
void graph_to_csr(sparse_graph_t* graph, int** pia, int** pja);


/* Given element connectivity, generate a sparse adjacency matrix 
 *   ix(numelt,nen) - Element connectivity array
 *   lo             - Index of first node for this processor
 *   hi             - Index of last node for this processor
 *   pia            - Newly-allocated row pointer array for CSR format
 *   pja            - Newly-allocated column index array for CSR format
 */
void make_csr(int* ix, int numelt, int lo, int hi, int nen,
              int** pia, int** pja);


#endif /* SPARSE_GRAPH_H */
