#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <mpi.h>
#include <parmetis.h>
#include "sparse-graph.h"


/* --- Arena list ---
 *
 * We initially use a dynamic data structure to contain the node graph:
 * we keep an adjacency list for each node, stored as a linked list.
 * All the nodes are the same size, and we never free a node unless we're
 * freeing the entire list, so we use an arena allocator system.
 */


/* Number of nodes to be allocated in one go */
#define NODE_POOL_SIZE 4000


typedef struct node_t {
    int data;
    struct node_t* next;
} node_t;


typedef struct node_pool_t {
    node_t nodes[NODE_POOL_SIZE];
    struct node_pool_t* next;
    int free;
} node_pool_t;


/*
 * Fast allocator for list nodes
 */
static node_t* node_alloc(node_pool_t** poolp)
{
    node_pool_t* pool = *poolp;
    if (!pool || pool->free == 0) {
        pool = (node_pool_t*) malloc(sizeof(node_pool_t));
        pool->next = NULL;
        pool->free = NODE_POOL_SIZE;
        *poolp = pool;
    } 
    return &(pool->nodes[--(pool->free)]);
}


/*
 * Reset the node pool
 */
static void free_all_nodes(node_pool_t** poolp)
{
    node_pool_t* pool = *poolp;
    while (pool) {
        node_pool_t* dead = pool;
        pool = pool->next;
        free(dead);
    }
    *poolp = NULL;
}


/*
 * Add a node to an ordered linked list
 */
static void add_node(node_pool_t** poolp, node_t** node, int data)
{
    node_t* new_node;
    while (*node != NULL && (*node)->data < data)
        node = &((*node)->next);
    if (*node && (*node)->data == data)
        return;
    new_node = node_alloc(poolp);
    new_node->data = data;
    new_node->next = *node;
    *node = new_node;
}


/*
 * Print list (debugging utility)
 */
void print_node_list(node_t* node)
{
    while (node) {
        printf(" %d", node->data);
        node = node->next;
    }
}


/* --- Graph representation ---
 *
 * We represent a horizontal slice of the adjacency matrix with each
 * processor, using linked lists for the adjacency lists.
 */


struct sparse_graph_t {
    node_pool_t* pool;    /* Pool for node memory         */
    node_t** adj;         /* Array of adjacency lists     */
    int lo;               /* Index of first node in slice */
    int hi;               /* Index of last node in slice  */
};


/*
 * Create a new dynamic sparse graph object
 */
sparse_graph_t* make_sparse_graph(int lo, int hi)
{
    int mlocal = (hi-lo)+1;
    sparse_graph_t* result = (sparse_graph_t*) malloc(sizeof(sparse_graph_t));
    result->pool = NULL;
    result->adj  = (node_t**) malloc(mlocal * sizeof(node_t*));
    result->lo   = lo;
    result->hi   = hi;
    memset(result->adj, 0, mlocal * sizeof(node_t*));
    return result;
}


/*
 * Free an existing sparse graph object
 */
void destroy_sparse_graph(sparse_graph_t* graph)
{
    free(graph->adj);
    free_all_nodes(&(graph->pool));
    free(graph);
}


/*
 * Add an edge to a sparse graph (no self-loops, no negative indices)
 */
void sparse_graph_add(sparse_graph_t* graph, int src, int dest)
{
    int lo = graph->lo;
    int hi = graph->hi;
    if (src < lo || src > hi || src == dest || dest < 0) 
        return;
    add_node(&(graph->pool), &(graph->adj[src-lo]), dest);
}


/*
 * Add an element contribution to a sparse graph
 */
void sparse_graph_assemble(sparse_graph_t* graph, int* ix, int numelt, int nen)
{
    int e, i, j;
    for (e = 0; e < numelt; ++e) {
        for (i = 0; i < nen; ++i) {
            for (j = 0; j < nen; ++j)
                sparse_graph_add(graph, ix[i], ix[j]);
        }
        ix += nen;
    }
}


/*
 * Convert the linked-list adjacency representation to a CSR index structure.
 */
void graph_to_csr(sparse_graph_t* graph, int** pia, int** pja)
{
    int i, p;
    node_t** adj = graph->adj;
    node_t* n;
    int* ia;
    int* ja;
    int mlocal = (graph->hi-graph->lo)+1;

    /* Build the CSR pointer list */
    ia = (int*) malloc((mlocal+1)*sizeof(int));
    p = 0;
    ia[0] = 0;
    for (i = 0; i < mlocal; ++i) {
        for (n = adj[i]; n != NULL; n = n->next)
            ++p;
        ia[i+1] = p;
    }

    /* Build the column index list */
    ja = (int*) malloc(ia[mlocal]*sizeof(int));
    p = 0;
    for (i = 0; i < mlocal; ++i) {
        for (n = adj[i]; n != NULL; n = n->next)
            ja[p++] = n->data;
    }

    *pia = ia;
    *pja = ja;
}


/*
 * Given element connectivity, generate a sparse adjacency matrix
 */
void make_csr(int* ix, int numelt, int lo, int hi, int nen,
              int** ia, int** ja)
{
    sparse_graph_t* graph = make_sparse_graph(lo, hi);
    sparse_graph_assemble(graph, ix, numelt, nen);
    graph_to_csr(graph, ia, ja);
    destroy_sparse_graph(graph);
}
