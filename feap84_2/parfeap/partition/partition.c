/*
c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    06/05/2012
c       1. Update call to ParMetis_V3_PartKway to           10/01/2013
c          work with parMETIS 4.0.2   
c-----[--.----+----.----+----.-----------------------------------------]
c     Purpose:  Standalone partioning program using parMETIS

c     Usage: mpirun -n npm partition feap_doms Ifileflat
c            npm -- number of processors to use when running parMETIS
c      feap_doms -- number of domains for which you are partitioning
c      Ifileflat -- flat input file that you wish to partition
c-----[--.----+----.----+----.-----------------------------------------]
*/



#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <mpi.h>
#include <parmetis.h>
#include "sparse-graph.h"

#define VERBOSE  0
#define fprintf0 if (id == 0) fprintf
#define dprintf0 if (id == 0 && VERBOSE) fprintf


/* --- Streaming operations ---
 *
 * MPI doesn't have a built-in stream primitive, let alone buffered
 * streams.  Here's a simple stream for integer data.
 */


typedef struct mpi_stream_t {
    int src;     /* Source node              */
    int count;   /* Remaining data           */
    int nchunk;  /* Size of buffer blocks    */
    int nbuf;    /* Number of buffered items */
    int* buf;    /* Buffer                   */
} mpi_stream_t;


/*
 * Stream a large block of data to dest in blocks of nchunk.
 */
void mpi_send_block(int* buf, int count, int dest, int nchunk)
{
    while (count > 0) {
        if (count < nchunk)
            nchunk = count;
        MPI_Send(buf, nchunk, MPI_INT, dest, 0, MPI_COMM_WORLD);
        count -= nchunk;
        buf += nchunk;
    }
}


/*
 * Receive a stream of data from the src in blocks of nchunk.
 */
mpi_stream_t* mpi_open_stream(int src, int count, int nchunk)
{
    mpi_stream_t* stream = (mpi_stream_t*) malloc(sizeof(mpi_stream_t));
    stream->src = src;
    stream->count = count;
    stream->nchunk = nchunk;
    stream->nbuf = 0;
    stream->buf = (int*) malloc(nchunk * sizeof(int));
    return stream;
}


/*
 * Clean up after streaming data from another processor
 */
void mpi_close_stream(mpi_stream_t* stream)
{
    free(stream->buf);
    free(stream);
}


/*
 * Fetch data from a remote processor stream
 */
int* mpi_fetch_stream(mpi_stream_t* stream, int nfetch)
{
    int* data;
    MPI_Status status;
    if (stream->nbuf == 0) {
        if (stream->count < stream->nchunk)
            stream->nchunk = stream->count;
        MPI_Recv(stream->buf, stream->nchunk, MPI_INT, stream->src, 
                 0, MPI_COMM_WORLD, &status);
        stream->count -= stream->nchunk;
        stream->nbuf  =  stream->nchunk;
    }
    data = stream->buf + (stream->nchunk - stream->nbuf);
    stream->nbuf -= nfetch;
    return data;
}


/* --- Collect and print ---
 *
 * Print the distributed connectivity information at processor 0
 */


void print_part(FILE* fp, int id, int ntasks, int* vtxdist, int* part)
{
    const int block_size = 1024;
    if (id != 0) {
        int mlocal = vtxdist[id+1]-vtxdist[id];
        mpi_send_block(part, mlocal, 0, block_size);
    } else {
        int p, i;
        for (i = vtxdist[id]; i < vtxdist[id+1]; ++i)
            fprintf(fp, "%d\n", part[i]+1);
        for (p = 1; p < ntasks; ++p) {
            int pmlocal = vtxdist[p+1]-vtxdist[p];
            mpi_stream_t* stream = mpi_open_stream(p, pmlocal, block_size);
            for (i = vtxdist[p]; i < vtxdist[p+1]; ++i)
                fprintf(fp, "%d\n", *mpi_fetch_stream(stream, 1)+1);
        }
    }
}


void print_csr(FILE* fp, int id, int ntasks, int* vtxdist, int* ia, int* ja)
{
    const int block_size = 1024;
    if (id != 0) {

        int mlocal = vtxdist[id+1]-vtxdist[id];
        MPI_Send(ia, mlocal+1, MPI_INT, 0, 0, MPI_COMM_WORLD);
        MPI_Send(&(ia[mlocal]), 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
        mpi_send_block(ja, ia[mlocal], 0, block_size);

    } else {
        int mlocal = vtxdist[1]-vtxdist[0];
        int i, p, start_row, nentries;
        mpi_stream_t* stream;
        MPI_Status status;

        /* -- Print pointer array -- */

        start_row = 1;
        for (i = 0; i < mlocal; ++i)
            fprintf(fp, "%d\n", ia[i]+start_row);

        start_row += ia[mlocal];
        for (p = 1; p < ntasks; ++p) {
            int pmlocal = vtxdist[p+1]-vtxdist[p];
            int* pia = (int*) malloc((pmlocal+1) * sizeof(int));
            MPI_Recv(pia, pmlocal+1, MPI_INT, p, 0, MPI_COMM_WORLD, &status);
            for (i = 0; i < pmlocal; ++i)
                fprintf(fp, "%d\n", pia[i]+start_row);
            start_row += pia[pmlocal];
            free(pia);
        }

	fprintf(fp, "%d\n", start_row);

        /* -- Print connectivity data -- */

        nentries = ia[mlocal];
        for (i = 0; i < nentries; ++i)
            fprintf(fp, "%d\n", ja[i]+1);

        for (p = 1; p < ntasks; ++p) {
            MPI_Recv(&nentries, 1, MPI_INT, p, 0, MPI_COMM_WORLD, &status);
            stream = mpi_open_stream(p, nentries, block_size);
            for (i = 0; i < nentries; ++i)
                fprintf(fp, "%d\n", *mpi_fetch_stream(stream, 1)+1);
            mpi_close_stream(stream);
        }
    }
}

 
/* --- Partition and distribute map ---
 */

/*
 * Call ParMETIS to partition the unweighted nodal graph
 */
int* partition(int id, int* vtxdist, int* ia, int* ja, int nparts)
{
    int wgtflag = 0; /* No weights        */
    int numflag = 0; /* C-style           */
    int ncon = 1;    /* Number of weights */
    int options[3] = {1, 0, 0};
    int edgecut;
    int* vwgt     = NULL;
    int* adjwgt   = NULL;
    int i;
    real_t tpwgts[nparts];
    real_t ubvec = 1.05;

    for(i=0;i<nparts;i++) tpwgts[i]=1.0/nparts;

    MPI_Comm comm = MPI_COMM_WORLD;

    int mlocal = vtxdist[id+1]-vtxdist[id];
    int* part = (int*) malloc(mlocal * sizeof(int));

    ParMETIS_V3_PartKway(vtxdist, ia, ja, vwgt, adjwgt, &wgtflag,
                         &numflag, &ncon, &nparts, tpwgts, &ubvec,
                         options, &edgecut, part, &comm);

    return part;
}


/* --- Parallel test driver ---
 *
 * Make node zero handle the file input and broadcast elements to
 * everyone else in chunks.  Each processor builds its own slice
 * of the adjacency matrix.  I assume for the moment that all the
 * input will be Fortran style (one-based indices).
 */


/*
 * Convert a string to all uppercase
 */
char* str_upper(char* s)
{
    char* t = s;
    if (s) {
        while (*t) {
            if (*t >= 'a' && *t <= 'z')
                *t += ('A'-'a');
            ++t;
        }
    }
    return s;
}


/*
 * Read a FEAP input file up through the header info, return header params
 */
void read_feap_header(FILE* fp, int* numnp, int* numel, int* nummat,
                      int* ndm, int* ndf, int* nen)
{
    char buf[256];
    while (str_upper(fgets(buf, sizeof(buf), fp)) &&
           strncmp(buf, "FEAP", 4) != 0);
    fscanf(fp, "%d %d %d %d %d %d",
           numnp, numel, nummat, ndm, ndf, nen);
}


/*
 * Skip to the start of the element cards
 */
void read_to_feap_elts(FILE* fp)
{
    char buf[256];
    while (str_upper(fgets(buf, sizeof(buf), fp)) && 
           strncmp(buf, "ELEM", 4) != 0);
}


/*
 * Read a block of max(num_left, max_chunk) element cards into elt_data.
 * Only connectivity is stored.  Return the number of elements read.
 */
int read_feap_elts(FILE* fp, int* elt_data, int nen, 
                   int* num_left, int max_chunk)
{
    int i, j;
    int nchunk = max_chunk;
    if (*num_left < nchunk)
        nchunk = *num_left;

    for (i = 0; i < nchunk; ++i, --*num_left) {
        int eln, inc, mat;
        fscanf(fp, "%d %d %d", &eln, &inc, &mat);
        for (j = 0; j < nen; ++j) {
            fscanf(fp, "%d", &(elt_data[i*nen+j]));
            elt_data[i*nen+j]--;
        }
    }
    return nchunk;
}


/*
 * Process the FEAP header and broadcast to all processors
 */
void process_feap_header(FILE* fp, int id, int* numnp, int* numel, int* nen)
{
    int data[3];
    if (id == 0) {
        int nummat, ndm, ndf;
        read_feap_header(fp, numnp, numel, &nummat, &ndm, &ndf, nen);
        data[0] = *numnp;
        data[1] = *numel;
        data[2] = *nen;
    }
    MPI_Bcast(data, 3, MPI_INT, 0, MPI_COMM_WORLD);
    *numnp = data[0];
    *numel = data[1];
    *nen   = data[2];
}


/*
 * Process the element cards on all processors.
 */
void process_elt_cards(FILE* fp, int id, sparse_graph_t* graph,
                       int numel, int nen, int max_chunk)
{
    int* elt_data = malloc((1 + max_chunk * nen) * sizeof(int));
    int nchunk;

    if (id == 0)
        read_to_feap_elts(fp);

    do {
        if (id == 0) {
            nchunk = read_feap_elts(fp, elt_data+1, nen, &numel, max_chunk);
            elt_data[0] = nchunk;
        }
        MPI_Bcast(elt_data, 1+max_chunk*nen, MPI_INT, 0, MPI_COMM_WORLD);
        nchunk = elt_data[0];
        sparse_graph_assemble(graph, elt_data+1, nchunk, nen);

    } while (nchunk == max_chunk);
    free(elt_data);
}


/*
 * Decide on the node partitioning between processors
 */
int* build_node_dist(int ntasks, int numnp)
{
    /* Everyone decide on the range of nodal points (zero-based) */
    int* vtxdist = (int*) malloc((ntasks+1) * sizeof(int));
    int  mlocal = ((numnp+ntasks-1) / ntasks);
    int  i;

    for (i = 0; i < ntasks; ++i)
        vtxdist[i] = mlocal*i;
    vtxdist[ntasks] = numnp;
    return vtxdist;
}


/*
 * Open a file on processor 0, and exit gracefully on failure.
 */
FILE* fopen0(int id, char* name, char* mode)
{
    FILE* fp = NULL;
    int success;
    if (id == 0) {
        fp = fopen(name, mode);
        success = (fp != NULL);
        if (!success)
            fprintf(stderr, "Could not open %s\n", name);
    }
    MPI_Bcast(&success, 1, MPI_INT, 0, MPI_COMM_WORLD);
    if (!success) {
        MPI_Finalize();
        exit(-1);        
    }
    return fp;
}


/*
 * Close a file on processor 0
 */
void fclose0(int id, FILE* fp)
{
    if (id == 0)
        fclose(fp);
}


/*
 * Build the graph file name using the following rule:
 *  - Strip off the path to get path + Iname
 *  - Then use the name path + graph.name
 */
char* make_graph_name(char* iname)
{
    int   len   = strlen(iname);
    char* gname = (char*) malloc(len + 6);
    int   i;

    for (i = len; i >= 0 && iname[i] != '\\' && iname[i] != '/'; --i);
    memcpy(gname, iname, i+1);
    gname[i+1] = 0;
    strcat(gname, "graph.");
    strcat(gname, iname+i+2);

    return gname;
}


int main(int argc, char** argv)
{
    const int IX_CHUNK = 10;
    int id, ntasks, nparts;
    int numnp, numel, nen;
    int* vtxdist;
    sparse_graph_t* graph;
    int* ia;
    int* ja;
    int* part;
    char* graphname;
    FILE* fp;

    if (MPI_Init(&argc, &argv) != MPI_SUCCESS) {
        printf("MPI initialization failed!\n");
        exit(1);
    }
    MPI_Comm_size(MPI_COMM_WORLD, &ntasks);
    MPI_Comm_rank(MPI_COMM_WORLD, &id);

    if (argc < 3) {
        fprintf0(stderr, "Syntax: partition nparts infile [outfile]\n");
        MPI_Finalize();
        exit(-1);
    }

    nparts = atoi(argv[1]);
    if (nparts <= 1)
        nparts = ntasks;

    if (argc > 3)
        graphname = argv[3];
    else
        graphname = make_graph_name(argv[2]);
    dprintf0(stderr, "Filename: %s\n", graphname);

    fp = fopen0(id, argv[2], "r+");
    dprintf0(stderr, "Partitioning %d ways\n", nparts);
            
    dprintf0(stderr, "Reading file data\n");
    process_feap_header(fp, id, &numnp, &numel, &nen);
    vtxdist = build_node_dist(ntasks, numnp);
    graph = make_sparse_graph(vtxdist[id],vtxdist[id+1]-1);
    process_elt_cards(fp, id, graph, numel, nen, IX_CHUNK);
    fclose0(id, fp);

    dprintf0(stderr, "Converting to CSR\n");
    graph_to_csr(graph, &ia, &ja);

    dprintf0(stderr, "Partitioning\n");
    part = partition(id, vtxdist, ia, ja, nparts);

    dprintf0(stderr, "Write partitions\n");
    fp = fopen0(id, graphname, "w");
    print_part(fp, id, ntasks, vtxdist, part);
    free(part);

    dprintf0(stderr, "Writing nodal graph\n");
    print_csr(fp, id, ntasks, vtxdist, ia, ja);
    fclose0(id, fp);

    dprintf0(stderr, "Finishing\n");
    destroy_sparse_graph(graph);
    free(ia);
    free(ja);

    MPI_Finalize();
    exit(0);
}
