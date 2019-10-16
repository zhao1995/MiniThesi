/*
c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    06/05/2012
c-----[--.----+----.----+----.-----------------------------------------]
c     Purpose:  C cover function to partition mesh using METIS ver 5.x

c     Inputs:
c        numnp   -- number of nodes
c        xadj    -- pointers into adjaceny matrix
c        adjncy  -- adjaceny matrix for node graph
c        domains -- number of domains to partition graph into

c     Outputs:
c        part    -- partitioned graph
c-----[--.----+----.----+----.-----------------------------------------]
*/

#include <stdio.h>
#include "metis.h"

int smetis_(int *numnp,int *xadj,int *adjncy,int *domains, int *part) {

	int rval;
 	int ncon = 1;
        int edgecut = 0;
        int options[METIS_NOPTIONS];
        
        // Default options 
        METIS_SetDefaultOptions(options);

        options[METIS_OPTION_OBJTYPE]   = 0;
        options[METIS_OPTION_CTYPE]     = 1;
        options[METIS_OPTION_DBGLVL]    = 0;
        options[METIS_OPTION_NITER]     = 10;
        options[METIS_OPTION_NCUTS]     = 1;
        options[METIS_OPTION_MINCONN]   = 0;
        options[METIS_OPTION_CONTIG]    = 0;
        options[METIS_OPTION_NUMBERING] = 1;

        if (*domains < 8) {
	//  Recursive
       		options[METIS_OPTION_IPTYPE]    = 0;
        	options[METIS_OPTION_RTYPE]     = 0;
        	options[METIS_OPTION_UFACTOR]   = 1;

        	rval = METIS_PartGraphRecursive(numnp,&ncon,xadj,adjncy,
                       NULL,NULL,NULL,
                       domains,NULL,NULL,options,&edgecut,part);
	}
	else {
	// Kway
        	options[METIS_OPTION_IPTYPE]    = 4;
       		options[METIS_OPTION_RTYPE]     = 1;
        	options[METIS_OPTION_UFACTOR]   = 30;

        	rval = METIS_PartGraphKway(numnp,&ncon,xadj,adjncy,
                       NULL,NULL,NULL,
                       domains,NULL,NULL,options,&edgecut,part);
	}

     	return(rval);

}
