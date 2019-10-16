c$Id:$
      function pdnrm2(n, v1, n1)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Change to use input argument n                   07/01/2013
c          Remove common block pfeapb.h
c-----[--+---------+---------+---------+---------+---------+---------+-]
c      Purpose: Compute Euclidian-norm of v1

c      Inputs:
c         n    - Length of vectors in current domain (numpeq)
c         v1(*) - Vector 1
c         n1    - Step size

c      Output:
c         pdnrm2 - real double precision
c-----[--+---------+---------+---------+---------+---------+---------+-]
       implicit   none

       integer    n,  n1
       real*8     pdnrm2,v1(*), pddot

       pdnrm2 = pddot(n, v1,n1, v1,n1 )
       pdnrm2 = sqrt(pdnrm2)

       end
