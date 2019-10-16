c$Id:$
      function pddot(n,v1,n1,v2,n2)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Added rdum(1) to fix compile error               05/01/2013
c          Change input to ddot to us passed in value
c          Remove common block pfeapb.h
c-----[--+---------+---------+---------+---------+---------+---------+-]
c      Purpose: Compute parallel vector dot product

c      Inputs:
c         n     - Length of vectors in current domain (numpeq)
c         v1(*) - Vector 1
c         n1    - Step size
c         v2(*) - Vector 1
c         n2    - Step size

c      Output:
c         pddot - real double precision
c-----[--+---------+---------+---------+---------+---------+---------+-]
       implicit   none

       integer    n, n1,n2
       real*8     pddot,v1(*),v2(*), ddot, tdatabuf(2), rdum(1)

       rdum(1)  = ddot(n,v1,n1,v2,n2)

       call pfeapsr(rdum,tdatabuf,1,.true.) ! Sum from all processors
       pddot = rdum(1)

       end
