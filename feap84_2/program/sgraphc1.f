c$Id:$
      subroutine sgraphc1(ixc,ncen,ncen1,numcels,xadj)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--+---------+---------+---------+---------+---------+---------+-]
c     Purpose: Form adjacency list for graph corresponding to a finite
c              element contact mesh

c     Input:
c     -----
c     ixc    - list of node numbers for each element
c     ncen    - maximum number of nodes on an element
c     ncen1   - dimension for ix(*,*) array with node connections
c     numcels - Number of contact elements
c     xadj   - undefined

c     Output:
c     -------
c     ixc     - unchanged
c     nce     - unchanged
c     ncen    - unchanged
c     ncen1   - unchanged
c     numcels - unchanged
c     xadj    - index vector for adj
c             - nodes adjacent to node i are found in adj(j), where
c               j = xadj(i), xadj(i)+1, ..., xadj(i+1)-1
c             - degree of node i given by xadj(i+1)-xadj(i)

c     Notes:
c     ------

c     This routine typically requires about 25 percent elbow room for
c     assembling the adj list (i.e. iadj/2e is typically around 1.25).
c     In some cases, the elbow room may be larger (iadj/2e is slightly
c     less than 2 for the 3-noded triangle) and in other cases it may be
c     zero (iadj/2e = 1 for bar elements)

c     Programmer:             Scott Sloan
c     -----------

c     Last modified:          1 March 1991        Scott Sloan
c     --------------

c     Copyright 1989:         Scott Sloan
c     ---------------         Department of Civil Engineering
c                             University of Newcastle
c                             NSW 2308
c-----[--+---------+---------+---------+---------+---------+---------+-]
      implicit   none

      integer    ncen,ncen1,numcels, i,j,nel,mel
      integer    ixc(ncen1,*), xadj(*)

      save

c     Estimate degree of each node (always an overestimate)

      do i = 1,numcels
        nel = 0
        mel = 0
        do j = 1,ncen
          if(ixc(j,i).gt.0) then
            nel = nel + 1
            mel = j
          endif
        end do ! j

        do j = 1,mel
          if(ixc(j,i).gt.0) then
            xadj(ixc(j,i)) = xadj(ixc(j,i)) + nel
          endif
        end do ! j
      end do ! i

      end
