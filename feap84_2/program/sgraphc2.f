c$Id:$
      subroutine sgraphc2(ne,ncen,ncen1,ixc,adj,xadj)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--+---------+---------+---------+---------+---------+---------+-]
c     Purpose:
c     --------

c     Form adjacency list for a graph corresponding to a finite element
c     mesh

c     Input:
c     -----

c     ne     - number of contact elements in finite element mesh
c     ncen   - maximum number of nodes on an element
c     ncen1  - dimension for ixc(*,*) array with node connections
c     ix     - list of node numbers for each element
c     adj    - undefined
c     xadj   - undefined

c     Output:
c     -------

c     n      - unchanged
c     ne     - unchanged
c     ncen   - unchanged
c     ncen1  - unchanged
c     ix     - unchanged
c     adj    - adjacency list for all nodes in graph
c            - list of length 2e where e is number of edges in graph
c              (note that 2e = xadj(n+1)-1 )
c     xadj   - index vector for adj
c            - nodes adjacent to node i are found in adj(j), where
c              j = xadj(i), xadj(i)+1, ..., xadj(i+1)-1
c            - degree of node i given by xadj(i+1)-xadj(i)

c     Notes:
c     ------

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

      include   'iofile.h'

      integer    i,j,k,l,m
      integer    ne
      integer    ncen,ncen1, mel
      integer    nodej,nodek
      integer    adj(*),ixc(ncen1,ne)
      integer    xadj(*)

      save

c     Form adjacency list (which may contain zeros)

      do i = 1,ne
        mel = 0
        do j = 1,ncen
          if(ixc(j,i).gt.0) then
            mel = j
          endif
        end do ! j
        do j = 1,mel-1
          if(ixc(j,i).gt.0) then
            nodej = ixc(j,i)
            do k = j+1,mel
              if(ixc(k,i).gt.0) then
                nodek = ixc(k,i)
                do l = xadj(nodej),xadj(nodej+1)-1
                  if(adj(l).eq.nodek) goto 70
                  if(adj(l).eq.    0) goto 55
                end do ! l

c               Error

                write(iow,1000)
                call plstop()

c               Add to adjacency graph

   55           continue
                adj(l) = nodek
                do m = xadj(nodek),xadj(nodek+1)-1
                  if(adj(m).eq.0)  goto 65
                end do ! m

c               Error

                write(iow,1000)
                call plstop()
   65           continue
                adj(m) = nodej
              endif
   70         continue
            end do ! k
          endif
        end do ! j
      end do ! i

c     Error message

 1000 format(//,1x,'***ERROR in Sgraph***',
     &       //,1x,'cannot assemble node adjacency list',
     &       //,1x,'check ix array')

      end
