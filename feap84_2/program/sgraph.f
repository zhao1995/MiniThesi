c$Id:$
      subroutine sgraph(n,ne,nen,nen1,ix,iadj,adj,xadj, e2)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--+---------+---------+---------+---------+---------+---------+-]
c     Purpose: Form adjacency list for graph corresponding to a finite
c              element mesh

c     Input:
c     -----

c     n      - number of nodes in graph (finite element mesh)
c     ne     - number of elements in finite element mesh
c     nen    - maximum number of nodes on an element
c     nen1   - dimension for ix(*,*) array with node connections
c     ix     - list of node numbers for each element
c     iadj   - length of vector adj
c            - set iadj=ne*nen*(nen-1) for a mesh of elements
c              with nen nodes
c     adj    - undefined
c     xadj   - undefined

c     Output:
c     -------

c     n      - unchanged
c     ne     - unchanged
c     nen    - unchanged
c     nen1   - unchanged
c     ix     - unchanged
c     iadj   - unchanged
c     adj    - adjacency list for all nodes in graph
c            - list of length 2e where e is number of edges in graph
c              (note that 2e = xadj(n+1)-1 )
c     xadj   - index vector for adj
c            - nodes adjacent to node i are found in adj(j), where
c              j = xadj(i), xadj(i)+1, ..., xadj(i+1)-1
c            - degree of node i given by xadj(i+1)-xadj(i)

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

      include   'comblk.h'
      include   'compac.h'
      include   'iofile.h'
      include   'pointer.h'

      integer    i,j,k,l,m,n
      integer    ne, e2
      integer    iadj,nen,nen1, nel,mel
      integer    jstop,jstrt,nodej,nodek
      integer    adj(iadj),ix(nen1,ne)
      integer    xadj(n+1)

      save

c     Initialise adjacency list and its index vector

      do i = 1,iadj
        adj(i) = 0
      end do ! i
      do i = 1,n
        xadj(i) = 0
      end do ! i

c     Estimate degree of each node (always an overestimate)

      do i = 1,ne
        nel = 0
        mel = 0
        do j = 1,nen
          if(ix(j,i).gt.0) then
            nel = nel + 1
            mel = j
          endif
        end do ! j

        do j = 1,mel
          if(ix(j,i).gt.0) then
            nodej       = ix(j,i)
            xadj(nodej) = xadj(nodej) + nel
          endif
        end do ! j
      end do ! i

c     Check the contact elements

      if(numcels.gt.0) then
        call sgraphc1(mr(np(168)),ncen,ncen1,numcels,xadj)
      endif

c     Reconstruct xadj to point to start of each set of neighbours

      e2 = 1
      do i = 1,n
        e2      = e2 + xadj(i)
        xadj(i) = e2 - xadj(i)
      end do ! i
      xadj(n+1) = e2

c     Form adjacency list (which may contain zeros)

      do i = 1,ne
        mel = 0
        do j = 1,nen
          if(ix(j,i).gt.0) then
            mel = j
          endif
        end do ! j
        do j = 1,mel-1
          if(ix(j,i).gt.0) then
            nodej = ix(j,i)
            do k = j+1,mel
              if(ix(k,i).gt.0) then
                nodek = ix(k,i)
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

      if(numcels.gt.0) then
        call sgraphc2(numcels,ncen,ncen1,mr(np(168)),adj,xadj)
      endif

c     Strip any zeros from adjacency list

      k     = 0
      jstrt = 1
      do i = 1,n
        jstop = xadj(i+1) - 1
        do j = jstrt,jstop
          if(adj(j).eq.0) goto 105
          k      = k + 1
          adj(k) = adj(j)
        end do ! j
  105   continue
        xadj(i+1) = k + 1
        jstrt     = jstop + 1
      end do ! i

      e2 = xadj(n+1) - 1

c     Error message

 1000 format(//,1x,'***ERROR in Sgraph***',
     &       //,1x,'cannot assemble node adjacency list',
     &       //,1x,'check ix array')

      end
