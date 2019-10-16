c$Id:$
      subroutine sdiamtr(n,e2,adj,xadj,mask,ls,xls,hlevel,snode,nc)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--+---------+---------+---------+---------+---------+---------+-]
c     Purpose: Find nodes which define psuedo-diameter of graph and
c              store distances from end node

c     Input:
c     ------

c     n      - the total number of nodes in the graph
c     e2     - twice the number of edges in the graph  = xadj(n+1)-1
c     adj    - adjacency list for all nodes in the graph
c            - list of length 2e where e is the number of edges in
c              the graph and 2e = xadj(n+1)-1
c     xadj   - index vector for adj
c            - nodes adjacent to node i are found in adj(j), where
c              j = xadj(i), xadj(i)+1, ...,xadj(i+1)-1
c            - degree of node i given by xadj(i+1)-xadj(i)
c     mask   - masking vector for graph
c            - visible nodes have mask = 0, node invisible otherwise
c     ls     - undefined
c     xls    - undefined
c     hlevel - undefined
c     snode  - undefined
c     nc     - undefined

c     Output:
c     ------

c     n      - unchanged
c     e2     - unchanged
c     adj    - unchanged
c     xadj   - unchanged
c     mask   - list of distances of nodes from the end node
c     ls     - list of nodes in the component
c     xls    - not used
c     hlevel - not used
c     snode  - starting node for numbering
c     nc     - the number of nodes in this component of graph

c     subroutines called:  srootls, sisorti
c     -------------------

c     Note:      snode and enode define a pseudo-diameter
c     -----

c     Programmer:             Scott Sloan
c     -----------

c     Last modified:          1 March 1991      Scott Sloan
c     --------------

c     Copyright 1989:         Scott Sloan
c     ---------------         Department of Civil Engineering
c                             University of Newcastle
c                             NSW 2308
c-----[--+---------+---------+---------+---------+---------+---------+-]
      implicit   none

      integer    i,j,n
      integer    e2,nc
      integer    node
      integer    depth,enode,hsize,istop,istrt,jstop,jstrt,snode,width
      integer    degree,ewidth,mindeg,sdepth
      integer    ls(n)
      integer    adj(e2),xls(n+1)
      integer    mask(n),xadj(n+1)
      integer    hlevel(n)

      save

c     Choose first guess for starting node by min degree
c     Ignore nodes that are invisible (mask ne 0)

      mindeg = n
      do i = 1,n
        if(mask(i).eq.0)then
          degree = xadj(i+1)-xadj(i)
          if(degree.lt.mindeg)then
            snode  = i
            mindeg = degree
          end if
        end if
      end do ! i

c     Generate level structure for node with min degree

      call srootls(n,snode,n+1,e2,adj,xadj,mask,ls,xls,sdepth,width)

c     Store number of nodes in this component

      nc = xls(sdepth+1)-1

c     Iterate to find start and end nodes

   15 continue

c     Store list of nodes that are at max distance from starting node
c     Store their degrees in xls

      hsize = 0
      istrt = xls(sdepth)
      istop = xls(sdepth+1)-1
      do i = istrt,istop
        node          =  ls(i)
        hsize         =  hsize+1
        hlevel(hsize) =  node
        xls(node)     =  xadj(node+1)-xadj(node)
      end do ! i

c     Sort list of nodes in ascending sequence of their degree
c     Use insertion sort algorithm

      if(hsize.gt.1) call sisorti(hsize,hlevel,n,xls)

c     Remove nodes with duplicate degrees

      istop  = hsize
      hsize  = 1
      degree = xls(hlevel(1))
      do i = 2,istop
        node = hlevel(i)
        if(xls(node).ne.degree)then
          degree        = xls(node)
          hsize         = hsize+1
          hlevel(hsize) = node
        endif
      end do ! i

c     Loop over nodes in shrunken level

      ewidth = nc+1
      do i = 1,hsize
        node = hlevel(i)

c       Form rooted level structures for each node in shrunken level

        call srootls(n,node,ewidth,e2,adj,xadj,mask,ls,xls,depth,width)
        if(width.lt.ewidth)then

c         Level structure was not aborted during assembly

          if(depth.gt.sdepth)then

c           Level structure of greater depth found
c           Store new starting node, new max depth, and begin
c           a new iteration

            snode  = node
            sdepth = depth
            goto 15
          endif

c         Level structure width for this end node is smallest so far
c         Store end node and new min width

          enode  = node
          ewidth = width
        end if
      end do ! i

c     Generate level structure rooted at end node if necessary

      if(node.ne.enode)then
        call srootls(n,enode,nc+1,e2,adj,xadj,mask,ls,xls,depth,width)
      endif

c     Store distances of each node from end node

      do i = 1,depth
        jstrt = xls(i)
        jstop = xls(i+1)-1
        do j = jstrt,jstop
          mask(ls(j)) = i-1
        end do ! j
      end do ! i

      end
