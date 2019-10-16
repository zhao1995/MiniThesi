c$Id:$
      subroutine srootls(n,root,maxwid,e2,adj,xadj,mask,ls,xls,depth,
     &                   width)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--+---------+---------+---------+---------+---------+---------+-]
c     Purpose: Generate rooted level structure using a fortran 77
c              implementation of the algorithm given by George and Liu

c     Input:
c     ------

c     n      - number of nodes
c     root   - root node for level structure
c     maxwid - max permissible width of rooted level structure
c            - abort assembly of level structure if width is ge maxwid
c            - assembly ensured by setting maxwid = n+1
c     e2     - twice the number of edges in the graph = xadj(n+1)-1
c     adj    - adjacency list for all nodes in graph
c            - list of length 2e where e is the number of edges in
c              the graph and 2e = xadj(n+1)-1
c     xadj   - index vector for adj
c            - nodes adjacent to node i are found in adj(j), where
c              j = xadj(i), xadj(i)+1, ..., xadj(i+1)-1
c            - degree of node i is xadj(i+1)-xadj(i)
c     mask   - masking vector for graph
c            - visible nodes have mask = 0
c     ls     - undefined
c     xls    - undefined
c     depth  - undefined
c     width  - undefined

c     Output:
c     -------

c     n      - unchanged
c     root   - unchanged
c     maxwid - unchanged
c     e2     - unchanged
c     adj    - unchanged
c     xadj   - unchanged
c     mask   - unchanged
c     ls     - list containing a rooted level structure
c            - list of length nc
c     xls    - index vector for ls
c            - nodes in level i are found in ls(j), where
c              j = xls(i), xls(i)+1, ..., xls(i+1)-1
c            - list of max length nc+1
c     depth  - number of levels in rooted level structure
c     width  - width of rooted level structure

c     Note:  If width ge maxwid then assembly has been aborted
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
      integer    nbr
      integer    node,root
      integer    depth,jstop,jstrt,lstop,lstrt,lwdth,width
      integer    maxwid
      integer    ls(n)
      integer    adj(e2),xls(n+1)
      integer    mask(n),xadj(n+1)

      save

c     Initialisation

      mask(root) = 1
      ls(1)      = root
      nc         = 1
      width      = 1
      depth      = 0
      lstop      = 0
      lwdth      = 1
      do while (lwdth.gt.0)

c       lwdth is the width of the current level
c       lstrt points to start of current level
c       lstop points to end of current level
c       nc counts the nodes in component

        lstrt      = lstop+1
        lstop      = nc
        depth      = depth+1
        xls(depth) = lstrt

c       Generate next level by finding all visible neighbours
c       of nodes in current level

        do i = lstrt,lstop
          node  = ls(i)
          jstrt = xadj(node)
          jstop = xadj(node+1)-1
          do j = jstrt,jstop
            nbr = adj(j)
            if(mask(nbr).eq.0)then
              nc        = nc+1
              ls(nc)    = nbr
              mask(nbr) = 1
            end if
          end do ! j
        end do ! i

c       Compute width of level just assembled and the width of the
c       level structure so far

        lwdth = nc-lstop
        width = max(lwdth,width)

c       Abort assembly if level structure is too wide

        if(width.ge.maxwid) goto 35

      end do ! while
      xls(depth+1) = lstop+1

c     Reset mask=0 for nodes in the level structure

   35 continue
      do i = 1,nc
        mask(ls(i)) = 0
      end do ! i

      end
