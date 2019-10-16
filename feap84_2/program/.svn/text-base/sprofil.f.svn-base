c$Id:$
      subroutine sprofil(n,nnn,e2,adj,xadj,oldpro,newpro)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--+---------+---------+---------+---------+---------+---------+-]
c     Purpose: Compute the profiles using both original and new node
c              numbers

c     Input:
c     ------

c     n      - total number of nodes in graph
c     nnn    - list of new node numbers for graph
c            - new node number for node i is given by nnn(i)
c     e2     - twice the number of edges in the graph = xadj(n+1)-1
c     adj    - adjacency list for all nodes in graph
c            - list of length 2e where e is the number of edges in
c              the graph and 2e = xadj(n+1)-1
c     xadj   - index vector for adj
c            - nodes adjacent to node i are found in adj(j), where
c              j = xadj(i), xadj(i)+1, ..., xadj(i+1)-1
c            - degree of node i given by xadj(i+1)-xadj(i)
c     oldpro - undefined
c     newpro - undefined

c     Output:
c     -------

c     n      - unchanged
c     nnn    - unchanged
c     e2     - unchanged
c     adj    - unchanged
c     xadj   - unchanged
c     oldpro - profile with original node numbering
c     newpro - profile with new node numbering

c     Note:      profiles include diagonal terms
c     -----

c     Programmer:             Scott Sloan
c     -----------

c     Last modified:          13 August 1991     Scott Sloan
c     --------------

c     Copyright 1989:         Scott Sloan
c     ---------------         Department of Civil Engineering
c                             University of Newcastle
c                             NSW 2308
c-----[--+---------+---------+---------+---------+---------+---------+-]
      implicit   none

      integer    i,j,n
      integer    e2
      integer    newmin,newpro,oldmin,oldpro
      integer    adj(e2),nnn(n)
      integer    xadj(n+1)

      save

c     Set profiles and loop over each node in graph

      oldpro = 0
      newpro = 0
      do i = 1,n

c       Find lowest numbered neighbour of node i
c       (using both old and new node numbers)

        oldmin = i
        newmin = nnn(i)
        do j = xadj(i),xadj(i+1)-1
          oldmin = min(oldmin,adj(j))
          newmin = min(newmin,nnn(adj(j)))
        end do ! j

c       Update profiles

        oldpro = oldpro + dim(    i ,oldmin)
        newpro = newpro + dim(nnn(i),newmin)
      end do ! i

      end
