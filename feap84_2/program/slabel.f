c$Id:$
      subroutine slabel(n,e2,adj,xadj,nnn,iw,oldpro,newpro)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--+---------+---------+---------+---------+---------+---------+-]
c     Purpose: Label a graph for small profile and rms wavefront

c     Input:
c     ------

c     n      - total number of nodes in graph
c     e2     - twice the number of edges in the graph = xadj(n+1)-1
c     adj    - adjacency list for all nodes in graph
c            - list of length 2e where e is the number of edges in
c              the graph and 2e = xadj(n+1)-1
c     xadj   - index vector for adj
c            - nodes adjacent to node i are found in adj(j), where
c              j = xadj(i), xadj(i)+1, ..., xadj(i+1)-1
c            - degree of node i given by xadj(i+1)-xadj(i)
c     nnn    - undefined
c     iw     - undefined
c     oldpro - undefined
c     newpro - undefined

c     Output:
c     -------

c     n      - unchanged
c     e2     - unchanged
c     adj    - unchanged
c     xadj   - unchanged
c     nnn    - list of new node numbers
c            - new number for node i given by nnn(i)
c            - if original node numbers give a smaller profile then
c              nnn is set so that nnn(i)=i for i=1,n
c     iw     - not used
c     oldpro - profile using original node numbering
c     newpro - profile for new node numbering
c            - if original profile is smaller than new profile, then
c              original node numbers are used and newpro=oldpro

c     subroutines called:  sdiamtr, snumber, sprofil
c     -------------------

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

      integer    i,n
      integer    e2,i1,i2,i3,nc
      integer    snode
      integer    lstnum,newpro,oldpro
      integer    iw(3*n+1)
      integer    adj(e2),nnn(n)
      integer    xadj(n+1)

      save

c     Set all new node numbers = 0
c     This is used to denote all visible nodes

      do i = 1,n
        nnn(i) = 0
      end do ! i

c     Define offsets

      i1 = 1
      i2 = i1+n
      i3 = i2+n+1

c     Loop while some nodes remain unnumbered

      lstnum = 0
      do while (lstnum.lt.n)

c       Find end points of p-diameter for nodes in this component
c       Compute distances of nodes from end node

        call sdiamtr(n,e2,adj,xadj,nnn,iw(i1),iw(i2),iw(i3),snode,nc)

c       Number nodes in this component

        call snumber(n,nc,snode,lstnum,e2,adj,xadj,nnn,iw(i1),iw(i2))

      end do ! while

c     Compute profiles for old and new node numbers

      call sprofil(n,nnn,e2,adj,xadj,oldpro,newpro)

c     Use original numbering if it gives a smaller profile

      if(oldpro.lt.newpro)then
        do i = 1,n
          nnn(i) = i
        end do ! i
        newpro = oldpro
      end if

      end
