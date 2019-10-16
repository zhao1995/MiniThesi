c$Id:$
      subroutine zzpro1(ix,ib,ip,
     &                  ma,ndm,nen,nen1,numnp,numel, ipmax)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Remove extra statements around counting of ib(j) 05/07/2011
c          Increase 'ip' to numnp+1 for efficiency
c-----[--+---------+---------+---------+---------+---------+---------+-]
c      Purpose: Zienkiewicz-Zhu patch determination

c      Inputs:
c        ix(nen1,numel)   - Element connection array
c        ma               - Material number to project (ma = 0 for all)
c        ndm              - Mesh space dimension
c        nen              - Maximum number of nodes/element
c        nen1             - Dimension of 'ix' array
c        numnp            - Number of nodal points
c        numel            - Number of elements

c      Outputs:
c        ib(numnp)        - Number of element connected to each node
c        ip(numnp+1)      - Pointers for element patches
c        ipmax            - Size of element patch array
c-----[--+---------+---------+---------+---------+---------+---------+-]
      implicit   none

      include   'pointer.h'
      include   'comblk.h'

      integer    ndm, ma, nen,nen1,numnp,numel, ipmax
      integer    ij, j, n, nnc

      integer    ix(nen1,numel), ib(numnp), ip(numnp+1)

      save

c     1.) Identify active vertex nodes

c     Set boundary node indicator

c     call iprint(ix,nen1,numel,nen1,'IX_B')

      do n = 1,numnp
        ib(n) =  0
      end do ! n

c     Loop through elements to set up list

      do n = 1,numel
        if(ma.eq.0 .or. ix(nen1,n).eq.ma) then

c         Compute number of nodes on element

          do ij = 1,nen
            j = ix(ij,n)
            if(j.gt.0) then
              ib(j) = ib(j) + 1
            endif
          end do ! ij

        endif
      end do ! n

c     2.) Set pointers for patch elements

c     Set number of vertex nodes for each element

      if(ndm.eq.1) then                      ! Line/pt
        nnc = 2
      elseif(ndm.eq.2) then                  ! Line/pt
        nnc = 3
      elseif(ndm.eq.3) then                  ! Line/pt
        nnc = 3
      endif

      ip(1) = 1
      do n = 1,numnp
        if(ib(n).lt.nnc) then
          ib(n)   = 0
        endif
        ip(n+1) = ip(n) + ib(n)
      end do ! n

      ipmax = ip(numnp+1)

      end
