c$Id:$
      subroutine shl2bd(d,xl,r,ndm,ndf,nel)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1.  Remove test on axisymmetric - shell is axisym   17/04/2008
c-----[--+---------+---------+---------+---------+---------+---------+-]
c      Purpose:  Compute body forces for 2-d shells

c      Inputs:
c        d(*)        - Body loading in 10-11
c        xl(ndm,nel) - Element coordinates
c        ndm         - Mesh dimenstion
c        ndf         - Degree of freedoms/node
c        nel         - Number of element nodes

c      Outputs:
c        r(ndf,nel)  - Body forces added
c-----[--+---------+---------+---------+---------+---------+---------+-]
      implicit   none

      integer    ndm,ndf,nel, i,j,l,lint
      real*8     rr,xjac
      real*8     d(*),xl(ndm,nel),r(ndf,nel), sw(2,3),shp(2,3)

      lint = nel
      call int1d(lint,sw)

      do l = 1,lint
        call shp1d(sw(1,l),xl,shp,ndm,nel,xjac)
        xjac = xjac*sw(2,l)
        rr = 0.0d0
        do i = 1,nel
          rr = rr + shp(2,i)*xl(1,i)
        end do ! i
        xjac = xjac*rr
        do i = 1,ndm
          rr = d(10+i)*xjac
          do j = 1,nel
            r(i,j) = r(i,j) + rr*shp(2,j)
          end do ! j
        end do ! i
      end do ! l

      end ! shl2db
