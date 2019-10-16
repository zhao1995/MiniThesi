c$Id:$
      subroutine hlcn2z(xx,jac,nshp,linr,lint,nps)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    04/01/2014
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: SPR projection of element history variables to nodes

c      Inputs:
c        xx(10,*)     - Coordinate for patch
c        jac(*)       - Jacobian
c        nshp(nel)    - Shape functions
c        linr         - Number of projection value points
c        lint         - Number quadrature points
c        nps          - Number of coordinate functions
c        nel          - Number nodes on element

c      Outputs:
c        ehp(30,*)    - Integral of variables
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'eldatp.h'
      include  'zzcom1.h'

      integer   linr,lint,nps
      real*8    xx(10,*), jac(*), nshp(16,*)

      integer   i, l
      real*8    plhm(30), cons

      save

c     Compute projections for SPR

      do l = 1,linr

c       Project history points to reduced points

        plhm(1:plhmax) = 0.0d0
        do i = 1,lint
          plhm(1:plhmax) = plhm(1:plhmax) + nshp(i,l)*plhis(1:plhmax,i)
        end do ! i

c       Compute projection matrix

        do i = 1,nps
          cons            = xx(i,l)*jac(l)
          ehp(1:plhmax,i) = ehp(1:plhmax,i) + cons*plhm(1:plhmax)
        end do ! i

      end do ! l

      end
