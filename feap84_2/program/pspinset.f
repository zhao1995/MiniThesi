c$Id:$
      subroutine pspinset(pspinv,spnum, nprop,theta,nn,xc,v0,edge)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    19/03/2009
c       1. Add edge storage to spinv()                      04/01/2011
c-----[--.----+----.----+----.-----------------------------------------]
c     Purpose: Set parameters for applying spin boundary conditions

c     Input:
c        spinv(12,*)   - Spin data
c        spnum         - Spin set number

c     Output:
c        nprop         - Translation proportional load number
c        theta         - Spin velocity
c        nn(3)         - Spin axis normal
c        xc(3)         - Spin center
c        v0(3)         - Translational velocity
c        edge          - Edge coordinate
c-----[--.----+----.----+----.-----------------------------------------]
      implicit   none

      integer    spnum, nprop, j
      real*8     theta,edge
      real*8     pspinv(12,*), nn(3),xc(3),v0(3)

      do j = 1,3
        xc(j) = pspinv(j  ,spnum)
        nn(j) = pspinv(j+3,spnum)
        v0(j) = pspinv(j+8,spnum)
      end do ! j
      theta = pspinv(7,spnum)
      nprop = nint(pspinv(8,spnum))
      edge  = pspinv(12,spnum)

      end
