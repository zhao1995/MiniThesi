c$Id:$
      subroutine bdivwt(xbm,sv,wt,numnp,jmax)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Divide plot values by weight

c      Inputs:
c         wt(*)      - Weights
c         numnp      - Number of nodes
c         jmax       -

c      Outputs:
c         xbm(3,*,*) - Beam surface coordinates
c         sv(jmax,*) - Plot values
c-----[--.----+----.----+----.-----------------------------------------]
      implicit none

      include 'iofile.h'

      integer  j,jmax,n,numnp
      real*8   rwt
      real*8   xbm(3,jmax,numnp),sv(jmax,numnp),wt(numnp)

      save

      do n = 1,numnp
        if(wt(n).gt.0.0d0) then
          rwt = 1.d0/wt(n)
        else
          rwt = 0.d0
        endif
        do j = 1,jmax
          xbm(1,j,n) = xbm(1,j,n)*rwt
          xbm(2,j,n) = xbm(2,j,n)*rwt
          xbm(3,j,n) = xbm(3,j,n)*rwt
          sv(j,n)    = sv(j,n)*rwt
        end do ! j
      end do ! n

      end
