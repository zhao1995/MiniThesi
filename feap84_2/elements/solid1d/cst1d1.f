c$Id:$
      subroutine cst1d1(shp,shpr,xsj,xsj0,xx,epsr,ul,dr,di,
     &                  si,pr,pi)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Complex 1-d strains

c      Inputs:

c      Outputs:
c-----[--.----+----.----+----.-----------------------------------------]
      implicit   none

      include   'cdata.h'
      include   'eldata.h'
      include   'eltran.h'
      include   'pmod2d.h'
      include   'sdata.h'

      integer    j,j1,k,k1
      real*8     shp(2,*),shpr(*),epsr(*),epsi(6),ul(ndf,*)
      real*8     dr(6,6),di(6,6), sigr(6),sigi(6),si(nst,*),pr(*),pi(*)
      real*8     xsj,xsj0,xx, aj1,aj2, bd11,bd13

c     Compute imaginary strains

      do j = 1,6
        epsi(j) = 0.0d0
      end do ! j
      do j = 1,nel
        epsi(1) = epsi(1) + shp(1,j)*ul(1,j)
        epsi(3) = epsi(3) + shp(2,j)*ul(1,j)
      end do ! j

c     Multiply jacobian by radius for axisymmetry

      if(stype.eq.3) then
        epsi(3) = epsi(3)/xx
      else
        epsi(3) = 0.0d0
      end if

c     compute imaginary stress

      do j = 1,4
        sigr(j) = 0.0d0
        sigi(j) = 0.0d0
        do k = 1,4
          sigr(j) = sigr(j) + di(j,k)*epsi(k)
          sigi(j) = sigi(j) + dr(j,k)*epsi(k) + di(j,k)*epsr(k)
        end do ! k
      end do ! j

      j1 = 1
      do j = 1,nel

        aj1 = shp(1,j)*xsj
        aj2 = shp(2,j)*xsj0

c       Compute B_trans * D * j * w

        bd11 = aj1*di(1,1) + aj2*di(3,1)
        bd13 = aj1*di(1,3) + aj2*di(3,3)

c       Loop over columns (symmetry noted)

        k1 = 1
        do k = 1,nel
          si(j1  ,k1  ) = si(j1  ,k1  ) + (bd11*shp(1,k)
     &                                  +  bd13*shpr(k))*ctan(1)
          k1 = k1 + ndf
        end do ! k

c       Residual for added real part

        pr(j1  ) = pr(j1  ) + aj1*sigr(1) + aj2*sigr(3)

c       Residual for imaginary part

        pi(j1  ) = pi(j1  ) - aj1*sigi(1) - aj2*sigi(3)

        j1 = j1 + ndf
      end do ! j

      end
