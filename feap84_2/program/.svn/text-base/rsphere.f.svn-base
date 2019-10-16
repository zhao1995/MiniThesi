c$Id:$
      subroutine rsphere(n1,n2,xx,rcg,rlam,inert,rresid, acc)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Spherical joint for explicit integrations

c      Inputs:
c         iln(2)    - Line style: 1 = type; 2 = width

c      Outputs:
c         none      - Set output data through commons
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'iofile.h'
      include  'ptdat6.h'
      include  'rigid1.h'
      include  'tdata.h'
      include  'tdato.h'

      integer   n1,n2, i,j
      real*8    m1,ratiog, gap(3),dgap(3)
      real*8    rcg(3,11,*), rlam(3,3,6,*), y1(3),y2(3), x1(3),x2(3)
      real*8    inert(6,6,*), rresid(6,*), acc(6,*), dda(9), xx(3,2)
      real*8    rgc(9),tgc(9,9), tt(3), c11(3,3), c21(3,3), c22(3,3)

      save

c     Compute y1,y2 at time t_n+1

      do i = 1,3
        y1(i)  = rlam(i,1,3,n1)*(xx(1,1) - rcg(1,1,n1))
     &         + rlam(i,2,3,n1)*(xx(2,1) - rcg(2,1,n1))
     &         + rlam(i,3,3,n1)*(xx(3,1) - rcg(3,1,n1))
        y2(i)  = rlam(i,1,3,n2)*(xx(1,1) - rcg(1,1,n2))
     &         + rlam(i,2,3,n2)*(xx(2,1) - rcg(2,1,n2))
     &         + rlam(i,3,3,n2)*(xx(3,1) - rcg(3,1,n2))
        gap(i) = (rcg(i,2,n1) - rcg(i,2,n2)) + (y1(i) - y2(i))
      end do ! i
      dgap(1) = (rcg(1,5,n1) - rcg(1,5,n2))
     &        + rlam(2,2,5,n1)*y1(3) - rlam(3,2,5,n1)*y1(2)
     &        - rlam(2,2,5,n2)*y2(3) + rlam(3,2,5,n2)*y2(2)
      dgap(2) = (rcg(2,5,n1) - rcg(2,5,n2))
     &        + rlam(3,2,5,n1)*y1(1) - rlam(1,2,5,n1)*y1(3)
     &        - rlam(3,2,5,n2)*y2(1) + rlam(1,2,5,n2)*y2(3)
      dgap(3) = (rcg(3,5,n1) - rcg(3,5,n2))
     &        + rlam(1,2,5,n1)*y1(2) - rlam(2,2,5,n1)*y1(1)
     &        - rlam(1,2,5,n2)*y2(2) + rlam(2,2,5,n2)*y2(1)

      epl(20) = sqrt(gap(1)**2 +gap(2)**2 +gap(3)**2)
      ratiog  = min(2.d0,4.d0/((dt+dtold)*dt))
      do i = 1,3
        gap (i) = gap (i)*ratiog
      end do ! i

c     Form tangent and residual

      do i = 1,9
        do j = 1,9
          tgc(j,i) = 0.0d0
        end do ! j
      end do ! i

c     Compute residual

      call sasbtdc(rlam(1,2,5,n1),rlam(1,2,5,n1),y1, x1)
      call sasbtdc(rlam(1,2,5,n2),rlam(1,2,5,n2),y2, x2)

      m1 = inert(1,1,n1)
      do i = 1,3
        x1(i)    = x1(i) - x2(i) - xx(1,2)*gap(i) - xx(2,2)*dgap(i)
        tt(i)    = rresid(i  ,1) - m1*x1(i)
        rgc(i  ) = rresid(i  ,2) + tt(i)
        rgc(i+3) = rresid(i+3,2)
        rgc(i+6) = rresid(i+3,1)
        dda(i  ) = rcg(i,6,n2)
        dda(i+3) = rlam(i,3,5,n2)
        dda(i+6) = rlam(i,3,5,n1)
      end do ! i

      rgc(4)   = rgc(4) + (y2(2)*tt(3) - y2(3)*tt(2))
      rgc(5)   = rgc(5) + (y2(3)*tt(1) - y2(1)*tt(3))
      rgc(6)   = rgc(6) + (y2(1)*tt(2) - y2(2)*tt(1))

      rgc(7)   = rgc(7) - (y1(2)*tt(3) - y1(3)*tt(2))
      rgc(8)   = rgc(8) - (y1(3)*tt(1) - y1(1)*tt(3))
      rgc(9)   = rgc(9) - (y1(1)*tt(2) - y1(2)*tt(1))

c     Constrained tangent part

      call  sasbt(y2,y2, c22)
      call  sasbt(y2,y1, c21)
      call  sasbt(y1,y1, c11)
      do i = 1,3
        do j = 1,3
          tgc(j  ,i  ) = tgc(j  ,i  ) + inert(j  ,i  ,n1)
          tgc(j+3,i+3) = tgc(j+3,i+3) + m1*c22(j,i)
          tgc(j+3,i+6) = tgc(j+3,i+6) - m1*c21(j,i)
          tgc(j+6,i+3) = tgc(j+6,i+3) - m1*c21(i,j)
          tgc(j+6,i+6) = tgc(j+6,i+6) + m1*c11(j,i)
        end do ! j
      end do ! i

      tgc(1,5) =  m1*y2(3)
      tgc(5,1) =  tgc(1,5)

      tgc(1,6) = -m1*y2(2)
      tgc(6,1) =  tgc(1,6)

      tgc(1,8) = -m1*y1(3)
      tgc(8,1) =  tgc(1,8)

      tgc(1,9) =  m1*y1(2)
      tgc(9,1) =  tgc(1,9)

      tgc(2,4) = -m1*y2(3)
      tgc(4,2) =  tgc(2,4)

      tgc(2,6) =  m1*y2(1)
      tgc(6,2) =  tgc(2,6)

      tgc(2,7) =  m1*y1(3)
      tgc(7,2) =  tgc(2,7)

      tgc(2,9) = -m1*y1(1)
      tgc(9,2) =  tgc(2,9)

      tgc(3,4) =  m1*y2(2)
      tgc(4,3) =  tgc(3,4)

      tgc(3,5) = -m1*y2(1)
      tgc(5,3) =  tgc(3,5)

      tgc(3,7) = -m1*y1(2)
      tgc(7,3) =  tgc(3,7)

      tgc(3,8) =  m1*y1(1)
      tgc(8,3) =  tgc(3,8)

c     Remainder of residual

      do i = 1,9
        do j = 1,9
          rgc(i) = rgc(i) - tgc(i,j)*dda(j)
        end do ! j
      end do ! i

c     Remainder of tangent

      do i = 1,3
        do j = 1,3
          tgc(j  ,i  ) = tgc(j  ,i  ) + inert(j  ,i  ,n2)
          tgc(j+3,i+3) = tgc(j+3,i+3) + inert(j+3,i+3,n2)
          tgc(j+6,i+6) = tgc(j+6,i+6) + inert(j+3,i+3,n1)
        end do ! j
      end do ! i

      call invert(tgc,9,9)

c     Compute constrained acceleration

      do i = 1,9
        dda(i) = 0.0d0
        do j = 1,9
          dda(i) = dda(i) + tgc(i,j)*rgc(j)
        end do ! j
      end do ! i

c     Store accelerations for n1-n2 constraint

      do i = 1,3
        acc(i  ,2) = dda(i)
        acc(i+3,2) = dda(i+3)
        acc(i  ,1) = dda(i) + x1(i)
        acc(i+3,1) = dda(i+6)
      end do ! i

c     Recover rest of n1 acceleration r-dot-dot

      acc(1,1) = acc(1,1) - (y2(2)*dda(6) - y2(3)*dda(5))
     &                    + (y1(2)*dda(9) - y1(3)*dda(8))
      acc(2,1) = acc(2,1) - (y2(3)*dda(4) - y2(1)*dda(6))
     &                    + (y1(3)*dda(7) - y1(1)*dda(9))
      acc(3,1) = acc(3,1) - (y2(1)*dda(5) - y2(2)*dda(4))
     &                    + (y1(1)*dda(8) - y1(2)*dda(7))

      end
