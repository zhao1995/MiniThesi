c$Id:$
      subroutine gausspw (nn,sw)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006

c-----[--.----+-!--.----+----.----+------------------------------------]
c     Input:
c     nn       = Number of Gauss Points to compute

c     Outputs:
c     sw(1,*)  = Gauss Point Coordinates
c     sw(2,*)  = Gauss Point Weights
c-----[--.----+-!--.----+----.----+------------------------------------]
      implicit   none

      integer    nn, n
      real*8     sw(2,*)
      real*8     fn, beta, cc, xt, dpn,pn1, flgama

      fn   = dble(nn)
      beta = exp(2.d0*flgama(1.d0) - flgama(2.d0))
      cc   = 2.d0*beta
      do n = 2,nn
        cc = cc*4.d0*dble(n-1)**4
     &       / (dble(2*n-1)*dble(2*n-3)*dble(2*n-2)**2)
      end do ! n

      do n = 1,nn

c       Largest zero

        if(n.eq.1) then
          xt = 1.d0 - 2.78d0/(4.d0 + fn*fn)

c       Second zero

        elseif(n.eq.2) then
          xt = xt - (4.1d0 + 0.246d0*(fn - 8.d0)/fn)*(1.d0 - xt)

c       Third zero

        elseif(n.eq.3) then
          xt = xt - (1.67d0 + 0.3674d0*(fn - 8.d0)/fn)*(sw(1,1) - xt)

c       Second last zero

        elseif(n.eq.nn-1) then
          xt = xt + (xt - sw(1,n-2))/0.766d0/(1.d0 + 0.639d0*(fn-4.d0)
     &                                      /(1.d0 + 0.710d0*(fn-4.d0)))

c       Last zero

        elseif(n.eq.nn) then
          xt = xt + (xt - sw(1,n-2))/1.67d0/(1.d0 + 0.22d0*(fn-8.d0)/fn)

c       Intermediate roots

        else
          xt    = 3.d0*sw(1,n-1) - 3.d0*sw(1,n-2) + sw(1,n-3)
        endif

c       Find root using xt-value

        call root (xt,nn,dpn,pn1)
        sw(1,n) = xt
        sw(2,n) = cc/(dpn*pn1)

      end do ! n

c     Reverse order of points

      do n = 1, nn
        sw(1,n) = - sw(1,n)
      end do ! n

      end
