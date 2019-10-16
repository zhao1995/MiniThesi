c$Id:$
      function flgama (w)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c     1. Use 'pi' from 'pconstant.h'                        14/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose:  Gamma function routine

c      Inputs:
c         w     -   Initial estimate of root

c      Outputs:
c         w     -   Gamma function
c-----[--.----+----.----+----.-----------------------------------------]
      implicit none

      include 'pconstant.h'

      integer  m, i
      real*8   flgama, w, x, p, fk, y, z,zz

      x  =  w
      fk = -1.d0

c     w less eq 0.5

      if (x .lt. 0.5d0) then
        m = 1
        p = pi/sin(x*pi)
        x = 1.d0 - x
      else
        m = 0
        p = 0.0d0
      endif

      do while(x + fk - 6.d0 .le.0.0d0)
        fk = fk + 1.d0
      end do ! while

      z  = x + fk
      zz = z*z

      y  = (z - 0.5d0)*log(z) - z + 0.9189385332047d0 + (((((-4146.d0/zz
     &        + 1820.d0)/zz - 1287.d0)/zz + 1716.d0)/zz - 6006d0)/zz
     &        + 180180.d0)/z/2162160.d0

      if(fk.gt.0.0d0) then
        do i = 1,int(fk)
          fk = fk - 1.d0
          y  = y - log(x + fk)
        end do ! i
      endif

      if(m.ne.0) then
        if(p.le.0.0d0) then
          write (*,2000) w
          y = 0.d0
        else
          y = log(p) - y
        endif
      endif
      flgama = y

2000  format (2x,'gamma(',e11.4,') is negative')

      end
