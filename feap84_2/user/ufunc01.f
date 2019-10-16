c$Id:$
      subroutine ufunc01(fu,xi, f, fsw)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: User function: f(xi)
c      Inputs:
c        fu(14)   - Parameters of function
c        xi       - Parameter value:  -1 < xi < 1
c        fsw      - Switch value: 0 - output type; 1 - compute f(xi)

c      Outputs:
c        f        - Function value
c-----[--.----+----.----+----.-----------------------------------------]
      implicit   none

      include   'iofile.h'
      include   'print.h'
      include   'pconstant.h'

      character  sgnc*17
      integer    fsw
      real*8     xi, f, fu(14)

      if(fsw.eq.0) then

        sgnc = '   '
        if(prt) then
          write(iow,2000) fu(1)
          if(fu(2).ne.0.0d0) then
            if(fu(2).lt.0.0d0) sgnc(16:16) = '-'
            if(fu(2).gt.0.0d0) sgnc(16:16) = '+'
            if(    fu(4).lt.0.0d0) then
              write(iow,2001) sgnc,abs(fu(2)),fu(3),' - ',abs(fu(4))
            elseif(fu(4).gt.0.0d0) then
              write(iow,2001) sgnc,abs(fu(2)),fu(3),' + ',abs(fu(4))
            else
              write(iow,2001) sgnc,abs(fu(2)),fu(3),' ) '
            endif
          endif
          if(fu(5).ne.0.0d0) then
            if(fu(5).lt.0.0d0) sgnc(16:16) = '-'
            if(fu(5).gt.0.0d0) sgnc(16:16) = '+'
            if(    fu(7).lt.0.0d0) then
              write(iow,2002) sgnc,abs(fu(5)),fu(6),' - ',abs(fu(7))
            elseif(fu(7).gt.0.0d0) then
              write(iow,2002) sgnc,abs(fu(5)),fu(6),' + ',abs(fu(7))
            else
              write(iow,2002) sgnc,abs(fu(5)),fu(6),' ) '
            endif
          endif
        endif

      else
        f = fu(1) + fu(2)*sin(fu(3)*pi*xi + fu(4))
     &            + fu(5)*cos(fu(6)*pi*xi + fu(7))
      endif

c     Format

2000  format(/10x,'f(s) = ',1p,1e10.3)

2001  format(a17,1p,1e10.3,'*sin(',1p,1e10.3,'*pi*s',a3:,1p,1e10.3,')')

2002  format(a17,1p,1e10.3,'*cos(',1p,1e10.3,'*pi*s',a3:,1p,1e10.3,')')

      end
