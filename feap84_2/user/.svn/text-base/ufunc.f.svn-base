c$Id:$
      function ufunc(fnum,fu, xi, fsw)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Dimension fu(14)                                 11/06/2011
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: User function driver
c      Inputs:
c        fnum     - Function number
c        fu(14)   - Parameters of function
c        xi       - Parameter value:  -1 < xi < 1
c        fsw      - Function switch: 0 - Output type; 1 - compute ufunc

c      Outputs:
c        ufunc    - Function value
c-----[--.----+----.----+----.-----------------------------------------]
      implicit   none

      integer    fnum, fsw
      real*8     ufunc, xi, fu(14)

      ufunc = 0.0d0
      if(fnum.eq.1) then
        call ufunc01(fu, xi, ufunc, fsw)
      elseif(fnum.eq.2) then
        call ufunc02(fu, xi, ufunc, fsw)
      elseif(fnum.eq.3) then
        call ufunc03(fu, xi, ufunc, fsw)
      elseif(fnum.eq.4) then
        call ufunc04(fu, xi, ufunc, fsw)
      elseif(fnum.eq.5) then
        call ufunc05(fu, xi, ufunc, fsw)
      elseif(fnum.eq.6) then
        call ufunc06(fu, xi, ufunc, fsw)
      elseif(fnum.eq.7) then
        call ufunc07(fu, xi, ufunc, fsw)
      elseif(fnum.eq.8) then
        call ufunc08(fu, xi, ufunc, fsw)
      elseif(fnum.eq.9) then
        call ufunc09(fu, xi, ufunc, fsw)
      elseif(fnum.eq.10) then
        call ufunc10(fu, xi, ufunc, fsw)
      endif

      end
