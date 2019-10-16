c$Id:$
      subroutine penares (ch2,ch3,resv)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Remove 'cp0' from routine -- unused              21/04/2007
c       2. Remove 'resfl' from argument; fn = ch3(p3(5))    31/01/2013
c          Add computation of resv(7) for lagrange mult soln
c-----[--.----+----.----+----.-----------------------------------------]
c      Rel. 1.0 - G.Z. - March 15, 1996

c      Acronym: Contact PENAlty RESidual

c      Purpose: Compute residual vector

c      Inputs :
c         ch2(*)  - Contact history variables (current)
c         ch3(*)  - Contact history variables (static)

c      Outputs:
c         resv(*) - RESidual Vector
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'c_0.h'
      include  'c_comnd.h'
      include  'c_keyh.h'
      include  'c_pair.h'
      include  'c_tanfl.h'
      include  'counts.h'
      include  'eltran.h'

      integer   kc,istgt
      real*8    ch2(*),ch3(*),resv(*)
      real*8    vns(6),vn0(6),vts(6)
      real*8    s21,c21,d21,d21o,csi,csi2,csibb,gn,fn,ft
      real*8    s21c,c21c,csic,csi2c,c1,c2

      save

      call cdebug0 ('      penares',-1)

c     Get data for normal stiffness

      istgt = nint(ch2(p1(3)))
      s21   = ch2(p1(5))
      c21   = ch2(p1(6))
      d21   = ch2(p1(7))
      csi   = ch2(p1(8))
      csi2  = 1.d0 - csi

c     Form NS, N0, TS vectors

      call stifv1 (s21,c21,csi,csi2,vns,'NS')

c     Form residual vector

      if (istgt.gt.1) then
        s21c  = ch2(p1(13))
        c21c  = ch2(p1(14))
        csic  = ch2(p1(15))
        csi2c = 1.d0 - csic
        call stifv1 (s21c,c21c,csic,csi2c,vns,'NS')
      endif

      fn = ch3(p3(5))

      do kc = 1, 6
        resv(kc) = -fn*vns(kc)
      end do

c     Extract gap

      gn    = ch2(p1(9))

      if( ifsolm.eq.2 ) then
        resv(7) = -gn
      endif

c     Friction contribution

      if (iffric.eq.1) then
        d21   = ch2(p1(7))
        csibb = ch2(p1(17))
        d21o  = ch3(p3(1))
        ft = 1.d0
        call stifv1 (s21,c21,0.d0,0.d0,vn0,'N0')
        call stifv1 (s21,c21,csi,csi2,vts,'TS')

        c1 = ft * d21o/d21
        c2 = c1 * gn/d21
        do kc =1,6
          resv(kc) = resv(kc) - c1*vts(kc) - c2*vn0(kc)
        end do
      endif

      end
