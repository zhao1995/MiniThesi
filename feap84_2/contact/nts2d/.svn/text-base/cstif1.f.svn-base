c$Id:$
      subroutine cstif1 (ch2,ch3,nsiz,tanm,resv)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Remove 'cp0' from argument list                  21/04/2007
c-----[--.----+----.----+----.-----------------------------------------]
c      Rel. 1.0 - G.Z. - March 15, 1996

c      Acronym: Contact STIFness 1

c      Purpose: Compute tangent stiffness and residual vector

c      Inputs :
c         ch2(*)  - Contact history variables (current)
c         ch3(*)  - Contact history variables (static)

c      Outputs:
c         tanm(*) - TANgent Matrix
c         resv(*) - RESidual Vector
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'c_0.h'
      include  'c_comnd.h'
      include  'c_contac.h'
      include  'c_keyh.h'
      include  'c_pair.h'
      include  'c_tanfl.h'
      include  'counts.h'
      include  'eltran.h'

      integer   kr,kc,istgt,nsiz
      real*8    s21,c21,d21,s21c,c21c,csic,csi2c,csi,csi2
      real*8    fn,ft,gn,dfngn,dftgn,dfttd,d21o,c1,c2,c3
      real*8    gw1,gw2,gw3,gw4,gw5,gw6
      real*8    vns(6),vn0(6),vts(6),vtc(6),vnc(6),vt0(6)
      real*8    ch2(*),ch3(*),tanm(nsiz,*),resv(*)

      save

      call cdebug0 ('      cstif1',-1)

c     Get data for normal stiffness

      istgt = nint(ch2(p1(3)))
      s21   = ch2(p1(5))
      c21   = ch2(p1(6))
      d21   = ch2(p1(7))
      csi   = ch2(p1(8))
      csi2  = 1.d0 - csi
      gn    = ch2(p1(9))
      fn    = ch2(p1(51))
      dfngn = ch2(p1(52))

c     Form NS, N0, TS vectors

      call stifv1 (s21,c21,csi,csi2,vns,'NS')
      call stifv1 (s21,c21,0.d0,0.d0,vn0,'N0')
      call stifv1 (s21,c21,csi,csi2,vts,'TS')

c     Clean vectors

      do kc = 1,8
        resv(kc) = 0.0d0
        do kr = 1,8
          tanm(kr,kc) = 0.0d0
        end do ! kr
      end do ! kc

c     Check on sign of contact force if > 0 set to zero

      if(fn.gt.0.0d0 .and. ifistgn) then
c     if(fn.gt.0.0d0) then
        ch2(p1(51)) = 0.0d0
        if(ifaugm.ne.1) then
           ch2(p1(151)) = 0.0d0 ! Augmented force
        endif
        if(ifsolm.ne.2) then
           ch2(p1( 21)) = 0.0d0 ! Lagrange multiplier force
        endif
        return                  ! N.B. Assembles zero tangent/residual
      endif

c     Compute  -fn/d21 * [TSxN0 + N0xTS + (gn/l)*N0xN0]   ([1])

      if (uafl) then
        if (istgt.eq.1) then
          c1 = -fn/d21
          c2 = c1 * gn/d21
          do kr = 1, 6
            gw1 = c1 * vts(kr)
            gw2 = c1 * vn0(kr)
            gw3 = c2 * vn0(kr)

            do kc = kr, 6
              tanm(kr,kc) = gw1*vn0(kc) + gw2*vts(kc) + gw3*vn0(kc)
            end do
          end do

c         Compute (d(fn)/d(gn))* [NSxNS]

          do kr = 1, 6
            gw1 = dfngn * vns(kr)
            do kc = kr, 6
              tanm(kr,kc) = tanm(kr,kc) + gw1*vns(kc)
            end do
          end do

c         Lagrange multiplier formation

          if( ifsolm.eq.2 ) then
            do kr = 1,6
              tanm(kr,7) = tanm(kr,7) + vns(kr)
            end do ! kr
          endif

c       Contribution for special corner cases

        else
          s21c = ch2(p1(13))
          c21c = ch2(p1(14))
          csic = ch2(p1(15))
          csi2c = 1.d0 - csic
          call stifv1 (s21c,c21c,csic,csi2c,vtc,'TC')
          call stifv1 (s21c,c21c,csic,csi2c,vnc,'NC')

          s21   = ch2(p1(5))
          c21   = ch2(p1(6))
          csi   = ch2(p1(8))

c         Compute (fn/gn)*[TCxTC] + (d(fn)/d(gn))*[NCxNC]

          if (gn.ne.0.d0) then
            c1 = fn/gn
          else
            c1 = 0
          endif
          c2 = dfngn
          do kr = 1, 6
            gw1 = c1 * vtc(kr)
            gw2 = c2 * vnc(kr)
            do kc = kr, 6
              tanm(kr,kc) = gw1*vtc(kc) + gw2*vnc(kc)
            end do
          end do

c         Lagrange multiplier formation

          if( ifsolm.eq.2 ) then
            do kr = 1,6
              tanm(kr,7) = tanm(kr,7) + vnc(kr)
            end do ! kr
          endif

        endif

c       Account for symmetry
c       REMARK - stiffness matrix passed completely full

        do kr = 1, nsiz
          do kc = kr, nsiz
            tanm(kc,kr) = tanm(kr,kc)
          end do
        end do

c       Get data for tangential stiffness

        if (iffric.eq.1) then
          call stifv1 (s21,c21,0.d0,0.d0,vt0,'T0')
          ft     = ch2(p1(53))
          dfttd  = ch2(p1(54))
          dftgn  = ch2(p1(55))

c         Check tangential solution mode

          d21o  = ch3(p3(1))

c         Compute (d(ft)/d(gt)) * d21o**2/d21**2 *
c         [TSxTS + gn/d21*[N0xTS + TSxN0] + (gn**2/d21**2)*[N0xN0]]

          if (dfttd.ne.0.d0) then
            c1 = dfttd * d21o**2/d21**2
            c2 = c1 * gn/d21
            c3 = c1 * (gn/d21)**2
            do kr = 1,6
              gw1 = c1 * vts(kr)
              gw2 = c2 * vts(kr)
              gw3 = c2 * vn0(kr)
              gw4 = c3 * vn0(kr)
              do kc = 1,6
                tanm(kr,kc) = tanm(kr,kc) + gw1*vts(kc) + gw2*vn0(kc)+
     &                        gw3*vts(kc) + gw4*vn0(kc)
              end do
            end do
          endif

c         Compute d(ft)/d(gn)) * d21o/d21 * [TSxNS + gn/d21 N0xN0]

          if (dftgn.ne.0.d0) then
            c1 = dftgn * d21o/d21
            c2 = c1 * gn/d21
            do kr = 1,6
              gw1 = c1 * vts(kr)
              gw2 = c2 * vn0(kr)
              do kc = 1,6
                tanm(kr,kc) = tanm(kr,kc) + gw1*vns(kc) + gw2*vns(kc)
              end do
            end do
          endif

c         Compute ft * d21o/d21**2 * [-[T0xTs - TSxT0]
c                                    + [NSxN0 + N0xNS]
c                         - 2*gn/d21 * [T0xN0 + N0xT0]]

          if (ft .ne. 0.d0) then
            c1 = ft * d21o/d21**2
            c2 = c1 * 2.d0*gn/d21
            do kr = 1,6
              gw1 = c1 * vns(kr)
              gw2 = c1 * vn0(kr)
              gw3 = c1 * vts(kr)
              gw4 = c1 * vt0(kr)
              gw5 = c2 * vn0(kr)
              gw6 = c2 * vt0(kr)
              do kc = 1,6
                tanm(kr,kc) = tanm(kr,kc) + gw1*vn0(kc) + gw2*vns(kc)-
     &                        gw3*vt0(kc) - gw4*vts(kc) - gw5*vt0(kc)-
     &                        gw6*vn0(kc)
              end do
            end do
          endif
        endif
      endif

c     fn residual

      fn = ch3(p3(5))

c     Form residual vector

      if (dbfl) then
        if (istgt.gt.1) then
          call stifv1 (s21c,c21c,csic,csi2c,vns,'NS')
        endif

        do kc = 1, 6
          resv(kc) = -fn*vns(kc)
        end do

c       Lagrange multiplier formation

        if( ifsolm.eq.2 ) then
          resv(7) = -gn
        endif

c       Friction contribution

        if (iffric.eq.1) then

          c1 = ft * d21o/d21  ! Could be just ft
          c2 = c1 * gn/d21
          do kc =1,6
            resv(kc) = resv(kc) - c1*vts(kc) - c2*vn0(kc)
          end do
        endif
      endif

c     Multiply tangent by integration parameter

      if(ctan(1).ne.1.d0) then
        do kc = 1,nsiz
          do kr = 1,nsiz
            tanm(kr,kc) = tanm(kr,kc)*ctan(1)
          end do
        end do
      endif

      end
