c$Id:$
      subroutine ref_cyl(d,eps,sig,dd, xa, iopt)

      implicit   none

      include   'refnd.h'

      integer    iopt
      real*8     d(*)
      real*8     eps(6),sig(6),dd(6,6), xa(3)

      real*8     vlen
      real*8     tt(3,3),v1(3),v2(3),v3(3)

      real*8     ee(6),ss(6), dm(6,6), tm(6,6)


      aref = nint(d(250))
      if(aref.eq.1) then

        refa(:,1) = d(251:253)
        refa(:,2) = d(254:256)
        v3(:) = refa(:,2) - refa(:,1)
        vlen  = 1.0d0/sqrt(v3(1)**2 + v3(2)**2 + v3(3)**2)
        v3(:) = v3(:)*vlen

        v1(:) = xa(:) - refa(:,1)

        v2(1) = v3(2)*v1(3) - v3(3)*v1(2)
        v2(2) = v3(3)*v1(1) - v3(1)*v1(3)
        v2(3) = v3(1)*v1(2) - v3(2)*v1(1)
        vlen  = 1.0d0/sqrt(v2(1)**2 + v2(2)**2 + v2(3)**2)
        v2(:) = v2(:)*vlen

        v1(1) = v2(2)*v3(3) - v2(3)*v3(2)
        v1(2) = v2(3)*v3(1) - v2(1)*v3(3)
        v1(3) = v2(1)*v3(2) - v2(2)*v3(1)

        tt(:,1) = v1(:)
        tt(:,2) = v2(:)
        tt(:,3) = v3(:)
      endif
c     Transform strain

      if(iopt.eq.1) then

        ee(:)    = eps(:)
        ee(4:6)  = eps(4:6)*0.5d0
        call pushr2(tt,ee,eps,1.0d0)
        eps(4:6) = eps(4:6)*2.d0

      elseif(iopt.eq.2) then

c       Transform stress

        ss(:) = sig(:)
        call pushr2(tt,ss,sig,1.0d0)

c       Transform moduli

        dm    = dd
        call tranr4(tt,tt,tm,.false.)
        call pushr4(tm,tm,dm,dd,1.0d0)

      endif

      end

