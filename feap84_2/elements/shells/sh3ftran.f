c$Id:$
      subroutine sh3ftran ( dt   , shp  , aa   ,
     &                      dr   , r1   , vrn  , vr1  , arn  , ar1,
     &                      dira , dir0 , rhoa , rhoi , p    , s  )

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c        Description:    SH3FTRAN is the subroutine which computes the
c                        transient (inertial) terms of the residual
c                        vector and tangent stiffness matrix
c                        for the Energy-Momentum method

c        Authors:        N.Tarnow & J.C. Simo

c        Date:           March 1993.

c        Caution:        Must use the conservative time-stepping
c                        algorithms in conjunction with this routine,
c                        i.e., beta,cons, etc... (nop=5).
c-----[--.----+----.----+----.-----------------------------------------]
c        Routine Input:
c        --------------
c        dt ............ Time step.
c        shp ........... Nodal shape function values.
c        aa ............ Gauss point acceleration at T-n+a.
c        dr ............ Localized nodal incremental rotation vector.
c        r1 ............ Localized nodal total rotation vector & T-n+1.
c        vrn,vr1 ....... Gauss point rot. velocities at T-n and T-n+1.
c        arn,ar1 ....... Gauss point rot accelerations at T-n and T-n+1.
c-----[--.----+----.----+----.-----------------------------------------]
c       Remark: In case of shell intersections dr,r1, vrn, vr1, arn and
c               ar1 are used to store d/dt Lambda at tn and tn+1
c-----[--.----+----.----+----.-----------------------------------------]
c        dira .......... Localized nodal director
c                        vectors, dir(3,ElementNode) at time T-n+a.
c        dir0 .......... Localized nodal director
c                        vectors, dir(3,ElementNode) at time T-0.
c        rhoa .......... Trans. thickness weighted density multiplyied
c                        by the jacobian and weight at the Gauss point.
c        rhoi .......... Rot. inertia weighted density multiplyied
c                        by the jacobian and weight at the Gauss point.

c        Routine Output:
c        ---------------
c        s ............. Element stiffness.
c        p ............. Element residual.
c-----[--.----+----.----+----.-----------------------------------------]

      implicit  none

      include  'crotas.h'

      integer   j,   k, l, j1, k1    , m
      real*8    dt    , dtr  , dtr4
      real*8    rhoa  , rhoi , fac   , fac1 , fac2

      real*8    w1(3,3)      , w0(3,3)      , shp  (4)     , aa   (3)
      real*8    dr   (3,9)   , r1   (3,9)   , vrn  (3,9)   , vr1  (3,9)
      real*8    arn  (3,9)   , ar1  (3,9)   , wd   (3)
      real*8    dira (3,9)   , dir0(3,9)
      real*8    s    (24,24) , p    (24)

      save

c     Compute director velocities at shell intersections

      dtr  = 1.d0/dt
      dtr4 = 4.d0*dtr*dtr

      do j=1,4
        if(ndof(j).eq.6)then
          do k=1,3
            w1(k,1) = dr(k,j)
            w1(k,2) = vr1(k,j)
            w1(k,3) = ar1(k,j)
            w0(k,1) = r1(k,j)
            w0(k,2) = vrn(k,j)
            w0(k,3) = arn(k,j)
          end do ! k
          do k=1,3
            vr1(k,j) = 0.0d0
            vrn(k,j) = 0.0d0
            do l=1,3
               vr1(k,j) = vr1(k,j) + w1(k,l)*dir0(l,j)
               vrn(k,j) = vrn(k,j) + w0(k,l)*dir0(l,j)
            end do ! l
          end do ! k
        endif
      end do ! j

c     Node Loop 1

      j1 = 0
      do j = 1 , 4

c       Multiply Delta.w * t-n+a

        do k  = 1 , 3
          wd(k) = (vr1(k,j)-vrn(k,j))*dtr
        end do ! k

c       Set Up Residual

        do k = 1 , 3
          p(j1+k  ) = p(j1+k  ) - shp(j) * rhoa * aa(k)
          p(j1+k+3) = p(j1+k+3) - shp(j) * rhoi * wd(k)
        end do ! k

        j1 = j1 + 6

c     End Node Loop 1

      end do ! j

c     Node Loop 2

      j1 = 0
      do j = 1 , 4
        fac = rhoi*shp(j)*dtr
        do l = 1 , 3
          wd(l) = fac*(vr1(l,j) - vrn(l,j))
        end do ! l
        fac1 = wd(1)*dira(1,j) + wd(2)*dira(2,j) + wd(3)*dira(3,j)

c       Final Displacement Stiffness

        fac = rhoa*shp(j) * dtr4
        k1 = 0
        do k = 1 , 4
          fac2 = fac*shp(k) - fac1
          do l = 1 , 3
            s(j1+l,k1+l) = s(j1+l,k1+l) + fac2
          end do ! l
          k1 = k1 + 6
        end do ! k

c       Final Rotational Stiffness

        fac = rhoi*shp(j) * dtr4
        do m = 4 , 6
          s(j1+m,j1+m) = s(j1+m,j1+m) + fac
        end do ! m

        j1 = j1 + 6

c     End Node Loop 2

      end do ! j

      end
