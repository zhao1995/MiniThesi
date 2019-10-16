c$Id:$
      subroutine viscfd(d,fi,finv,xj,xin,sn,hn,ntm, tau,aa)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c     1. Use constants from 'pconstant.h'                   14/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose:  Viscoelastic/damage for 3-d finite deformation

c      Inputs:

c      Outputs:
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'elengy.h'
      include  'pconstant.h'
      include  'tdata.h'

      logical   load
      integer   i, j, n, nv, ntm
      real*8    xj,xin,ddam,gg,xji23,trsmh,gammai,gexpi,wengy,damg
      real*8    d(*), fi(3,3),finv(3,3), sn(*),hn(ntm,*),tau(*),aa(6,*)
      real*8    tau0(6),PI2(6),smh(6),bigh(6),one(6),snold(6),ht(6)

      save

      data      one/3*1.d0,3*0.d0/

c     Save elastic stress for updates

      tau0(5) = 0.0d0
      tau0(6) = 0.0d0

      do i = 1,ntm
        tau0(i) = tau(i)
      end do ! i
      damg = 1.0d0
      ddam = 1.0d0

c     Compute damage parameters

      load = .false.
      if(d(59).ne.0.0d0) then

        wengy = estore
        call edam3f(d(58),wengy,xin,damg,ddam,load)

c       Set elastic tangent parameter and update tau for damage

        gg     = 1.0d0
        do i = 1,ntm
          tau0(i) = damg*tau0(i)
          tau(i)  = tau0(i)
        end do ! i

c       Modify elastic tangent by damage parameter

        do i = 1,ntm
          do j = 1,ntm
            aa(i,j) = damg*aa(i,j)
          end do ! j
        end do ! i
      endif

c     Compute viscoelastic parameters

      nv   = nint(d(57))
      if(nv.gt.0) then

c       Jn**(-2/3)

        xji23 = xj**(-two3)

c       Compute PI-n+1 and update for any damage

        call pushr2(finv, tau0, PI2, 1.0d0)

c       Update history parameters for 2nd Piola-Kirchhoff stress

        do i = 1,ntm
          snold(i) = sn(i)
          sn(i)    = PI2(i)/xji23
          bigh(i)  = 0.0d0
        end do ! i

c       Sum over series of terms in viscoelasticity

        gg     = 1.0d0
        do n = 1,nv
          gammai = d(2*n+49)
          gexpi  = exp(- 0.5d0*dt/d(2*n+50) )

c         Update history variables for viscoelasticity

          do i = 1,ntm
            ht(i)   = gexpi*(gexpi*hn(i,n) - snold(i))
            hn(i,n) = ht(i) + gexpi*sn(i)

c           Aaccumulate: gammai*bigh-tilde and tangent factor

            bigh(i) = bigh(i) + gammai*ht(i)
          end do ! i
          gg     = gg - gammai*(1.d0 - gexpi)
        end do ! n

c       Push forward to compute stress deviator

        call pushr2(fi, bigh, smh, 1.0d0)
        do i = 1,ntm
          smh(i) = xji23*smh(i)
        end do ! i
        if(ntm.eq.1) then
          trsmh = 0.0d0
        else
          trsmh = (smh(1) + smh(2) + smh(3))*one3
          do i = 1,3
            smh(i) = smh(i) - trsmh
          end do ! i
        endif

c       Update Kirchhoff stress

        do i = 1,ntm
          tau(i) = gg*tau0(i) + smh(i)
        end do ! i

c       Update deviatoric tangent

        do i = 1,ntm
          do j = 1,ntm
            aa(i,j) = gg*aa(i,j)
     1      - two3*((smh(i)+trsmh*one(i))*one(j) + one(i)*smh(j))
          end do ! j
          aa(i,i) = aa(i,i) + trsmh*(1.d0 + one(i))
        end do ! i
      endif

c     If damage loading add final elastic/viscoelastic tangent update

      if(load) then

c       Normalize damage function derivative by current damage parameter
c       and damage function squared (to compensate for stress below)

        ddam = ddam/(xin*damg**2)
        ddam = gg*ddam
        do i = 1,ntm
          gexpi = ddam*tau0(i)
          do j = 1,ntm
            aa(i,j) = aa(i,j) + gexpi*tau0(j)
          end do ! j
        end do ! i
      endif

      end
