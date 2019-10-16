c$Id:$
      subroutine modlsd(ii,d,ta,eps,h1,h2,nh,istrt, dd,sig,alam,ha,isw)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved
c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Add augmented parameters 'alam' and 'ha';        20/07/2007
c          Update stress for augments (except for plane stress)
c       2. Add axisymmetric with thickness plane stress     18/01/2011
c       3. Add 'ii' to modlsd argument                      05/01/2012
c       4. Add g33 to fpstrs call                           10/01/2012
c       5. Add 'd' to call of srvemat                       10/05/2012
c       6. Subtract volumetric thermal strain value to 'ha' 10/07/2012
c       7. Add elaugm.h and insert alam and ha              23/10/2012
c       8. Add anisotropic plasticity module                05/02/2013
c-----[--.----+----.----+----.-----------------------------------------]
c     Purpose: Small Deformation Constitutive Equation Driver

c     Input parameters
c          ii        -  Point being processed
c          d(*)      -  up to ndd-nud-1 material parameters
c          ta        -  Temperature change form stress free state
c          eps(9,3)  -  Current strains and fluxes at point
c          h(nh)     -  History terms at point
c          nh        -  Number of history terms
c          istrt     -  Start state: 0 = elastic; 1 = last solution
c          alam      -  Augmented multiplier
c          im        -  Material type
c     Ouput parameters
c          dd(6,6,5) -  Current material tangent moduli
c                       Coupled problems:    | dd_1   dd_2 |
c                                            | dd_3   dd_4 |
c                       Rayleigh damping: dd_5
c          sig(6)    -  Stresses at point.
c          ha        -  Augmented function (trace-eps)
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'cdat1.h'
      include  'complx.h'
      include  'elaugm.h'
      include  'eldata.h'
      include  'eldatp.h'
      include  'elcount.h'
      include  'pmod2d.h'
      include  'sdata.h'

      include  'pointer.h'
      include  'comblk.h'

      logical   state, conv, strfl, pflps
      integer   ii,it,i,j,nh,istrt,isw, c33,ict, ntm, umat,uprm
      real*8    alam, ha, ta,e33,g33, epsth
      real*8    d(*),eps(9,*),h1(*),h2(*),dd(6,6,5),sig(*),theta(3)

      save

c     Extract analysis type: 1=plane stress; 2=plane strain; 3=axisy;
c                            8=axisy+torsion

      stype = nint(d(16))
      pflps = d(199).ne.0.0d0
      if(stype.eq.1) then
        e33 = eps(3,1)
        c33 = 3
        ict = 0
      elseif(pflps) then
        e33 = eps(2,1)
        c33 = 2
        ict = 0
      endif
      g33 = e33
      strfl = .true.

c     Check for user model

      uprm  = ndd-nud
      umat  = int(d(uprm)) - 100

c     Set number stress/strain components

      if(ndm.eq.3 .or. stype.eq.8) then
        ntm = 6
      elseif(ndm.eq.2) then
        ntm = 4
      else
        ntm = 1
      end if

c     Plane stress return point

100   continue

c     Zero stress and dd arrays

      do i = 1,6
        sig(i) = 0.0d0
        do j = 1,6
          dd(j,i,1) = 0.0d0
          dd(j,i,2) = 0.0d0
          dd(j,i,3) = 0.0d0
          dd(j,i,4) = 0.0d0
          dd(j,i,5) = 0.0d0
        end do ! j
      end do ! i
      epsth = 0.0d0 ! Thermal strain for augmenting

c     Set constant initial stresses

      if(nint(d(160)).eq.1) then
        do i = 1,6
          sig(i) = d(160+i)
        end do !
      end if

c     Program material models

      if(umat.lt.0) then

c       Set model type

        plasfl = nint(d(40)).eq.1 .or. nint(d(40)).eq.3
     &                            .or. nint(d(40)).eq.6
        viscfl = nint(d(40)).eq.2

c       P l a s t i c i t y

        if(plasfl) then

          if(nint(d(40)).eq.1) then

c           Plane stress plasticity

            if (stype.eq.1 .or.stype.eq.4) then

c             Move h1 to h2

              do i = 1,nh
                h2(i) = h1(i)
              end do ! i

              call epps2d(d,eps,h2,h2(4),h2(7),istrt, sig,dd,dd(1,1,5))
              state = h2(8).eq.0.0d0
              it    = 2
              strfl = .false.

c           Plane strain, axisymmetric: 3D plastic or viscoplastic

            else

c             Move h1 to h2

              do i = 1,nh
                h2(i) = h1(i)
              end do ! i

              call mises(d,eps,h2(3),h2,ntm,istrt, sig,dd,dd(1,1,5))
              state = h2(2).eq.0.0d0
              it    = 2
              if(.not.state .and. isw.eq.8) sig(10) = h2(1)

            endif

c         G e n e r a l i z e d    P l a s t i c i t y

          elseif(nint(d(40)).eq.3) then

c           Move h1 to h2

            do i = 1,nh
              h2(i) = h1(i)
            end do ! i

            call gplas3d(d,eps,h2(4),h2,ntm,istrt, sig,dd,dd(1,1,5))
            state = h2(2).eq.0.0d0
            it    = 3

c         A n i s o t r o p i c    P l a s t i c i t y

          elseif(nint(d(40)).eq.6) then

c           Move h1 to h2

            do i = 1,nh
              h2(i) = h1(i)
            end do ! i

            call plast_apl(d,eps,h2,nh, istrt, sig,dd,isw)

          endif

c       Representative volume element model (two-scale solutions)

        elseif(nint(d(40)).eq.4) then

          call srvemat(d,eps,ta,h1,h2,nh, sig,dd, isw)

c       E l a s t i c i t y

        else

c         Piezoelectric case

          if(d(150).gt.0.0d0) then

            call pzstrs(d,eps(7,1),eps,sig,dd)

c         Thermoelastic case

          else

            call estrsd(d,ta,eps,sig,dd,dd(1,1,5))
            epsth       = d(3)*ta
            nomats(1,1) = nomats(1,1) + 1

          endif

        end if

c       V i s c o e l a s t i c i t y

        if(viscfl) then

c         Complex modulus form

          if(cplxfl) then
            call cvisco(d,eps,sig,dd,dd(1,1,3))

c         Time moduli integration

          else

c           Move h1 to h2

            do i = 1,nh
              h2(i) = h1(i)
            end do ! i

            call viscoe(d,ta,eps,h2(1),h2(ntm+1),ntm,sig,dd,dd(1,1,5))
            epsth = d(3)*ta
          endif

        end if

c     U s e r    M o d e l    I n t e r f a c e

      else

c       Compute trace to pass to user module

        theta(1) = eps(1,1) + eps(2,1) + eps(3,1)
        theta(2) = eps(1,2) + eps(2,2) + eps(3,2)
        theta(3) = eps(1,3) + eps(2,3) + eps(3,3)

        elamd = alam
        eha   = ha
        call umodel(umat,eps,theta,ta,d(1),d(uprm+1),h1(1),h2(1),nh,
     &              ii,istrt,sig,dd,isw)
        alam  = elamd
        ha    = eha

      end if

c     Plane stress modification

      if((stype.eq.1 .or. pflps) .and. strfl) then
        ict        = ict + 1
        call fpstrs(sig,dd,g33,e33,c33,ict,6, conv,.false.)
        eps(c33,1) = e33
        if(.not.conv .and. ict.lt.6) go to 100
      endif

c     User model

      if(umat.gt.0) then ! Currently does nothing

c     Plastic update

      elseif(plasfl) then
        if(state) then
          nomats(2,it) = nomats(2,it) + 1
        else
          nomats(1,it) = nomats(1,it) + 1
        endif

c     Viscoelastic update

      elseif(viscfl) then

        nomats(1,4) = nomats(1,4) + 1

c     Piezoelectric case

      elseif(d(150).gt.0.0d0) then
        nomats(1,5) = nomats(1,5) + 1

c     Elastic update

      else
        nomats(1,1) = nomats(1,1) + 1
      endif

c     Compute augmented function (trace-epsilon = 0)

      if(stype.gt.1) then
        ha     = eps(1,1) + eps(2,1) + eps(3,1) - epsth
        sig(1) = sig(1)   + alam
        sig(2) = sig(2)   + alam
        sig(3) = sig(3)   + alam
      endif

c     Set history nodal data if requested

      if(histpltfl) then
        call elhispl(mr(np(303)),h2, ma,nh, ii)
      endif

      end
