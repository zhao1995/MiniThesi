c$Id:$
      subroutine modlfd(ii,d,f,finv,df,detf,ta,hn,hn1,nh,istrt,
     &                  dd,sig,bb,xlamd,ha, bbar, isw)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Add variables to permit viscoelastic behavior    07/06/2007
c          for modmnrv (modified Mooney-Rivlin model)
c       2. Add rvemat for multi-scale analysis              13/06/2009
c       3. Put move hn -> hn1 into only the program models  14/06/2009
c       4. Add argument to detf calls                       06/10/2010
c       5. Add axisymmetric thickness plane stress          18/01/2011
c       6. Modify for cases of small deformation            21/11/2011
c       7. Send displacement gradient to revmat             10/12/2011
c       8. Add 'ii' to argument list (quadrature point)     05/01/2012
c          Add call to elhispl to set history plots.
c       9. Modify plane stress computation                  10/01/2012
c      10. Set 'finflg' true for Hill-Mandel evaluations    03/04/2012
c      11. Add 'd' to call of revmat                        10/05/2012
c      12. Increase ictmax to 10 for plane stress inters    30/08/2012
c      13. Add elaugm.h and insert xlamd and ha             23/10/2012
c      14. Add anisotropic plasticity module                05/02/2013
c      15. Pass f(1,1,1) not f(1,1,3) to neoh3f             29/06/2013
c      16. Add 'fibe'r model based on I_4                   03/08/2013
c-----[--.----+----.----+----.-----------------------------------------]
c     Purpose: Driver for finite deformation constitutive models

c     Inputs:
c       ii           - Point being processed
c       d(*)         - Material parameters array
c       f(3,3,4)     - Deformation and Displacement gradient
c       finv(3,3)    - Inverse deformation gradient
c       df(3,3)      - Incremental deformation gradient
c       detf(4)      - Determinant of deformation gradient
c       ta           - Temperature at point
c       hn(*)        - History parameters at t_n
c       nh           - Number history parameters/stress point
c       istrt        - Start state: 0 = elastic; 1 = last solution
c       xlamd        - Augmentation "penalty" value
c       bbar         - Flag (true for B-bar, false for others)
c       isw          - Solution option from elements

c     Outputs
c       hn1(*)       - History parameters at t_n+1
c       dd(*,*,5)    - Expanded material moduli for mixed computation
c                      N.B. Computed from spatial tangent, aa(6,6,5);
c                      Otherwise copy of aa.
c       sig(10)      - Cauchy stress values at t_n+1
c       bb(6)        - Left Cauchy-Green deformation tensor
c       ha           - Augmentation function
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'cdat1.h'
      include  'elaugm.h'
      include  'elcount.h'
      include  'eldata.h'
      include  'eldatp.h'
      include  'elengy.h'
      include  'elpers.h'
      include  'sdata.h'

      include  'pointer.h'
      include  'comblk.h'

      include  'iofile.h'

      logical   bbar, plasfl,viscfl, state, conv, fplps
      integer   nh,istrt,isw,jsw, ii,i,j, imat, ntm, uprm,umat, stype
      integer   c22,c33,ict,ictmax,   nstv, nprop,nvect

      real*8    ta,xlamd,ha, detf22,f33,g33,detg22,trg22
      real*8    d(*),  bb(6),  detf(*),f(3,3,*),finv(3,3),df(3,3)
      real*8    hn(*), hn1(*), dd(*), sig(*), aa(6,6,5), b2(6,2)

      save

      data       ictmax / 10 /

c     Set number of active stress components

      stype = nint(d(16))
      fplps = d(199).ne.0.0d0
      finflg = .true.

c     Initialize plane stress iterations

      if((stype.eq.1 .or. stype.eq.4 .or. fplps) .and. isw.ne.14) then
        if(fplps) then
          c33 = 2
          c22 = 3
        else
          c33 = 3
          c22 = 2
        endif
        detf22  = f(1,1,1)*f(c22,c22,1) - f(1,c22,1)*f(c22,1,1)
        detg22  = f(1,1,3)*f(c22,c22,3) - f(1,c22,3)*f(c22,1,3)
        trg22   = f(1,1,3) + f(2,2,3)
        f(c33,c33,3) = hn1(nh)
        f(c33,c33,1) = hn1(nh) + 1.0d0
c       f33     = g33 + 1.0d0
        f33     = f(c33,c33,1)
        g33     = f(c33,c33,3)
c       g33     = f33 - 1.0d0
        ict     = 0
      endif

c     Set active stress numbers

      if(ndm.eq.3 .or. stype.eq.8) then
        ntm = 6
      elseif(ndm.eq.2) then
        ntm = 4
      elseif(ndm.eq.1) then
        ntm = 3
      endif

c     Set material model type and user material pointers

      uprm  = ndd-nud
      umat  = int(d(uprm)) - 100

c     Compute Left Cauchy-Green deformation tensor

      bb(1) = f(1,1,1)*f(1,1,1) + f(1,2,1)*f(1,2,1) + f(1,3,1)*f(1,3,1)
      bb(2) = f(2,1,1)*f(2,1,1) + f(2,2,1)*f(2,2,1) + f(2,3,1)*f(2,3,1)
      bb(3) = f(3,1,1)*f(3,1,1) + f(3,2,1)*f(3,2,1) + f(3,3,1)*f(3,3,1)
      bb(4) = f(1,1,1)*f(2,1,1) + f(1,2,1)*f(2,2,1) + f(1,3,1)*f(2,3,1)
      bb(5) = f(2,1,1)*f(3,1,1) + f(2,2,1)*f(3,2,1) + f(2,3,1)*f(3,3,1)
      bb(6) = f(1,1,1)*f(3,1,1) + f(1,2,1)*f(3,2,1) + f(1,3,1)*f(3,3,1)

c     Copy for use in constitution

      do i = 1,6
        b2(i,1) = bb(i)
      end do ! i

c     Compute Left Cauchy-Green deformation tensor - 1.0d0

      b2(1,2) = f(1,1,3)*f(1,1,3)+f(1,2,3)*f(1,2,3)+f(1,3,3)*f(1,3,3)
     &        + f(1,1,3)*2.0d0
      b2(2,2) = f(2,1,3)*f(2,1,3)+f(2,2,3)*f(2,2,3)+f(2,3,3)*f(2,3,3)
     &        + f(2,2,3)*2.0d0
      b2(3,2) = f(3,1,3)*f(3,1,3)+f(3,2,3)*f(3,2,3)+f(3,3,3)*f(3,3,3)
     &        + f(3,3,3)*2.0d0
      b2(4,2) = f(1,1,3)*f(2,1,3)+f(1,2,3)*f(2,2,3)+f(1,3,3)*f(2,3,3)
     &        + f(1,2,3) + f(2,1,3)
      b2(5,2) = f(2,1,3)*f(3,1,3)+f(2,2,3)*f(3,2,3)+f(2,3,3)*f(3,3,3)
     &        + f(2,3,3) + f(3,2,3)
      b2(6,2) = f(1,1,3)*f(3,1,3)+f(1,2,3)*f(3,2,3)+f(1,3,3)*f(3,3,3)
     &        + f(3,1,3) + f(1,3,3)

c     Plane stress return point

100   continue

c     Zero stress and moduli

      do i = 1,6
        sig(i) = 0.0d0
        do j = 1,6
          aa(j,i,1) = 0.0d0
          aa(j,i,2) = 0.0d0
          aa(j,i,3) = 0.0d0
          aa(j,i,4) = 0.0d0
          aa(j,i,5) = 0.0d0
        end do ! j
      end do ! i

c     Set constant initial stress state

      if(nint(d(160)).eq.1) then
        do i = 1,6
          sig(i) = d(160+i)
        end do ! i
      endif

c     Program material models

      if(umat.lt.0) then

c       Set model type

        plasfl = nint(d(40)).eq.1 .or. nint(d(40)).eq.3
     &                            .or. nint(d(40)).eq.6
        viscfl = nint(d(40)).eq.2

        imat = nint(d(20))

c       COMPUTE STRESS AND TANGENTS

c       Standard neo-Hookean elastic model

        if    (imat.eq.1) then

          call stnh3f(d,detf(1), b2(1,2), sig,aa,xlamd,ha,estore)

c       Modified neo-Hookean model

        elseif(imat.eq.2) then

          do i = 1,nh
            hn1(i) = hn(i)
          end do ! i

          call neoh3f(d,f,finv,detf(1), b2(1,2), hn1,ntm,
     &                 sig,aa, xlamd,ha,estore)

c       Ogden model

        elseif(imat.eq.3) then

          do i = 1,nh
            hn1(i) = hn(i)
          end do ! i
          jsw = nint(d(170))
          call nalp3f(d,f,finv,detf(1),b2(1,2), hn1,ntm, sig,aa,
     &                estore, 1,jsw)

c       Finite stretch plasticity model

        elseif(imat.eq.4) then

          do i = 1,nh
            hn1(i) = hn(i)
          end do ! i

          if(plasfl) then

c           Isotropic plasticity

            call plasfd(d,detf(1),f, hn1(1),hn1(3),hn1(9),
     &                  ntm,istrt, aa,sig,isw,state)
          else
            jsw = nint(d(170))
            call nalp3f(d,f,finv,detf(1),b2(1,2), hn1,ntm, sig,aa,
     &                  estore, 2,jsw)
          endif

c       Orthotropic models

        elseif(imat.eq.5 .or. imat.eq.6) then

          if(plasfl) then

c           Orthotropic plasticity

            do i = 1,nh
              hn1(i) = hn(i)
            end do ! i

            call modlf_apl(d,detf(1),f, hn1, istrt, aa,sig,isw)

          else

c           Saint-Venant-Kirchhoff (has energy conserving capability)

            call stvk(d,detf(1),f(1,1,3),df,sig,aa, estore)
          endif

c       Fung Pseudo-exponential model

        elseif(imat.eq.7) then

          call pfung(d,f(1,1,3),detf(1),sig,aa, estore)

c       Mooney-Rivlin stress and tangents

        elseif(imat.eq.9) then

          call mnrv3f(d,detf(1),b2(1,1),b2(1,2), sig,aa,xlamd,ha,estore)

c       Modified Mooney-Rivlin stress and tangents

        elseif(imat.eq.10) then

          do i = 1,nh
            hn1(i) = hn(i)
          end do ! i
          call modmnrv(d,f,finv,detf(1),b2, hn1,ntm, sig,aa,xlamd,ha,
     &                 estore)

c       Arruda-Boyce model

        elseif(imat.eq.11) then

          call arruda(d,detf(1),b2(1,2), sig,aa,xlamd,ha,estore)

c       Yeoh model

        elseif(imat.eq.12) then

          call yeoh3f(d,detf(1),b2(1,2), sig,aa,xlamd,ha,estore)

c       RVE model from multi-scale analysis

        elseif(imat.eq.13) then

          call rvemat(d,f(1,1,3),detf(1),ta,hn,hn1,nh, sig,aa, isw)

c       Fiber models based on I_4

        elseif(imat.eq.14) then

          sig(1:6)  = 0.0d0
          aa(:,:,1) = 0.0d0
          nstv = nint(d(260))

          do i = 1,nstv
            nprop = 268 + 3*nstv
            nvect = 258 + 3*nstv
            call pfiber(d(nprop),d(nvect),f(1,1,3),detf(1), sig, aa, 3)
          end do ! i

        endif

c     U s e r    M o d e l    I n t e r f a c e

      elseif(umat.le.100) then

        elamd = xlamd
        eha   = ha
        call umodel(umat,f,detf(1),ta,d(1),d(uprm+1),hn(1),hn1(1),nh,
     &              ii,istrt, sig,aa, isw)
        xlamd = elamd
        ha    = eha

c     U s e r    F i n i t e    M o d e l    I n t e r f a c e

      else

        umat = umat - 100
        call umodelf(umat,f,finv,df,detf(1),b2,ta,d(1),d(uprm+1),
     &               hn(1),hn1(1),nh,ntm,ii,istrt, sig,aa,xlamd,ha,isw)

      endif

c     Plane stress modification

      if((stype.eq.1 .or. stype.eq.4 .or. fplps) .and. isw.ne.14) then
        ict       = ict + 1
        call fpstrs(sig,aa,g33,f33,c33,ict,ictmax, conv,.true.)
        f(c33,c33,3)  = g33
        f(c33,c33,1)  = f33
        bb(c33)       = f33*f33
        b2(c33,1)     = f33*f33
        b2(c33,2)     = g33*g33 + 2.0d0*g33
        finv(c33,c33) = 1.d0/f33
        detf(1)       = detf22*f33
        detf(3)       = g33 + trg22 + g33*trg22 + f33*detg22
        if(.not.conv .and. ict.lt.ictmax) go to 100
        hn1(nh)       = g33
      endif

c     Count computation for models

      if(umat.lt.0) then

c       Plastic update

        if(plasfl) then
          if(state) then
            nomats(2,2) = nomats(2,2) + 1
          else
            nomats(1,2) = nomats(1,2) + 1
          endif

c       Viscoelastic update

        elseif(viscfl) then

          nomats(1,4) = nomats(1,4) + 1

c       Elastic update

        else
          nomats(1,1) = nomats(1,1) + 1
        endif
      endif

c     Project aa to D-matrix for B-bar

      if(bbar) then

        call dmatmx ( aa, dd )

      else

        call pmove  ( aa, dd, 36 )

      end if

c     Set history nodal data if requested

      if(histpltfl .and. isw.ne.14) then
        call elhispl(mr(np(303)),hn1, ma,nh, ii)
      endif

      end
