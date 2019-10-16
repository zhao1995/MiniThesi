c$Id:$
      subroutine fld3d3(d,ul,xl,s,p,ndf,ndm,nst,isw, finc)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Correct augmented update                         26/01/2007
c       2. Revise augment: augfp -> augf                    14/04/2009
c       3. Add prints for Green and Almansi strains         19/10/2009
c-----[--.----+----.----+----.-----------------------------------------]
c     Mechanical element:
c          8-Node Enhaced Strain Trilinear Element (idinc = 0)
c          8-Node Displacement Formulation         (idinc = 1)

c     Thermal element:
c          8-Node Displacement Formulation

c     NOTE: Recommended adiabatic split is (ARMERO & SIMO [1991,1992])
c           iops = 3
c           1991-1996 F. ARMERO, JCS & RLT

c                   D A T A   B A S E   M A P

c      Transfer variables:
c      ------------------
c          hr(nh1   )       :   9 Jacobians
c          hr(nh1+9 )       :   9 deformation gradient.
c                              9 deformation gradient inverse at nh2.

c          total # of words:   nhf = 9
c                              nhl = 9 + nhf*9

c      Thermal variables:
c      -----------------
c          hr(nhg   )       :   8 Heat flux residual vector
c                                (for Crank-Nickolson integration)

c        Each gauss point 1
c          hr(nhg+8 )       :   1 Temperature
c          hr(nhg+9 )       :   1 source term
c                                (for general. integration)

c          total # of words:   nht = 2
c                              nhj = 8 + nht*9

c      Mechanical variables:
c      --------------------
c        Enhanced strain parameters

c          hr(nhg   )       :  12 Internal parameters       (ui)
c          total # of words:   nhi = 12

c        Each gauss point 1 (imat = 2)
c          hr(nhg   )       :   1  Augmented parameter - NOT IMPLEMENTED
c          hr(nhg+1 )       :   1  Damage accumulation variable (XI)
c          total # of words:   nhv = 2

c        Each gauss point 1 (imat = 3)

c           hr(nn+1)           :  Not used
c           hr(nn+2)           :  Jacobian
c           hr(nn+3)           :  Damage accumulation variable   (XI)
c           hr(nn+4)-hr(nh+9)   :  2nd Piola-Kirchhoff stresses   (Sn)
c           hr(nn+10)-hr(nn+15) :  Viscoelastic-term 1 (if used)  (Hn)
c           hr(nn+16)-hr(nn+21) :  Viscoelastic-term 2 (if used)  (Hn)
c           hr(nn+22)-hr(nn+27) :  Viscoelastic-term 3 (if used)  (Hn)

c           m   = 0 for nv = 0 and m=1 for nv > 0
c           nhv = 2 + 1 + m * (6 + nv * 6)

c        Each gauss point 1 (imat = 1 or 4)

c          hr(nhg   )       :   1  Augmented parameter - NOT IMPLEMENTED
c          hr(nhg+1 )       :   1  Damage accumulation variable (XI)
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'augdat.h'
      include  'hdata.h'
      include  'iofile.h'
      include  'prstrs.h'
      include  'part0.h'
      include  'part1.h'
      include  'comblk.h'

      logical   ther,mech,finc
      integer   ndf,ndm,nst,isw
      real*8    d(*),ul(ndf,*),xl(ndm,*),s(nst,*),p(*)

      save

c     Augmented lagrangian update

      if(isw.eq.10) then

        hr(nh2) = hr(nh2) + augf*d(185)*hr(nh2+1)

c     Set initial values for plastic strains

      elseif(isw.eq.12) then

c     Tangent Stiffness and Residual Vector

      else

c       DECIDE ACTUAL SOLUTION PHASE

        mech = npart.eq.ndfp(1)
        ther = npart.eq.ndfp(4)

        if(isw.eq.4 .or. isw.eq.8) then
          mech = .true.
          ther =  ndf.gt.ndm
        endif

c       HEAT CONDUCTION SUB-ELEMENT

        if(ther) then
c         call  ther3f (d,ul,xl,tl,s,p,ndf,ndm,nst,isw,
c    &                     hr(nh1),hr(nh2),hr(nh1+9),hr(nh2+9))
        endif

c       MECHANICAL SUB-ELEMENT

        if( mech ) then
          call mech3f (d,ul,xl,s,p,ndf,ndm,nst,isw, finc)
        endif

      endif

      end
