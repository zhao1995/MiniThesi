c$Id:$
      subroutine inmate(d,tdof,ndv,eltype)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Use 'pi' from 'pconstant.h'                      14/11/2006
c       2. Add option for Model 4 for finite volume fn.     20/12/2006
c       3. Add 'augm,<on,off>' option                       14/03/2007
c       4. Add 'load radial' option                         18/05/2007
c       5. Set bulk modulus to zero for incompressible      20/07/2007
c          non-linear isotropic materials
c       6. Allow: 'start <elas,inel,plas>' options          12/10/2007
c       7. Change 'text(1)' to 'text(2)' for 'acti'         06/11/2007
c       8. Add options 'unif' or 'stab' to 'mixed' for      17/01/2008
c          etype = 8
c       9. Add options 'unif' or 'stab' to 'disp' for       23/01/2008
c          etype = 7
c      10. Set user material d(15) for history variables    20/07/2008
c      11. Set d(188) to indicate multidimensional frame    30/10/2008
c          constitution.
c      12. Add 'kappa' option for frame, plate and shells   31/10/2008
c      13. Add history storage for multi-dimensional        04/11/2008
c          reduction case.
c      14. Add d(189) for nurb indicator                    14/11/2008
c      15. Set section history for multi-dimensional case   15/11/2008
c          Increase maximum quadrature: 10x10 for rectangle
c          Revise prints on quadrature, etc.
c      16. Add check on positive rectangle area for SECT    12/12/2008
c          Revise format 2013 to 25x
c      17. Allow for penalty to be set on 'unif' option     02/02/2009
c      18. Add thermal softening values and print thermal   22/02/2009
c          expansion for finite deformations.
c      19. Add 'nu31' to write of 2078 format               19/04/2009
c      20. Add 'rve' (imat = 13) for multi-scale solutions  13/06/2009
c      21. Set history storage for NURB type elements       04/08/2009
c      22. Multiply finite deformation yield and hardening  10/08/2009
c          by sqrt(2/3) and 2/3 factors.
c      23. Remove text on frame inputs of elastic material  08/09/2009
c      24. Add 'plot <on,off>' for projections              28/09/2009
c      25. Add 'trve' thermal multiscale set; 'four heat c' 28/11/2009
c          or 'four spec  c'
c      26. Add 'tmat' for material model numbers in pure    14/01/2010
c          thermal problems
c      27. Return 'ietype' in d(194) for coupled problems   15/01/2010
c      28. Add use of global grou(p or nd) factors and      06/03/2010
c          proportional loads
c      29. Add finite deformation generalize plasticity     04/04/2010
c      30. Add 'B-spline' option: d(189) = 3                21/04/2010
c      31. Correct spelling in 4005; add compliance prints  14/08/2010
c      32. Remove print of d(37) for shear beam, reorder    24/09/2010
c          output of beam element descriptors.
c      33. Output normal moduli for orthotropic materials   08/10/2010
c      34. Add 't-sp'line to 'tspl'ine for quadratures      28/12/2010
c      35. Add axisymmetric 1-d thickness plane stress      18/01/2011
c      36. Add solution type to Ogden and Princ. Stretch    05/04/2011
c      37. Use 'nh2' to add extra quad pt for enhanced      17/04/2011
c      38. Set 'imat = 0' when user model is specified      02/07/2011
c      39. Add 'Bezier' option: d(189) = 5                  13/07/2011
c      40. Add set of reference vector for shells           19/07/2011
c      41. Add d(241) to set max number dofs/node           08/10/2011
c      42. Add d(240) to set element or nodal based loops   30/10/2011
c      43. Add input of history plot tables                 05/01/2012
c      44. Add storage for plane stress thickness strain    10/01/2012
c      45. Set d(170) to 5 for incompressible solutions     16/02/2012
c      46. Set prtyp(2,ma) for FE2 problems                 06/05/2012
c      47. Set ietype15 for both solids and shells          14/05/2012
c          Set default finite volume type for all.
c      48. Remove 'trve' option (only 'rve')                18/05/2012
c      49. Add 'H-Spline' option: d(189) = 7                08/11/2012
c      50. Add 'mixed u-p' nodal form                       04/02/2013
c      51. Add anisotropic plasticity modules               05/02/2013
c      52. Allow for isotropic + anisotropic plasticity     13/03/2013
c      53. For 'fini'te command leave shell linear model    29/03/2013
c      54. Set d(242) = 1 for undefined plastic vectors     25/04/2013
c          Set nh1 using nhd not 5 for small strain
c      55. Change 'axial' to 'cylindrical' on format 2079   06/05/2013
c          Add aref, refa(3,2) for cylindrical transforms
c      56. Allow for quadrature set of thermal problems     08/05/2013
c      57. Add element formulation to finite shell          28/06/2013
c      58. Add nstv - no. struct vects: d(260) = 0 to 3     26/07/2013
c          vects: d(261:263) to d(267:269)
c      58. Add fiber model based on I_4 as individual one   03/08/2013
c      59. Move 'g12' before computation of d(21) for       03/09/2013
c          Output Regular Mooney-Rivlin material.
c      60. Add 'global equations' to d(283)                 14/12/2013
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Input material parameters for material models

c      Inputs:
c         ietype    - Element type indicator
c                     1: Mechanical solid
c                     2: Truss
c                     3: Frame
c                     4: Plate
c                     5: Shell/Membrane
c                     6: Thermal solid
c                     7: Three-dimensional
c                     8: Unspecified

c      Outputs:
c         d(*)      - Material set parameters
c         tdof      - Dof to specify temperature location in
c                     coupled thermo-mechanical problems
c         ndv       - Number of element history variables
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'cdata.h'
      include  'cdat1.h'
      include  'chdata.h'
      include  'eldata.h'
      include  'eldatp.h'
      include  'elplot.h'
      include  'hdata.h'
      include  'iofile.h'
      include  'pconstant.h'
      include  'pglob1.h'
      include  'pmod2d.h'
      include  'qudshp.h'
      include  'refnd.h'
      include  'refng.h'
      include  'rigid1.h'
      include  'sdata.h'
      include  'setups.h'
      include  'strnum.h'
      include  'umac1.h'

      include  'pointer.h'
      include  'comblk.h'

      character text(2)*15,wd(9)*13,qd(4)*21,ty(2)*12

      logical   pcomp,erflag,errck,pinput,tinput,palloc,setvar
      logical   cflag,eflag,fflag,iflag,sflag,tflag,uflag,efl,flg,oned
      logical   uparfl,qdflg,qdinp,softfl, muld, inputs,incomp,damgfl
      logical   ietype15, fiberfl

      integer   eltype,ietype,imat,ndv,nvis,tdof,i,ii,j,jj,nn,umat,uprm
      integer   n1,n3,nlay,ntub,nlob,nseg,nrct,nqud,nhd,ntm,naug,nsiz,ti
      integer   tmat,nstv

      real*8    bulk,e1,e2,e3,nu12,nu23,nu31,g12,g23,g31,rad,sigc,sigt
      real*8    rrmin, t1, enorm
      real*8    d(*)
      real*8    ev(10),cc(3,3),dd(6,6),alp(3), rr(6)

      save

      data      wd/'Plane Stress' ,'Plane Strain' ,'Axisymmetric',
     &             'Plate & Shell','Truss & Frame','Thermal',
     &             '3-Dimensional','Axisy+Torsion','Unspecified' /

      data      qd/'G a u s s' ,'N o d a l','T r i a n g l e',
     &             'T e t r a h e d r o n'/

      data      ty/'Displacement' ,'Mixed       '/

c     PART 1: Set default values

      tdof     = gtdof
      etype    = 1
      dtype    = gdtype
      ietype   = abs(eltype)
      ietype15 = ietype.eq.1 .or. ietype.eq.5

      if(dtype.gt.0 .or. ietype.ne.1) then
        fflag  = .false.
        sflag  = .true.
      else
        fflag  = .true.
        sflag  = .false.
      endif

      cflag  = .false.
      eflag  = .false.
      erflag = .false.
      hflag  = .false.
      iflag  = .false.
      incomp = .false.
      oned   = .false.
      muld   = .false.
      plasfl = .false.
      qdflg  = .true.
      qdinp  = .false.
      softfl = .false.
      tflag  = .false.
      uflag  = .false.
      damgfl = .false.
      viscfl = .false.
      uparfl = .false.
      fiberfl= .false.

      if(nen.gt.11.and. ndm.eq.2) then
        ii = 4
      elseif(nen.gt.4 .and. ndm.eq.2 .or.
     &   nen.gt.8 .and. ndm.eq.3 ) then
        ii = 3
      else
        ii = 2
      endif
      jj   = ii - 1
      imat = 1
      tmat = 1
      lref = 0
      aref = 0
      naug = 0
      nvis = 0
      nlay = 0
      nqud = 0
      nrct = 0
      nseg = 0
      nstv = 0
      ntub = 0
      nlob = 2

c     Default rigid body number (denotes flexible)

      nrmatl= 0

c     Default Mass: Consistent

      d( 7) = 1.d0

c     Default angle in degrees

      rad   = pi/180.0d0
      d(31) = d(31)/rad

c     Default no thermal split

      d(136) = -1.0d0  !  No split

c     Default ground factors

      if(ietype.ne.6)then
        if(groufl) then
          do i = 1,ndf
            d(70+i) = gfac(i)
          end do ! i
        endif
        if(groupl) then
          do i = 1,ndf
            d(73+i) = gprop(i)
          end do ! i
        endif
      end if

c     Set element type

      d(170) = 1.d0  ! Default volume model

c     Solid type

      if(ietype.eq.1 .or. ietype.eq.6) then
        stype  = g2type

        if(tsplfl) then
          d(189) = 4.0d0
          do i = 1,ndm
            d(189+i) = min(4,nint(sqrt(dble(nen))))
          end do ! i
        endif ! tsplfl

c     Truss/frame type

      elseif(ietype.eq.2 .or.ietype.eq.3) then
        if(ietype.eq.2) oned   = .true.
        stype = 5
        lref  = gref
        sref  = 0
        if(lref.eq.1) then
          do i = 1,ndm
            refx(i) = grefx(i)
            tref(i) = 0.0d0
          end do ! i
        else
          do i = 1,ndm
            refx(i) = gtref(i)
            tref(i) = 0.0d0
          end do ! i
        end if

c     Plate type

      elseif(ietype.eq.4) then
        stype = 4
        ii    = 3

c     Shell type: Set reference vector for stress orientations

      elseif(ietype.eq.5) then
        stype = 4
        ii    = 2
        lref  = gref
        sref  = 0
        if(lref.eq.1) then
          do i = 1,ndm
            refx(i) = grefx(i)
            tref(i) = 0.0d0
          end do ! i
        else
          do i = 1,ndm
            refx(i) = gtref(i)
            tref(i) = 0.0d0
          end do ! i
        end if

c     Three dimensional model

      elseif(ietype.eq.7) then
        stype = 7
      endif

c     Angular velocity and Rayleigh damping set

      if(ietype.ne.6) then
        d(65) = gomega(1)
        d(77) = gray(1)
        d(78) = gray(2)
      endif

c     Quadrature type

      if(gquadn.gt.0.0d0) then
        d(182) = 1.d0                         ! Nodal Quadrature
      endif

      d(14)  = 1.d0
      alp(1) = 0.d0

c     Augment element flag

      d(185) = gaugm

c     PART 2: Poll for legal input records: Stop on blank record

      inputs = .true.

      do while(inputs)

c       Input record

        errck = .true.
        do while(errck)
          errck = tinput(text,2,ev,10)
        end do ! while

c       Reset element type

        if    (pcomp(text(1),'soli',4)) then
          ietype = 1

        elseif(pcomp(text(1),'trus',4)) then
          ietype = 2

        elseif(pcomp(text(1),'fram',4)) then
          ietype = 3

        elseif(pcomp(text(1),'plat',4)) then
          ietype = 4

        elseif(pcomp(text(1),'shel',4)) then
          ietype = 5

        elseif(pcomp(text(1),'memb',4)) then
          ietype = 5

c       Plate/Shell Thickness

        elseif(ietype.ne.2.and.ietype.ne.3
     &                    .and. pcomp(text(1),'thic',4)) then

          if(pcomp(text(2),'func',4)) then
            d( 14) = ev(1)
            d(242) = ev(2)
          else
            d( 14) = ev(1)
            if(ev(2).eq.0.0d0) then
              d(37) = 5.d0/6.d0
            else
              d(37) = ev(2)
            endif
          endif
          d(102) = ev(3)

c       Frame/Plate/Shell Kappa

        elseif(ietype.ge.3 .and. ietype.le.5
     &                     .and. pcomp(text(1),'kapp',4)) then

          if(ev(1).eq.0.0) then
            d(37) = 5.d0/6.d0
          else
            d(37) = ev(1)
          endif

c       Mass density

        elseif(pcomp(text(1),'dens',4)) then
          d(4)  = ev(1)
          if(ev(2).eq.0.0d0) then
            d(8) = 1.d0
          else
            d(8) = max(0.0d0,ev(2))
          endif

c       Damping factor

        elseif(pcomp(text(1),'damp',4)) then
          if(pcomp(text(2),'rayl',4)) then
            d(77) = ev(1)
            d(78) = ev(2)
          else
            d(70) = ev(1)
          endif

c       Mass matrix type

        elseif(pcomp(text(1),'mass',4)) then
          if(pcomp(text(2),'lump',4) .or.
     &       pcomp(text(2),'diag',4)) then
            d(  7) =  0.0d0
            d(183) =  ev(2)
          elseif(pcomp(text(2),'cons',4)) then
            d(  7) =  1.0d0
            d(183) =  ev(2)
          elseif(pcomp(text(2),'off',3)) then
            d(  7) = -1.0d0
            d(183) =  0.0d0
          else
            d(  7) =  ev(1)
            d(183) =  ev(2)
          endif

c       Method number

        elseif(pcomp(text(1),'meth',4)) then

          d(80) = ev(1)
          d(81) = ev(2)

c       Rigid body specifier

        elseif(pcomp(text(1),'rigi',4)) then
          nrmatl = nint(ev(1))
          nrbody = max(nrbody,nrmatl)

c       Hierarchical flag

        elseif(pcomp(text(1),'hier',4)) then

          if(pcomp(text(2),'off',3)) then
            d(120) = 0.0d0   ! Set formulation to normal
          else
            d(120) = 1.0d0   ! Set formulation to hierarchical
          endif

c       Normal surface loading

        elseif(((ietype.eq.4.or.ietype.eq.5) .or.
     &          (ietype.eq.3.and.ndm.eq.2)) .and.
     &           pcomp(text(1),'load',4)  ) then

          d(10) = ev(1)
          if(pcomp(text(2),'axia',4) .or. pcomp(text(2),'tang',4)) then
            d(83)  = ev(2)
            d(196) = ev(3)
          elseif(pcomp(text(2),'foll',4)) then
            d(68)  = 1.0d0
            d(196) = ev(2)
          else
            d(68)  = 0.0d0
            d(196) = ev(2)
          endif

c       Solid material body forces

c       elseif(ietype.ne.4.and.ietype.ne.6
c    &                    .and. pcomp(text(1),'body',4)) then
        elseif(pcomp(text(1),'body',4)) then

          if(pcomp(text(2),'heat',4)) then
            d(66) = ev(1)
          elseif(pcomp(text(2),'patc',4)) then
            d(197) = ev(1)
            d(198) = ev(2)
          else
            d(11) = ev(1)
            d(12) = ev(2)
            d(13) = ev(3)
            if(pcomp(text(2),'loca',4)) then
              d(69) = 2.0d0
            elseif(pcomp(text(2),'gfol',4)) then
              d(69) = 3.0d0
            elseif(pcomp(text(2),'lfol',4)) then
              d(69) = 4.0d0
            elseif(pcomp(text(2),'radi',4)) then
              d(69) = 5.0d0
            else
              d(69) = 1.0d0
            endif
          endif

c       Orthotropic material principal direction angle

        elseif((ietype.eq.1.or.ietype.eq.4.or.ietype.eq.5)
     &                    .and. pcomp(text(1),'angl',4)) then

          if(pcomp(text(2),'pola',4)) then
            d(85) = 1.0d0
            d(86) = ev(1)
            d(87) = ev(2)
            d(88) = ev(3)
          else
            d(31) = ev(1)
          endif

c       Transient indicator for element use

        elseif(pcomp(text(1),'tran',4)) then

          if(pcomp(text(2),'back',4)) then
            d(89) = 1
          elseif(pcomp(text(2),'newm',4)) then
            d(89) = 2
          elseif(pcomp(text(2),'user',4)) then
            d(89) = 3
          elseif(pcomp(text(2),'expl',4)) then
            d(187) = 1.0d0
          elseif(pcomp(text(2),'impl',4)) then
            d(187) = 0.0d0
          endif

c       Frame/Truss: Cross section properties

        elseif((ietype.eq.2.or.ietype.eq.3)
     &                    .and. pcomp(text(1),'cros',4)) then

          cflag = .true.

          do i = 1,7
            d(31+i) = ev(i)
          end do ! i

          if(ev(5).eq.0.0d0) then
            d(36) = d(33) + d(34)
          endif
          if(ev(6).eq.0.0d0) then
            d(37) = 5.d0/6.d0
          endif
          if(ev(7).eq.0.0d0) then
            d(38) = 5.d0/6.d0
          endif

c       Frame/Truss: Layer section properties

        elseif((ietype.eq.2 .or. ietype.eq.3)
     &         .and. pcomp(text(1),'laye',4)) then

          if(ndm.eq.2 .or. ietype.eq.2) then
            if(nint(d(188)).eq.0) then
              oned            = .true.
            endif
            nlay            =  nlay + 1
            if(nlay.gt.13) then
              write(ilg,4010)
              write(iow,4010)
              call plstop()
            endif
            d(101+2*nlay) =  ev(1)
            d(102+2*nlay) =  ev(2)
            nlob          =  min(6,max(nlob,nint(ev(3))))

c           Set area/inertia to prevent error message

            d(32)           = -999.d0
            d(33)           = -999.d0
          else
            write(ilg,4016)
            write(iow,4016)
            if(ior.lt.0) write(*,4016)
            erflag = .true.
          endif

c       Truss/Frame: Shaped section properties

        elseif((ietype.eq.2 .or. ietype.eq.3) .and.
     &               pcomp(text(1),'sect',4)) then

          if(nint(d(188)).eq.0) then
            oned            = .true.
          endif

c         Layer section (2-d frame elements)

          if    (pcomp(text(2),'laye',4)) then
            if(ndm.eq.2 .or. ietype.eq.2) then
              nlay            =  nlay + 1
              if(nlay.gt.13) then
                write(ilg,4010)
                write(iow,4010)
                call plstop()
              endif
              d(101+2*nlay) =  ev(1)
              d(102+2*nlay) =  ev(2)
              nlob          =  min(6,max(nlob,nint(ev(3))))
            else
              write(ilg,4016)
              write(iow,4016)
              if(ior.lt.0) write(*,4016)
              erflag = .true.
            endif

c         Tube section (thin wall)

          elseif(pcomp(text(2),'tube',4)) then

            ntub   =  max(4,nint(ev(3)))
            nlob   =  min(6,max(nlob,nint(ev(4))))
            d(100) =  1
            d(101) =  ntub         ! number of sector
            d(102) =  nlob         ! quadrature/sector
            d(103) =  ev(1)        ! radius
            d(104) =  ev(2)        ! thickness

c         Rectangular section

          elseif(pcomp(text(2),'rect',4)) then

            nrct            =  nrct + 1
            if(nrct.gt.5) then
              write(ilg,4011)
              write(iow,4011)
              call plstop()
            endif
            n1              =  max(3,min(10,nint(ev(5))))
            n3              =  max(3,min(10,nint(ev(6))))
            nqud            =  nqud + n3*n1
            d(100)          =  2
            d(101+5*nrct-4) =  ev(1)
            d(101+5*nrct-3) =  ev(2)
            d(101+5*nrct-2) =  ev(3)
            d(101+5*nrct-1) =  ev(4)

            d(101+5*nrct  ) =  100*n1 + n3

c           Check area of rectangle

            if(ev(1).ge.ev(3) .or. ev(2).ge.ev(4)) then
              write(ilg,4017) (ev(i),i=1,4)
              write(iow,4017) (ev(i),i=1,4)
              call plstop()
            endif

c         Wide-flange section

          elseif(pcomp(text(2),'wide',4)) then
            d(100)  = 3
            nqud    = 12
            do i = 1,6
              d(100+i) = ev(i)
            end do ! i

c         Channel section

          elseif(pcomp(text(2),'chan',4)) then
            d(100)  = 4
            nqud    = 12
            do i = 1,6
              d(100+i) = ev(i)
            end do ! i

c         Angle section

          elseif(pcomp(text(2),'angl',4)) then
            d(100)  = 5
            nqud    = 8
            do i = 1,6
              d(100+i) = ev(i)
            end do ! i

c         Circular solid section

          elseif(pcomp(text(2),'circ',4)) then
            d(100)  = 6
            i       = nint(ev(2))
            call int2dc(i,nqud,tt)
            do i = 1,6
              d(100+i) = ev(i)
            end do ! i

          else
            write(iow,*) ' *ERROR* - no section ',text(2)
            write(ilg,*) ' *ERROR* - no section ',text(2)
            call plstop()
          endif

c         Set area/inertia to prevent error message

          d(32)           = -999.d0
          d(33)           = -999.d0

c       Plot projection controls

        elseif(pcomp(text(1),'plot',4)) then

          if(pcomp(text(2),'off',3)) then
            d(171) = 1.0d0
          else
            d(171) = 0.0d0
          endif

c       Piezo-electric model inputs

        elseif(pcomp(text(1),'piez',4)) then

          if(pcomp(text(2),'perm',4)) then
            d(151) = ev(1)
            d(152) = ev(2)
            d(153) = ev(3)
          elseif(pcomp(text(2),'line',4) .or.
     &           pcomp(text(2),'lin+',4)) then
            d(150) = max(d(150),1.d0)
            d(154) = ev(1)
            d(155) = ev(2)
            d(156) = ev(3)
          elseif(pcomp(text(2),'lin-',4)) then
            d(150) = 4.d0
            d(157) = ev(1)
            d(158) = ev(2)
            d(159) = ev(3)
          endif

c       Frame: Reference node/vector set

        elseif(pcomp(text(1),'refe',4)) then

          if    (pcomp(text(2),'node',4)) then
            lref = 1
            do i = 1,ndm
              refx(i) = ev(i)
            end do ! i
          elseif(pcomp(text(2),'vect',4)) then
            lref = 2
            do i = 1,ndm
              refx(i) = ev(i)
            end do ! i
          elseif(pcomp(text(2),'pola',4)) then
            lref = 3
            do i = 1,ndm
              refx(i) = ev(i)
            end do ! i
          elseif(pcomp(text(2),'axia',4)) then
            lref = 4
            do i = 1,ndm
              refx(i) = ev(i)
            end do ! i
          elseif(pcomp(text(2),'shea',4)) then
            sref = 1
            do i = 1,ndm
              tref(i) = ev(i)
            end do ! i
          else
            lref = 0
            write(ilg,4008)
            write(iow,4008)
            if(ior.lt.0) then
              write(*,4008)
            end if
            erflag = .true.
          endif

c       Frame: No Shear Option

        elseif((ietype.eq.3 .or. ietype.eq.4)
     &          .and. pcomp(text(1),'shea',4)) then

          if(pcomp(text(2),'off',3)) then
            d(79) = 1.0d0
          else
            d(79) = 0.0d0
          endif

c       Truss/Frame: Nonlinear flag

        elseif((ietype.eq.2.or.ietype.eq.3)
     &                    .and. pcomp(text(1),'nonl',4)) then

          if(ev(1).eq.0.0d0) then
            d(39) = 1.0d0
          else
            d(39) = 0.0d0
          endif

c       Kinematics type: Small deformation

        elseif(pcomp(text(1),'smal',4)) then

          dtype = 1
          if(ietype.lt.6) then
            sflag = .true.
            fflag = .false.
          endif

c       Kinematics type: Finite deformation

        elseif(pcomp(text(1),'fini',4)) then

          dtype = -1
          if(ietype.lt.5) then   ! Leave shell model elastic
            fflag = .true.
            sflag = .false.
            if(pcomp(text(2),'volu',4)) then
              d(170) = max(1,min(4,nint(ev(1))))
            endif
          endif

c       Element type: Displacement

        elseif(pcomp(text(1),'disp',4)) then

          etype = 1
          if(pcomp(text(2),'unif',4)   .or.
     &       pcomp(text(2),'stab',4) ) then
            etype = 7
            if(d(60).eq.0.0d0) then
              d(60) = ev(1)
            endif
          endif

c       Element type: Mixed (B-Bar)

        elseif(pcomp(text(1),'mixe',4)) then

          etype = 2

          if(pcomp(text(2),'enha',4)) then
            etype = 5
          elseif(pcomp(text(2),'u-p',3)) then
            etype = 9
          elseif(pcomp(text(2),'unif',4)   .or.
     &           pcomp(text(2),'stab',4) ) then
            etype = 8
            if(d(60).eq.0.0d0) then
              d(60) = ev(1)
            endif
          endif

c       Element type: Enhanced Strain

        elseif(pcomp(text(1),'enha',4)) then

          etype = 3

          if(pcomp(text(2),'mixe',4)) then
            etype = 5
          endif

c       Element type: Energy-Momentum Conserving

        elseif(pcomp(text(1),'cons',4)) then

          etype = 4

c       Element type: Co-rotational Formulation

        elseif(pcomp(text(1),'coro',4)) then

          etype = 6

c       Interpolation type: Nurbs

        elseif(pcomp(text(1),'nurb',4)) then

          if(    pcomp(text(2),'tspl',4) .or.
     &           pcomp(text(2),'t-sp',4)) then
            d(189) = 4.0d0
          elseif(pcomp(text(2),'hspl',4) .or.
     &           pcomp(text(2),'h-sp',4)) then
            d(189) = 7.0d0
          elseif(pcomp(text(2),'mixe',4)) then
            d(189) = 6.0d0
          elseif(pcomp(text(2),'bezi',4)) then
            d(189) = 5.0d0
          elseif(pcomp(text(2),'bspl',4)) then
            d(189) = 3.0d0
          elseif(pcomp(text(2),'loca',4)) then
            d(189) = 2.0d0
          else
            d(189) = 1.0d0
          endif
          do i = 1,3
            if(ev(i).gt.0.0d0) then
              d(189+i) = ev(i)
            elseif(ev(i).lt.0.0d0) then
              d(189+i) = 0.0d0
            else
              d(189+i) = 5.d0
            endif
          end do ! i

c       Interpolation type: B-Spline

        elseif(pcomp(text(1),'bspl',4)) then

          d(189) = 3.0d0
          do i = 1,3
            if(ev(i).gt.0.0d0) then
              d(189+i) = ev(i)
            else
              d(189+i) = 5.d0
            endif
          end do ! i

c       Interpolation type: T-Spline

        elseif(pcomp(text(1),'tspl',4)  .or.
     &         pcomp(text(1),'t-sp',4)) then

          d(189) = 4.0d0
          do i = 1,3
            if(ev(i).gt.0.0d0) then
              d(189+i) = ev(i)
            elseif(ev(i).lt.0.0d0) then
              d(189+i) = 0.0d0
            else
              d(189+i) = 5.d0
            endif
          end do ! i

c       Interpolation type: H-Spline

        elseif(pcomp(text(1),'hspl',4)  .or.
     &         pcomp(text(1),'h-sp',4)) then

          d(189) = 7.0d0
          do i = 1,3
            if(ev(i).gt.0.0d0) then
              d(189+i) = ev(i)
            elseif(ev(i).lt.0.0d0) then
              d(189+i) = 0.0d0
            else
              d(189+i) = 5.d0
            endif
          end do ! i

c       Initial data

        elseif(pcomp(text(1),'init',4)) then

c         Constant initial stress/strain state

          if( pcomp(text(2),'stre',4)) then
            d(160) = 1
            do i = 1,6
              d(160+i) = ev(i)
            end do ! i
          elseif( pcomp(text(2),'stra',4)) then
            d(160) = 3
            do i = 1,6
              d(160+i) = ev(i)
            end do ! i
          elseif( pcomp(text(2),'augm',4)) then
            d(160) = 2
          endif

c       Quadrature data

        elseif(pcomp(text(1),'quad',4)) then

c         Set nodal quadrature flag

          qdinp = .true.
          qdflg = .true.
          if(pcomp(text(2),'noda',4) .or. pcomp(text(2),'node',4)) then
            d(182) = 1.0d0
          elseif(pcomp(text(2),'loba',4)) then
            d(182) = 0.0d0
            jj     = -1
          elseif(pcomp(text(2),'tria',4)) then
            qdflg = .false.
            ti    = 3
          elseif(pcomp(text(2),'tetr',4)) then
            qdflg = .false.
            ti    = 4
          elseif(pcomp(text(2),'gaus',4)) then
            d(182) = 0.0d0
          endif

c         Limit quadrature for built-in elements

          if(ev(1).ne.0.0d0 .or. ev(2).ne.0.0d0) then
            if(eltype.gt.0 .and. qdflg) then
              ii = max(-1,min(5,nint(ev(1))))
              jj = max( 1,min(5,nint(ev(2))))
            else
              ii = nint(ev(1))
              jj = nint(ev(2))
            endif
          endif

c       Temperature dof

        elseif(ietype.ne.6 .and.
     &     (pcomp(text(1),'temp',4) .or. pcomp(text(1),'volt',4))) then

          tdof = nint(ev(1))

c       Input solution type

        elseif((ietype15 .or. ietype.eq.6)  .and.
     &                     pcomp(text(1),'plan',4)) then
          if( pcomp(text(2),'stre',4)) then
            stype = 1
          else
            stype = 2
          endif

        elseif((ietype.eq.1 .or. ietype.eq.6)  .and.
     &                     pcomp(text(1),'axis',4)) then
          if( pcomp(text(2),'tors',4)) then
            stype = 8
          else
            stype = 3
          endif
          if(ndm.eq.1 .and. pcomp(text(2),'plan',4)) then
            d(199) = 1.0d0
          endif

c       Penalty parameter

        elseif(pcomp(text(1),'pena',4)) then

          d(60)  = ev(1)
          d(121) = ev(2)

c       Augmentation switch

        elseif(pcomp(text(1),'augm',4)) then

          if(pcomp(text(2),'off',3)) then
            d(185) = 0.0d0
          else
            d(185) = 1.0d0
          endif

          if(pcomp(text(2),'expl',4)) then
            d(186) = 1.0d0
          else
            d(186) = 0.0d0
          endif

c       Multidimensional frame constitution indicator

        elseif(pcomp(text(1),'mult',4)) then

          if(ev(1).eq.0.0d0) then
            d(188) = 3.0d0
          else
            d(188) = ev(1)
          endif
          oned = .false.
          muld = .true.

c       One-dimensional frame constitution indicator

        elseif(pcomp(text(1),'one',3)) then

          d(188) = 0.0d0
          oned   = .true.
          muld   = .false.

c       Incompressible flag

        elseif(pcomp(text(1),'inco',4)) then

          if(pcomp(text(2),'off',3)) then
            incomp = .false.
          else
            incomp = .true.
            d(170) = 5.0d0
          endif

c       Thermal properties

        elseif(ietype.ne.6 .and. pcomp(text(1),'ther',4)) then

c         Reference level from zero temperature

c         Orthotropic input

          if(pcomp(text(2),'orth',4)) then

            tflag  = .true.
            alp(1) = ev(1)
            alp(2) = ev(2)
            alp(3) = ev(3)
            d(9)   = ev(4)

c         Transversely isotropic input

          elseif(pcomp(text(2),'tran',4)) then

            tflag  = .true.
            alp(1) = ev(1)
            alp(2) = ev(2)
            alp(3) = ev(1)
            d(9)   = ev(3)

c         Thermal coupling dissipation parameters

          elseif(pcomp(text(2),'diss',4)) then

            softfl = .true.
            d(130) = ev(1)
            d(131) = ev(2)

c         Thermal coupling softening parameters

          elseif(pcomp(text(2),'soft',4)) then

            softfl = .true.
            d(132) = ev(1)
            d(133) = ev(2)
            d(134) = ev(3)
            d(135) = 0.0d0

c         Thermal coupling softening parameters

          elseif(pcomp(text(2),'spli',4)) then

            d(136) = ev(1)

c         Isotropic inputs

          else

            tflag  = .true.
            alp(1) = ev(1)
            alp(2) = ev(1)
            alp(3) = ev(1)
            d(9)   = ev(2)

          endif

c       Input error estimator value

        elseif(pcomp(text(1),'adap',4)) then

          if(pcomp(text(2),'erro',4)) then
            if(ev(1).eq.0.0d0) then
              d(50) = 0.05
            else
              d(50) = ev(1)
            endif
          endif

c       Elastic properties

        elseif(ietype.ne.6 .and. pcomp(text(1),'elas',4)) then

          eflag  = .true.

c         Transverse isotropy inputs

          if(pcomp(text(2),'tran',4)) then

            iflag = .false.
            imat  =  5

            e1   = ev(1)
            e2   = ev(2)
            e3   = ev(1)
            nu12 = ev(3)
            nu23 = ev(3)
            nu31 = ev(4)
            g12  = ev(5)
            g23  = ev(5)
            g31  = ev(1)/(2.d0 + nu31 + nu31)

c         Orthotropic inputs

          elseif(pcomp(text(2),'orth',4)) then

            iflag = .false.
            imat  =  5

            e1   = ev(1)
            e2   = ev(2)
            e3   = ev(3)
            nu12 = ev(4)
            nu23 = ev(5)
            nu31 = ev(6)
            g12  = ev(7)
            g23  = ev(8)
            g31  = ev(9)

c         Finite Elastic Models

c         Regular compressible Neo-Hookean

          elseif(ietype15 .and. pcomp(text(2),'neoh',4)) then

            imat  =  1
            dtype = -1
            fflag = .true.
            sflag = .false.

            e1    = ev(1)
            nu12  = ev(2)

c         Modified compressible Neo-Hookean

          elseif(ietype15 .and. pcomp(text(2),'mneo',4)) then

            imat  =  2
            dtype = -1
            fflag = .true.
            sflag = .false.

            e1    = ev(1)
            nu12  = ev(2)

c         Ogden compressible model

          elseif(ietype15 .and. pcomp(text(2),'ogde',4)) then

            imat  =  3
            dtype = -1
            fflag = .true.
            sflag = .false.

            bulk  = ev(1)
            g12   = ev(2)
            nu12  = ev(3)
            g23   = ev(4)
            nu23  = ev(5)
            g31   = ev(6)
            nu31  = ev(7)
            e1    = max(bulk,abs(g12),abs(g23),abs(g31))

c         Regular compressible Mooney-Rivlin

          elseif(ietype15 .and. pcomp(text(2),'moon',4) .or.
     &                          pcomp(text(2),'mmoo',4)) then

            if(pcomp(text(2),'moon',4)) then
              imat  =  9
            else
              imat  =  10
            endif
            dtype = -1
            fflag = .true.
            sflag = .false.

            e1    = ev(1)
            nu12  = ev(2)
            nu23  = ev(3)

c         Fung type exponential model

          elseif(pcomp(text(2),'fung',4)) then

            imat   =  7
            dtype  = -1
            fflag  = .true.
            sflag  = .false.

            d(30) = ev(1)
            do i = 1,9
              d(20+i) = ev(i+1)
            end do ! i
            e1     = d(30)

c         Full inputs

          elseif(pcomp(text(2),'modu',4) .or.
     &           pcomp(text(2),'comp',4)) then

            imat  =  8

c           Input by rows (store transpose)

            nsiz   = nint(ev(1))
            if(nsiz.eq.0 .or. nsiz.gt.6) then
              write(ilg,4006) text(2),nsiz
              write(iow,4006) text(2),nsiz
              call plstop()
            endif
            d(200) = nsiz
            do i =  1,nsiz
              errck = pinput(dd(1,i),nsiz)
            end do ! i

c           Convert to moduli

            if(pcomp(text(2),'comp',4)) then
              do i = 1,nsiz
                do j = 1,i
                  dd(i,j) = dd(j,i)
                end do ! j
              end do ! i
              call invert(dd,nsiz,6)
            endif

c           Save triangular part

            n1 = 200
            do i = 1,nsiz
              do j = 1,i
                n1 = n1 + 1
                d(n1) = dd(j,i)
              end do ! j
            end do ! i
            nsiz = n1
            e1   = d(201)

c         Compressible Arruda-Boyce model

          elseif(ietype15 .and. pcomp(text(2),'arru',4)) then

            imat  = 11
            fflag = .true.
            sflag = .false.
            dtype = -1

            e1    = ev(1)
            nu12  = ev(2)
            nu23  = ev(3)

c         Compressible Yeoh model

          elseif(ietype15 .and. pcomp(text(2),'yeoh',4)) then

            imat  = 12
            fflag = .true.
            sflag = .false.
            dtype = -1

            e1    = ev(1)
            nu12  = ev(2)
            nu23  = ev(3)
            nu31  = ev(4)

c         Fiber models based on I_4

          elseif(ietype15 .and. pcomp(text(2),'holz',4)
     &                     .or. pcomp(text(2),'weis',4)) then

            imat  = 14
            fflag = .true.
            sflag = .false.
            dtype = -1
            fiberfl = .true.

c           Set fiber moduli and direction values

            enorm         = sqrt(ev(3)**2 + ev(4)**2 + ev(5)**2)
            if(enorm.gt.0.0d0) then
              ev(3:5) = ev(3:5)/enorm
            else
              write(iow,4019)
              write(ilg,4019)
              call plstop()
            endif
            if(pcomp(text(2),'holz',4)) then       ! Holzapfel-Gasser model
              d(271+3*nstv) = 1.0d0
            elseif(pcomp(text(2),'weis',4)) then   ! Weiss model
              d(271+3*nstv) = 2.0d0
            endif
            d(272+3*nstv) = ev(1)
            d(273+3*nstv) = ev(2)
            d(261+3*nstv) = ev(3)
            d(262+3*nstv) = ev(4)
            d(263+3*nstv) = ev(5)
            nstv          = nstv + 1
            d(260)        = nstv

c         Isotropic inputs

          else

            imat  = 4
            if(pcomp(text(2),'stve',4).or.pcomp(text(2),'stvk',4)) then
              imat  = 5
              dtype = -1
              fflag = .true.
              sflag = .false.
            elseif(pcomp(text(2),'cons',4)) then
              imat  = 6
              dtype = -1
              etype =  4
              fflag = .true.
              sflag = .false.
            endif
            iflag = .true.

            e1    = ev(1)
            e2    = e1
            e3    = e1
            nu12  = ev(2)

            nu23  = nu12
            nu31  = nu12
            g12   = 0.5d0*e1/(1.d0 + nu12)
            g23   = g12
            g31   = g12
            bulk  = e1/(1.d0 - 2.d0*min(0.49999999d0,nu12))/3.d0

          endif

c       Tension/Compression only properties

        elseif(ietype.eq.2 .and. pcomp(text(1),'tens',4)) then
          d(167) = 1.d0
        elseif(ietype.eq.2 .and. pcomp(text(1),'comp',4)) then
          d(167) = -1.d0

c       Representative Volume Material behavior

        elseif(ietype.ne.6 .and. (pcomp(text(1),'srve',4) .or.
     &                            pcomp(text(1),'rvem',4) .or.
     &                            pcomp(text(1),'rve',3))) then

          if(ntasks.le.1) then
            write(ilg,*) ' *ERROR* Use RVE material only for np > 1'
            write(iow,*) ' *ERROR* Use RVE material only for np > 1'
            call plstop()

          elseif(rank.eq.0) then

            prtyp(1,ma) = 2
            if(pcomp(text(1),'srve',4)) then  ! Small deformation
              d(40) = 4.0d0 ! Stres  RVE
              prtyp(2,ma) =  1
            else                              ! Finite deformation
              dtype = -1
              prtyp(2,ma) = -1
            endif
            imat = 13
            e1   = 1.0d0 ! Prevents error message
            d(1) = e1
            nh1  = nh1 + 6 + 1 ! Store stresses plus one
            if(pcomp(text(1),'rvem',4)) then
              nh1 = nh1 + 4
            endif
            write(iow,2091) nh1

            if(.not.pcomp(text(2),'    ',4)) then
              i = index(xxx,',')
              if(i.eq.0) then
                i = index(xxx,' ')
              endif

              do j = i+1,i+20
                if( xxx(j:j).ne.' ' ) then
                    exit
                endif
              end do ! j
              i = index(xxx(j:j+70),' ')
              matfile(ma)   = xxx(j:j+i)
              write(iow,2092) xxx(j:j+i)
            endif

          endif

c       Representative Volume Material behavior

        elseif(ietype.eq.6 .and. pcomp(text(1),'rve',3)) then

          if(ntasks.le.1) then
            write(ilg,*) ' *ERROR* Use RVE material only for np > 1'
            write(iow,*) ' *ERROR* Use RVE material only for np > 1'
            call plstop()

          elseif(rank.eq.0) then

            prtyp(1,ma) = 1
            if(pcomp(text(1),'rve',3)) then
              d(40)       = 5.0d0 ! Thermal RVE
              prtyp(2,ma) =  1
            else
              dtype       = -1
              prtyp(2,ma) = -1
            endif
            tmat = 1
            t1   = 1.0d0 ! Prevents error message
            nh1  = nh1 + 4 + 1 ! Store temperature and gradient
            write(iow,2091) nh1

            if(.not.pcomp(text(2),'    ',4)) then
              i = index(xxx,',')
              if(i.eq.0) then
                i = index(xxx,' ')
              endif

              do j = i+1,i+20
                if( xxx(j:j).ne.' ' ) then
                    exit
                endif
              end do ! j
              i = index(xxx(j:j+70),' ')
              matfile(ma)   = xxx(j:j+i)
              write(iow,2092) xxx(j:j+i)
            endif

          endif

c       Activation indicators

        elseif(pcomp(text(1),'acti',4)) then

          if(pcomp(text(2),'ther',4)) then
            if(ev(1).eq.0.0d0) then
              d(168) =  ev(1)
            else
              d(168) = -1.0d0
            endif
          else ! mechanical
            if(ev(1).eq.0.0d0) then
              d(169) =  ev(1)
            else
              d(169) = -1.0d0
            endif
          endif

c       Viscoelastic-damgage properties

        elseif(ietype.ne.6 .and. pcomp(text(1),'dama',4)) then

          damgfl = .true.
          d(58)  = ev(1)
          d(59)  = ev(2)

c       Viscoelastic properties

        elseif(ietype.ne.6 .and. pcomp(text(1),'visc',4)) then

          viscfl = .true.
          plasfl = .false.
          d(40)  = 2.d0  ! Constitution Viscoelastic
          nvis   = nvis + 1
          if(nvis.le.3) then
            d(50+2*nvis-1) = ev(1)
            d(50+2*nvis  ) = ev(2)
            d(57)          = nvis
          else
            write(ilg,4015)
            write(iow,4015)
            call plstop()
          endif

c       Plasticity properties

        elseif(ietype.ne.6 .and. pcomp(text(1),'plas',4)) then

          plasfl = .true.
          viscfl = .false.

c         Mises yield/loading function

          if(pcomp(text(2),'mise',4)) then

            d(40)  = 1.d0  ! Constitution Isotropic Plastic
            d(41) = ev(1)
            d(42) = ev(2)
            d(43) = ev(3)
            if(d(41).eq.0.0d0) d(41) = 1.d+20
            if(d(42).eq.0.0d0) d(42) = d(41)
            if(d(43).eq.0.0d0) d(43) = 1.d0
            d(46) = 1.d0  ! Mises flag

c         Drucker-Prager yield/loading function

          elseif(pcomp(text(2),'druc',4)) then

            d(40)  = 1.d0  ! Constitution Isotropic Plastic
            d(41) = ev(1)
            d(42) = ev(2)
            if(d(42).eq.0.0d0) d(42) = d(41)
            d(46) = 2.d0  ! Drucker-Prager flag

c         Prager-Lode yield/loading function

          elseif(pcomp(text(2),'lode',4)) then

            d(40)  = 1.d0  ! Constitution Isotropic Plastic
            d(41) = ev(1)
            d(42) = ev(2)
            if(d(42).eq.0.0d0) d(42) = d(41)
            d(46) = 3.d0  ! Prager-Lode flag

c         PLAStic YIELd: sigma_y,sig_inf,beta,H_i

          elseif(pcomp(text(2),'yiel',4)) then
            d(40) = 6.0d0   ! Anisotropic plasticity indicator
            d(41) = ev(1)   ! sigma_y0
            d(42) = ev(2)   ! sigma_yinf
            d(43) = ev(3)   ! beta
            d(44) = ev(4)   ! Hiso
            d(46) = 5.0d0   ! Saturation orthotropic

c         PLAStic SWIFt: Hiso,K,eps0,n
c         N.B. sigma_y = K*eps0**n

          elseif(pcomp(text(2),'swif',4)) then
            d(40) = 6.0d0   ! Anisotropic plasticity indicator
            d(41) = ev(1)*ev(2)**ev(3) ! sigma_y (for output, not used)
            d(42) = ev(1)   ! K
            d(43) = ev(2)   ! eps0
            d(44) = ev(3)   ! n
            d(45) = ev(4)   ! Hiso
            d(46) = 6.0d0   ! Swift orthotropic

c         PLAStic HILL functional form: R_11, R_22, R_33, R_12, R_23, R_31

          elseif(pcomp(text(2),'hill',4)) then
            rrmin  = 100.0d0
            do i = 1,6
              rr(i) = ev(i)
              rrmin = min(rrmin,rr(i))
            end do ! i

            if(rrmin.le.0.0d0) then
              write(iow,4018)
              call mprint(rr,1,6,1,'Hill Ratios')
              call plstop()
            endif

c           Hill parameters: Reciprocal squares (for printing)

            d(54) = 1.0d0/rr(1)**2
            d(55) = 1.0d0/rr(2)**2
            d(56) = 1.0d0/rr(3)**2
            d(51) = 0.5d0*(d(55) + d(56) - d(54))  ! F
            d(52) = 0.5d0*(d(56) + d(54) - d(55))  ! G
            d(53) = 0.5d0*(d(54) + d(55) - d(56))  ! H
            d(54) = 0.5d0/rr(5)**2                 ! L
            d(55) = 0.5d0/rr(6)**2                 ! M
            d(56) = 0.5d0/rr(4)**2                 ! N

c         Hardening parameters

          elseif(pcomp(text(2),'hard',4)) then

            d(44)  = ev(1)
            d(45)  = ev(2)

c         PLAStic KINEmatic hardening:   h2,h3,h4,h5,h6,j1

          elseif(pcomp(text(2),'kine',4)) then
              d(150) = 1.0d0
              do i = 1,6
                d(i+150) = ev(i)
              end do ! i

c         PLAStic orientation VECTors: V_1 & V_2

          elseif(pcomp(text(2),'vect',4)) then
            d(242) = 1.0d0
            do i = 1,6
              d(i+242) = ev(i)
            end do ! i

c         Linear segment hardening parameters

          elseif(pcomp(text(2),'segm',4)) then

            nseg          = nseg + 1
            if(nseg.gt.6) then
              write(ilg,4012)
              write(iow,4012)
              call plstop()
            endif
            d(130)        = nseg
            d(128+3*nseg) = ev(1)
            d(129+3*nseg) = ev(2)
            d(130+3*nseg) = ev(3)

c         Viscoplastic rate parameter

          elseif(pcomp(text(2),'rate',4) .or.
     &           pcomp(text(2),'time',4)) then

            d(180) = ev(1)
            d(181) = max(1.d0,ev(2))

c        Generalized model parameters

          elseif(pcomp(text(2),'gene',4)) then

            d(40)  = 3.d0  ! Constitution Generalized Plastic
            d(41)  = ev(1)
            d(42)  = ev(2)
            d(43)  = ev(3)
            d(46)  = 4.d0  ! Finite deformation indicator for G-plast

          endif

c       VECTor ORTHotropic

        elseif(ietype.ne.6 .and. pcomp(text(1),'vect',4)) then

c         Orthotropic orientation VECTors: V_1 & V_2

          if(pcomp(text(2),'orth',4)) then
            d(242) = 1.0d0
            do i = 1,6
              d(i+242) = ev(i)
            end do ! i
         endif

c       Orientation for anisotropic models

        elseif((ietype.eq.1 .or. ietype.ge.6)  .and.
     &                     pcomp(text(1),'refe',4)) then

          if(pcomp(text(2),'cyli',4)) then
            aref = 1
            do i = 1,3
              refa(i,1) = ev(i)
              refa(i,2) = ev(i+3)
            end do ! i
          endif

c       Anisotropic structure vectors

        elseif(pcomp(text(1),'stru',4)) then

          enorm         = sqrt(ev(1)**2 + ev(2)**2 + ev(3)**3)
          if(enorm.gt.0.0d0) then
            ev(1:3) = ev(1:3)/enorm
          else
            write(iow,4019)
            write(ilg,4019)
            call plstop()
          endif
          d(261+3*nstv) = ev(1)
          d(262+3*nstv) = ev(2)
          d(263+3*nstv) = ev(3)
          nstv          = nstv + 1
          d(260)        = nstv

c       Anisotropic fiber parameters

        elseif(pcomp(text(1),'fibe',4)) then

c         Fiber reference structure vector

          enorm         = sqrt(ev(3)**2 + ev(4)**2 + ev(5)**2)
          if(enorm.gt.0.0d0) then
            ev(3:5) = ev(3:5)/enorm
          else
            write(iow,4019)
            write(ilg,4019)
            call plstop()
          endif
          d(261+3*nstv) = ev(3)
          d(262+3*nstv) = ev(4)
          d(263+3*nstv) = ev(5)

c         Holzapfel-Gasser model

          if(pcomp(text(2),'holz',4)) then
            d(271+3*nstv) = 1.0d0
            d(272+3*nstv) = ev(1)
            d(273+3*nstv) = ev(2)
            fiberfl       = .true.

c         Weiss model

          elseif(pcomp(text(2),'weis',4)) then
            d(271+3*nstv) = 2.0d0
            d(272+3*nstv) = ev(1)
            d(273+3*nstv) = ev(2)
            fiberfl       = .true.
          endif

c         Set number of fibers

          nstv          = nstv + 1
          d(260)        = nstv

c       History plot desinations

        elseif(pcomp(text(1),'hist',4)) then

c         Allocate the history table

          if(np(303).eq.0) then
            setvar = palloc(303,'HPLTB',10*nummat, 1)
            histpltfl = .true.
          endif
          call sethisplt(mr(np(303)),ev)

c       Global equations

        elseif(pcomp(text(1),'glob',4)) then

          if(pcomp(text(2),'    ',4) .or. pcomp(text(2),'equa',4)) then
            d(283) = ev(1)
          endif

c       Angular velocity: rad/sec

        elseif(ietype.ne.6 .and. pcomp(text(1),'omeg',4)) then

          if(pcomp(text(2),'cycl',4)) then
            d(65) = ev(1)*2.d0*pi
          else
            d(65) = ev(1)
          endif

c       Fourier Heat Conduction properties

        elseif(pcomp(text(1),'four',4)) then

          hflag = .true.
          d(67) = 1.0d0 ! Heat constitution added
          tmat  = 2    ! Fourier model

c         Surface convection

          if(pcomp(text(2),'conv',4)) then

            d(127) = ev(1)    ! Surface convection      (h)
            d(128) = ev(2)    ! Free stream temperature (T_inf)

c         Reference absolute temperature

          elseif(pcomp(text(2),'temp',4)  .or.
     &           pcomp(text(2),'refe',4)) then

            d(129) = ev(1)    ! Reference absolute temperature (T_ref)
            d(195) = ev(2)    ! Fraction mechanical work going to heat

c         Specific heat

          elseif(pcomp(text(2),'spec',4) .or.
     &           pcomp(text(2),'heat',4)) then

            d(64) = ev(1)     ! Specific heat (c)

c         Orthotropic inputs

          elseif(pcomp(text(2),'orth',4)) then

            d(61) = ev(1)     ! k_1 - conductivity
            d(62) = ev(2)     ! k_2 - conductivity
            d(63) = ev(3)     ! k_3 - conductivity
            d(64) = ev(4)     ! c   - specific heat
            t1    = d(61)

c         Isotropic inputs

          else

            d(61) = ev(1)
            d(62) = ev(1)
            d(63) = ev(1)
            d(64) = ev(2)
            t1    = d(61)

          endif

c         Thermal dof for mechanical solid and truss problems

          if(ietype.eq.1 .and. tdof.eq.0 .or. ietype.eq.2) then

            tdof = ndm + 1

          endif

c       Ground motion acceleration factors/proportional load numbers

        elseif(ietype.ne.6 .and. pcomp(text(1),'grou',4)) then

          do i = 1,ndm
            d(70+i) = ev(2*i-1)
            d(73+i) = ev(2*i)
          end do ! i

c       Element (Lagrange multiplier) unknowns

        elseif(pcomp(text(1),'elem',4).or.pcomp(text(1),'lagr',4)) then

          d(230) = ev(1)

c       Formulation type:

        elseif(pcomp(text(1),'form',4)) then

          if(pcomp(text(2),'elem',4)) then
            d(240) = 0.0d0
          elseif(pcomp(text(2),'noda',4)) then
            d(240) = 1.0d0
          endif

c       Element degree of freedoms

        elseif(pcomp(text(1),'dofs',4)) then

          d(241) = ev(1)

c       Consititutive solution start state: Default = elastic (0)

        elseif(pcomp(text(1),'star',4)) then

          if(pcomp(text(2),'elas',4)) then
            d(84)  = 0.0d0
          elseif(pcomp(text(2),'inel',4).or.pcomp(text(2),'plas',4))then
            d(84)  = 1.0d0
          else
            d(84)  = ev(1)
          endif

c       User Material Parameters

        elseif(pcomp(text(1),'upar',4)) then

          uparfl = .true.
          do i = 1,10
            d(230+i) = ev(i)
          end do ! i

c       User Material Model interface

        elseif(pcomp(text(1),'ucon',4).or.pcomp(text(1),'fcon',4)) then

c         Default user constitutive equation number

          imat    = 0    ! Clear internal material number for user one
          umat    = 0
          uprm    = ndd - nud
          n1      = 0
          n3      = 0

          uct     = 'read'
          call uconst(text(2),ev, d(1), d(uprm+1),n1,n3, umat)

          if(pcomp(text(1),'ucon',4)) then
            d(uprm) = umat + 100
          else
            d(uprm) = umat + 200
            dtype   = -1
          endif

c         Activate user program models

          uflag  = .true.
          if(e1.eq.0.0d0) then
            e1 = 1.0d0
          endif

c         Increase number of history terms/quadrature point

          d(15)  = n1
          nh1    = nh1 + n1
          nh3    = nh3 + n3

c       Check end of data

        elseif(pcomp(text(1),'    ',4)) then

c         Transfer to sets and checks

          inputs = .false.

        endif

      end do ! while

c     Number of stress/strain history terms/pt

      if(ndm.eq.3 .or. stype.eq.8) then
        ntm = 6
      elseif(ndm.eq.2) then
        ntm = 4
      else
        ntm = 1
      endif

c     PART 3: Set final parameters and output material state

c     Set moduli

      if(ietype.ne.6) then


c       Small deformation options

        if(sflag) then

          if(imat.ne.8 .and. imat.ne.13) then
            d(1)    = e1
            d(2)    = nu12
c           d(3)    = alp(1)

            dd(1,1) =  1.d0/e1
            dd(2,2) =  1.d0/e2
            dd(3,3) =  1.d0/e3

            dd(1,2) = -dd(1,1)*nu12
            dd(1,3) = -dd(3,3)*nu31
            dd(2,3) = -dd(2,2)*nu23

c           1-Dimensional Models

            if(stype.eq.5 .and. oned) then
              dd(2,2) =  1.d0
              dd(3,3) =  1.d0
              dd(1,2) =  0.d0
              dd(1,3) =  0.d0
              dd(2,3) =  0.d0

c           Plane Stress Models

            elseif(stype .eq.4 .or. (stype.eq.1 .and. plasfl)) then

              d(90)   =  dd(1,3)
              d(91)   =  dd(2,3)

              dd(3,3) =  1.d0
              dd(2,3) =  0.d0
              dd(1,3) =  0.d0

            endif

            dd(2,1) = dd(1,2)
            dd(3,1) = dd(1,3)
            dd(3,2) = dd(2,3)

c           Mechanical modulus properties

            if(.not.incomp) then
              do j = 1,3
                do i = 1,3
                  cc(i,j) = dd(i,j)
                end do ! i
              end do ! j
              call invert(dd,3,6)
              if(min(dd(1,1),dd(2,2),dd(3,3)).lt.0.0d0) then
                write(iow,4005) ((cc(i,j),j=1,3),i=1,3),
     &                          ((dd(i,j),j=1,3),i=1,3)
                write(ilg,4005) ((cc(i,j),j=1,3),i=1,3),
     &                          ((dd(i,j),j=1,3),i=1,3)
                if(ior.lt.0) then
                  write(*,4005) ((cc(i,j),j=1,3),i=1,3),
     &                          ((dd(i,j),j=1,3),i=1,3)
                endif
              endif
            elseif(iflag) then
              dd(1,1) =  4.d0*g12/3.d0
              dd(2,2) =  dd(1,1)
              dd(3,3) =  dd(1,1)
              dd(1,2) = -dd(1,1)*0.5d0
              dd(2,3) =  dd(1,2)
              dd(1,3) =  dd(1,2)
            endif

c           Set for plane stress

            if(stype.eq.5 .and. ietype.ne.3) then
              dd(3,3) = 0.0d0
            endif

c           Save moduli

            d(21) = dd(1,1)
            d(22) = dd(2,2)
            d(23) = dd(3,3)
            d(24) = dd(1,2)
            d(25) = dd(2,3)
            d(26) = dd(1,3)
            d(27) = g12
            d(28) = g23
            d(29) = g31

c           Thermal properties

            d(47) = dd(1,1)*alp(1) + dd(1,2)*alp(2) + dd(1,3)*alp(3)
            d(48) = dd(2,1)*alp(1) + dd(2,2)*alp(2) + dd(2,3)*alp(3)
            d(49) = dd(3,1)*alp(1) + dd(3,2)*alp(2) + dd(3,3)*alp(3)

c           Set for plane stress problems

            if(stype .eq.4 .or. (stype.eq.1 .and. plasfl)) then
              d(23) = 0.0d0
              d(49) = 0.0d0
              d(92) = alp(3)
            endif
          endif ! imat.ne.8

c         Output parameters for element

          if(plasfl) then
            jj   = ii
            d(6) = jj
          endif

c         Output elastic properties

          if(eflag) then
            if(stype.eq.5 .or. iflag) then
              write(iow,2000) wd(stype),e1,nu12
              if(ior.lt.0) then
                write(*,2000) wd(stype),e1,nu12
              endif
            elseif(imat.eq.8) then
              write(iow,2063) (d(i),i=201,nsiz)
              if(ior.lt.0) then
                write(*,2063) (d(i),i=201,nsiz)
              endif
            else
              write(iow,2001) wd(stype),e1,e2,e3,nu12,nu23,nu31,
     &                        g12,g23,g31,d(31)
              write(iow,2096) ((dd(i,j),j=1,3),i=1,3)
              if(ior.lt.0) then
                write(*,2001) wd(stype),e1,e2,e3,nu12,nu23,nu31,
     &                        g12,g23,g31,d(31)
                write(*,2096) ((dd(i,j),j=1,3),i=1,3)
              endif
              if(nint(d(85)).gt.0) then
                write(iow,2065) (d(i),i=86,85+ndm)
                if(ior.lt.0) then
                  write(*,2065) (d(i),i=86,85+ndm)
                endif
              endif
            endif
            if(ietype15 .and. qdinp) then
              if(d(182).gt.0.0d0) then
                i = 2
              else
                i = 1
              endif
              write(iow,2019) qd(i),ii,jj
              if(ior.lt.0) then
                write(*,2019) qd(i),ii,jj
              endif
            endif
          elseif(fiberfl) then
            continue
          elseif(.not.uflag .and. imat.ne.13) then
            write(ilg,4000)
            write(iow,4000)
            if(ior.lt.0) write(*,4000)
            erflag = .true.
          endif

c       Finite deformation options

        elseif(fflag)then

c         Set default thermal expansion coefficients

          d(47) = alp(1)
          d(48) = alp(2)
          d(49) = alp(3)

c         Check for isotropic elastic/orthotropic plastic

          if(imat.eq.4 .and. nint(d(40)).eq.6) then
            imat  = 5
            iflag = .true.
          endif

c         Output Regular NeoHookean

          if(imat.eq.1) then

            if(incomp) then
              bulk = 0.0d0
            else
              bulk = e1/(1.d0 - nu12*2.d0)/3.d0
            endif
            g12  = e1/(1.d0 + nu12)/2.d0
            if(ior.lt.0) then
              write(*,2010) ' ',wd(stype),e1,nu12,bulk,g12
            endif
            write(iow,2010) ' ',wd(stype),e1,nu12,bulk,g12
            if(nint(d(170)).lt.5) then
              if(ior.lt.0) then
                write(*,2062) nint(d(170))
              endif
              write(iow,2062) nint(d(170))
            endif

c           Compute Lame' parameters

            d(1)  = e1
            d(2)  = nu12
c           d(3)  = alp(1)

            if(incomp) then
              d(21) = 0.0d0
            else
              d(21) = bulk - 2.d0/3.d0*g12
            endif
            d(22) = g12

c         Output Modified NeoHookean

          elseif(imat.eq.2) then

            if(incomp) then
              bulk = 0.0d0
            else
              bulk = e1/(1.d0 - nu12*2.d0)/3.d0
            endif
            g12  = e1/(1.d0 + nu12)/2.d0
            if(ior.lt.0) then
              write(*,2010) ' Modified ',wd(stype),e1,nu12,bulk,g12
            endif
            write(iow,2010) ' Modified ',wd(stype),e1,nu12,bulk,g12
            if(nint(d(170)).lt.5) then
              if(ior.lt.0) then
                write(*,2062) nint(d(170))
              endif
              write(iow,2062) nint(d(170))
            endif

c           Set material parameters

            d(1)  = e1
            d(2)  = nu12
c           d(3)  = alp(1)

            d(21) = bulk
            d(22) = g12

c         Output Ogden

          elseif(imat.eq.3) then

            nn = 1
            if(nu23.ne.0.0d0) nn = 2
            if(nu31.ne.0.0d0) nn = 3
            if(nn.eq.1) then
              if(ior.lt.0) then
                write(*,2024) wd(stype),bulk,g12,nu12
              endif
              write(iow,2024) wd(stype),bulk,g12,nu12
            elseif(nn.eq.2) then
              if(ior.lt.0) then
                write(*,2024) wd(stype),bulk,g12,nu12,g23,nu23
              endif
              write(iow,2024) wd(stype),bulk,g12,nu12,g23,nu23
            elseif(nn.eq.3) then
              if(ior.lt.0) then
                write(*,2024) wd(stype),bulk,g12,nu12,g23,nu23,g31,nu31
              endif
              write(iow,2024) wd(stype),bulk,g12,nu12,g23,nu23,g31,nu31
            endif
            if(nint(d(170)).lt.5) then
              if(ior.lt.0) then
                write(*,2062) nint(d(170))
              endif
              write(iow,2062) nint(d(170))
            endif

c           Set material parameters

            d(21) = bulk
            d(22) = g12
            d(23) = nu12
            d(24) = g23
            d(25) = nu23
            d(26) = g31
            d(27) = nu31
            d(28) = nn

c         Logarithmic stretch model: Like linear elastic

          elseif(imat.eq.4) then

            if(ior.lt.0) then
              write(*,2027) wd(stype),e1,nu12,bulk,g12
            endif
            write(iow,2027) wd(stype),e1,nu12,bulk,g12

c         Set material parameters

            d(1)  = e1
            d(2)  = nu12
c           d(3)  = alp(1)

            d(21) = bulk
            d(22) = g12

c         St. Venant-Kirchhoff and Orthotropic elastic

          elseif(imat.eq.5 .or. imat.eq.6) then

            write(iow,2094)
            if(ior.lt.0) write(*,2094)
            if(iflag) then
              write(iow,2000) wd(stype),e1,nu12
              if(ior.lt.0) then
                write(*,2000) wd(stype),e1,nu12
              endif
            else
              write(iow,2001) wd(stype),e1,e2,e3,nu12,nu23,nu31,
     &                        g12,g23,g31,d(31)
              if(ior.lt.0) then
                write(*,2001) wd(stype),e1,e2,e3,nu12,nu23,nu31,
     &                        g12,g23,g31,d(31)
              endif
              if(d(85).gt.0) then
                write(iow,2065) (d(i),i=1,86,85+ndm)
                if(ior.lt.0) then
                  write(*,2065) (d(i),i=1,86,85+ndm)
                endif
              endif
            endif

c           Set material parameters

            d(1)    = e1
            d(2)    = nu12
c           d(3)    = alp(1)

            dd(1,1) =  1.d0/e1
            dd(2,2) =  1.d0/e2

            dd(1,2) = -dd(1,1)*nu12
            dd(2,1) =  dd(1,2)

            if(stype.ne.1 .or. nint(d(40)).eq.6) then
              dd(3,3) =  1.d0/e3

              dd(1,3) = -dd(3,3)*nu31
              dd(2,3) = -dd(2,2)*nu23

              dd(3,1) =  dd(1,3)
              dd(3,2) =  dd(2,3)
            else
              dd(3,3) =  1.d0

              dd(1,3) =  0.d0
              dd(2,3) =  0.d0

              dd(3,1) =  0.d0
              dd(3,2) =  0.d0
            endif

c           Mechanical modulus properties

            if(.not.incomp .or. stype.eq.4) then
              call invert(dd,3,6)
              write(iow,2096) ((dd(i,j),j=1,3),i=1,3)
              if(ior.lt.0) then
                write(*,2096) (( dd(i,j),j=1,3),i=1,3)
              endif
            endif

c           Save moduli

            d(21) = dd(1,1)
            d(22) = dd(2,2)
            d(23) = dd(3,3)
            d(24) = dd(1,2)
            d(25) = dd(2,3)
            d(26) = dd(3,1)
            d(27) = g12
            d(28) = g23
            d(29) = g31

c         Fung model

          elseif(imat.eq.7) then
            if(ietype.eq.5) then
              write(iow,2054) wd(stype),(d(20+i),i=0,4)
              if(ior.lt.0) then
                write(*,2054) wd(stype),(d(20+i),i=0,4)
              endif
            else
              write(iow,2055) wd(stype),(d(20+i),i=0,9)
              if(ior.lt.0) then
                write(*,2055) wd(stype),(d(20+i),i=0,9)
              endif
            endif

c         Output Regular Mooney-Rivlin

          elseif(imat.eq.9) then

            g12  = e1/(1.d0 + nu12)/2.d0
            if(incomp) then
              bulk  = 0.0d0
              d(21) = 0.0d0
            else
              bulk  = e1/(1.d0 - nu12*2.d0)/3.d0
              d(21) = bulk - 2.d0/3.d0*g12
            endif
            if(ior.lt.0) then
              write(*,2074) ' ',wd(stype),e1,nu12,bulk,g12,nu23
            endif
            write(iow,2074) ' ',wd(stype),e1,nu12,bulk,g12,nu23
            if(nint(d(170)).lt.5) then
              if(ior.lt.0) then
                write(*,2062) nint(d(170))
              endif
              write(iow,2062) nint(d(170))
            endif

c           Compute Lame' parameters

            d(1)  = e1
            d(2)  = nu12
c           d(3)  = alp(1)

            d(22) = g12
            d(23) = nu23

c         Output Modified Mooney-Rivlin

          elseif(imat.eq.10) then

            d(21) = e1/(1.d0 - nu12*2.d0)/3.d0
            d(22) = e1/(1.d0 + nu12)/2.d0
            d(23) = nu23
            if(ior.lt.0) then
              write(*,2074) ' Modified ',wd(stype),e1,nu12,
     &                      (d(i),i=21,23)
            endif
            write(iow,2074) ' Modified ',wd(stype),e1,nu12,
     &                      (d(i),i=21,23)
            if(nint(d(170)).lt.5) then
              if(ior.lt.0) then
                write(*,2062) nint(d(170))
              endif
              write(iow,2062) nint(d(170))
            endif

            d(1)  = e1
            d(2)  = nu12
c           d(3)  = alp(1)

c         Output compressible Arruda-Boyce

          elseif(imat.eq.11) then

            if(incomp) then
              bulk = 0.0d0
            else
              bulk = e1/(1.d0 - nu12*2.d0)/3.d0
            endif
            g12  = e1/(1.d0 + nu12)/2.d0
            if(ior.lt.0) then
              write(*,2077) ' Modified ',wd(stype),e1,nu12,nu23,bulk,g12
            endif
            write(iow,2077) ' Modified ',wd(stype),e1,nu12,nu23,bulk,g12
            if(nint(d(170)).lt.5) then
              if(ior.lt.0) then
                write(*,2062) nint(d(170))
              endif
              write(iow,2062) nint(d(170))
            endif

c           Set material parameters

            d(1)  = e1
            d(2)  = nu12
c           d(3)  = alp(1)

            d(21) = bulk
            d(22) = g12/(1+ 0.6d0*nu23 + 99.d0*nu23*nu23/175.d0)
            d(23) = nu23

c         Output compressible Yeoh

          elseif(imat.eq.12) then

            if(incomp) then
              bulk = 0.0d0
            else
              bulk = e1/(1.d0 - nu12*2.d0)/3.d0
            endif
            g12  = e1/(1.d0 + nu12)/2.d0
            if(ior.lt.0) then
              write(*,2078) ' Modified ',wd(stype),e1,nu12,nu23,nu31,
     &                                   bulk,g12
            endif
            write(iow,2078) ' Modified ',wd(stype),e1,nu12,nu23,nu31,
     &                                   bulk,g12
            if(nint(d(170)).lt.5) then
              if(ior.lt.0) then
                write(*,2062) nint(d(170))
              endif
              write(iow,2062) nint(d(170))
            endif

c           Set material parameters

            d(1)  = e1
            d(2)  = nu12
c           d(3)  = alp(1)

            d(21) = bulk
            d(22) = g12
            d(23) = nu23
            d(24) = nu31

          endif
        endif

c       Axisymmetric Plane Stress option output

        if(d(199).gt.0.0d0) then
          write(iow,2097)
          if(ior.lt.0) then
            write(*,2097)
          endif
        endif

c       Output thermal expansions

        if(tflag) then
          write(iow,2002) alp,d(9),tdof
          if(ior.lt.0) then
            write(*,2002) alp,d(9),tdof
          endif
          if(d(129).ne.0.0d0) then
            write(iow,2076) d(129),d(195)
            if(ior.lt.0) then
              write(*,2076) d(129),d(195)
            endif
          endif
          if(softfl) then
            write(iow,2089) (d(i),i=130,135),nint(d(136))
            if(ior.lt.0) then
              write(*,2089) (d(i),i=130,135),nint(d(136))
            endif
          endif
        endif

c       Output Fourier heat conduction properties

        if(hflag) then
          write(iow,2020) wd(stype),(d(i),i=61,64),d(66)
          if(d(127).gt.0.0d0) then
            write(iow,2073) d(127),d(128)
          endif
          if(ior.lt.0) then
            write(*,2020) wd(stype),(d(i),i=61,64),d(66)
            if(d(127).gt.0.0d0) then
              write(*,2073) d(127),d(128)
            endif
          endif
        endif

c       Output tension/compression only states

        if(d(167).gt.0.0d0) then
          if(ior.lt.0) then
            write(*,2052)
          endif
          write(iow,2052)
        elseif(d(167).lt.0) then
          if(ior.lt.0) then
            write(*,2053)
          endif
          write(iow,2053)
        endif

c       Output activation indicators

        if(d(168).ne.0.0d0) then
          if(ior.lt.0) then
            write(*,2070)
          endif
          write(iow,2070)
        endif
        if(d(169).ne.0) then
          if(ior.lt.0) then
            write(*,2071)
          endif
          write(iow,2071)
        endif

c       Output Transient type

        if(d(89).gt.0.0d0) then
          if(ior.lt.0) then
            write(*,2072) nint(d(89))
          endif
          write(iow,2072) nint(d(89))
        endif

c       Output shell/plate thickness

        if(stype.ne.5 .and. stype.ne.7) then
          write(iow,2018) d(14)
          if(nint(d(102)).gt.1) write(iow,2051) nint(d(102))
          if(nint(d(242)).ge.1) write(iow,3000) nint(d(242))
        endif
        if(stype.eq.4) then
          write(iow,2021) d(10),nint(d(196))
          if(dtype.lt.0 .and. d(68).gt.0.0d0) then
            write(iow,2033)
          endif
        endif
        if(ior.lt.0) then
          if(stype.ne.5 .and. stype.ne.7) then
            write(*,2018) d(14)
            if(nint(d(102)).gt.1) write(*,2051) nint(d(102))
            if(nint(d(242)).ge.1) write(*,3000) nint(d(242))
          endif
          if(stype.eq.4) then
            write(*,2021) d(10),nint(d(196))
            if(dtype.lt.0 .and. d(68).gt.0.0d0) then
              write(iow,2033)
            endif
          endif
        endif

c       Output density and body loading

        if(ior.lt.0) then
          write(*,2029) d(4),d(11),d(12),d(13)
        endif
        write(iow,2029) d(4),d(11),d(12),d(13)

c       Patch loading

        if(max(abs(d(197)),abs(d(198))).gt.0.0d0) then
          write(iow,2095) d(197),d(198)
        endif

c       Maximum wave speed estimate for isotropy

        if(e1.gt.0.0d0 .and. nu12.lt.0.5d0 .and. d(4).gt.0.0d0) then
          d(184) = e1*(1.d0-nu12)/((1.d0+nu12)*(1.d0-2.d0*nu12)*d(4))
          d(184) = sqrt(d(184))
          write(iow,2083) d(184)
          if(d(187).gt.0.0d0) then
            write(iow,2085)
          else
            write(iow,2086)
          endif
          if(ior.lt.0) then
            write(*,2083) d(184)
            if(d(187).gt.0.0d0) then
              write(*,2085)
            else
              write(*,2086)
            endif
          endif
        endif

c       Output constant initial stresses

        if(nint(d(160)).eq.1) then
          if(ietype.eq.2 .or. ietype.eq.3) then
            if(oned) then
              j = 1
            else
              j = 6
            endif
          else
            j = 6
          endif
          write(iow,2044) (d(160+i),i=1,j)
          if(ior.lt.0) then
            write(*,2044) (d(160+i),i=1,j)
          endif
        elseif(nint(d(160)).eq.3) then
          if(ietype.eq.2 .or. ietype.eq.3) then
            if(oned) then
              j = 1
            else
              j = 6
            endif
          else
            j = 6
          endif
          write(iow,2045) (d(160+i),i=1,j)
          if(ior.lt.0) then
            write(*,2045) (d(160+i),i=1,j)
          endif
        elseif(nint(d(160)).eq.2) then
          if(ietype.eq.2 .or. ietype.eq.3) then
            naug = 1
            write(iow,2014)
            if(ior.lt.0) then
              write(*,2014)
            endif
          endif
        endif

c       Output frame distributed loading type

        if(d(69).eq.1.0d0) then
          write(iow,2034)
          if(ior.lt.0) write(*,2034)
        elseif(d(69).eq.2.0d0) then
          write(iow,2035)
          if(ior.lt.0) write(*,2035)
        elseif(d(69).eq.3.0d0) then
          write(iow,2036)
          if(ior.lt.0) write(*,2036)
        elseif(d(69).eq.4.0d0) then
          write(iow,2037)
          if(ior.lt.0) write(*,2037)
        elseif(d(69).eq.5.0d0) then
          write(iow,2087)
          if(ior.lt.0) write(*,2087)
        endif

c       Output angular velocity

        if(d(65).ne.0.0d0) then
          if(ior.lt.0) then
            write(*,2030) d(65)
          endif
          write(iow,2030) d(65)
        endif

c       Output ground acceleration factors

        flg = .false.
        efl = .false.
        do i = 1,ndm
          if(nint(d(73+i)).gt.0 ) flg = .true.
          if(nint(d(73+i)).gt.10) then
            efl = .true.
          endif
        end do ! i
        if(flg) then
          write(iow,2023) (d(70+i),nint(d(73+i)),i=1,ndm)
          if(ior.lt.0) then
            write(*,2023) (d(70+i),nint(d(73+i)),i=1,ndm)
          endif
        endif

        if(efl) then
          write(ilg,4007)
          write(iow,4007)
          if(ior.lt.0) write(*,4007)
          erflag = .true.
        endif

c       Output section properties

        if(cflag) then
          if(ietype.eq.2) then
            j = 32
          else
            j = 38
          endif
          write(iow,2015) (d(i),i=32,j)
          if(ior.lt.0) then
            write(*,2015) (d(i),i=32,j)
          endif

          if(ndm.eq.2 .and. ietype.eq.3) then
            if(d(33).lt.0.0d0) then
              write(ilg,4009) d(33)
              write(iow,4009) d(33)
              if(ior.lt.0) then
                write(*,4009) d(33)
              endif
              erflag = .true.
            endif
          elseif(ndm.eq.3 .and. ietype.eq.3) then
            ev(1) = d(33)*d(34) - d(35)*d(35)
            ev(2) = 1.d-8*max(abs(d(33)),abs(d(34)))
            if((ev(1).lt.ev(2) .and. d(35).ne.0.0d0)
     &                          .or. ev(1).le.0.0d0) then
              write(ilg,4009) ev(1)
              write(iow,4009) ev(1)
              if(ior.lt.0) then
                write(*,4009) ev(1)
              endif
              erflag = .true.
            endif
          endif
        endif

c       Output layer data

        if(nlay.gt.0) then
          d(101) = nlay
          d(102) = nlob
          nqud   = nlay*nlob
          if(ior.lt.0) then
            write(*,2038) (i,d(101+2*i),d(101+2*i+1),i=1,nlay)
            write(*,2039) nlob
          endif
          write(iow,2038) (i,d(101+2*i),d(101+2*i+1),i=1,nlay)
          write(iow,2039) nlob
          if(d(37).eq.0.0d0) d(37) = 5.d0/6.d0

c         Compute location of centroid and shift input data

c         call bmcent(nlay,d(102))

c         Output shifted layer data

c         if(ior.lt.0) then
c           write(*,2038) (i,d(101+2*i),d(101+2*i+1),i=1,nlay)
c           write(*,2039) nlob
c         endif
c         write(iow,2038) (i,d(101+2*i),d(101+2*i+1),i=1,nlay)
c         write(iow,2039) nlob
        endif

c       Set default shear factor for 3-d beams

        if(nint(d(100)) .gt. 0) then
          if(d(37).eq.0.0d0) d(37) = 5.d0/6.d0
          if(d(38).eq.0.0d0) d(38) = 5.d0/6.d0
        endif

c       Output tube data

        if(ntub.gt.0) then
          if(ior.lt.0) then
            write(*,2040) d(103),d(104),ntub,nlob
          endif
          write(iow,2040) d(103),d(104),ntub,nlob
        endif

c       Output rectangle data

        if(nrct.gt.0) then
          d(101) = nrct
          if(ior.lt.0) then
            write(*,2042) ((d(96+5*j+i),i=1,4),nint(d(101+5*j))/100,
     &                      mod(nint(d(101+5*j)),100),j=1,nrct)
          endif
          write(iow,2042) ((d(96+5*j+i),i=1,4),nint(d(101+5*j))/100,
     &                      mod(nint(d(101+5*j)),100),j=1,nrct)
        endif

c       Output shaped section data

        if(nint(d(100)).eq.3) then   ! Wide flange
          write(iow,2048) (d(100+i),i=1,6)
          if(ior.lt.0) then
            write(*,2048) (d(100+i),i=1,6)
          endif
        elseif(nint(d(100)).eq.4) then   ! Channel
          write(iow,2049) (d(100+i),i=1,6)
          if(ior.lt.0) then
            write(*,2049) (d(100+i),i=1,6)
          endif
        elseif(nint(d(100)).eq.5) then   ! Angle
          write(iow,2050) (d(100+i),i=1,4)
          if(ior.lt.0) then
            write(*,2050) (d(100+i),i=1,4)
          endif
        elseif(nint(d(100)).eq.6) then   ! Solid circle
          d(102) = max(1.d0,min(4.d0,d(102)))  ! quadr: 1 le nqd le 4
          write(iow,2059) d(101),nint(d(102))
          if(ior.lt.0) then
            write(*,2059) d(101),nint(d(102))
          endif
        endif

c       Output multidimensional constitution

        if(muld) then
          write(iow,2088) nint(d(188))
          if(ior.lt.0) then
            write(*,2088) nint(d(188))
          endif
        endif

c       Output plasticity parameters

        if(plasfl) then

c         Small deformation plasticity (also for one-d problems)

          if(sflag .or. oned) then

            if(nint(d(40)).eq.1) then
              write(iow,2003) (d(i),i=41,45)
              if(ior.lt.0) then
                write(*,2003) (d(i),i=41,45)
              endif
              if(d(180).gt.0.0d0) then
                write(iow,2075) d(180),nint(d(181))
                if(ior.lt.0) then
                  write(*,2075) d(180),nint(d(181))
                endif
              endif
              nn  = 1
              nhd = ntm + 1  ! epp, ep(i),i=1,ntm
            elseif(nint(d(40)).eq.3) then
              write(iow,2047) (d(i),i=41,45)
              if(ior.lt.0) then
                write(*,2047) (d(i),i=41,45)
              endif
              nn  = 2
              nhd = ntm + 1  ! epp, ep(i),i=1,ntm

c           Anisotropic plasticity

            elseif(nint(d(40)).eq.6) then

c             Plastic parameters

              if(nint(d(46)).eq.5) then ! Plastic saturation yield
                if(ior.lt.0) then
                  write(*,2102) (d(i),i=41,44)
                endif
                write(iow,2102) (d(i),i=41,44)
              elseif(nint(d(46)).eq.6) then ! Swift yield
                if(ior.lt.0) then
                  write(*,2103) (d(i),i=41,45)
                endif
                write(iow,2103) (d(i),i=41,45)
              endif

c             OUTPUT PLASTIC PARAMETERS

c             Orthogonalize plastic vectors

              if(nint(d(242)).eq.1) then
                call orthf_apl(d(243))
              else
                d(242) = 1.0d0
                d(243) = 1.0d0
                d(247) = 1.0d0
              endif

              if(ior.lt.0) then
                write(*,2104) (rr(i),i=1,6),(d(i),i=51,56),
     &                        (d(i),i=243,248)
              endif
              write(iow,2104) (rr(i),i=1,6),(d(i),i=51,56),
     &                        (d(i),i=243,248)

              nn  = 1
              nhd = 6 + 1 +1 ! ep(i),i=1,ntm, epp, cock
            endif

c         Plane stress reductions

          if(stype.eq.1) then

          endif

c         Finite deformation plasticity

          elseif(fflag) then

c           Isotropic Plasticity

            if(imat.eq.4) then
              nhd = 2 + 1 + 6 + 3 ! epp(2), dete, be(i),i=1,6,
                                  ! e_pl(i),i=1,3
              nn  = 1
              if(nint(d(46)).eq.1) then
                write(iow,2056) 'von Mises Yield Function',
     &                          (d(i),i=41,45)
                d(41) = sqt23*d(41)
                d(42) = sqt23*d(42)
              elseif(nint(d(46)).eq.2) then
                sigt  = d(41)
                sigc  = d(42)
                d(41) = sqrt( 8.d0/3.d0)*(sigc * sigt)/(sigt + sigc)
                d(42) = sqrt(18.d0/3.d0)*(sigc - sigt)/(sigt + sigc)
                write(iow,2057) sigt,sigc,d(41),d(42),d(44),d(45)
              elseif(nint(d(46)).eq.3) then
                sigt  = d(41)
                sigc  = d(42)
                d(41) = sqrt( 8.d0/3.d0)*(sigc * sigt)/(sigt + sigc)
                d(42) =                  (sigc - sigt)/(sigt + sigc)
                write(iow,2058) sigt,sigc,d(41),d(42),d(44),d(45)
                if(abs(d(42)).gt. 0.125d0) then
                  write(ilg,4013)
                  write(iow,4013)
                  if(ior.lt.0) then
                    write(*,4013)
                  endif
                endif
              elseif(nint(d(46)).eq.4) then
                write(iow,2056) 'Generalized Plasticity Mises Function',
     &                          (d(i),i=41,45)
                d(41) = sqt23*d(41)
                d(42) = sqt23*d(42)
              endif

c             Multiply hardening moduli by sqrt(2/3)

              d(44) = sqt23*d(44)
              d(45) = two3 *d(45)

c           Anisotropic Plasticity

            elseif(imat.eq.5) then

c             Modify and output alternative elastic parameters

              if(nint(d(46)).eq.5 .or. nint(d(46)).eq.6) then

c               Plastic parameters

                if(nint(d(46)).eq.5) then ! Plastic saturation yield
                  if(ior.lt.0) then
                    write(*,2102) (d(i),i=41,44)
                  endif
                  write(iow,2102) (d(i),i=41,44)
                elseif(nint(d(46)).eq.6) then ! Swift yield
                  if(ior.lt.0) then
                    write(*,2103) (d(i),i=41,45)
                  endif
                  write(iow,2103) (d(i),i=41,45)
                endif

c               OUTPUT PLASTIC PARAMETERS

c               Orthogonalize plastic vectors

                if(nint(d(242)).eq.1) then
                  call orthf_apl(d(243))
                else
                  d(243) = 1.0d0
                  d(247) = 1.0d0
                endif

                if(ior.lt.0) then
                  write(*,2104) (rr(i),i=1,6),(d(i),i=51,56),
     &                          (d(i),i=243,248)
                endif
                write(iow,2104) (rr(i),i=1,6),(d(i),i=51,56),
     &                          (d(i),i=243,248)

c               Set memory for history data at each point

                nn  = 1
                nhd = 14  ! Cockroft-Lathem in slot 14

              endif

c           Error on inputs

            else

              write(ilg,4014)
              write(iow,4014)
              call plstop()

            endif

          endif ! fflag/sflag

c         One dimensional model history storage

          if(oned) then
            if(nseg.eq.0) then
              nh1 = nh1 + nn + 3   ! ep(1); epp; state each layer
            else
              write(iow,2041) (d(i),i=131,130+3*nseg)
              if(ior.lt.0) then
                write(*,2041) (d(i),i=131,130+3*nseg)
              endif
              nh1 = nh1 + 4   ! ep(1); alp(1); epp; state each layer
            endif

c         Multi dimensional model history storage

          else
            if(ndm.eq.3 .or. stype.eq.8) then
              nh1 = nh1 + nhd*nn + 1
            else
              if(stype.eq.1 .and. sflag) then  ! Plane stress
                nh1 = nh1 + nn + 7   ! ep(3), beta(3), ep; state
              else
                if(fflag) then
                  nh1 = nh1 + nhd*nn + 1
                else                           ! Plane strain
!                 nh1 = nh1 + 5*nn + 1 ! ep(4); ep; state
                  nh1 = nh1 + nhd*nn + 1 ! ep(ntm); ep; state
                endif
              endif
            endif
          endif

        endif ! plasfl

c       Output visco-elastic properties

        if(viscfl) then
          write(iow,2004) (d(i),i=51,50+2*nvis)
          if(ior.lt.0) then
            write(*,2004) (d(i),i=51,50+2*nvis)
          endif
          if(fflag) then
            i = 1
          else
            i = 0
          endif
          if(oned) then
            nh1 = nh1 + i + 1 + nvis   ! xi; eps_n; hh(*) each point
          else
            if(ndm.eq.3 .or. stype.eq.8) then
              nh1 = nh1 + i + 6 + 6*nvis  ! xi; eps_n(6); hh(6) @ point
            else
              nh1 = nh1 + i + 4 + 4*nvis  ! xi; eps_n(4); hh(4) @ point
            endif
          endif
        endif ! viscfl

c       Output damage properties

        if(damgfl) then
          write(iow,2028) d(58),d(59)
          if(ior.lt.0) then
            write(*,2028) d(58),d(59)
          endif
        endif ! damgfl

c       Output anisotropic structure vectors

        if(nstv.gt.0) then

          if(imat.eq.14 .or. fiberfl) then

            do j = 0,nstv-1
              if(nint(d(271+3*j)).eq.1) then
                write(iow,2100) 'Holzapfel-Gasser',
     &                          j+1,(d(271+3*j+i),i=1,2)
                if(ior.lt.0) then
                  write(*,2100) 'Holzapfel-Gasser',
     &                          j+1,(d(271+3*j+i),i=1,2)
                endif
              elseif(nint(d(271+3*j)).eq.2) then
                write(iow,2100) 'Weiss',
     &                          j+1,(d(271+3*j+i),i=1,2)
                if(ior.lt.0) then
                  write(*,2100) 'Weiss',
     &                          j+1,(d(271+3*j+i),i=1,2)
                endif
              endif
            end do ! j

          endif

          write(iow,2101) (j+1,(d(260+3*j+i),i=1,3),j=0,nstv-1)
          if(ior.lt.0) then
            write(*,2101) (j+1,(d(260+3*j+i),i=1,3),j=0,nstv-1)
          endif

        endif

c       Output Piezo-electric properties

        if(d(150).gt.0.0d0) then
          if(tdof.eq.0) tdof = ndm + 1
          write(iow,2043) tdof,(d(i),i=151,155 + nint(d(150)))
          if(ior.lt.0) then
            write(*,2043) tdof,(d(i),i=151,155 + nint(d(150)))
          endif
        endif

c       Output quadrature data

        if(.not.qdflg .and. qdinp) then
          write(iow,2019) qd(ti),ii,jj
          if(ior.lt.0) then
            write(*,2019) qd(ti),ii,jj
          endif
        endif

c       Constitutive start (.ne.0 for nonclassical elastic)

        if(d(84).eq.0.0d0) then
          write(iow,2066) 'Elastic'
          if(ior.lt.0) then
            write(*,2066) 'Elastic'
          endif
        else
          write(iow,2066) 'Inelastic'
          if(ior.lt.0) then
            write(*,2066) 'Inelastic'
          endif
        endif

c       Kinematics type

        if(ietype.ne.6) then
          if(dtype.gt.0) then
            write(iow,2005)
            if(ior.lt.0) then
              write(*,2005)
            endif
          else
            i = min(etype,2)
            write(iow,2006) ty(i)
            if(ior.lt.0) then
              write(*,2006) ty(i)
            endif
          endif
        endif

c       Hierarchical formulation

        if(d(120).gt.0) then

          write(iow,2081)
          if(ior.lt.0) then
            write(*,2081)
          endif

        endif

c       Element type

        if(ietype.eq.1) then
          if(etype.eq.1) then
            write(iow,2012) 'Displacement'
            if(ior.lt.0) then
              write(*,2012) 'Displacement'
            endif
          elseif(etype.eq.2) then
            write(iow,2012) 'Mixed'
            if(ior.lt.0) then
              write(*,2012) 'Mixed'
            endif
          elseif(etype.eq.3) then
            write(iow,2012) 'Enhanced'
            if(ior.lt.0) then
              write(*,2012) 'Enhanced'
            endif
          elseif(etype.eq.4) then
            write(iow,2012) 'Energy Conserving'
            if(ior.lt.0) then
              write(*,2012) 'Energy Conserving'
            endif
          elseif(etype.eq.5) then
            write(iow,2012) 'Mixed-Enhanced'
            if(ior.lt.0) then
              write(*,2012) 'Mixed-Enhanced'
            endif
          elseif(etype.eq.6) then
            write(iow,2012) 'Co-Rotational'
            if(ior.lt.0) then
              write(*,2012) 'Co-Rotational'
            endif
          elseif(etype.eq.7) then
            write(iow,2012) 'Displacement Uniform Deformation Gradient'
            if(ior.lt.0) then
              write(*,2012) 'Displacement Uniform Deformation Gradient'
            endif
          elseif(etype.eq.8) then
            write(iow,2012) 'Mixed Uniform Deformation Gradient'
            if(ior.lt.0) then
              write(*,2012) 'Mixed Uniform Deformation Gradient'
            endif
          elseif(etype.eq.9) then
            write(iow,2012) 'Mixed U-p'
            if(ior.lt.0) then
              write(*,2012) 'Mixed U-p'
            endif
          endif
          if(incomp) then
            write(iow,2067)
            if(ior.lt.0) then
              write(*,2067)
            endif
          endif
        elseif(ietype.eq.3) then
          if(dtype.gt.0) then
            if(nint(d(79)).eq.0.0d0 .and. etype.ne.3) then
              write(iow,2012) 'Linear Displacement'
              if(ior.lt.0) then
                write(*,2012) 'Linear Displacement'
              endif
            else
              write(iow,2012) 'Cubic Displacement'
              if(ior.lt.0) then
                write(*,2012) 'Cubic Displacement'
              endif
            endif
          else
            if(nint(d(79)).eq.0.0d0) then
              write(iow,2012) 'Linear Displacement'
              if(ior.lt.0) then
                write(*,2012) 'Linear Displacement'
              endif
              if(ndm.eq.3) then
                if(etype.eq.2 .or. etype.eq.3) then
                  write(iow,2013) 'Incremental Rotation'
                  if(ior.lt.0) then
                    write(*,2013) 'Incremental Rotation'
                  endif
                elseif(etype.eq.4) then
                  write(iow,2013) 'Energy Conserving'
                  if(ior.lt.0) then
                    write(*,2013) 'Energy Conserving'
                  endif
                elseif(etype.eq.6) then
                  write(iow,2013) 'Co-Rotational'
                  if(ior.lt.0) then
                    write(*,2013) 'Co-Rotational'
                  endif
                else
                  write(iow,2013) 'Simo-Vu Quoc Displacement'
                  if(ior.lt.0) then
                    write(*,2013) 'Simo-Vu Quoc Displacement'
                  endif
                endif
              endif
            else
              write(iow,2012) 'Cubic Displacement, 2nd Order'
              if(ior.lt.0) then
                write(*,2012) 'Cubic Displacement, 2nd Order'
              endif
            endif
          endif
        endif

c       Shear deformation on/off

        if(ietype.eq.3 .or. ietype.eq.4) then
          if(d(79).gt.0.0d0) then
            write(iow,2012) 'No Shear Deformation'
            if(ior.lt.0) then
              write(*,2012) 'No Shear Deformation'
            endif
          else
            write(iow,2012) 'Includes Shear Deformation'
            if(ior.lt.0) then
              write(*,2012) 'Includes Shear Deformation'
            endif
          endif
        endif

c       Output plot indicator

        if(nint(d(171)).eq.1) then
          if(ior.lt.0) then
            write(*,2093) 'Off'
          endif
          write(iow,2093) 'Off'
        else
          if(ior.lt.0) then
            write(*,2093) 'On'
          endif
          write(iow,2093) 'On'
        endif

        if(histpltfl) then
          call outhisplt(mr(np(303)))
        endif

c       Output augmentation option

        if(etype.gt.1) then
          if(d(185).gt.0.0d0) then
            write(iow,2084) 'On'
            if(ior.lt.0) then
              write(*,2084) 'On'
            endif
          else
            write(iow,2084) 'Off'
            if(ior.lt.0) then
              write(*,2084) 'Off'
            endif
          endif
        endif

c     Thermal element only

      elseif(ietype.eq.6) then

        write(iow,2020) wd(stype),(d(i),i=61,64),d(66),d(4)
        if(d(127).gt.0.0d0) then
          write(iow,2073) d(127),d(128)
        endif
        if(ior.lt.0) then
          write(*,2020) wd(stype),(d(i),i=61,64),d(66),d(4)
          if(d(127).gt.0.0d0) then
            write(*,2073) d(127),d(128)
          endif
        endif
        if(qdinp) then
          if(d(182).gt.0.0d0) then
            j = 2
          else
            j = 1
          endif
          write(iow,2019) qd(j),ii,jj
          if(ior.lt.0) then
            write(*,2019) qd(j),ii,jj
          endif
        endif

      endif

c     Interpolation type

      if(nint(d(189)).eq.1) then
        write(iow,2090) 'Global NURBS',(i,nint(d(189+i)),i=1,ndm)
        if(ior.lt.0) then
          write(*,2090) 'Global NURBS',(i,nint(d(189+i)),i=1,ndm)
        endif
      elseif(nint(d(189)).eq.2) then
        write(iow,2090) 'Local NURBS',(i,nint(d(189+i)),i=1,ndm)
        if(ior.lt.0) then
          write(*,2090) 'Local NURBS',(i,nint(d(189+i)),i=1,ndm)
        endif
      elseif(nint(d(189)).eq.3) then
        write(iow,2090) 'B-Spline',(i,nint(d(189+i)),i=1,ndm)
        if(ior.lt.0) then
          write(*,2090) 'B-Spline',(i,nint(d(189+i)),i=1,ndm)
        endif
      elseif(nint(d(189)).eq.4) then
        write(iow,2090) 'T-Spline',(i,nint(d(189+i)),i=1,ndm)
        if(ior.lt.0) then
          write(*,2090) 'T-Spline',(i,nint(d(189+i)),i=1,ndm)
        endif
      elseif(nint(d(189)).eq.5) then
        write(iow,2090) 'Bezier  ',(i,nint(d(189+i)),i=1,ndm)
        if(ior.lt.0) then
          write(*,2090) 'Bezier  ',(i,nint(d(189+i)),i=1,ndm)
        endif
      elseif(nint(d(189)).eq.6) then
        write(iow,2090) 'NURBS Mixed',(i,nint(d(189+i)),i=1,ndm)
        if(ior.lt.0) then
          write(*,2090) 'NURBS Mixed',(i,nint(d(189+i)),i=1,ndm)
        endif
      elseif(nint(d(189)).eq.7) then
        write(iow,2090) 'H-Spline',(i,nint(d(189+i)),i=1,ndm)
        if(ior.lt.0) then
          write(*,2090) 'H-Spline',(i,nint(d(189+i)),i=1,ndm)
        endif
      endif

c     Mass type

      if(d(4).gt.0.0d0) then
        if(d(7).eq.0.0d0) then
          write(iow,2007)
          if(ior.lt.0) then
            write(*,2007)
          endif
        elseif(d(7).eq.1.0d0) then
          write(iow,2008)
          if(ior.lt.0) then
            write(*,2008)
          endif
        else
          write(iow,2009) d(7)
          if(ior.lt.0) then
            write(*,2009) d(7)
          endif
        endif
        if(ietype.eq.3 .or. ietype.eq.5.and.dtype.lt.0) then
          write(iow,2032) d(8)
          if(ior.lt.0) then
            write(*,2032) d(8)
          endif
        endif
        if(d(183).ne.0.0d0) then
          write(iow,2082) d(183)
          if(ior.lt.0) then
            write(*,2082) d(183)
          endif
        endif
      endif

c     Damping factor

      if(d(70).gt.0.0d0) then
        write(iow,2046) d(70)
        if(ior.lt.0) then
          write(*,2046) d(70)
        endif
      endif

c     Rayleigh Damping factors

      if(max(abs(d(77)),abs(d(78))).gt.0.0d0) then
        write(iow,2060) d(77),d(78)
        if(ior.lt.0) then
          write(*,2060) d(77),d(78)
        endif
      endif

c     Number of internal element (Lagrange multiplier) variables

      if(nint(d(230)).ne.0) then
        write(iow,2069) nint(d(230))
        if(ior.lt.0) then
          write(iow,2069) nint(d(230))
        endif
      endif

c     Number of active element degrees of freedom per node

      if(nint(d(241)).gt.0) then
        write(iow,2098) nint(d(241))
        if(ior.lt.0) then
          write(iow,2098) nint(d(241))
        endif
      endif

c     Output nodal based formulation

      if(nint(d(240)).ne.0) then
        write(iow,2099)
        if(ior.lt.0) then
          write(iow,2099)
        endif
      endif

c     Method type factors

      if(max(abs(d(80)),abs(d(81))).gt.0.0d0) then
        write(iow,2061) d(80),d(81)
        if(ior.lt.0) then
          write(*,2061) d(80),d(81)
        endif
      endif

c     Error value

      if(d(50).ne.0.0d0) then
        write(iow,2011) d(50)
        if(ior.lt.0) then
          write(*,2011) d(50)
        endif
      endif

      if(d(39).ne.0.0d0) then

        write(iow,2016)
        if(ior.lt.0) then
          write(*,2016)
        endif
      endif

c     Output penalty value

      if(d(60).ne.0.0d0) then
        write(iow,2022) d(60)
        if(ior.lt.0) then
          write(*,2022) d(60)
        endif
      endif

c     Set beam/shell kappa

      if((ietype.eq.4 .or. ietype.eq.5)      .and.
     &   (nint(d(79)).eq.0) .and. dtype.lt.0) then
        write(iow,2017) d(37)
        if(ior.lt.0) then
          write(*,2017) d(37)
        endif
      endif

c     Output reference coordinate/vector

      if    (sref.eq.1) then

        d(93) = sref
        do i = 1,2
          d(93+i) = tref(i)
        end do ! i

        write(iow,2031) (i,tref(i),i=1,2)
        if(ior.lt.0) then
          write(*,2031) (i,tref(i),i=1,2)
        endif

      endif

c     Output reference coordinate/vector

      if(lref.gt.0) then
        d(96) = lref
        do i = 1,3
          d(96+i) = refx(i)
        end do ! i
      endif
      if    (lref.eq.1) then

        write(iow,2025) (i,refx(i),i=1,ndm)
        if(ior.lt.0) then
          write(*,2025) (i,refx(i),i=1,ndm)
        endif

      elseif(lref.eq.2) then

        write(iow,2026) (i,refx(i),i=1,ndm)
        if(ior.lt.0) then
          write(*,2026) (i,refx(i),i=1,ndm)
        endif

      elseif(lref.eq.3) then

        write(iow,2079)
        if(ior.lt.0) then
          write(*,2079)
        endif

      elseif(lref.eq.4) then

        write(iow,2064)
        if(ior.lt.0) then
          write(*,2064)
        endif

      endif

c     Output orientation vectors

      if(aref.gt.0) then
        d(249)     = aref
        d(250:252) = refa(:,1)
        d(253:255) = refa(:,2)
        if(aref.eq.1) then        ! cylindrical

          write(iow,2079)
          write(iow,2105) 1,(refa(i,1),i=1,ndm)
          write(iow,2105) 2,(refa(i,2),i=1,ndm)
          if(ior.lt.0) then
            write(*,2079)
            write(*,2105) 1,(refa(i,1),i=1,ndm)
            write(*,2105) 2,(refa(i,2),i=1,ndm)
          endif

        endif
      endif

c     Output number global equations

      if(nint(d(283)).gt.0) then
        write(iow,2106) nint(d(283))
        if(ior.lt.0) then
          write(*,2106) nint(d(283))
        endif
      endif

c     Output rigid body indicator

      if(nrmatl.gt.0) then
        write(iow,2080) nrmatl
        if(ior.lt.0) then
          write(*,2080) nrmatl
        endif
      endif

c     Output user parameters

      if(uparfl) then

        j = 0
        do i = 1,10
          if(d(230+i).ne.0.0d0) then
            j = i
          endif
        end do ! i

        write(iow,2068) (i,d(230+i),i=1,j)
        if(ior.lt.0) then
          write(*,2068) (i,d(230+i),i=1,j)
        endif

      endif

c     Add history storage for multi-dimensional reduction case

      if(muld) then
        nh1 = nh1 + nint(d(188))
      endif

c     Save quadrature

      if(qdinp .or. (ietype.ne.1 .and. ietype.ne.6)) then
        d(5) = ii
        d(6) = jj
      endif

c     Save types

      istp   = stype
      d(3)   = alp(1) + alp(2) + alp(3)
      d(15)  = nh1 + naug ! nh
      d(16)  = stype
      d(17)  = etype
      d(18)  = dtype
      d(19)  = tdof
      d(20)  = imat
      d(31)  = d(31)*rad
      d(193) = tmat
      d(194) = ietype

c     Solid (1) or shell/membrane (5)

      if(ietype15) then
        nh2   = nh1

c       Set 2 and 3 dimensional quadrature order for FE elements

        jj    = ii
        j     = ii

c       Set 2 and 3 dimensional quadrature order for NURB elements

        if(nint(d(189)).gt.0) then
          ii = nint(d(190))
          jj = nint(d(191))
          j  = nint(d(192))
        endif

        if(ndm.eq.2) then
          nh1   = nh1*ii*jj
          if(etype.eq.5) nh1 = nh1 + 2 ! mixed-enhance unknown
        elseif(ndm.eq.3) then
          nh1   = nh1*ii*jj*j
          if(etype.eq.5) nh1 = nh1 + 9 ! mixed-enhance unknown
        endif
        if(etype.eq.3 .and. sflag) nh1 = nh1 + 5     ! linear element
        if(etype.eq.3 .and. fflag) nh1 = nh1 + nh2   ! extra quadrature
        if(stype.eq.1) then
          nh1   = nh1 + ii*jj ! store F_33/q_pt
          d(15) = d(15) + 1.0d0
        endif

c     Truss (2) or Frame (3)

      elseif(ietype.eq.2 .or. ietype.eq.3) then
        d(82) = max(1,nqud)
        if(nint(d(79)).eq.1) then
          nh1 = nh1*3
        endif
        if((muld.or.oned).and.nlay.gt.1) then
          nh1 = nh1*((nlay-1)*(nlob-1) + 1)  ! Total data/cross section
        elseif((muld.or.oned).and.ntub.gt.1) then
          nh1 = nh1*ntub*(nlob-1)            ! Total data/cross section
        elseif((muld.or.oned).and.nqud.gt.1) then
          nh1 = nh1*nqud                     ! Total data/cross section
        endif

c     Plate (4)

      elseif(ietype.eq.4) then
        nh1   = nh1*ii

c     Thermal (6)

      elseif(ietype.eq.6) then
        nh1   = nh1*ii**ndm

c     Three dimensional (7)

      elseif(ietype.eq.7) then
        nh1   = nh1*ii**3
      endif

c     Set history for saving element variables

      d(149) = nh1 ! Total number of variables on frame cross section
      nh3    = nh3 + ndv

c     Set augmenting base value

      d(185) = d(185)*d(21)

c     Set number of element variables

      nlm  = nint(d(230))

c     Check for warnings

      if(d(4).eq.0.0d0) then
        write(iow,3100)
        if(ior.lt.0) then
          write(*,3100)
        endif
      endif

c     Check for errors

      if     (ietype.eq.1) then
        if(e1.eq.0.0d0) then
          write(ilg,4001)
          write(iow,4001)
          if(ior.lt.0) then
            write(*,4001)
          endif
        endif
      elseif (ietype.eq.2) then
        if(e1.eq.0.0d0 .or. d(32).eq.0.0d0) then
          write(ilg,4002) 'Truss',e1,d(32)
          write(iow,4002) 'Truss',e1,d(32)
          if(ior.lt.0) then
            write(*,4002) 'Truss',e1,d(32)
          endif
        endif
      elseif (ietype.eq.3) then
        if(e1.eq.0.0d0 .or. d(32).eq.0.0d0 .or. d(33).eq.0.0d0) then
          write(ilg,4002) 'Frame',e1,d(32),d(33)
          write(iow,4002) 'Frame',e1,d(32),d(33)
          if(ior.lt.0) then
            write(*,4002) 'Frame',e1,d(32),d(33)
          endif
        endif
      elseif (ietype.eq.4) then
        if(e1.eq.0.0d0 .or. d(37).eq.0.0d0 .or. d(14).eq.0.0d0) then
          write(ilg,4003) 'Plate',e1,d(14),d(37)
          write(iow,4003) 'Plate',e1,d(14),d(37)
          if(ior.lt.0) then
            write(*,4003) 'Plate',e1,d(14),d(37)
          endif
        endif
      elseif (ietype.eq.5) then
        if(e1.eq.0.0d0 .or. d(37).eq.0.0d0 .or. d(14).eq.0.0d0) then
          write(ilg,4003) 'Shell',e1,d(14),d(37)
          write(iow,4003) 'Shell',e1,d(14),d(37)
          if(ior.lt.0) then
            write(*,4003) 'Shell',e1,d(14),d(37)
          endif
        endif
      elseif (ietype.eq.6) then
        if(t1.eq.0.0d0) then
          write(ilg,4004)
          write(iow,4004)
          if(ior.lt.0) then
            write(*,4004)
          endif
        endif
      endif

c     Output user history information

      if(uflag) then
        write(iow,5000) nint(d(15)),nh1,nh3
      endif

c     Errors detected in inputs

      if(erflag) call plstop()

c     I/O Formats

2000  format( 5x,'M e c h a n i c a l   P r o p e r t i e s'//
     & 10x,a,' Analysis'//
     & 10x,'Modulus E       ',1p,1e12.5/
     & 10x,'Poisson ratio   ',0p,1f8.5/)

2001  format( 5x,'M e c h a n i c a l   P r o p e r t i e s'//
     & 10x,a,' Analysis'//
     & 10x,'Modulus E-1     ',1p,1e12.5/
     & 10x,'Modulus E-2     ',1p,1e12.5/
     & 10x,'Modulus E-3     ',1p,1e12.5/
     & 10x,'Poisson ratio 12',0p,1f8.5 /
     & 10x,'Poisson ratio 23',0p,1f8.5 /
     & 10x,'Poisson ratio 31',0p,1f8.5 /
     & 10x,'Modulus G-12    ',1p,1e12.5/
     & 10x,'Modulus G-23    ',1p,1e12.5 /
     & 10x,'Modulus G-31    ',1p,1e12.5/
     & 10x,'Angle (psi)   ',0p,1f10.5/)

2002  format(/5x,'T h e r m a l   E x p a n s i o n s'//
     & 10x,'Th. Alpha-1',1p,1e17.5/10x,'Th. Alpha-2',1p,1e17.5/
     & 10x,'Th. Alpha-3',1p,1e17.5/10x,'T_0        ',1p,1e17.5/
     & 10x,'Th. D.O.F. ',i7/)

2003  format(/5x,'M i s e s   P l a s t i c   P a r a m e t e r s'//
     &  10x,'Yield stress short  ',1p,1e15.5/
     &  10x,'Yield stress infin. ',1p,1e15.5/
     &  10x,'Hardening exponent  ',1p,1e15.5/
     &  8x,'Linear hardening parts'/
     &  10x,'Isotropic hardening ',1p,1e15.5/
     &  10x,'Kinematic hardening ',1p,1e15.5/)

2004  format(/5x,'V i s c o e l a s t i c   P a r a m e t e r s'//
     &   (10x,'nu-i :ratio',1p,1e15.5/10x,'lam-i:time ',1p,1e15.5:))

2005  format(/5x,'E l e m e n t    M o d e l'//
     &       10x,'Formulation  : Small Deformation.')
2006  format(/5x,'E l e m e n t    M o d e l'//
     &       10x,'Formulation  : Finite Deformation ',a)

2007  format( 10x,'Mass type    : Lumped.')
2008  format( 10x,'Mass type    : Consistent.')
2009  format( 10x,'Mass type    : Interpolated:',1p,1e11.4)

2010  format( 5x,'M e c h a n i c a l   P r o p e r t i e s'//
     &        4x,a,'NeoHookean Stored Energy Function '//
     &       10x,a,' Analysis'//
     &       10x,'Modulus E       ',1p,1e12.5/
     &       10x,'Poisson ratio   ',0p,1f8.5/
     &       10x,'Bulk Modulus    ',1p,1e12.5/
     &       10x,'Shear Modulus   ',1p,1e12.5/1x)

2011  format(10x,'Error estimator ',1p,1e12.5 )

2012  format(10x,'Element type : ',a)

2013  format( 25x,a:/,25x,'Shear factor (k) =',1p,1e12.5)

2014  format(/10x,'Augmented Solution for Inextensible Behavior')

2015  format(/5x,'C r o s s   S e c t i o n   P a r a m e t e r s'//
     &    10x,'Area      ',1p,1e15.5:/10x,'I_xx      ',1p,1e15.5/
     &    10x,'I_yy      ',1p,1e15.5 /10x,'I_xy      ',1p,1e15.5/
     &    10x,'J_zz      ',1p,1e15.5 /10x,'Kappa_x   ',1p,1e15.5/
     &    10x,'Kappa_y   ',1p,1e15.5)

2016  format( 10x,'Non-linear analysis')

2017  format( 10x,'Plate/Shell  : Kappa ',1p,1e15.5)

2018  format( 10x,'Thickness       ',1p,1e12.5)
2019  format(/5x,a,'   Q u a d r a t u r e'/
     & 10x,'Quadrature: Arrays',i3    /10x,'Quadrature: Output',i3/)

2020  format(/5x,'T h e r m a l   P r o p e r t i e s'//
     & 10x,a,' Analysis'//
     & 10x,'Cond. K-1    ',1p,1e15.5/ 10x,'Cond. K-2    ',1p,1e15.5/
     & 10x,'Cond. K-3    ',1p,1e15.5/ 10x,'Specific Heat',1p,1e15.5/
     & 10x,'Heat Source  ',1p,1e15.5/:10x,'Density      ',1p,1e15.5)

2021  format( 10x,'Loading - q ',1p,1e16.5/
     &        10x,'q Prop. load',i6)

2022  format( 10x,'Penalty - k ',1p,1e16.5)

2023  format(
     &    /5x,'P r o p o r t i o n a l   B o d y   L o a d i n g s',//
     &    10x,'1-Dir. Factor',1p,1e15.5,': Proportional Load No.',i3/:
     &    10x,'2-Dir. Factor',1p,1e15.5,': Proportional Load No.',i3/:
     &    10x,'3-Dir. Factor',1p,1e15.5,': Proportional Load No.',i3/)

2024  format( 5x,'M e c h a n i c a l   P r o p e r t i e s'//
     &        5x,'Ogden Stored Energy Function '//
     &       10x,a,' Analysis'//
     &       10x,'Bulk  Modulus   ',1p,1e12.5//
     &       10x,'Modulus  - C1   ',1p,1e12.5/
     &       10x,'Exponent - n1   ',1p,1e12.5/:/
     &       10x,'Modulus  - C2   ',1p,1e12.5/
     &       10x,'Exponent - n2   ',1p,1e12.5/:/
     &       10x,'Modulus  - C3   ',1p,1e12.5/
     &       10x,'Exponent - n3   ',1p,1e12.5/1x)

2025  format(/5x,'R e f e r e n c e    C o o r d i n a t e s'//
     &        (10x,'X-',i1,' = ',1p,1e12.5:))

2026  format(/5x,'R e f e r e n c e    V e c t o r'//
     &        (10x,'v-',i1,' = ',1p,1e12.5:))

2027  format( 5x,'M e c h a n i c a l   P r o p e r t i e s'//
     &        5x,'Logarithmic Stretch Stored Energy Function '//
     &       10x,a,' Analysis'//
     &       10x,'Modulus E       ',1p,1e12.5/
     &       10x,'Poisson ratio   ',0p,1f8.5/
     &       10x,'Bulk Modulus    ',1p,1e12.5/
     &       10x,'Shear Modulus   ',1p,1e12.5/1x)

2028  format(/5x,'D a m a g e    P a r a m e t e r s'//
     &    10x,'Damage Limit',1p,1e13.5/10x,'Decay rate',1p,1e15.5)

2029  format(10x,'Density         ',1p,1e12.5//
     &       10x,'1-Gravity Load  ',1p,1e12.5/
     &       10x,'2-Gravity Load  ',1p,1e12.5/
     &       10x,'3-Gravity Load  ',1p,1e12.5)

2030  format(/10x,'Angular Velocity (radians/time)',1p,1e15.5/)

2031  format(/5x,'S h e a r   C e n t e r   C o o r d i n a t e s'//
     &        (10x,'X-',i1,' = ',1p,1e12.5:))

2032  format( 10x,'Rotational Mass Factor:',1p,1e12.5)

2033  format( 15x,'Follower loading.')
2034  format( 10x,'Global axis loading.')
2035  format( 10x,'Local axis loading.')
2036  format( 10x,'Global axis follower loading.')
2037  format( 10x,'Local axis follower loading.')

2038  format(/5x,'B e a m   L a y e r   D a t a'/
     &       10x,'Layer  Z-location        Width'/(i14,1p2e13.4))

2039  format(10x,'Number of Lobatto points/layer =',i3)

2040  format(/5x,'T u b u l a r   B e a m   D a t a'/
     &       10x,'Radius            =',1p,e13.4/
     &       10x,'Thickness         =',1p,e13.4/
     &       10x,'No. Sectors       =',i6/
     &       10x,'Quad. Pts/Sectors =',i6/)

2041  format(/5x,'H a r d e n i n g    P a r a m e t e r s'//
     &    12x,'Strain e_p  Isotropic Yield-Y  Kinematic Hardening-H'/
     &     (11x,1p,1e12.5,1p,1e19.5,1p,1e22.5))

2042  format(/5x,'R e c t a n g u l a r   B e a m   D a t a'/
     &       10x,'  y-low left  z-low left  y-up right  z-up right',
     &           '  y-quadr  z-quadr'/(10x,1p,4e12.4,i7,i9))

2043  format(/5x,'P i e z o - E l e c t r i c   D a t a'//
     &       10x,'Voltage D.O.F.    =',i3//
     &       10x,'Permeability -  1 =',1p,e13.4/
     &       10x,'Permeability -  2 =',1p,e13.4/
     &       10x,'Permeability -  3 =',1p,e13.4//
     &       10x,'Coupling(+)  e_12 =',1p,e13.4/
     &       10x,'Coupling(+)  e_22 =',1p,e13.4/
     &       10x,'Coupling(+)  e_45 =',1p,e13.4/:
     &       10x,'Coupling(-)  e_12 =',1p,e13.4/
     &       10x,'Coupling(-)  e_22 =',1p,e13.4/
     &       10x,'Coupling(-)  e_45 =',1p,e13.4/)

2044  format(/5x,'I n i t i a l   S t r e s s   D a t a'//
     &       10x,'11-Stress         =',1p,e13.4/:
     &       10x,'22-Stress         =',1p,e13.4/
     &       10x,'33-Stress         =',1p,e13.4/
     &       10x,'12-Stress         =',1p,e13.4/
     &       10x,'23-Stress         =',1p,e13.4/
     &       10x,'31-Stress         =',1p,e13.4/)

2045  format(/5x,'I n i t i a l   S t r a i n   D a t a'//
     &       10x,'11-Strain         =',1p,e13.4/:
     &       10x,'22-Strain         =',1p,e13.4/
     &       10x,'33-Strain         =',1p,e13.4/
     &       10x,'12-Strain         =',1p,e13.4/
     &       10x,'23-Strain         =',1p,e13.4/
     &       10x,'31-Strain         =',1p,e13.4/)

2046  format( 10x,'Damping     ',1p,1e16.5)

2047  format(/5x,'M i s e s   G e n e r a l i z e d   ',
     &           'P l a s t i c   P a r a m e t e r s'//
     &       10x,'Yield stress        ',1p,1e15.5,' (Sigma_0)'/
     &       10x,'Yield stress infin. ',1p,1e15.5,' (Sigma_inf)'/
     &       10x,'Transition parameter',1p,1e15.5,' (Delta)'/
     &        8x,'Linear hardening parts'/
     &       10x,'Isotropic hardening ',1p,1e15.5,' (H_iso)':/
     &       10x,'Kinematic hardening ',1p,1e15.5,' (H_kin)'/)

2048  format(/5x,'W i d e f l a n g e   B e a m   D a t a'/
     &       10x,'Height             ',1p,1e15.5/
     &       10x,'Top    flange width',1p,1e15.5/
     &       10x,'Bottom flange width',1p,1e15.5/
     &       10x,'Top    flange thick',1p,1e15.5/
     &       10x,'Bottom flange thick',1p,1e15.5/
     &       10x,'Web thickness      ',1p,1e15.5/)

2049  format(/5x,'C h a n n e l   B e a m   D a t a'/
     &       10x,'Height             ',1p,1e15.5/
     &       10x,'Top    flange width',1p,1e15.5/
     &       10x,'Bottom flange width',1p,1e15.5/
     &       10x,'Top    flange thick',1p,1e15.5/
     &       10x,'Bottom flange thick',1p,1e15.5/
     &       10x,'Web thickness      ',1p,1e15.5/)

2050  format(/5x,'A n g l e   B e a m   D a t a'/
     &       10x,'Height             ',1p,1e15.5/
     &       10x,'Width              ',1p,1e15.5/
     &       10x,'Height thickness   ',1p,1e15.5/
     &       10x,'Width  thickness   ',1p,1e15.5/)

2051  format( 10x,'Thickness pts',i12)

2052  format( 10x,'Tension only material')

2053  format( 10x,'Compression only material')

2054  format( 5x,'M e c h a n i c a l   P r o p e r t i e s'//
     & 10x,a,' Analysis'//5x,'Fung Exponential Model parameters:'/
     & 10x,'Modulus C   ',    1p,1e16.5/
     & 10x,'A-11        ',    1p,1e16.5/
     & 10x,'A-22        ',    1p,1e16.5/
     & 10x,'A-12        ',    1p,1e16.5/
     & 10x,'A-44        ',    1p,1e16.5/)

2055  format( 5x,'M e c h a n i c a l   P r o p e r t i e s'//
     & 10x,a,' Analysis'//5x,'Fung Exponential Model parameters:'/
     & 10x,'Modulus C   ',    1p,1e16.5/
     & 10x,'A-11        ',    1p,1e16.5/
     & 10x,'A-22        ',    1p,1e16.5/
     & 10x,'A-33        ',    1p,1e16.5/
     & 10x,'A-12        ',    1p,1e16.5/
     & 10x,'A-23        ',    1p,1e16.5/
     & 10x,'A-31        ',    1p,1e16.5/
     & 10x,'A-44        ',    1p,1e16.5/
     & 10x,'A-55        ',    1p,1e16.5/
     & 10x,'A-66        ',    1p,1e16.5/)

2056  format(/8x,a/
     &       10x,'Yield stress      (initial) :  ',1p,e12.5/
     &       10x,'Yield stress     (infinity) :  ',1p,e12.5/
     &       10x,'Hardening exponent  (delta) :  ',1p,e12.5/
     &        8x,'Linear hardening parts'/
     &       10x,'Isotropic hardening (H_iso) :  ',1p,1e12.5/
     &       10x,'Kinematic hardening (H_kin) :  ',1p,1e12.5/)

2057  format(/8x,'Drucker-Prager Yield function'/
     &       10x,'Yield stress in tension     :  ',1p,e12.5/
     &       10x,'Yield stress in compression :  ',1p,e12.5/
     &       10x,'Yield function radius       :  ',1p,e12.5/
     &       10x,'Alpha parameter             :  ',1p,e12.5/
     &        8x,'Linear hardening parts'/
     &       10x,'Isotropic hardening (H_iso) :  ',1p,1e12.5/
     &       10x,'Kinematic hardening (H_kin) :  ',1p,1e12.5/)

2058  format(/8x,'Prager-Lode Yield function'/
     &       10x,'Yield stress in tension     :  ',1p,e12.5/
     &       10x,'Yield stress in compression :  ',1p,e12.5/
     &       10x,'Yield function radius       :  ',1p,e12.5/
     &       10x,'Lode angle parameter        :  ',1p,e12.5/
     &        8x,'Linear hardening parts'/
     &       10x,'Isotropic hardening (H_iso) :  ',1p,1e12.5/
     &       10x,'Kinematic hardening (H_kin) :  ',1p,1e12.5/)

2059  format(/5x,'C i r c u l a r   B e a m   D a t a'/
     &       10x,'Radius             ',1p,1e15.5/
     &       10x,'Quadrature order   ',i9/
     &       15x,'1 =  4 point (4-interior)'/
     &       15x,'2 =  5 point (4-perimeter, center)'/
     &       15x,'3 =  9 point (4-perimeter, 4-interior, center)'/
     &       15x,'4 = 17 point (8-perimeter, 8-interior, center)')

2060  format(/8x,'Rayleigh Damping Ratios'/
     &       10x,'Mass  value: a0',1p,1e14.5/
     &       10x,'Stiff value: a1',1p,1e14.5)

2061  format(/8x,'Method Type Values'/
     &       10x,'Type 1: Value  ',1p,1e14.5/
     &       10x,'Type 2: Value  ',1p,1e14.5)

2062  format(10x,'Volume Model    ',i2/
     &       15x,'1: U(J) = 0.25*(J**2 - 1 - 2 * ln J)'/
     &       15x,'2: U(J) = 0.50*(J - 1)**2'/
     &       15x,'3: U(J) = 0.50*(ln J)**2'/
     &       15x,'4: U(J) = 2.00*(J - 1 - ln J)'/1x)

2063  format( 5x,'M e c h a n i c a l   P r o p e r t i e s'//
     &       10x,' Moduli From Full Inputs'//
     &       10x,1p,1e11.3/10x,1p,2e11.3/10x,1p,3e11.3/
     &       10x,1p,4e11.3/10x,1p,5e11.3/10x,1p,6e11.3/)

2064  format(/10x,'Axial Reference Vector')

2065  format(/5x,'P o l a r   A n g l e   O r i g i n'/
     &        10x,'1-Coord   = ',1p,1e12.5/
     &        10x,'2-Coord   = ',1p,1e12.5:/
     &        10x,'3-Coord   = ',1p,1e12.5:)

2066  format(/5x,'C o n s t i t u t i v e    S t a r t'//
     &       10x,'Start state  : ',a)

2067  format( 10x,'Incompressible Formulation')

2068  format( 5x,'U s e r   P a r a m e t e r s'//
     &      (10x,'Parameter',i3,'   =',1p,e13.4:))

2069  format( 5x,'E l e m e n t   V a r i a b l e s'//
     &       10x,'Number/element    =',i3/)

2070  format( 10x,'Activation thermal')
2071  format( 10x,'Activation mechanical')

2072  format( 10x,'Constitutive equation transient type =',i5/
     &        15x,' 1 = Euler   Method'/
     &        15x,' 2 = Newmark Method'/
     &        15x,' 3 = User    Method')

2073  format(/5x,'S u r f a c e   C o n v e c t i o n'/
     &       10x,'Convection  (h)    ',1p,1e15.5/
     &       10x,'Temperature (T_inf)',1p,1e15.5)

2074  format( 5x,'M e c h a n i c a l   P r o p e r t i e s'//
     &        4x,a,'Mooney-Rivlin Stored Energy Function '//
     &       10x,a,' Analysis'//
     &       10x,'Modulus E       ',1p,1e12.5/
     &       10x,'Poisson ratio   ',0p,1f8.5/
     &       10x,'Bulk Modulus    ',1p,1e12.5/
     &       10x,'Shear Modulus   ',1p,1e12.5/
     &       10x,'2nd Invariant c ',1p,1e12.5/1x)

2075  format( 8x,'Viscoplastic parameters'/
     &       10x,'Rate time parameter ',1p,1e15.5/
     &       10x,'Yield function power',i8/)

2076  format( 8x,'Absolute Reference Temperature'//
     &       10x,'Temperature (T_ref)',1p,1e15.5/
     &       10x,'Heat fraction (r_H)',1p,1e15.5)

2077  format( 5x,'M e c h a n i c a l   P r o p e r t i e s'//
     &        4x,a,'Arruda-Boyce Stored Energy Function '//
     &       10x,a,' Analysis'//
     &       10x,'Modulus E     - E ',1p,1e12.5/
     &       10x,'Poisson ratio - v ',0p,1f8.5/
     &       10x,'Reciprocal n  - m ',1p,1e12.5/
     &       10x,'Bulk Modulus  - K ',1p,1e12.5/
     &       10x,'Shear Modulus - G ',1p,1e12.5/1x)

2078  format( 5x,'M e c h a n i c a l   P r o p e r t i e s'//
     &        4x,a,'Yeoh Stored Energy Function '//
     &       10x,a,' Analysis'//
     &       10x,'Modulus E     - E ',1p,1e12.5/
     &       10x,'Poisson ratio - v ',0p,1f8.5/
     &       10x,'Factor 1      - k1',0p,1f8.5/
     &       10x,'Factor 2      - k2',0p,1f8.5/
     &       10x,'Bulk Modulus  - K ',1p,1e12.5/
     &       10x,'Shear Modulus - G ',1p,1e12.5/1x)

2079  format(/5x,'C y l i n d r i c a l    R e f e r e n c e   ',
     &           'V e c t o r'/)

2080  format(/10x,'Material Set is Rigid Body Number',i5)

2081  format( 8x,'Hierarchical formulation')

2082  format( 10x,'Mass scaling factor  - b_m:',1p,1e12.5)

2083  format( 10x,'Wave speed - c_o',1p,1e12.5)

2084  format( 10x,'Augmenting   : ',a)

2085  format( 10x,'Explicit time integration')

2086  format( 10x,'Implicit time integration')

2087  format( 10x,'Cylindrical coordinate loading.')

2088  format(/8x,'Multi-dimensional constitution'/
     &       10x,'Number components ',i8)

2089  format(/5x,'Plastic Dissipation Parameters'/
     &       10x,'Easy glide dissipation (chi_g) ',1p,1e12.5/
     &       10x,'Hardening dissipation  (chi_h) ',1p,1e12.5/
     &       /5x,'Linear Thermal Softening Parameters'/
     &       10x,'Initial  yield softening (w_0) ',1p,1e12.5/
     &       10x,'Infinite yield softening (w_i) ',1p,1e12.5/
     &       10x,'Isotropic hardening soft (w_h) ',1p,1e12.5/
     &       10x,'Kinematic hardening soft (=0)  ',1p,1e12.5/
     &       /5x,'Thermal Split Option'/
     &       10x,'Isothermal (0); Adiabatic (>0) ',i8/)

2090  format(10x,'Interpolation: ',a/(15x,'Quadrature ',i1,' = ',i5:))

2091  format(10x,'Material behavior from RVE Analysis'/
     &       15x,'History/pt =',i5)

2092  format(10x,'RVE Model Filename : ',a)

2093  format(10x,'Plot Output  : ',a)

2094  format( 5x,'St. Venant-Kirchhoff   E l a s t i c   M o d e l')

2095  format(/10x,'Body Patch Loading'/
     &        10x,'1-Body Load     ',1p,1e12.5/
     &        10x,'2-Body Load     ',1p,1e12.5/)

2096  format(10x,'Normal Direction Moduli'/
     &       10x,'    |',1p,1e12.5,1p,2e13.5,' |'/
     &       10x,'D = |',1p,1e12.5,1p,2e13.5,' |'/
     &       10x,'    |',1p,1e12.5,1p,2e13.5,' |'/)

2097  format(10x,'Plane Stress in Thickness')

2098  format( 5x,'E l e m e n t   A c t i v e   V a r i a b l e s'//
     &       10x,'Number dof/node   =',i3/)

2099  format(10x,'Nodal Based Element Groups.')

2100  format(10x,a,' Parameters'//
     &      (10x,'Fiber',i2,': H1 =',1p,1e13.5,' H2 =',1p,1e13.5:))

2101  format(/5x,'Anistropic Structure vectors'//
     &      (10x,'Vector',i2,':',1p,3e13.5:))

2102  format(/5x,'Plastic Parameters          '//
     &       10x,'Yield Stress Y0 ',1p,1e13.5/
     &       10x,'Yield Stress Yi ',1p,1e13.5/
     &       10x,'Saturate,beta   ',1p,1e13.5/
     &       10x,'Iso-hardening   ',1p,1e13.5)

2103  format(/5x,'Swift Power Law Plastic Parameters'//
     &       10x,'Yield Stress    ',1p,1e13.5/
     &       10x,'Strength Coeff. ',1p,1e13.5/
     &       10x,'Init. Strain    ',1p,1e13.5/
     &       10x,'Hard Exponent   ',1p,1e13.5/
     &       10x,'Iso-hardening   ',1p,1e13.5)

2104  format(/5x,'Yield Stress Ratios'//
     &       10x,'Yield Ratio R_11',1p,1e13.5/
     &       10x,'Yield Ratio R_22',1p,1e13.5/
     &       10x,'Yield Ratio R_33',1p,1e13.5/
     &       10x,'Yield Ratio R_12',1p,1e13.5/
     &       10x,'Yield Ratio R_23',1p,1e13.5/
     &       10x,'Yield Ratio R_31',1p,1e13.5/
     &       /5x,'Hill Yield Parameters'//
     &       10x,'Yield Stress   F',1p,1e13.5/
     &       10x,'Yield Stress   G',1p,1e13.5/
     &       10x,'Yield Stress   H',1p,1e13.5/
     &       10x,'Yield Stress   L',1p,1e13.5/
     &       10x,'Yield Stress   M',1p,1e13.5/
     &       10x,'Yield Stress   N',1p,1e13.5/
     &       /5x,'Orthotropic Material Axes'//
     &        8x,'Vector 1'/
     &       10x,'Component 1     ',1p,1e13.5/
     &       10x,'Component 2     ',1p,1e13.5/
     &       10x,'Component 3     ',1p,1e13.5/
     &        8x,'Vector 2'/
     &       10x,'Component 1     ',1p,1e13.5/
     &       10x,'Component 2     ',1p,1e13.5/
     &       10x,'Component 3     ',1p,1e13.5)

2105  format(10x,'Point',i2,':',3(1p,1e13.5:))

2106  format(10x,'Global Equations ',i8)

3000  format(10x,'Thickness function =',i3)

3100  format(/5x,'W a r n i n g s   &   E r r o r s'/
     &       10x,'Material density is zero.')

4000  format(/' *ERROR* INPT: No elastic properties input.'/)

4001  format(/' *ERROR* INPT: Solid Element: Modulus zero.')

4002  format(/' *ERROR* INPT: ',a,': Modulus  = ',1p,1e12.5/
     &        '                      Area     = ',1p,1e12.5/:
     &        '                      Intertia = ',1p,1e12.5/)

4003  format(/' *ERROR* INPT: ',a,': Modulus  = ',1p,1e12.5/
     &        '                      Thickness= ',1p,1e12.5/
     &        '                      Kappa    = ',1p,1e12.5/)

4004  format(/' *ERROR* INPT: Thermal: Conductivity zero.')

4005  format(/' *ERROR* INPT: Elastic properties not positive definite'/
     &        '         Computed Compliances'/3(10x,1p,3e15.5/),
     &        '         Computed Moduli'/(10x,1p,3e15.5))

4006  format(/' *ERROR* INPT: Incorrect size for ',a,' =',i5/)

4007  format(/' *ERROR* INPT: Incorrect proportional load number'/)

4008  format(/' *ERROR* INPT: Incorrect reference vector or node',
     &        ' specification.'/)

4009  format(/' *ERROR* INPT: Inertia determinant zero or negative'/
     &       15x,'Det = I_xx*I_yy - I_xy*I_xy =',1p,1e12.5/)

4010  format(/' *ERROR* INPT: Too many layer sections: Limit = 13.'/)

4011  format(/' *ERROR* INPT: Too many rectangular sections:',
     &        ' Limit = 5.'/)

4012  format(/' *ERROR* INPT: Too many hardening segments: Limit = 6.'/)

4013  format(/' *WARNING* Non-convex yield function'/)

4014  format(/' *ERROR* INPT: Incorrect elastic model for plasticity',
     &       ' solution: Use ISOTropic'/)

4015  format(/' *ERROR* INPT: Too many viscoelastic terms: Limit = 3.'/)

4016  format(/' *ERROR* INPT: Layers are only for 2-d problems.'/)

4017  format(/' *ERROR* INPT: Area of SECTion RECTangle zero or',
     &        ' negative'/
     &        8x,' y_left  = ',1p,1e12.4/8x,' z_left  = ',1p,1e12.4/
     &        8x,' y_right = ',1p,1e12.4/8x,' z_right = ',1p,1e12.4/)

4018  format(/' *ERROR* INPT: Hill plasticity ratios must be postitive')

4019  format(/' *ERROR* INPT: Structure vector has zero length')

c     User function Formats

5000  format(//5x,'U s e r   H i s t o r y   I n f o r m a t i o n'//
     &        10x,'Number of history variables/quad point    ',i8/
     &        10x,'Number of history variables total         ',i8/
     &        10x,'Number of element history variables total ',i8/)

      end

      subroutine sethisplt(hpltab,ev)

      implicit   none
      include   'eldata.h'
      integer    hpltab(10,*)
      real*8     ev(10)

      if(nint(ev(2)).gt.0 .and. nint(ev(2)).le.10) then
        hpltab(nint(ev(2)),ma) = nint(ev(1))
      else
        write(*,*) ' PLOT HISTORY TABLE ERROR: GLOBAL =',nint(ev(2))
      endif

      end

      subroutine outhisplt(hpltab)

      implicit   none

      include   'iofile.h'
      include   'eldata.h'
      include   'eldatp.h'

      integer    i, hpltab(10,*)

      do i = 1,10
        if(hpltab(i,ma).gt.0) then
          write(iow,2002) i,hpltab(i,ma)
          hplmax = max(hplmax,i)
        endif
      end do ! i

c     Format

2002  format(10x,'History Plot : Contour =',i5,' for variable =',i5)

      end
