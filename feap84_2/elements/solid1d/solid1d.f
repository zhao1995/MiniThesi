c$Id:$
      subroutine solid1d(d,ul,xl,ix,tl,s,p,ndf,ndm,nst,isw)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Add implicit/explicit integration option         04/04/2007
c       2. Add mther and pdof to control computation of     29/03/2011
c          thermal parts
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Plane and axisymmetric linear thermo-elastic element
c               driver

c      Inputs:
c         d(*)  - Element parameters
c         ul(ndf,*) - Current nodal solution parameters
c         xl(ndm,*) - Nodal coordinates
c         ix(*)     - Global nodal connections
c         tl(*)     - Nodal temp vector
c         ndf       - Degree of freedoms/node
c         ndm       - Mesh coordinate dimension
c         nst       - Element array dimension
c         isw       - Solution option switch

c      Outputs:
c         s(nst,*)  - Element array
c         p(*)      - Element vector
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'cdata.h'
      include  'eldata.h'
      include  'eltran.h'
      include  'eqsym.h'
      include  'evdata.h'
      include  'hdata.h'
      include  'iofile.h'
      include  'part0.h'
      include  'part1.h'
      include  'pmod2d.h'
      include  'strnum.h'
      include  'comblk.h'

      logical   mech,ther,mther
      integer   ndf,ndm,nst,isw, i,tdof,pdof, ix(*)
      real*8    d(*),ul(ndf,4),xl(ndm,*),tl(*),s(nst,nst),p(nst),th(4)
      real*8    ctan1

      save

c     Extract type data

      stype = nint(d(16))
      etype = nint(d(17))
      dtype = nint(d(18))
      hflag = d(67).eq.1.0d0

c     Set nodal temperatures: Can be specified or computed

      if(isw.gt.1) then
        tdof = nint(d(19))
        if(tdof.le.0) then
          do i = 1,nel ! {
            th(i) = tl(i)
          end do ! i     }
        else
          do i = 1,nel ! {
            th(i) = ul(tdof,i)
          end do ! i     }
        endif
      endif

c     Input material properties

      if(isw.eq.1) then

        write(iow,2001)
        if(ior.lt.0) write(*,2001)
        call inmate(d,tdof,0,1)
        if(etype.eq.2) then
          nh1 = nh1 + 2
        endif

c       Set tdof to zero if 1 or larger than ndf

        if(tdof.gt.ndf) then
          write(iow,3003)
          if(ior.lt.0) write(*,3003)
          tdof = 0
        elseif(tdof.eq.1) then
          write(iow,3004)
          if(ior.lt.0) write(*,3004)
          tdof = 0
        endif

c       Deactivate dof in element for dof > 2

        if(nint(d(16)).eq.8) then
          do i = 3,ndf
            ix(i) = 0
          end do ! i
        else
          do i = 2,ndf
            ix(i) = 0
          end do ! i
        endif

c       If temperature dof is specified activate dof

        if(tdof.gt.0) then
          ix(tdof) = 1
        endif

c       Set plot sequence for triangles

        pstyp = 1

c       Set number of projected stresses

        istv = max(istv,6)

c     Check element for errors in input data

      elseif(isw.eq.2) then

c     Compute mass or geometric stiffness matrix

      elseif(isw.eq.5) then

        mech =  npart.eq.ndfp(1)
        if(tdof.gt.ndm .and. hflag) then
          ther = npart.eq.ndfp(tdof)
        else
          ther = .false.
        endif

        if(imtyp.eq.1) then
          if(mech) then
            call mass1d(d,xl,s,p,ndf,ndm,nst)
          endif
          if(ther) then
            call therm1d(d,ul(tdof,1),xl,ix,s(tdof,tdof),p(tdof),
     &                   ndf,ndm,nst,isw)
          endif
        else
c         put call to geometric stiffness routine here
        endif

c     Compute damping matrix

      elseif(isw.eq.9) then

        mech =  npart.eq.ndfp(1)
        if(tdof.gt.ndm .and. hflag) then
          ther = npart.eq.ndfp(tdof)
        else
          ther = .false.
        endif

        if(mech) then
          call damp1d(d,xl,s,ndf,ndm,nst)
        endif
        if(ther) then
          call therm1d(d,ul(tdof,1),xl,ix,s(tdof,tdof),p(tdof),
     &                 ndf,ndm,nst,isw)
        endif

c     History manipulation: None required

      elseif(isw.eq.12) then
        return

c     Critical time step computation

      elseif(isw.eq.21) then
        call tcritnd(d,xl,ul,ndm,ndf,nel)

c     Body Force computation

      elseif(isw.eq.15 .or. isw.eq.23) then
        call sbody1d(d,xl, p,ndm,ndf ,isw)

c     Normal to surface computation for projections

      elseif(isw.eq.24) then
c       call pnorm1d(p,xl,ndm,nel)

c     External node check

      elseif(isw.eq.26) then

c     Compute residuals and tangents for parts

      else
        mech =  npart.eq.ndfp(1) .or. isw.eq.4 .or. isw.eq.8
        if(tdof.gt.ndm .and. hflag) then
          ther = npart.eq.ndfp(tdof) .or. isw.eq.4 .or. isw.eq.8
        else
          ther = .false.
        endif

      endif

c     Check if symmetric or unsymmetric tangent

      if(neqs.eq.1) then
        mther = ther
      else
        mther = .false.
      endif


c     Compute stress-divergence vector (p) and stiffness matrix (s)

      if(mech) then

c     Explicit/Implicit element solutions

        if(isw.eq.3) then
          ctan1 = ctan(1)
          if(d(187).gt.0.0d0 .and.
     &       min(ctan(2),ctan(3)).gt.0.0d0) then
            ctan(1) = 0.0d0
          endif
        endif

c       Displacement Model

        if(etype.eq.1) then

          if(dtype.gt.0) then
            call sld1d1(d,ul,xl,th,s,p,ndf,ndm,nst,isw,mther)
          else
            call fld1d1(d,ul,xl,s,p,ndf,ndm,nst,isw)
          endif

c       Mixed Model (B-Bar)

        elseif(etype.eq.2) then

          if(dtype.gt.0) then
            call sld1d2(d,ul,xl,th,s,p,ndf,ndm,nst,isw,mther)
          else
            call fld1d2(d,ul,xl,s,p,ndf,ndm,nst,isw)
          endif

c       Enhanced Strain Model (B-Bar)

        elseif(etype.eq.3) then

          if(dtype.gt.0) then
            call sld1d3(d,ul,xl,th,s,p,ndf,ndm,nst,isw)
          else
c           call fld1d3(d,ul,xl,s,p,ndf,ndm,nst,isw)
          endif

c       Energy Conserving Model

        elseif(etype.eq.4) then

          if(dtype.gt.0) then
            call sld1d1(d,ul,xl,th,s,p,ndf,ndm,nst,isw,mther)
          else
            call fld1d4(d,ul,xl,s,p,ndf,ndm,nst,isw)
          endif

c       Mixed-Enhanced Strain Model

        elseif(etype.eq.5) then

          if(dtype.gt.0) then
c           call sld1d5(d,ul,xl,th,s,p,ndf,ndm,nst,isw)
          else
c           call fld1d5(d,ul,xl,s,p,ndf,ndm,nst,isw)
            write(ilg,4000)
            write(iow,4000)
            call plstop()
          endif

        endif
        if(isw.eq.3) then
          ctan(1) = ctan1
        endif
      endif

c     Compute thermal vector (p) and matrix (s)

      if(ther) then

c       Prevent accumulation of weight for projections

        if(isw.eq.8) then
          pdof = 1
          do i = 1,nen
            p(i) = 0.0d0
          end do ! i
        else
          pdof = tdof
        endif

        call therm1d(d,ul(tdof,1),xl,ix,s(pdof,pdof),p(pdof),
     &               ndf,ndm,nst,isw)
      endif

c     Formats for input-output

2001  format(
     & /5x,'O n e   D i m e n s i o n a l   S o l i d   E l e m e n t'/)

3003  format(' *WARNING* Thermal d.o.f. > active d.o.f.s : Set to 0')

3004  format(' *WARNING* Thermal d.o.f. can not be 1: Set to 0')

4000  format(' *ERROR* Mixed-Enhanced element type does not currently',
     &       ' exist.')

      end
