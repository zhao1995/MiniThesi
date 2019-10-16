c$Id:$
      subroutine solid3dtm(d,ul,xl,ix,tl,s,p,ndf,ndm,nst,isw)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    06/02/2010
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Three Dimensional Coupled Thermo-Mechanical Solid
c               Element Driver

c      Inputs:
c         d(*)  - Element parameters
c         ul(ndf,*) - Current nodal solution parameters
c         xl(ndm,*) - Nodal coordinates
c         ix(*)     - Global nodal connections
c         ndf       - Degree of freedoms/node
c         ndm       - Mesh coordinate dimension
c         nst       - Element array dimension
c         isw       - Solution option switch

c      Outputs:
c         s(nst,*)  - Element array
c         p(ndf,*)  - Element vector
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

      logical   errck
      integer   ndf,ndm,nst,isw, i,tdof, ix(*)
      real*8    d(*),ul(ndf,*),xl(ndm,*),tl(*),s(nst,nst),p(nst)
      real*8    shp(4,64), ctan1

      save

c     Extract type data

      stype = nint(d(16))
      etype = nint(d(17))
      dtype = nint(d(18))
      hflag = d(67).eq.1.d0

c     Set zero state

c     Output element type

      if(isw.eq.0 .and. ior.lt.0) then
        write(*,*) '   TM_solid: 3-d Solid Thermo-Mechanical Element.'

c     Input material properties

      elseif(isw.eq.1) then

        write(iow,2000)
        if(ior.lt.0) write(*,2000)
        call inmate(d,tdof,   0   ,1)

c       Set to check in each element

        stype = nint(d(16))
        etype = nint(d(17))
        dtype = nint(d(18))

c       Deactivate dof in element for dof > 4: DOFS = u1,u2,u3,T

        do i = 5,ndf
          ix(i) = 0
        end do ! i

c       Set plot sequence

        pstyp = 3

c       Check for errors in problem size

        errck = .false.
        if(ndf.lt.4 .or. nen.lt.4) then
          write(ilg,4000) ndf,nen
          write(iow,4000) ndf,nen
          errck = .true.
        endif

        if(etype.eq.5) then
          write(ilg,4001)
          write(iow,4001)
          errck = .true.
        endif

        if(errck) then
          call plstop()
        endif

c     Check element for errors

      elseif(isw.eq.2) then

        if(nel.eq.4 .or. nel.eq.10) then
          call cktets ( n, ix, xl, ndm, nel, shp )
        else
          call ckbrk8 ( n, ix, xl, ndm, nel, shp )
        endif
        return

c     Compute mass or geometric stiffness matrix

      elseif(isw.eq.5) then

c       Mass contribution

        if(imtyp.eq.1) then

          call mass3d(d,xl,s,p,ndf,ndm,nst)

c       Geometric stiffness here

        else
c         Put call to geometric stiffness routine here
        endif
        return

c     Compute damping matrix

      elseif(isw.eq.9) then

c       Damping contribution

        call damp3d(d,xl,s,ndf,ndm,nst)
        return

c     History manipulation: None currently required

c     elseif(isw.eq.12) then

c       return

c     Critical time step computation

      elseif(isw.eq.21) then

        call tcritnd(d,xl,ul,ndm,ndf,nel)

c     Body Force computation

      elseif(isw.eq.15 .or.isw.eq.23) then

        call sbody3d(d,xl, p,ndm,ndf, isw)
        return

c     Compute face normals

      elseif(isw.eq.24) then

        call pnorm3d(p,xl,ndm,nel)
        return

c     Compute boundary nodes

      elseif(isw.eq.26) then

        call pcorner3d()
        return

      endif

c     Compute stress-divergence vector (p) and stiffness matrix (s)

c     Explicit/Implicit element solutions

      if(isw.eq.3) then
        ctan1 = ctan(1)
        if(d(187).gt.0.0d0 .and.
     &     min(ctan(2),ctan(3)).gt.0.0d0) then
          ctan(1) = 0.0d0
        endif
      endif

c     Displacement Model

      if(etype.eq.1) then

        if(dtype.gt.0) then
          call sld3d1tm(d,ul,xl,tl,s,p,ndf,ndm,nst,isw)
        else
          call fld3d1tm(d,ul,xl,s,p,ndf,ndm,nst,isw)
        endif

c     Mixed Model (B-Bar)

      elseif(etype.eq.2) then

        if(dtype.gt.0) then
          call sld3d2tm(d,ul,xl,s,p,ndf,ndm,nst,isw)
        else
          call fld3d2tm(d,ul,xl,s,p,ndf,ndm,nst,isw)
        endif

c     Enhanced Strain Model

      elseif(etype.eq.3) then

        if(dtype.gt.0) then
          write(iow,4003)
          write(ilg,4003)
          call plstop()
        else
          write(iow,4003)
          write(ilg,4003)
          call plstop()
        endif

c     Energy Conserving Model

      elseif(etype.eq.4) then

        if(dtype.gt.0) then
          call sld3d1tm(d,ul,xl,tl,s,p,ndf,ndm,nst,isw)
        else
          write(iow,4005)
          write(ilg,4005)
          call plstop()
        endif

c     Mixed-Enhanced Model

      elseif(etype.eq.5) then

        write(iow,4003)
        write(ilg,4003)
        call plstop()

c     Co-rotational formulation

      elseif(etype.eq.6) then

        write(iow,4004)
        write(ilg,4004)
        call plstop()

c     Displacement Uniform Deformation gradient Model (B-Bar)

      elseif(etype.eq.7) then

        if(dtype.gt.0) then
          write(iow,4002)
          write(ilg,4002)
          call plstop()
        else
c         call fld3d1utm(d,ul,xl,s,p,ndf,ndm,nst,isw)
        endif

c     Mixed Uniform Deformation gradient Model (B-Bar)

      elseif(etype.eq.8) then

        if(dtype.gt.0) then
          write(iow,4002)
          write(ilg,4002)
          call plstop()
        else
c         call fld3d2utm(d,ul,xl,s,p,ndf,ndm,nst,isw)
        endif

      endif
      if(isw.eq.3) then
        ctan(1) = ctan1
      endif

c     Formats for input-output

2000  format(
     & /3x,'T h r e e   D i m e n s i o n a l   C o u p l e d',
     &     '   S o l i d   E l e m e n t'/)

4000  format(' *ERROR* Problem control record incorrect:'/
     &       '         DOFs/Node (ndf) = ',i4,' - Should be 4 or more'/
     &       '         Nodes/Elm (nen) = ',i4,' - Should be 4 or more')

4001  format(' *ERROR* Mixed/Enhanced element type does not currently',
     &       ' exist.')

4002  format(' *ERROR* Small deformation, uniform deformation',
     &       ' gradient element type'/9x,' does not currently exist.')

4003  format(' *ERROR* No Mixed/Enhanced 3D coupled element')

4004  format(' *ERROR* No Co-rotational 3D coupled element')

4005  format(' *ERROR* No Conserving 3D coupled element')

      end
