c$Id:$
      subroutine shell2d(d,ul,xl,ix,s,r,ndf,ndm,nst,isw)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Add implicit/explicit integration option         04/04/2007
c-----[--.----+----.----+----.-----------------------------------------]
c     Shell elements: 2-D Axisymmetric
c         (a) Finite motions with trigonometric functions
c         (a) Large displacement  with small rotations

c     Three integration points are used.

c     Arguments:
c        d(*)      - specified parameter array
c        ul(ndf,*) - local nodal solution values
c        xl(ndm,*) - local nodal coordinate values
c        ix(*)     - node numbers
c        s(nst,nst) - finite element array (stiffness, mass, geometric
c                                           stiffness)
c        r(nst)     - finite element array (residual, lumped mass)
c        ndf        - number of degree of freedoms at node ( > or = 3 )
c        ndm        - spatial dimension of element         (      = 2 )
c        nst        - size of finite element arrays        ( > or = 6 )
c        isw        - solution option
c                   = 1: Input values and store in d(*) array
c                   = 2: Check mesh coordinate and connection inputs
c                        for errors
c                   = 3: Compute element residual (p) and stiffness (s)
c                   = 4: Output element results
c                   = 5: Compute mass (p,s)/geometric stiffness array(s)
c                   = 6: Compute element residual (p)
c                   = 8: Compute nodal projections
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'cdata.h'
      include  'eldata.h'
      include  'eltran.h'
      include  'iofile.h'
      include  'mdata.h'
      include  'pmod2d.h'
      include  'strnum.h'
      include  'comblk.h'

      integer   ndf,ndm,nst,isw, i, tdof, ix(*)
      real*8    d(*),ul(ndf,nen,*),xl(ndm,*),s(nst,*),r(ndf,*), ctan1

      save

c     Go to correct array processor

      if(isw.eq.0 .and. ior.lt.0) then
        write(*,*) '   Shell: 2-d Axisymmetric Model'

c     Input material properties

      elseif(isw.eq.1) then
        write(iow,2000)
        if(ior.lt.0) write(*,2000)
        call inmate(d,tdof,0,4)

c       Set plot sequence for 2-node line

        pstyp = 1

c       Check control record dimensions

        if(ndf.lt.3 .or. nen.lt.2) then
          write(  *,2001) ndf,nen
          write(iow,2001) ndf,nen
          write(ilg,2001) ndf,nen
          call plstop()
        end if

c       Deactivate dof in element for dof > 3

        do i = 4,ndf
          ix(i) = 0
        end do ! i

c       Number of projections

        istv = 11

c     History variable manipulation: None currently required

      elseif(isw.eq.12) then

        return

c     Critical time step computation

      elseif(isw.eq.21) then

        call tcritnd(d,xl,ul,ndm,ndf,nel)

c     External node check

      elseif(isw.eq.26) then

      else

c       Explicit/Implicit element solutions

        if(isw.eq.3) then
          ctan1 = ctan(1)
          if(d(187).gt.0.0d0 .and.
     &       min(ctan(2),ctan(3)).gt.0.0d0) then
            ctan(1) = 0.0d0
          endif
        endif

        dtype = nint(d(18))

        if(nel.eq.2) then

c         Small displacement formulation

          if(dtype.gt.0) then
            call shl2ds(d,ul,xl,s,r,ndf,ndm,nst,isw)

c         Finite displacement formulation

          elseif(dtype.lt.0) then
            call shl2df(d,ul,xl,s,r,ndf,ndm,nst,isw)

c         Second order formulation

          else
            call shl2dn(d,ul,xl,s,r,ndf,ndm,nst,isw)
          endif
        else
          write(iow,3000) n
          if(ior.lt.0) write(*,3000) n
          call plstop()
        endif
        if(isw.eq.3) then
          ctan(1) = ctan1
        endif
      endif

c     Formats for input-output

2000  format(/5x,'E l a s t i c   S h e l l   E l e m e n t'/)

2001  format(/5x,'*ERROR* Some of following control record values',
     &       ' are incorrect:'
     &       /13x,'ndf (should be > or = 3): Specified =',i3
     &       /13x,'nen (should be > or = 2): Specified =',i3)

3000  format('*ERROR* Shell Element',i8,' has more than 2-nodes')

      end
