c$Id:$
      subroutine frame2d(d,ul,xl,ix,s,r,ndf,ndm,nst,isw)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Add implicit/explicit integration option         04/04/2007
c-----[--.----+----.----+----.-----------------------------------------]
c     Two dimensional frame element

c     Control Information:

c       ndm  - Spatial dimension of problem       = 2
c       ndf  - Number degree-of-freedoms at node >= 3
c              ( 1 = u_1 ; 2 = u_2 ; 3 = theta )
c       nen  - Two node element (nel = 2)        >= 2
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'bdata.h'
      include  'cdata.h'
      include  'eldata.h'
      include  'eltran.h'
      include  'iofile.h'
      include  'comblk.h'

      integer   ndf,ndm,nst,isw,tdof,i,ix(*)
      real*8    cs,sn,le,d(*),ul(ndf,nen,*),xl(ndm,*),s(nst,*),r(*)
      real*8    ctan1

      save

c     Set Element type

c     CHECK ELEMENT FOR ERRORS

      if(isw.eq.2) then

        cs = xl(1,2) - xl(1,1)
        sn = xl(2,2) - xl(2,1)
        le = sqrt(cs*cs+sn*sn)

        if(ix(1).le.0 .or. ix(2).le.0 .or. ix(1).eq.ix(2)) then
          write(iow,4000) n,ix(1),ix(2)
          if(ior.lt.0) write(*,4000) n,ix(1),ix(2)
        endif
        if(le.le.0.0d0) then
          write(iow,4001) n
          if(ior.lt.0) write(*,4001) n
        endif

c     Critical time step computation

      elseif(isw.eq.21) then

        call tcritnd(d,xl,ul,ndm,ndf,nel)

c     Body force computation

      elseif(isw.eq.15 .or. isw.eq.23) then

        call fbody2d(d,xl,ul, r,s, ndm,ndf,nst, isw)

c     External nodes

      elseif(isw.eq.26) then

c     Compute element arrays

      else

c       INPUT MATERIAL PARAMETERS

        if(isw.eq.1) then
          if(ior.lt.0) write(*,2000)
          write(iow,2000)
          call inmate(d,tdof,0,3)

c         Set plot sequence

          pstyp = 1

c         Check control data prameters

          if(ndf.lt.3 .or. nen.lt.2) then
            write(ilg,3000) ndf,nen
            write(iow,3000) ndf,nen
            call plstop()
          endif

c         Deactivate dof in element for dof > 3

          do i = 4,ndf
            ix(i) = 0
          end do ! i

        endif

c       Explicit/Implicit element solutions

        if(isw.eq.3) then
          ctan1 = ctan(1)
          if(d(187).gt.0.0d0 .and.
     &       min(ctan(2),ctan(3)).gt.0.0d0) then
            ctan(1) = 0.0d0
          endif
        endif

c       Small deformaion

        if(d(18).gt.0.0d0 ) then

c         Shear deformable

          if(d(79).eq.0.0d0) then

c           2-node: cubic/quadratic enhanced interpolations)

            if(nint(d(17)).eq.3) then

              call frams2e(d,ul,xl,s,r,ndf,ndm,nst,isw)

c           2-node: linear interpolations

            else

              call frams2d(d,ul,xl,s,r,ndf,ndm,nst,isw)

            endif

c         Euler-Bernoulli  (2-node: cubic interpolations)

          else

            call franf2d(d,ul,xl,s,r,ndf,ndm,nst,isw)

          endif

c       Finite deformation (2-node: linear interpolations)

        else

c         Shear deformable (2-node: linear, finite displacements)

          if(d(79).eq.0.0d0) then

            call framf2d(d,ul,xl,s,r,ndf,ndm,nst,isw)

c         No shear case    (2-node: cubic, 2-nd order displacements)

          else

            call franf2d(d,ul,xl,s,r,ndf,ndm,nst,isw)

          endif
        endif
        if(isw.eq.3) then
          ctan(1) = ctan1
        endif
      endif

c     Formats

2000  format(5x,'T w o    D i m e n s i o n a l    F r a m e'/)

3000  format(/' *ERROR* Control parameters set incorrectly: Must be:'/
     &        '         ndf = ',i3,': (Should be 3 or more)'/
     &        '         nen = ',i3,': (Should be 2 or more)'/)

4000  format(' *ERROR* Element',i7,' has unspecified node: ix =',2i7)
4001  format(' *ERROR* Element',i7,' has zero length')

      end
