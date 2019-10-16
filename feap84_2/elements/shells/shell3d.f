c$Id:$
      subroutine shell3d(d,ul,xl,ix,s,p,ndf,ndm,nst,isw)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Add implicit/explicit integration option         04/04/2007
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose:  Quadrilateral shell element for feap

c     Input parameters set as follows:

c       Small deformation
c         ndm = 3 (x,y,z cartesian coordinates at nodes)
c         ndf = 6 (u-x,u-y,u-z,r-x,r-y,r-z at nodes)
c         nen = 4 nodes (counterclockwise around element)

c       Note: 1-direction bisects diagonals between 2-3 element and
c             2-direction bisects diagonals between 3-4 element nodes.

c       Finite deformation
c         ndm = 3 (x,y,z cartesian coordinates at nodes)
c         ndf = 5 (u-x,u-y,u-z,r-a,r-b at nodes)
c         nen = 4 nodes (counterclockwise around element)
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      implicit  none

      include  'cdata.h'
      include  'eldata.h'
      include  'eltran.h'
      include  'iofile.h'
      include  'hdata.h'
      include  'mdata.h'
      include  'strnum.h'

      include  'comblk.h'

      integer   ndm,ndf,nst,isw, tdof, i
      integer   ix(*)

      real*8    d(*),xl(ndm,*),ul(ndf,*),s(nst,*),p(nst), ctan1

      save

c     Input material properties

      if(isw.eq.1) then

        if(ior.lt.0) write(*,2000)
        write(iow,2000)
        call inmate(d,tdof,0,5)

c       Check control record dimensions

        if(d(18).gt.0.0d0) then
          i = 6
        else
          i = 5
        endif

        if(ndf.lt.i .or. nen.lt.3) then
          write(  *,2001) i,ndf,nen
          write(iow,2001) i,ndf,nen
          write(ilg,2001) i,ndf,nen
          call plstop()
        end if

c       Deactivate dof in element for dof > 6

        do i = 7,ndf
          ix(i) = 0
        end do ! i

c       History for through thickness integrations

        if(nint(d(102)).gt.1) then
          nh1 = nh1*nint(d(102))
        endif

c       Construct rotation parameters: u-x=1; u-y=2 (same as defaults)

        ea(1,-iel) = 1
        ea(2,-iel) = 2

c       Construct rotation parameters: theta-x = 4; theta-y = 5

        er(1,-iel) = 4
        er(2,-iel) = 5

c       Set plot sequence

        pstyp = 2

c       Set maximum number of stress plots

        istv = 24

c       Initialize the finite deformation shell

        if(d(18).lt.0.0d0) then
          call shl3df(d,ul,xl,s,p,ndf,ndm,nst,isw)
        endif

c     History variable manipulation: None currently required

      elseif(isw.eq.12) then

        return

c     Critical time step computation

      elseif(isw.eq.21) then

        call tcritnd(d,xl,ul,ndm,ndf,nel)

c     External node check

      elseif(isw.eq.26) then

        call pcorner2d()

C     Remaining options

      else

c       Explicit/Implicit element solutions

        if(isw.eq.3) then
          ctan1 = ctan(1)
          if(d(187).gt.0.0d0 .and.
     &       min(ctan(2),ctan(3)).gt.0.0d0) then
            ctan(1) = 0.0d0
          endif
        endif

c       Small deformation

        if(d(18).gt.0.0d0) then
          if(nint(d(102)).gt.1) then
            call shl3di(d,ul,xl,ix,s,p,ndf,ndm,nst,isw)
          else
            call shl3ds(d,ul,xl,ix,s,p,ndf,ndm,nst,isw)
          endif

c       Finite deformation

        else
          call shl3df(d,ul,xl,s,p,ndf,ndm,nst,isw)
        endif
        if(isw.eq.3) then
          ctan1 = ctan(1)
        endif

      endif

c     Format

2000  format(5x,'T h r e e   D i m e n s i o n a l   S h e l l',
     &          '   E l e m e n t')

2001  format(/5x,'*ERROR* Some of following control record values',
     &       ' are incorrect:'
     &       /13x,'ndf (should be > or = ',i1,'): Specified =',i3
     &       /13x,'nen (should be > or = 3): Specified =',i3)

      end
