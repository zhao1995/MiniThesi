c$Id:$
      subroutine membr3d(d,ul,xl,ix,s,p,ndf,ndm,nst,isw)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Add implicit/explicit integration option         04/04/2007
c       2. Add 'ix' on call to elements                     15/09/2010
c       3. Set 'j=0' then add to 'nh3' for int8 use         17/04/2011
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose:  Quadrilateral membrane element for feap

c     Input parameters set as follows:

c       Small and finite deformation
c         ndm =  3 (x,y,z cartesian coordinates at nodes)
c         ndf =  3 (u-x,u-y,u-z at nodes)
c         nen = >3 nodes
c-----[--.----+----.----+----.-----------------------------------------]
      implicit   none

      include   'eldata.h'
      include   'eltran.h'
      include   'hdata.h'
      include   'iofile.h'
      include   'mdata.h'
      include   'strnum.h'
      include   'comblk.h'

      integer    ndm,ndf,nst,isw, tdof, i,j, ix(*)
      real*8     d(*),xl(ndm,*),ul(ndf,*),s(nst,*),p(nst), ctan1

      save

c     Input material properties

      if(isw.eq.1) then

        if(ior.lt.0) write(*,2000)
        write(iow,2000)
        call inmate(d,tdof, 12 ,5)

c       Deactivate dof in element for dof > 3

        do i = 4,ndf
          ix(i) = 0
        end do ! i

c       Set plot sequence

        pstyp = 2

c       Initialize the finite deformation shell

        if(d(18).lt.0.0d0) then
          call mem3df(d,ul,xl,s,p,ndf,ndm,nst,isw)
        endif

c       Set maximum plot variable number

        istv = 17

c     History variable manipulation: None currently required

      elseif(isw.eq.12) then

        return

c     Initialize element strains for activation

      elseif(isw.eq.17) then

        j = 0
        do i = 1,nel
          hr(nh3+j  ) = ul(1,i)
          hr(nh3+j+1) = ul(2,i)
          hr(nh3+j+2) = ul(3,i)
          j       = j + 3
        end do ! i

c     Initialize element strains for deactivation

      elseif(isw.eq.18) then

        do i = 1,3*nel
          hr(nh3+i-1) = 0.0d0
        end do ! i

c     Critical time step computation

      elseif(isw.eq.21) then

        call tcritnd(d,xl,ul,ndm,ndf,nel)

c     External node check

      elseif(isw.eq.26) then

        call pcorner2d()

c     Remaining options

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
          call mem3ds(d,ul,xl,ix,s,p,ndf,ndm,nst,isw)

c       Finite deformation

        else
          call mem3df(d,ul,xl,s,p,ndf,ndm,nst,isw)
        endif
        if(isw.eq.3) then
          ctan(1) = ctan1
        endif

      endif

c     Format

2000  format(5x,'T h r e e   D i m e n s i o n a l   M e m b r a n e',
     &          '   E l e m e n t')

      end
