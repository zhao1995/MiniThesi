c$Id:$
      subroutine point(d,ul,ix,s,r,ndf,ndm,nst,isw)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Set pstyp = -1 for no plots                      03/01/2007
c       2. Set pstyp =  0 for no plots                      31/08/2008
c       3. Add 'jj' to 'nh3' for isw.eq.17                  17/04/2011
c-----[--.----+----.----+----.-----------------------------------------]
c     Purpose:  Point Stiffness, Mass, Damping ELEMENT

c     Input data: (isw = 1)

c       Record ('mass',m)
c         m - Mass / dof

c       Record ('damp',d)
c         d - Damper / dof

c       Record ('spri',k)
c         k - Spring / dof

c       Record ('orie',k(i),i=1,2*ndm)
c         X_1,X_2 - Coordinates

c     Outputs: (isw = 4)
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'augdat.h'
      include  'bdata.h'
      include  'cdata.h'
      include  'eldata.h'
      include  'eltran.h'
      include  'hdata.h'
      include  'iofile.h'
      include  'part0.h'
      include  'rdata.h'
      include  'tdata.h'
      include  'comblk.h'

      logical   errck, tinput, pcomp, spring, masses, damper
      character type*15
      integer   ndf, ndm, nst, isw, i,j,jj, ix(*)
      real*8    tang(3,3),damp(3,3),mass(3,3),dirc(3),leng,ctan2,ctan3
      real*8    du(3),dv(3),da(3),td(6)
      real*8    d(*),ul(ndf,nen,*),s(nst,*),r(ndf,*)

      save

c     INPUT MATERIAL PROPERTIES

      if(isw.eq.1) then

c       Record:

        spring = .false.
        damper = .false.
        masses = .false.
        do i = 11,15
          d(i) = 0.0d0
        end do ! i
        d(11) = -1.d0

c       Disable dof's

        do i = ndm+1,ndf
          ix(i) = 0
        end do ! i

        type = 'start'
        do while(.not.pcomp(type,'    ',4))

          if(ior.lt.0) then
            write(*,3000)
            call pprint('  >')
          endif
          errck = tinput(type,1,td,6)

c         Mass parameter

          if    (pcomp(type,'mass',4)) then

            d(1) = td(1)
            masses = .true.

c         Damper parameter

          elseif(pcomp(type,'damp',4)) then

            d(2) = td(1)
            damper = .true.

c         Spring parameter

          elseif(pcomp(type,'spri',4)) then

            d(3) = td(1)
            spring = .true.

c         Orientation parameter

          elseif(pcomp(type,'orie',4)) then

            do i = 1,ndm
               d(10+i) = td(i)
            end do ! i

c         Inadmissible data

          elseif(.not.pcomp(type,'    ',4)) then
            write(  *,4000) type
            write(ilg,4000) type
            write(iow,4000) type

          end if

        end do ! while

c       Output material properties

        write(iow,2000)
        if(masses) then
          write(iow,2001) d(1)
        endif
        if(damper) then
          write(iow,2002) d(2)
        endif
        if(spring) then
          write(iow,2003) d(3)
        endif
        write(iow,2004) (d(i),i=11,10+ndm)

c       Set history terms

        nh1 = 1
        nh3 = 2*ndm

c       Set for no plots

        pstyp = 0

c     CHECK ELEMENTS

      elseif(isw.eq.2) then

c     COMPUTE ELEMENT STIFFNESS AND RESIDUAL ARRAYS

      elseif(isw.eq.3  .or. isw.eq.6 ) then

        leng = 0.0d0
        do i = 1,ndm
          dirc(i) =  d(10+i)
          leng    =  leng + dirc(i)**2
          du(i)   =  ul(i,2,1) - ul(i,1,1)
          dv(i)   =  ul(i,2,4) - ul(i,1,4)
          da(i)   =  ul(i,2,5) - ul(i,1,5)
        end do ! i
        if(leng.gt.0.0d0) then
          leng = 1.d0/sqrt(leng)
        else
          leng = 1.d0
        endif
        do i = 1,ndm
          dirc(i) = dirc(i)*leng
        end do ! i
        do i = 1,ndm
          do j = 1,ndm
            mass(j,i) =         d(1)*dirc(i)
            damp(j,i) = dirc(j)*d(2)*dirc(i)
            tang(j,i) = dirc(j)*d(3)*dirc(i)
          end do ! j
        end do ! i

c       Compute residual and tangent

        if(ndfo(1).gt.0 .or. shflg) then
          ctan2 = ctan(2)
          ctan3 = ctan(3)
        else
          ctan2 = 0.0d0
          ctan3 = 0.0d0
        endif
        do i = 1,ndm
          do j = 1,ndm
            r(i,1) =  r(i,1) + tang(i,j)*du(j)
     &                       + damp(i,j)*dv(j)
     &                       + mass(i,j)*da(j)
            s(i    ,j    ) =   tang(i,j)*ctan(1)
     &                       + damp(i,j)*ctan2
     &                       + mass(i,j)*ctan3
            s(i+ndf,j    ) = - s(i,j)
            s(i    ,j+ndf) =   s(i,j)
            s(i+ndf,j+ndf) =   s(i,j)
          end do ! j
          r(i,2) = -r(i,1)
        end do ! i

c     Initialize element strains for activation

      elseif(isw.eq.17) then

        jj = 0
        do i = 1,2
          do j = 1,ndm
            hr(nh3+jj) = ul(j,i,1)
            jj     = jj + 1
          end do ! j
        end do ! i

c     Initialize element strains for deactivation

      elseif(isw.eq.18) then

        do i = 0,2*ndm-1
          hr(nh3+i) = 0.0d0
        end do ! i

c     External node check

      elseif(isw.eq.26) then

      endif

c     I/O Formats

2000  format(9x,'P o i n t    E l e m e n t'//)
2001  format(14x,'MASS     :',1p,1e12.4)
2002  format(14x,'DAMPER   :',1p,1e12.4)
2003  format(14x,'SPRING   :',1p,1e12.4)
2004  format(14x,'DIRECTION:',1p,3e12.4)

3000  format(' Input: TYPE, VALUES')

4000  format(' *WARNING* Property: ',a,' is not valid')

      end
