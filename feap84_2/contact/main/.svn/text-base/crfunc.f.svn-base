c$Id:$
      subroutine crfunc (nsopt,td,emax)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Save transformation data in td(:) array which    15/05/2013
c          is moved to cs0(:,:) array on return
c-----[--.----+----.----+----.-----------------------------------------]
c      Coded/modified by:                           Date:          rel.:
c               Robert L. Taylor        10 Novermber 2001            1.0

c      Acronym: Contact Read FUNCtions

c      Purpose: Input of contact surface descriptions as functions

c      Inputs:
c         nsopt   - Function number
c         td(14)  - Parameters for function

c      Outputs:
c         emax    - Element max number found
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'iofile.h'
      include  'print.h'
      include  'trdata.h'

      integer   nsopt,emax, i,j, ii, imax,imin
      real*8    td(14), vr(3,3),vs(3), vlen,vmin,vmax

      save

c     Check sub-command option and subcommand data not needed

      call cdebug0 ('      crfunc',-1)

c     Set defaults

      emax  = 0

c     Cylinder
      if    (nsopt.eq.1) then
        write(iow,2001) nint(td(1)),td(2),td(3),nint(td(4))

c     Sphere
      elseif(nsopt.eq.2) then
        write(iow,2002) nint(td(1)),td(2),td(3),nint(td(4))

c     Cartesian plane
      elseif(nsopt.eq.3) then
        write(iow,2003) nint(td(1)),td(2),td(3),nint(td(4))

c     Normal plane
      elseif(nsopt.eq.4) then
        write(iow,2004) nint(td(1)),td(2),td(3),td(4),td(5),
     &                  td(6),nint(td(7))

c     Polynomial surface
      elseif(nsopt.eq.5) then
        do i = 7,14
          if(td(i).ne.0.0d0) ii = i
        end do ! i
        write(iow,2005) (td(i),i=1,3),(nint(td(i)),i=4,6),
     &                  (i-6,td(i),i=7,ii)
      endif

      if(nsopt.le.2) then
        if(nint(td(5)).eq.0) then
          write(iow,2010) tr,xr
          vr = tr
          vs = xr
        else
          write(iow,2011) (td(i),i=6,11)

c         Specified point case compute orientation

          vr(1,3) = td( 9) - td(6)
          vr(2,3) = td(10) - td(7)
          vr(3,3) = td(11) - td(8)
          vs(1)   = td( 5)
          vs(2)   = td( 7)
          vs(3)   = td( 8)

          vlen = vr(1,3)*vr(1,3) + vr(2,3)*vr(2,3) + vr(3,3)*vr(3,3)
          vlen  = 1.d0/sqrt(vlen)
          vr(:,3) = vr(:,3)*vlen
          vmax = vr(1,3)
          vmin = vr(1,3)
          do i = 2,3
            if(vr(i,3).gt.vmax) then
              vmax = vr(i,3)
              imax = i
            endif
            if(vr(i,3).lt.vmin) then
              vmin = vr(i,3)
              imin = i
            endif
          end do ! i
          vr(:,1) = 0.0d0
          vr(imin,1) =  vmax
          vr(imax,1) = -vmin
          vlen = 1.d0/sqrt(vr(imin,1)**2 + vr(imax,1)**2)
          vr(:,1) = vr(:,1)*vlen
          vr(1,2) = vr(2,3)*vr(3,1) - vr(3,3)*vr(2,1)
          vr(2,2) = vr(3,3)*vr(1,1) - vr(1,3)*vr(3,1)
          vr(3,2) = vr(1,3)*vr(2,1) - vr(2,3)*vr(1,1)
        endif

c       Save transformation data for surface

        ii = 4
        do j = 1,3
          do i = 1,3
            ii = ii + 1
            td(ii) = vr(i,j)
          end do ! i
        end do ! j
        do i = 1,3
          ii = ii + 1
          td(ii) = vs(i)
        end do ! i
      endif

c     Formats

2001  format(/7x,'Cylindrical Surface: Radial Expansion'/
     &       10x,'Motion Direction            =',i8/
     &       10x,'Initial Radius              =',1p,1e12.4/
     &       10x,'Displacement value (0: fix) =',1p,1e12.4/
     &       10x,'Expansion proportional load =',i8)

2002  format(/7x,'Spherical Surface: Radial Expansion'/
     &       10x,'Motion Direction            =',i8/
     &       10x,'Initial Radius              =',1p,1e12.4/
     &       10x,'Displacement value (0: fix) =',1p,1e12.4/
     &       10x,'Expansion proportional load =',i8)

2003  format(/7x,'Cartesian Surface: Radial Expansion'/
     &       10x,'Coordinate Direction        =',i8/
     &       10x,'Initial Coordinate          =',1p,1e12.4/
     &       10x,'Displacement value (0: fix) =',1p,1e12.4/
     &       10x,'Expansion proportional load =',i8/)

2004  format(/7x,'Plane Surface: Normal Motion'/
     &       10x,'Motion Direction            =',i8/
     &       10x,'Normal Vector (n_1)         =',1p,1e12.4/
     &       10x,'Normal Vector (n_2)         =',1p,1e12.4/
     &       10x,'Normal Vector (n_3)         =',1p,1e12.4/
     &       10x,'Initial Coordinate          =',1p,1e12.4/
     &       10x,'Displacement value (0: fix) =',1p,1e12.4/
     &       10x,'Expansion proportional load =',i8/)

2005  format(/7x,'Polynomial Surface:'/
     &       10x,'Displacement Vector (u_1)   =',1p,1e12.4/
     &       10x,'Displacement Vector (u_2)   =',1p,1e12.4/
     &       10x,'Displacement Vector (u_3)   =',1p,1e12.4/
     &       10x,'Proportional load   (p_1)   =',i8/
     &       10x,'Proportional load   (p_2)   =',i8/
     &       10x,'Proportional load   (p_3)   =',i8/
     &      (10x,'Polynomial term ',i3,'       =',1p,1e12.4))

2010  format(10x,'Transformation data'/
     &       12x,'Orientation Matrix',3(/15x,1p,3e12.4)/
     &       12x,'Translation'/15x,1p,3e12.4/)

2011  format(10x,'Point 1: Coordinates        =',1p,3e12.4/
     &       10x,'Point 2: Coordinates        =',1p,3e12.4/)

      end
