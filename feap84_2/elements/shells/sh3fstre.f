c$Id:$
      subroutine sh3fstre ( xl, shp, ce, cx, cr, sn, sq, sm, xjw,
     &                      lint, ndm, isw )

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Move set of istv to isw.eq.1                     31/05/2013
c       2. Correct print labels for xy components           01/06/2013
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose:  Stress output routine

c      Inputs:
c         xl(3,*)  - Nodal coordinates
c         shp(4,*) - Element shape functions
c         ce(3,*)  - Membrane strains
c         cx(2,*)  - Shear    strains
c         cr(3,*)  - Bending  strains
c         sn(3,*)  - Membrane stresses
c         sq(3,*)  - Shear    stresses
c         sm(3,*)  - Bending  stresses
c         xjw(*)   - Jacobian weight
c         lint     - Number of output points
c         ndm      - Mesh spatial dimension
c         isw      - Output switch: 4 = print; 8 = project.

c      Outputs:
c         Prints or projected values (in hr array)
c-----[--.----+----.----+----.-----------------------------------------]

      implicit  none

      include  'cdata.h'
      include  'eldata.h'
      include  'iofile.h'
      include  'pointer.h'
      include  'prstrs.h'
      include  'strnum.h'
      include  'comblk.h'

      integer   i, j, lint, ndm, isw

      real*8    xl(ndm,*), shp(4,4), zphg(3,4)
      real*8    ce(3,4)  , cx(2,4) , cr(3,4)
      real*8    sn(3,4)  , sq(2,4) , sm(3,4) , xjw(4)

      save

c     Output Element Stresses:

      if(isw.eq.4) then
        mct = mct - 1
        if (mct.le.0) then
          if (ior.lt.0) write (*,4000)
          write (iow,4000)
          mct = 12
        endif

c       Gauss Loop

        do i = 1, lint

c         Interpolate Position

          do j = 1, 3
            zphg(j,i) = shp(1,i)*xl(j,1) + shp(2,i)*xl(j,2)
     &                + shp(3,i)*xl(j,3) + shp(4,i)*xl(j,4)
          end do ! i

c         Output to Screen

          if (ior.lt.0) then
            write (*,4100)  n, i, (zphg(j,i),j=1,3),
     &                     (ce(j,i),j=1,3), (sn(j,i),j=1,3),
     &                     (cr(j,i),j=1,3), (sm(j,i),j=1,3),
     &                     (cx(j,i),j=1,2), (sq(j,i),j=1,2)
          endif

c         Output to File

          write (iow,4100)  n, i, (zphg(j,i),j=1,3),
     &                     (ce(j,i),j=1,3), (sn(j,i),j=1,3),
     &                     (cr(j,i),j=1,3), (sm(j,i),j=1,3),
     &                     (cx(j,i),j=1,2), (sq(j,i),j=1,2)

        end do ! i

c       Blank record between elements

        if (ior.lt.0) write (*,*) ' '
        write (iow,*) ' '


c     Compute Nodal Stresses:

      elseif(isw.eq.8) then

c       Project Stresses

        call sh3fplst ( sn(1,1)   , sq(1,1)   , sm(1,1),
     &                  hr(np(35)), hr(np(36)), shp(1,1), xjw(1), lint)

      endif

c     Element Stress Format:

4000  format(/
     & ' ----------------------------------------',
     & '------------------------------------- '/,
     & '  Dynamic Finite-Deformation Elastic ',
     & '4-Node Shell Element -- Gauss Point Info  '/,
     & ' ----------------------------------------',
     & '------------------------------------- '//,
     & '  Elment Num   GaussPnt #   Mid-surf-X ',
     & '  Mid-surf-Y   Mid-surf-Z'/,
     & '  MemStrn-xx   MemStrn-yy   MemStrn-xy ',
     & '  MemStrs-xx   MemStrs-yy   MemStrs-xy  '/,
     & '  BndStrn-xx   BndStrn-yy   BndStrn-xy ',
     & '  BndStrs-xx   BndStrs-yy   BndStrs-xy  '/,
     & '  ShrStrn-1    ShrStrn-2    ShrStrs-1  ',
     & '  ShrStrs-2                             '/,
     & ' ------------ ------------ ------------',
     & ' ------------ ------------ ------------ '/)

4100  format(i7,12x,i1,6x,1p,3e13.5/(1p,6e13.5))

      end
