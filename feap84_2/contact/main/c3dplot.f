c$Id:$
      subroutine c3dplot (x,ix1,ix2,pen1,pen2,ifsurf)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Coded/modified by:                           Date:          rel.:
c               Giorgio Zavarise           April 10, 1996            1.0

c      Acronym: Contact 3d GEOmetry PLot

c      Purpose: Plot geometry of 3D contact pair

c      Inputs:
c         x(*)    - Nodal coordinates in deformed position
c         ix1(*)  - Element nodal connection list for surface 1
c         ix2(*)  - Element nodal connection list for surface 2

c      Outputs:
c         none      - Plot outputs to screen/file
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'c_geom.h'
      include  'cdata.h'
      include  'sdata.h'
      include  'pdatay.h'

      integer   nsy,ke,kl,n1
      integer   ix1(dnope1,*),ix2(dnope2,*),pen1,pen2,ifsurf
      real*8    x(3,*)

      save

      call cdebug0 ('  c3dplot',-1)

c     Reflect coordinates for symmetry plots

      do nsy = 1,nsym
        lsym = isym(nsy)
        call pltsym(x,3,numnp,lsym)

c       Plot slave surfaces

        if (ifsurf.ne.2) then
          call pppcol(pen1,0)

          do ke = 1, neps1
c           call pppcol(pen1,0)
            n1 = ix1(nope1,ke)
            call plotl(x(1,n1),x(2,n1),x(3,n1),3)
            do kl = 1,nope1
              n1 = ix1(kl,ke)
              call plotl(x(1,n1),x(2,n1),x(3,n1),2)
            end do
          end do
        endif

c       Plot master surfaces

        if (ifsurf.ne.1) then
          call pppcol(pen2,0)

          do ke = 1, neps2
            n1 = ix2(nope2,ke)
            call plotl(x(1,n1),x(2,n1),x(3,n1),3)
            do kl = 1,nope2
              n1 = ix2(kl,ke)
              call plotl(x(1,n1),x(2,n1),x(3,n1),2)
            end do
          end do
        endif

c       Reset coordinates

        call pltsym(x,3,numnp,lsym)
      end do

      end
