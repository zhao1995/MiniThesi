c$Id:$
      subroutine resid3d(xsj,shp,sig,d,xl,vl,al,r,ndm,ndf, aflg)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Revise to have 'aflg' option for inertial effect 28/12/2007
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: 3-D residual routine

c      Inputs:
c        xsj       - Jacobian
c        shp(4,*)  - Shape functions
c        sig(*)    - Stress
c        d(*)      - Material parameters
c        xl(ndm,*) - Nodal coordinates
c        vl(ndf,*) - Nodal velocities
c        al(ndf,*) - Nodal accelerations
c        ndm       - Mesh dimentions
c        ndf       - Degree of freedoms at nodes
c        aflg      - Include acceleration/velocity effects if true

c      Outputs:
c        r(ndf,*)  - Residual
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'eldata.h'
      include  'eltran.h'

      logical   aflg
      integer   ndm,ndf, i,j
      real*8    aj0,dvm,lfac,cfac,xsj
      real*8    d(*),xl(ndm,*),vl(ndf,*),al(ndf,*),r(ndf,*)
      real*8    shp(4,*),sig(*),xx(3),ac(3),vc(3), bf(3),bt(3,3)

      save

c     Compute body force values

      do i = 1,ndm
        bf(i) = 0.0d0
      end do ! i
      call sbodyf(d, bf)

c     Angular velocity body force: d(4) = density, d(65) = omega

      if(d(4).gt.0.0d0 .and. d(65).gt.0.0d0) then
        do i = 1,3
          xx(i) = 0.0d0
          do j = 1,nel
            xx(i) = xx(i) + shp(4,j)*xl(i,j)
          end do ! j
        end do ! i
        call sbodyw(d(4),d(65),xx, bf,bt, .false.)
      endif

c     Compute gravity, thermal, inertia, and stress contributions

      do j = 1,nel
        r(1,j) = r(1,j) + (bf(1)*shp(4,j) - sig(1)*shp(1,j)
     &                                    - sig(4)*shp(2,j)
     &                                    - sig(6)*shp(3,j))*xsj
        r(2,j) = r(2,j) + (bf(2)*shp(4,j) - sig(4)*shp(1,j)
     &                                    - sig(2)*shp(2,j)
     &                                    - sig(5)*shp(3,j))*xsj
        r(3,j) = r(3,j) + (bf(3)*shp(4,j) - sig(6)*shp(1,j)
     &                                    - sig(5)*shp(2,j)
     &                                    - sig(3)*shp(3,j))*xsj
      end do ! j

c     Inertial/Rayleigh mass damping effect

      if(aflg .and. d(7).ge.0.0d0) then
        cfac = d(7)
        lfac = 1.d0 - cfac

c       Compute accelerations

        do i = 1,3
          ac(i) = 0.0d0
          do j = 1,nel
            ac(i) = ac(i) + shp(4,j)*al(i,j)
          end do ! j
        end do ! i

c       Compute inertia contributions

        dvm = d(4)*xsj
        do j = 1,nel
          aj0 = shp(4,j)*dvm
          do i = 1,3
            r(i,j) = r(i,j) - (cfac*ac(i) + lfac**al(i,j))*aj0
          end do ! i
        end do ! j

c       Compute Rayleigh mass damping

        if(d(77).ne.0.0d0) then
          do i = 1,3
            vc(i) = 0.0d0
            do j = 1,nel
              vc(i) = vc(i) + shp(4,j)*vl(i,j)
            end do ! j
          end do ! i
          dvm = dvm*d(77)
          do j = 1,nel
            aj0   = shp(4,j)*dvm
            do i = 1,3
              r(i,j) = r(i,j) - (cfac*vc(i) + lfac*vl(i,j))*aj0
            end do ! i
          end do ! j
        endif

      endif

      end
