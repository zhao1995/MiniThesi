c$Id:$
      subroutine pfrac3f(f,detf,sig,weng,dvol, r,ndf,ndm,ntyp)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Fix dimensions on f() and def()                  16/12/2013
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Material Force Calculations for finite deformation
c               elements.

c      Inputs:
c        f(3,3,2,*)  - Deformation gradient
c        detf(2,*)   - Determinant of deformation gradient
c        sig(10,*)   - Cauchy stress
c        weng(*)     - stored energy
c        dvol(*)     - Volume elements
c        nel         - Number nodes on element
c        ndf         - Number dof/node (max)
c        ndm         - Spatial dimension of mesh
c        ntyp        - Stress type: 1 = 1st P-K; 2 = 2nd P-K; 3 = Cauchy

c      Outputs:
c        r(ndf,*)    - Element material force
c-----[--.----+----.----+----.-----------------------------------------]
      implicit   none

      include   'eldata.h'
      include   'iofile.h'
      include   'prstrs.h'
      include   'qudshp.h'

      integer    ndf,ndm,ntyp, i,j,k,l
      real*8     thi,thj,dvol0
      real*8     f(3,3,4,*),detf(4,*), dvol(*), r(ndf,*)
      real*8     fadj(3,3), sig(10,*),weng(*),enmom(3,3),PK1(3,3)

c     Perform quadrature

      do l = 1,lint

c       Set 1st-PK stress

        dvol0 = dvol(l)
        if(ntyp.eq.1) then
          PK1(1,1) = sig(1,l)
          PK1(2,1) = sig(2,l)
          PK1(3,1) = sig(3,l)
          PK1(1,2) = sig(4,l)
          PK1(2,2) = sig(5,l)
          PK1(3,2) = sig(6,l)
          PK1(1,3) = sig(7,l)
          PK1(2,3) = sig(8,l)
          PK1(3,3) = sig(9,l)

c       Set 1st-PK stress from 2nd-PK stress

        elseif(ntyp.eq.2) then
          PK1(1,1) = f(1,1,1,l)*sig(1,l) + f(1,2,1,l)*sig(4,l)
     &             + f(1,3,1,l)*sig(6,l)
          PK1(1,2) = f(1,1,1,l)*sig(4,l) + f(1,2,1,l)*sig(2,l)
     &             + f(1,3,1,l)*sig(5,l)
          PK1(1,3) = f(1,1,1,l)*sig(6,l) + f(1,2,1,l)*sig(3,l)
     &             + f(1,3,1,l)*sig(3,l)
          PK1(2,1) = f(2,1,1,l)*sig(1,l) + f(2,2,1,l)*sig(4,l)
     &             + f(2,3,1,l)*sig(6,l)
          PK1(2,2) = f(2,1,1,l)*sig(4,l) + f(2,2,1,l)*sig(2,l)
     &             + f(2,3,1,l)*sig(5,l)
          PK1(2,3) = f(2,1,1,l)*sig(6,l) + f(2,2,1,l)*sig(5,l)
     &             + f(2,3,1,l)*sig(3,l)
          PK1(2,3) = f(3,1,1,l)*sig(1,l) + f(3,2,1,l)*sig(4,l)
     &             + f(3,3,1,l)*sig(6,l)
          PK1(2,3) = f(3,1,1,l)*sig(4,l) + f(3,2,1,l)*sig(2,l)
     &             + f(3,3,1,l)*sig(5,l)
          PK1(2,3) = f(3,1,1,l)*sig(6,l) + f(3,2,1,l)*sig(5,l)
     &             + f(3,3,1,l)*sig(3,l)

c       Set 1st-PK stress from Cauchy stress

        elseif(ntyp.eq.3) then

c         Compute adjoint

          fadj(1,1) = f(2,2,1,l)*f(3,3,1,l) - f(2,3,1,l)*f(3,2,1,l)
          fadj(2,1) = f(2,3,1,l)*f(3,1,1,l) - f(2,1,1,l)*f(3,3,1,l)
          fadj(3,1) = f(2,1,1,l)*f(3,2,1,l) - f(2,2,1,l)*f(3,1,1,l)

          fadj(1,2) = f(3,2,1,l)*f(1,3,1,l) - f(3,3,1,l)*f(1,2,1,l)
          fadj(2,2) = f(3,3,1,l)*f(1,1,1,l) - f(3,1,1,l)*f(1,3,1,l)
          fadj(3,2) = f(3,1,1,l)*f(1,2,1,l) - f(3,2,1,l)*f(1,1,1,l)

          fadj(1,3) = f(1,2,1,l)*f(2,3,1,l) - f(1,3,1,l)*f(2,2,1,l)
          fadj(2,3) = f(1,3,1,l)*f(2,1,1,l) - f(1,1,1,l)*f(2,3,1,l)
          fadj(3,3) = f(1,1,1,l)*f(2,2,1,l) - f(1,2,1,l)*f(2,1,1,l)

c         Transform: P = J*sigma*f^T
          PK1(1,1)  = sig(1,l)*fadj(1,1)
     &              + sig(4,l)*fadj(1,2)
     &              + sig(6,l)*fadj(1,3)
          PK1(2,1)  = sig(4,l)*fadj(1,1)
     &              + sig(2,l)*fadj(1,2)
     &              + sig(5,l)*fadj(1,3)
          PK1(3,1)  = sig(6,l)*fadj(1,1)
     &              + sig(5,l)*fadj(1,2)
     &              + sig(3,l)*fadj(1,3)

          PK1(1,2)  = sig(1,l)*fadj(2,1)
     &              + sig(4,l)*fadj(2,2)
     &              + sig(6,l)*fadj(2,3)
          PK1(2,2)  = sig(4,l)*fadj(2,1)
     &              + sig(2,l)*fadj(2,2)
     &              + sig(5,l)*fadj(2,3)
          PK1(3,2)  = sig(6,l)*fadj(2,1)
     &              + sig(5,l)*fadj(2,2)
     &              + sig(3,l)*fadj(2,3)

          PK1(1,3)  = sig(1,l)*fadj(3,1)
     &              + sig(4,l)*fadj(3,2)
     &              + sig(6,l)*fadj(3,3)
          PK1(2,3)  = sig(4,l)*fadj(3,1)
     &              + sig(2,l)*fadj(3,2)
     &              + sig(5,l)*fadj(3,3)
          PK1(3,3)  = sig(6,l)*fadj(3,1)
     &              + sig(5,l)*fadj(3,2)
     &              + sig(3,l)*fadj(3,3)
          dvol0    = dvol0/detf(1,l)
        else
          call pzero(PK1,9)
          write(  *,*) '  *ERROR* in PFRAC3D: NTYP =',ntyp
          write(iow,*) '  *ERROR* in PFRAC3D: NTYP =',ntyp
          write(ilg,*) '  *ERROR* in PFRAC3D: NTYP =',ntyp
          call plstop()
        endif

c       Compute Energy-Momentum tensor

        enmom(1,1) = weng(l) -(f(1,1,1,l) - jshft)*PK1(1,1)
     &                       - f(2,1,1,l)         *PK1(2,1)
     &                       - f(3,1,1,l)         *PK1(3,1)
        enmom(1,2) =         -(f(1,1,1,l) - jshft)*PK1(1,2)
     &                       - f(2,1,1,l)         *PK1(2,2)
     &                       - f(3,1,1,l)         *PK1(3,2)
        enmom(1,3) =         -(f(1,1,1,l) - jshft)*PK1(1,3)
     &                       - f(2,1,1,l)         *PK1(2,3)
     &                       - f(3,1,1,l)         *PK1(3,3)
        enmom(2,1) =         - f(1,2,1,l)         *PK1(1,1)
     &                       -(f(2,2,1,l) - jshft)*PK1(2,1)
     &                       - f(3,2,1,l)         *PK1(3,1)
        enmom(2,2) = weng(l) - f(1,2,1,l)         *PK1(1,2)
     &                       -(f(2,2,1,l) - jshft)*PK1(2,2)
     &                       - f(3,2,1,l)         *PK1(3,2)
        enmom(2,3) =         - f(1,2,1,l)         *PK1(1,3)
     &                       -(f(2,2,1,l) - jshft)*PK1(2,3)
     &                       - f(3,2,1,l)         *PK1(3,3)
        enmom(3,1) =         - f(1,3,1,l)         *PK1(1,1)
     &                       - f(2,3,1,l)         *PK1(2,1)
     &                       -(f(3,3,1,l) - jshft)*PK1(3,1)
        enmom(3,2) =         - f(1,3,1,l)         *PK1(1,2)
     &                       - f(2,3,1,l)         *PK1(2,2)
     &                       -(f(3,3,1,l) - jshft)*PK1(3,2)
        enmom(3,3) = weng(l) - f(1,3,1,l)         *PK1(1,3)
     &                       - f(2,3,1,l)         *PK1(2,3)
     &                       -(f(3,3,1,l) - jshft)*PK1(3,3)

c       Convert shape function derivatives to material form

        if(ntyp.eq.3) then
          do k = 1,nel
            thi         = shp3(1,k,l)*f(1,1,1,l)
     &                  + shp3(2,k,l)*f(2,1,1,l)
     &                  + shp3(3,k,l)*f(3,1,1,l)
            thj         = shp3(1,k,l)*f(1,2,1,l)
     &                  + shp3(2,k,l)*f(2,2,1,l)
     &                  + shp3(3,k,l)*f(3,2,1,l)
            shp3(3,k,l) = shp3(1,k,l)*f(1,3,1,l)
     &                  + shp3(2,k,l)*f(2,3,1,l)
     &                  + shp3(3,k,l)*f(3,3,1,l)
            shp3(1,k,l) = thi
            shp3(2,k,l) = thj
          end do ! k
        endif

c       Accumulate integral

        do k = 1,nel
          do i = 1,ndm
            thi = shp3(i,k,l)*dvol0
            do j = 1,ndm
              r(j,k)   = r(j,k) - enmom(j,i)*thi
            end do ! j
          end do ! i
        end do ! k

      end do ! l

      end
