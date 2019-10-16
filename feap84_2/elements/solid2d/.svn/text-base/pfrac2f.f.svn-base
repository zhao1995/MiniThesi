c$Id:$
      subroutine pfrac2f(f,detf,sig,weng,shp,dvol, r, lint,ndf,ndm,ntyp)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Updated dimensions of shp,f,detf,sig             13/12/2013
c       2. Removed computations of PK1(3,3)                 13/12/2013
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Material Force Calculations for finite deformation
c               elements.

c      Inputs:
c        f(3,3,2,*)  - Deformation gradient
c        detf(2,*)   - Determinant of deformation gradient
c        sig(9,*)    - Cauchy stress
c        weng(*)     - stored energy
c        shp(3,9,*)  - Shape function
c        dvol(*)     - Volume elements
c        lint        - Number quadrature points
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

      integer    lint,ndf,ndm,ntyp, i,j,k,l
      real*8     f(3,3,4,*),detf(4,*), shp(3,64,*),dvol(*), r(ndf,*)
      real*8     sig(16,*),weng(*),enmom(3,3),PK1(3,3), thi,dvol0

c     Perform quadrature

      do l = 1,lint

c       Set 1st-PK stress

        dvol0 = dvol(l)
        if(ntyp.eq.1) then
          PK1(1,1) = sig(1,l)
          PK1(2,1) = sig(2,l)
          PK1(1,2) = sig(3,l)
          PK1(2,2) = sig(4,l)

c       Set 1st-PK stress from 2nd-PK stress

        elseif(ntyp.eq.2) then
          PK1(1,1) = f(1,1,1,l)*sig(1,l) + f(1,2,1,l)*sig(4,l)
          PK1(1,2) = f(1,1,1,l)*sig(4,l) + f(1,2,1,l)*sig(2,l)
          PK1(2,1) = f(2,1,1,l)*sig(1,l) + f(2,2,1,l)*sig(4,l)
          PK1(2,2) = f(2,1,1,l)*sig(4,l) + f(2,2,1,l)*sig(2,l)

c       Set 1st-PK stress from Cauchy stress

        elseif(ntyp.eq.3) then
          PK1(1,1) =  f(2,2,1,l)*sig(1,l) - f(1,2,1,l)*sig(4,l)
          PK1(1,2) = -f(2,1,1,l)*sig(1,l) + f(1,1,1,l)*sig(4,l)
          PK1(2,1) =  f(2,2,1,l)*sig(4,l) - f(1,2,1,l)*sig(2,l)
          PK1(2,2) = -f(2,1,1,l)*sig(4,l) + f(1,1,1,l)*sig(2,l)
          dvol0    = dvol0/detf(1,l)
        else
          call pzero(PK1,9)
          write(  *,*) '  *ERROR* in PFRAC2F: NTYP -',ntyp
          write(iow,*) '  *ERROR* in PFRAC2F: NTYP -',ntyp
          write(ilg,*) '  *ERROR* in PFRAC2F: NTYP -',ntyp
          call plstop()
        endif

c       Compute Energy-Momentum tensor

        enmom(1,1) = weng(l) -(f(1,1,1,l) - jshft)*PK1(1,1)
     &                       - f(2,1,1,l)         *PK1(2,1)
        enmom(1,2) =         -(f(1,1,1,l) - jshft)*PK1(1,2)
     &                       - f(2,1,1,l)         *PK1(2,2)
        enmom(2,1) =         - f(1,2,1,l)         *PK1(1,1)
     &                       -(f(2,2,1,l) - jshft)*PK1(2,1)
        enmom(2,2) = weng(l) - f(1,2,1,l)         *PK1(1,2)
     &                       -(f(2,2,1,l) - jshft)*PK1(2,2)

c       Convert shape function derivatives to material form

        if(ntyp.eq.3) then
          do k = 1,nel
            thi        = shp(1,k,l)*f(1,1,1,l) + shp(2,k,l)*f(2,1,1,l)
            shp(2,k,l) = shp(1,k,l)*f(1,2,1,l) + shp(2,k,l)*f(2,2,1,l)
            shp(1,k,l) = thi
          end do ! k
        endif

c       Accumulate integral

        do k = 1,nel
          do i = 1,ndm
            thi = shp(i,k,l)*dvol0
            do j = 1,ndm
              r(j,k)   = r(j,k) - enmom(j,i)*thi
            end do ! j
          end do ! i
        end do ! k

      end do ! l

      end
