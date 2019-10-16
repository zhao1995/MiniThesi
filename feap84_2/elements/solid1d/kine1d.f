c$Id:$
      subroutine kine1d (shp,xl,ul,f,fi,df,detf,ndm,ndf,nel,nen)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Compute kinematic quantities for finite deformations

c      Inputs:
c         shp(2,nel)  - Reference configuration shape functions
c         xl(ndm,nel) - Nodal reference coordinates
c         ul(ndf,nel) - Nodal displacements
c         ndm         - Number mesh dimensions
c         ndf         - Number dof/node
c         nel         - Number nodes/element
c         nen         - Maximum number nodes/element

c      Outputs:
c         f(3,3,2)    - deformation gradient
c         fi(3,3)     - inverse deformation gradient
c         df(3,3)     - incremental deformation gradient
c         detf(2)     - determinant of deformation gradient
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'elcoor.h'
      include  'pmod2d.h'

      integer   ndm,ndf,nel,nen, i,j,k
      real*8    shp(2,*),xl(ndm,*),ul(ndf,nen,*)
      real*8    df(3,3),f(3,3,*),fi(3,3),detf(*),xx1

      save

c     Deformation gradient at t_n+1 : F_n+1 = I + GRAD u_n+1

      do i = 1,3
        do j = 1,3
          f(j,i,1) = 0.0d0
          f(j,i,2) = 0.0d0
          fi(j,i)  = 0.0d0
          df(j,i)  = 0.0d0
        end do ! j
        f(i,i,1) = f(i,i,1) + 1.0d0
        f(i,i,2) = f(i,i,2) + 1.0d0
      end do ! i
      xref(1) = 0.0d0
      xcur(1) = 0.0d0
      do j = 1,nel
        f(1,1,1) = f(1,1,1) + ul(1,j,1)*shp(1,j)
        df(1,1 ) = df(1,1 ) + ul(1,j,2)*shp(1,j)
        xref(1) = xref(1)   + xl(1,j)*  shp(2,j)
        xcur(1) = xcur(1)   + ul(1,j,1)*shp(2,j)
      end do ! j
      xref(2) = 0.0d0
      xref(3) = 0.0d0
      xcur(1) = xcur(1) + xref(1)
      xcur(2) = 0.0d0
      xcur(3) = 0.0d0

c     Deformation gradient at t_n: F_n

      f(1,1,2)  = f(1,1,1) - df(1,1)

      if(stype.eq.3) then
        f(3,3,1) = 0.0d0
        xx1      = 0.0d0
        df(3,3)  = 0.0d0
        do k = 1,nel
          xx1      = xx1      + xl(1,k  )*shp(2,k)
          f(3,3,1) = f(3,3,1) + ul(1,k,1)*shp(2,k)
          df(3,3)  = df(3,3)  + ul(1,k,2)*shp(2,k)
        end do ! k
        f(3,3,1) = 1.d0 + f(3,3,1)/xx1
        df(3,3)  = df(3,3)/xx1
        f(3,3,2) = f(3,3,1) - df(3,3)
      endif

c     Invert F

      detf(1) = f(1,1,1)
      detf(2) = f(1,1,2)

      fi(1,1) =  1.0d0/f(1,1,1)
      fi(2,2) =  1.0d0
      fi(3,3) =  1.0d0/f(3,3,1)

c     Determinants

      detf(1) = detf(1)*f(3,3,1)
      detf(2) = detf(2)*f(3,3,2)

c     Transform shape functions to current configuration

      do k = 1,nel
        shp(1,k) = fi(1,1)*shp(1,k)
      end do ! k

      end
