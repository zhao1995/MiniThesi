c$Id:$
      subroutine kine3d1u(shp,dvol0,ul,f,fi,df,fdet,detf,shp0,shp1,
     &                    V0,ndf,nel,nen,lint)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    11/27/2007
c       1. Increase arrays to store 125 node brick          20/12/2010
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Compute kinematic quantities of uniform deformation
c               form

c      Reference: J. Bonet and P. Bhargava, "A uniform deformation
c                 gradient hexahedron element with artificial hourglass
c                 control", IJNME, Vol. 38, pp 2809-2828 (1995)

c      Inputs:
c         shp(4,125,*)  - Reference configuration shape functions
c         dvol0(*)      - Differential reference volume
c         ul(ndf,nen,*) - Nodal displacements
c         ndf           - Number dof/node
c         nel           - Number nodes/element
c         nen           - Maximum number nodes/element
c         lint          - Number of quadrature points

c      Outputs:
c         f(3,3,2)      - Deformation gradient
c         fi(3,3)       - Inverse deformation gradient
c         df(3,3)       - Incremental deformation gradient
c         fdet(*)       - Determinant F at Gauss points
c         detf(2)       - Determinant of deformation gradient
c         shp0(3,8)     - Averaged reference shape functions
c         shp1(3,8)     - Averaged current shape functions
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'pconstant.h'

      integer   ndf,nel,nen,lint, i,j,k,l
      real*8    detfi, V0, V0i
      real*8    shp(4,125,*),dvol0(*),ul(ndf,nen,*)
      real*8    df(3,3),f(3,3,2),dfl(3,3),fl(3,3,2),fi(3,3)
      real*8    shp0(3,8), shp1(3,8), fdet(*),detf(2)

      save

c     Deformation gradient at t_n+1 : F_n+1 = I + GRAD u_n+1

      do i = 1,3
        do j = 1,3
          f(j,i,1) = 0.0d0
          df(j,i)  = 0.0d0
        end do ! j
      end do ! i

c     Shape functions at t_n

      do i = 1,8
        do j = 1,3
          shp0(j,i)   = 0.0d0
        end do ! j
      end do ! i
      V0       = 0.0d0

      do l = 1,lint

c       Compute deformation gradient at Gauss point

        do i = 1,3
          do j = 1,3
            fl(j,i,1) = 0.0d0
            dfl(j,i)  = 0.0d0
          end do ! j
        end do ! i
        do i = 1,3
          do j = 1,3
            do k = 1,nel
              fl(i,j,1) = fl(i,j,1) + ul(i,k,1)*shp(j,k,l)
              dfl(i,j)  = dfl(i,j)  + ul(i,k,2)*shp(j,k,l)
            end do ! k
          end do ! j
        end do ! i

c       Accumulate integral of deformation gradient

        do i = 1,3
          do j = 1,3
            f(j,i,1) = f(j,i,1) + fl(j,i,1)*dvol0(l)
            df(j,i)  = df(j,i)  +  dfl(j,i)*dvol0(l)
          end do ! j
          fl(i,i,1) = fl(i,i,1) + 1.0d0
          do j = 1,3
            fl(j,i,2) = fl(j,i,1) - dfl(j,i)
          end do ! j
        end do ! i

c       Compute integral of reference shape functions

        do i = 1,8
          do j = 1,3
            shp0(j,i) = shp0(j,i) + shp(j,i,l)*dvol0(l)
          end do ! j
        end do ! i

c       Accumulate integral of determinant of deformation gradient
c       i = 1: det F at t_n+1; i = 2: det F at t_n

        do i = 1,2
          detf(i) = fl(1,1,i)*(fl(2,2,i)*fl(3,3,i)-fl(3,2,i)*fl(2,3,i))
     &            + fl(1,2,i)*(fl(2,3,i)*fl(3,1,i)-fl(3,3,i)*fl(2,1,i))
     &            + fl(1,3,i)*(fl(2,1,i)*fl(3,2,i)-fl(3,1,i)*fl(2,2,i))
        end do ! i

c       Save Gauss point det F at t_n+1

        fdet(l) = detf(1)

c       Invert F_n+1

        detfi     = 1.d0/detf(1)
        fi(1,1) = (fl(2,2,1)*fl(3,3,1)-fl(3,2,1)*fl(2,3,1))*detfi
        fi(1,2) = (fl(3,2,1)*fl(1,3,1)-fl(1,2,1)*fl(3,3,1))*detfi
        fi(1,3) = (fl(1,2,1)*fl(2,3,1)-fl(2,2,1)*fl(1,3,1))*detfi

        fi(2,1) = (fl(2,3,1)*fl(3,1,1)-fl(3,3,1)*fl(2,1,1))*detfi
        fi(2,2) = (fl(3,3,1)*fl(1,1,1)-fl(1,3,1)*fl(3,1,1))*detfi
        fi(2,3) = (fl(1,3,1)*fl(2,1,1)-fl(2,3,1)*fl(1,1,1))*detfi

        fi(3,1) = (fl(2,1,1)*fl(3,2,1)-fl(3,1,1)*fl(2,2,1))*detfi
        fi(3,2) = (fl(3,1,1)*fl(1,2,1)-fl(1,1,1)*fl(3,2,1))*detfi
        fi(3,3) = (fl(1,1,1)*fl(2,2,1)-fl(2,1,1)*fl(1,2,1))*detfi

c       Transform shape functions to configuration at t_n+1

        do j = 1,nel
          do i = 1,3
            shp1(i,j) = fi(1,i)*shp(1,j,l)
     &                + fi(2,i)*shp(2,j,l)
     &                + fi(3,i)*shp(3,j,l)
          end do ! i
        end do ! j

c       Return Gauss point shape functions in current configuration

        do j = 1,nel
          do i = 1,3
            shp(i,j,l) = shp1(i,j)
          end do ! i
        end do ! j

c       Accumulate integral of element volume

        V0 = V0 + dvol0(l)
      end do ! l

c     Uniform determinant and deformation gradients

      V0i   = 1.d0/V0
      do i = 1,3
        do j = 1,3
          f(j,i,1) = f(j,i,1)*V0i
          df(j,i)  = df(j,i)*V0i
        end do ! j
        f(i,i,1) = f(i,i,1) + 1.0d0
      end do ! i

c     Uniform shape function derivatives

      do i = 1,8
        do j = 1,3
          shp0(j,i)   = shp0(j,i)*V0i
        end do ! j
      end do ! i

c     Deformation gradient at t_n: F_n

      do i = 1,3
        do j = 1,3
          f(j,i,2) = f(j,i,1) - df(j,i)
        end do ! j
      end do ! i

c     Determinant of F

      do i = 1,2
        detf(i) = f(1,1,i)*(f(2,2,i)*f(3,3,i)-f(3,2,i)*f(2,3,i))
     &          + f(1,2,i)*(f(2,3,i)*f(3,1,i)-f(3,3,i)*f(2,1,i))
     &          + f(1,3,i)*(f(2,1,i)*f(3,2,i)-f(3,1,i)*f(2,2,i))
      end do ! i

c     Invert F_n+1

      detfi     = 1.d0/detf(1)
      fi(1,1) = (f(2,2,1)*f(3,3,1)-f(3,2,1)*f(2,3,1))*detfi
      fi(1,2) = (f(3,2,1)*f(1,3,1)-f(1,2,1)*f(3,3,1))*detfi
      fi(1,3) = (f(1,2,1)*f(2,3,1)-f(2,2,1)*f(1,3,1))*detfi

      fi(2,1) = (f(2,3,1)*f(3,1,1)-f(3,3,1)*f(2,1,1))*detfi
      fi(2,2) = (f(3,3,1)*f(1,1,1)-f(1,3,1)*f(3,1,1))*detfi
      fi(2,3) = (f(1,3,1)*f(2,1,1)-f(2,3,1)*f(1,1,1))*detfi

      fi(3,1) = (f(2,1,1)*f(3,2,1)-f(3,1,1)*f(2,2,1))*detfi
      fi(3,2) = (f(3,1,1)*f(1,2,1)-f(1,1,1)*f(3,2,1))*detfi
      fi(3,3) = (f(1,1,1)*f(2,2,1)-f(2,1,1)*f(1,2,1))*detfi

c     Transform shape functions to t_n+1 configuration

      do j = 1,nel
        do i = 1,3
          shp1(i,j) = fi(1,i)*shp0(1,j)
     &              + fi(2,i)*shp0(2,j)
     &              + fi(3,i)*shp0(3,j)
        end do ! i
      end do ! j

      end
