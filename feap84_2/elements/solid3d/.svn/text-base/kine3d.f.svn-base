c$Id:$
      subroutine kine3d(shp,ul,f,fi,df,detf,ndf,nel,nen,lint)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Increase arrays to store 64 node brick           03/02/2009
c       2. Increase arrays to store 125 node brick          20/12/2010
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Compute kinematic quantities for finite deformations

c      Inputs:
c         shp(4,125,*)  - Reference configuration shape functions
c         ul(ndf,nen,*) - Nodal displacements
c         ndf           - Number dof/node
c         nel           - Number nodes/element
c         nen           - Maximum number nodes/element
c         lint          - Number of quadrature points

c      Outputs:
c         shp(4,125,*)  - Spatial configuration shape functions
c         f(3,3,4,*)    - Displacement gradient
c         fi(3,3,*)     - Inverse deformation gradient
c         df(3,3,*)     - Incremental deformation gradient
c         detf(4,*)     - Determinant of deformation gradient
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      integer   ndf,nel,nen,lint, i,j,k,l
      real*8    detfi
      real*8    shp(4,125,*),ul(ndf,nen,*)
      real*8    df(3,3,*),f(3,3,4,*),fi(3,3,*),detf(4,*), cc(3)

      save

c     Deformation gradient at t_n+1 : F_n+1 = I + GRAD u_n+1

      do l = 1,lint
        do i = 1,3
          do j = 1,3
            f(j,i,3,l) = 0.0d0
            df(j,i,l)  = 0.0d0
          end do ! j
        end do ! i
        do i = 1,3
          do j = 1,3
            do k = 1,nel
              f(i,j,3,l) = f(i,j,3,l) + ul(i,k,1)*shp(j,k,l)
              df(i,j,l ) = df(i,j,l ) + ul(i,k,2)*shp(j,k,l)
            end do ! k
            f(i,j,4,l) = f(i,j,3,l) - df(i,j,l)
          end do ! j
        end do ! i

c       Displacement gradient at t_n: F_n

        do i = 1,3
          do j = 1,3
            f(j,i,1,l)   = f(j,i,3,l)
            f(j,i,2,l)   = f(j,i,4,l)
          end do ! j
          f(i,i,1,l) = f(i,i,1,l) + 1.0d0
          f(i,i,2,l) = f(i,i,2,l) + 1.0d0
        end do ! i

c       Determinant of  F - 1

        do i = 3,4
          detf(i,l) = f(1,1,i,l) + f(2,2,i,l) + f(3,3,i,l)
     &              + f(1,1,i,l)*f(2,2,i,l) - f(1,2,i,l)*f(2,1,i,l)
     &              + f(2,2,i,l)*f(3,3,i,l) - f(2,3,i,l)*f(3,2,i,l)
     &              + f(3,3,i,l)*f(1,1,i,l) - f(3,1,i,l)*f(1,3,i,l)
     &              + f(1,1,i,l)*(f(2,2,i,l)*f(3,3,i,l)
     &                          - f(2,3,i,l)*f(3,2,i,l))
     &              + f(1,2,i,l)*(f(2,3,i,l)*f(3,1,i,l)
     &                          - f(2,1,i,l)*f(3,3,i,l))
     &              + f(1,3,i,l)*(f(2,1,i,l)*f(3,2,i,l)
     &                          - f(2,2,i,l)*f(3,1,i,l))
        end do ! i

c       Deformation gradients at t_n+1 and t_n

        detf(1,l) = detf(3,l) + 1.0d0
        detf(2,l) = detf(4,l) + 1.0d0

c       Invert F

        detfi     = 1.d0/detf(1,l)
        fi(1,1,l) = (f(2,2,1,l)*f(3,3,1,l)-f(3,2,1,l)*f(2,3,1,l))*detfi
        fi(1,2,l) = (f(3,2,1,l)*f(1,3,1,l)-f(1,2,1,l)*f(3,3,1,l))*detfi
        fi(1,3,l) = (f(1,2,1,l)*f(2,3,1,l)-f(2,2,1,l)*f(1,3,1,l))*detfi

        fi(2,1,l) = (f(2,3,1,l)*f(3,1,1,l)-f(3,3,1,l)*f(2,1,1,l))*detfi
        fi(2,2,l) = (f(3,3,1,l)*f(1,1,1,l)-f(1,3,1,l)*f(3,1,1,l))*detfi
        fi(2,3,l) = (f(1,3,1,l)*f(2,1,1,l)-f(2,3,1,l)*f(1,1,1,l))*detfi

        fi(3,1,l) = (f(2,1,1,l)*f(3,2,1,l)-f(3,1,1,l)*f(2,2,1,l))*detfi
        fi(3,2,l) = (f(3,1,1,l)*f(1,2,1,l)-f(1,1,1,l)*f(3,2,1,l))*detfi
        fi(3,3,l) = (f(1,1,1,l)*f(2,2,1,l)-f(2,1,1,l)*f(1,2,1,l))*detfi

c       Transform shape functions to current configuration

        do k = 1,nel
          do i = 1,3
            cc(i) = fi(1,i,l)*shp(1,k,l)
     &            + fi(2,i,l)*shp(2,k,l)
     &            + fi(3,i,l)*shp(3,k,l)
          end do ! i
          do i = 1,3
            shp(i,k,l) = cc(i)
          end do ! i
        end do ! k
      end do ! l

      end
