c$Id:$
      subroutine kine1m(shp,xl,ul,f,fi,df,detf,ndm,ndf,nel,nen,lint)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Compute kinematic quantities for finite deformations

c      Inputs:
c         shp(2,8,*)    - Reference configuration shape functions
c         xl(ndm,*)     - Nodal reference coordinates
c         ul(ndf,nen,*) - Nodal displacements
c         ndm           - Number mesh dimensions
c         ndf           - Number dof/node
c         nel           - Number nodes/element
c         nen           - Maximum number nodes/element
c         lint          - Number of quadrature points

c      Outputs:
c         f(3,3,2,*)    - Deformation gradient
c         fi(3,3,*)     - Inverse deformation gradient
c         df(3,3,*)     - Incremental deformation gradient
c         detf(2,*)     - Determinant of deformation gradient
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'pmod2d.h'

      integer   ndm,ndf,nel,nen,lint, i,j,k,l
      real*8    detfi,xx1,uu1,du1
      real*8    shp(2,8,*),xl(ndm,*),ul(ndf,nen,*)
      real*8    df(3,3,*),f(3,3,2,*),fi(3,3,*),detf(2,*)

      save

c     Deformation gradient at t_n+1 : F_n+1 = I + GRAD u_n+1

      do l = 1,lint
        do i = 1,3
          do j = 1,3
            f(j,i,1,l) = 0.0d0
            df(j,i,l)  = 0.0d0
          end do ! j
        end do ! i
        do k = 1,nel
          f(1,1,1,l) = f(1,1,1,l) + ul(1,k,1)*shp(1,k,l)
          df(1,1,l ) = df(1,1,l ) + ul(1,k,2)*shp(1,k,l)
        end do ! k
        f(1,1,1,l) = f(1,1,1,l) + 1.0d0
        f(2,2,1,l) = 1.0d0

c       Axisymmetry

        if(stype.eq.3) then
          xx1 = 0.0d0
          uu1 = 0.0d0
          du1 = 0.0d0
          do k = 1,nel
            xx1 = xx1 + xl(1,k  )*shp(2,k,l)
            uu1 = uu1 + ul(1,k,1)*shp(2,k,l)
            du1 = du1 + ul(1,k,2)*shp(2,k,l)
          end do ! k
          df(3,3,l)  = du1/xx1
          f(3,3,1,l) = uu1/xx1 + 1.d0
        else
          f(3,3,1,l) = 1.0d0
        endif

c       Deformation gradient at t_n: F_n

        do i = 1,3
          do j = 1,3
            f(j,i,2,l) = f(j,i,1,l) - df(j,i,l)
          end do ! j
        end do ! i

c       Invert F

        detf(1,l) = f(1,1,1,l)*f(2,2,1,l) - f(1,2,1,l)*f(2,1,1,l)
        detf(2,l) = f(1,1,2,l)*f(2,2,2,l) - f(1,2,2,l)*f(2,1,2,l)

        detfi   =  1.d0/detf(1,l)
        fi(1,1,l) =  f(2,2,1,l)*detfi
        fi(1,2,l) = -f(1,2,1,l)*detfi
        fi(1,3,l) =  0.0d0
        fi(2,1,l) = -f(2,1,1,l)*detfi
        fi(2,2,l) =  f(1,1,1,l)*detfi
        fi(2,3,l) =  0.0d0
        fi(3,3,l) =  1.0d0/f(3,3,1,l)

c       Determinants

        detf(1,l) = detf(1,l)*f(3,3,1,l)
        detf(2,l) = detf(2,l)*f(3,3,2,l)

        fi(3,1,l) =  0.0d0
        fi(3,2,l) =  0.0d0

c       Transform shape function to current configuration

        do k = 1,nel
          shp(1,k,l) = fi(1,1,l)*shp(1,k,l)
        end do ! k
      end do ! l

      end
