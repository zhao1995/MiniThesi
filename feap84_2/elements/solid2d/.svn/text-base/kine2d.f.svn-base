c$Id:$
      subroutine kine2d(shp,xl,ul,f,fi,df,detf,ndm,ndf,nel,nen,lint)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Dimension arrays for 36 nodes & quadrature pts.  20/03/2009
c       2. Dimension arrays for 64 nodes & quadrature pts.  04/05/2009
c       3. Add computation of displacement gradient and     21/11/2011
c          (det F - 1) at t_n and t_n+1
c       4. Initialize f(i,j,1,l) to zero                    07/01/2012
c       5. Correct computation of axisymmetric determinant  05/12/2012
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Compute kinematic quantities for finite deformations

c      Inputs:
c         shp(3,64,*)   - Reference configuration shape functions
c         xl(ndm,*)     - Nodal reference coordinates
c         ul(ndf,nen,*) - Nodal displacements
c         ndm           - Number mesh dimensions
c         ndf           - Number dof/node
c         nel           - Number nodes/element
c         nen           - Maximum number nodes/element
c         lint          - Number of quadrature points

c      Outputs:
c         f(3,3,4,*)    - Deformation gradient
c         fi(3,3,*)     - Inverse deformation gradient
c         df(3,3,*)     - Incremental deformation gradient
c         detf(4,*)     - Determinant of deformation gradient
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'pmod2d.h'

      integer   ndm,ndf,nel,nen,lint, i,j,k,l
      real*8    detfi,temp,xx1,uu1,du1
      real*8    shp(3,64,*),xl(ndm,*),ul(ndf,nen,*)
      real*8    df(3,3,*),f(3,3,4,*),fi(3,3,*),detf(4,*)

      save

c     Deformation gradient at t_n+1 : F_n+1 = I + GRAD u_n+1

      do l = 1,lint
        do i = 1,3
          do j = 1,3
            f(j,i,1,l) = 0.0d0
            f(j,i,3,l) = 0.0d0
            df(j,i,l)  = 0.0d0
          end do ! j
        end do ! i
        do i = 1,2
          do j = 1,2
            do k = 1,nel
              f(i,j,3,l) = f(i,j,3,l) + ul(i,k,1)*shp(j,k,l)
              df(i,j,l ) = df(i,j,l ) + ul(i,k,2)*shp(j,k,l)
            end do ! k
          end do ! j
          do j = 1,2
            f(i,j,1,l) = f(i,j,3,l)
          end do ! j
          f(i,i,1,l) = f(i,i,1,l) + 1.0d0
        end do ! i

c       Axisymmetry

        if(stype.eq.3 .or. stype.eq.8) then
          xx1 = 0.0d0
          uu1 = 0.0d0
          du1 = 0.0d0
          do k = 1,nel
            xx1 = xx1 + xl(1,k  )*shp(3,k,l)
            uu1 = uu1 + ul(1,k,1)*shp(3,k,l)
            du1 = du1 + ul(1,k,2)*shp(3,k,l)
          end do ! k
          df(3,3,l)  = du1/xx1
          f(3,3,3,l) = uu1/xx1
          f(3,3,1,l) = f(3,3,3,l) + 1.d0

c         Torsion terms

          if(stype.eq.8) then
            do k = 1,nel
              f(3,1,1,l) = f(3,1,1,l) + ul(3,k,1)*shp(1,k,l)
              f(3,2,1,l) = f(3,2,1,l) + ul(3,k,1)*shp(2,k,l)
              df(3,1,l)  = df(3,1,l)  + ul(3,k,2)*shp(1,k,l)
              df(3,2,l)  = df(3,2,l)  + ul(3,k,2)*shp(2,k,l)
            end do ! k
            uu1        = uu1 + xx1
            df(3,1,l)  = f(3,1,1,l)*du1 + df(3,1,l)*(uu1 - du1)
            df(3,2,l)  = f(3,2,1,l)*du1 + df(3,2,l)*(uu1 - du1)
            f(3,1,1,l) = f(3,1,1,l)*uu1
            f(3,2,1,l) = f(3,2,1,l)*uu1
            f(3,1,3,l) = f(3,1,1,l)
            f(3,2,3,l) = f(3,2,1,l)
          endif

c       Plane case

        else
          f(3,3,1,l) = 1.0d0
          f(3,3,2,l) = 1.0d0
          f(3,3,3,l) = 0.0d0
          f(3,3,4,l) = 0.0d0
          df(3,3,l)  = 0.0d0
        endif

c       Deformation gradient at t_n: F_n

        do i = 1,3
          do j = 1,3
            f(j,i,2,l) = f(j,i,1,l) - df(j,i,l)
            f(j,i,4,l) = f(j,i,3,l) - df(j,i,l)
          end do ! j
        end do ! i

c       Determinant of F (plane part)

        detf(1,l) = f(1,1,1,l)*f(2,2,1,l) - f(1,2,1,l)*f(2,1,1,l)
        detf(2,l) = f(1,1,2,l)*f(2,2,2,l) - f(1,2,2,l)*f(2,1,2,l)

c       Determinant of  F - 1

        detf(3,l) = f(1,1,3,l) + f(2,2,3,l) + f(3,3,3,l)
     &            + f(1,1,3,l) * f(2,2,3,l) - f(1,2,3,l) * f(2,1,3,l)
     &            + f(3,3,3,l) *(f(1,1,3,l) + f(2,2,3,l)
     &            + f(1,1,3,l) * f(2,2,3,l) - f(1,2,3,l) * f(2,1,3,l))
        detf(4,l) = f(1,1,4,l) + f(2,2,4,l) + f(3,3,4,l)
     &            + f(1,1,4,l) * f(2,2,4,l) - f(1,2,4,l) * f(2,1,4,l)
     &            + f(3,3,4,l) *(f(1,1,4,l) + f(2,2,4,l)
     &            + f(1,1,4,l) * f(2,2,4,l) - f(1,2,4,l) * f(2,1,4,l))

c       Invert F

        detfi   =  1.d0/detf(1,l)
        fi(1,1,l) =  f(2,2,1,l)*detfi
        fi(1,2,l) = -f(1,2,1,l)*detfi
        fi(1,3,l) =  0.0d0
        fi(2,1,l) = -f(2,1,1,l)*detfi
        fi(2,2,l) =  f(1,1,1,l)*detfi
        fi(2,3,l) =  0.0d0
        fi(3,3,l) =  1.0d0/f(3,3,1,l)

c       Axisymmetric determinants

        detf(1,l) = detf(1,l)*f(3,3,1,l)
        detf(2,l) = detf(2,l)*f(3,3,2,l)

        if(stype.eq.8) then
          detfi     = 1.d0/detf(1,l)
          fi(3,1,l) = (f(2,1,1,l)*f(3,2,1,l)
     &              -  f(2,2,1,l)*f(3,1,1,l))/detfi
          fi(3,2,l) = (f(3,1,1,l)*f(1,2,1,l)
     &              -  f(3,2,1,l)*f(1,1,1,l))/detfi
        else
          fi(3,1,l) =  0.0d0
          fi(3,2,l) =  0.0d0
        endif

c       Transform shape functions to current configuration

        do k = 1,nel
          temp       = fi(1,1,l)*shp(1,k,l) + fi(2,1,l)*shp(2,k,l)
          shp(2,k,l) = fi(1,2,l)*shp(1,k,l) + fi(2,2,l)*shp(2,k,l)
          shp(1,k,l) = temp
        end do ! k
      end do ! l

      end
