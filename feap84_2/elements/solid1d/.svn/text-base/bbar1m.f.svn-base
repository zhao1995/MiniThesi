c$Id:$
      subroutine bbar1m(phi,shp,shpr,dvol,xji,lint,nel,npm,
     &                  hh,theta,shpbar)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c     Purpose: Compute mixed formulation for the volumetric response

c     Inputs:
c        phi(3,*)       - Pressure/volume functions
c        shp(2,8,*)     - Shape functions & derivatives at gauss point
c        shpr(8,*)      - Shape functions divided by radius or zero
c        vol(*)         - Volume elements at gauss points at t_n+1
c        xji(2,*)       - Jacobian determinant at gauss points
c        lint           - Number of quadrature points
c        nel            - Number of nodes on element (should be 8)
c        npm            - Number of mixed modes

c     Outputs:
c        hh(3,3)        - Reference config. shape integrals (inverse)
c        theta(2,*)     - Mixed jacobian determinant for element
c        shpbar(8,*)    - Mixed derivatives of shape functions.
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'pmod2d.h'

      integer   lint,  nel,  npm,  i,  j,  k,  l
      real*8    phi(3,*),shp(2,8,*),shpr(8,*),dvol(*),theta(2,*)
      real*8    shpbar(8,*),gg(3,8),hh(3,3),ji(3,2),xji(2,*)
      real*8    hj(3,2), hg(3,8), voln, h0, h1, h2

      save

c     Constant pressure elements

      if(npm.eq.1) then

        do j = 1,nel
          shpbar(j,1) = 0.0d0
        end do ! j
        hh(1,1) = 0.d0
        h1      = 0.d0
        h2      = 0.d0

        do l = 1,lint

c         H-array and D-array

          voln    = dvol(l) / xji(1,l)
          hh(1,1) = hh(1,1) + voln
          h1      = h1      + dvol(l)
          h2      = h2      + voln*xji(2,l)

c         G-array

          do j = 1,nel
            shpbar(j,1) = shpbar(j,1)
     &                  + (shp(1,j,l) + shpr(j,l))* dvol(l)
          end do ! j
        end do ! l

c       Modify shpbar for B-bar type computations

        h0 = 1.d0/h1

        do j = 1,nel
          shpbar(j,1) = shpbar(j,1)*h0
        end do ! j

c       Average Jacobian

        hh(1,1)    = 1.d0 / hh(1,1)
        theta(1,1) = h1   * hh(1,1)
        theta(2,1) = h2   * hh(1,1)

        do l = 2,lint
          theta(1,l) = theta(1,1)
          theta(2,l) = theta(2,1)
          do j = 1,nel
            shpbar(j,l) = shpbar(j,1)
          end do ! j
        end do ! l

c     Higher order elements

      else

        do i = 1,npm
          do j = 1,nel
            gg(i,j) = 0.0d0
          end do ! j
          do j = 1,npm
            hh(j,i) = 0.0d0
          end do ! j
          ji(i,1) = 0.0d0
          ji(i,2) = 0.0d0
        end do ! i

c       Quadrature loop

        do l = 1,lint
          do j = 1,npm

            h1 = phi(j,l)*dvol(l)
            h0 = h1/xji(1,l)
            h2 = h1*xji(2,l)

c           Ji-array

            ji(j,1) = ji(j,1) + h1
            ji(j,2) = ji(j,2) + h1*xji(2,l)/xji(1,l)

c           H-array

            do i = 1,npm
              hh(i,j)    = hh(i,j)  + phi(i,l)*h0
            end do ! i

c           G-array

            do i = 1,nel
              gg(j,i) = gg(j,i) + (shp(1,i,l) + shpr(i,l))*h1
            end do ! i
          end do ! j

        end do ! l

c       Invert H-array

        call invert(hh,npm,3)

        do j = 1,2
          do i = 1,npm
            hj(i,j) = 0.0d0
            do k = 1,npm
              hj(i,j) = hj(i,j) + hh(i,k)*ji(k,j)
            end do ! k
          end do ! i
        end do ! j

        do j = 1,nel
          do i = 1,npm
            hg(i,j) = 0.0d0
            do k = 1,npm
              hg(i,j) = hg(i,j) + hh(i,k)*gg(k,j)
            end do ! k
          end do ! i
        end do ! j

        do l = 1,lint
          theta(1,l) = hj(1,1)
          theta(2,l) = hj(1,2)
          do k = 2,npm
            theta(1,l) = theta(1,l) + phi(k,l)*hj(k,1)
            theta(2,l) = theta(2,l) + phi(k,l)*hj(k,2)
          end do ! k
          h0         = 1.d0/theta(1,l)
          do j = 1,nel
            shpbar(j,l) = hg(1,j)
            do k = 2,npm
              shpbar(j,l) = shpbar(j,l) + phi(k,l)*hg(k,j)
            end do ! k
            shpbar(j,l) = h0*shpbar(j,l)
          end do ! j
        end do ! l

      endif

      end
