c$Id:$
      subroutine bbar3m(phi,shp,jac,xji,lint,nel,npm,
     &                  hh,theta,shpbar)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Increase arrays to store 64 node brick.          03/02/2009
c          Dimension phi(10,*) for quadratic elements.
c       2. Increase arrays to store 125 node brick.         20/12/2010
c       3. Increase dimension on hh, etc.to 10              22/12/2010
c       4. Convert computation of theta to relative form    21/11/2011
c          i.e. theta(1) --> 1 + theta(3), xji dimensioned to 4.
c       5. Correct invert from 4 to 10                      08/05/2013
c-----[--.----+----.----+----.-----------------------------------------]
c     Purpose: Compute mixed formulation for the volumetric response

c     Inputs:
c        shp(4,125,*)   - Shape function & derivative at gauss point
c        vol(*)         - Volume elements at gauss points at t_n+1
c        xji(4,*)       - Jacobian determinant at gauss points
c        lint           - Number of quadrature points
c        nel            - Number of nodes on element
c        npm            - Number of pressure modes/element

c     Outputs:
c        hh(10,10)      - Reference config shape integrals (inverse)
c        theta(4,*)     - Mixed jacobian determinant for element
c        shpbar(3,125,*)- Mixed derivatives of shape functions.
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      integer   lint, nel, npm, i, j, k, l
      real*8    phi(10,*),shp(4,125,*),jac(*),theta(4,*)
      real*8    shpbar(3,125,*),gg(10,3,125),hh(10,10),ji(10,2),xji(4,*)
      real*8    hj(10,2), hg(10,3,125), dvol, h0,h1,h2,h3

c     Constant pressure elements

      if(npm.eq.1) then

        do j = 1,nel
          do i = 1,3
            shpbar(i,j,1) = 0.0d0
          end do ! i
        end do ! j
        hh(1,1)  = 0.d0
        h1       = 0.0d0
        h2       = 0.0d0
        h3       = 0.0d0

        do l = 1,lint

c         H-array and D-array

          dvol    = jac(l) * xji(1,l)
          hh(1,1) = hh(1,1) + jac(l)
          h1      = h1      + jac(l) * xji(1,l)
          h2      = h2      + jac(l) * xji(3,l)
          h3      = h3      + jac(l) * xji(4,l)

c         G-array

          do j = 1,nel
            do i = 1,3
              shpbar(i,j,1) = shpbar(i,j,1) + shp(i,j,l) * dvol
            end do ! i
          end do ! j
        end do ! l

c       Modify shpbar for B-bar type computations

        h0 = 1.d0/h1

        do j = 1,nel
          do i = 1,3
            shpbar(i,j,1) = shpbar(i,j,1)*h0
          end do ! i
        end do ! j

c       Average Jacobian at t_n+1, t_n

        hh(1,1)     = 1.d0 / hh(1,1)
        theta(3,1)  = h2  *  hh(1,1)
        theta(4,1)  = h3  *  hh(1,1)

        theta(1,1) = theta(3,1) + 1.0d0
        theta(2,1) = theta(4,1) + 1.0d0

        do l = 2,lint
          do j = 1,4
            theta(j,l) = theta(j,1)
          end do ! j
          do j = 1,nel
            do i = 1,3
              shpbar(i,j,l) = shpbar(i,j,1)
            end do ! i
          end do ! j
        end do ! l

c     Higher order elements

      else

        do i = 1,npm
          do j = 1,nel
            gg(i,1,j) = 0.0d0
            gg(i,2,j) = 0.0d0
            gg(i,3,j) = 0.0d0
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

            h0 = phi(j,l) * jac(l)
            h1 = h0 * xji(3,l)
            h2 = h0 * xji(4,l)
            h3 = h0 * xji(1,l)

c           Ji-array

            ji(j,1) = ji(j,1) + h1
            ji(j,2) = ji(j,2) + h2

c           H-array

            do i = 1,npm
              hh(i,j)    = hh(i,j)  + phi(i,l)*h0
            end do ! i

c           G-array

            do i = 1,nel
              do k = 1,3
                gg(j,k,i) = gg(j,k,i) + shp(k,i,l)*h3
              end do ! k
            end do ! i
          end do ! j

        end do ! l

c       Invert H-array

        call invert(hh,npm,10)

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
            hg(i,1,j) = 0.0d0
            hg(i,2,j) = 0.0d0
            hg(i,3,j) = 0.0d0
            do k = 1,npm
              hg(i,1,j) = hg(i,1,j) + hh(i,k)*gg(k,1,j)
              hg(i,2,j) = hg(i,2,j) + hh(i,k)*gg(k,2,j)
              hg(i,3,j) = hg(i,3,j) + hh(i,k)*gg(k,3,j)
            end do ! k
          end do ! i
        end do ! j

        do l = 1,lint
          theta(3,l) = hj(1,1)
          theta(4,l) = hj(1,2)
          do k = 2,npm
            theta(3,l) = theta(3,l) + phi(k,l)*hj(k,1)
            theta(4,l) = theta(4,l) + phi(k,l)*hj(k,2)
          end do ! k
          theta(1,l) = theta(3,l) + 1.0d0
          theta(2,l) = theta(4,l) + 1.0d0
          h0         = 1.d0/theta(1,l)
          do j = 1,nel
            do i = 1,3
              shpbar(i,j,l) = hg(1,i,j)
              do k = 2,npm
                shpbar(i,j,l) = shpbar(i,j,l) + phi(k,l)*hg(k,i,j)
              end do ! k
              shpbar(i,j,l) = h0*shpbar(i,j,l)
            end do ! k
          end do ! j
        end do ! l

      endif

      end
