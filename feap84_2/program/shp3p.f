c$Id:$
      subroutine shp3p(sg,xl,shps,xsj,nef)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Add shape for 7-node surface triangle            06/09/2007
c       2. Compute jacobian determinant in xsj(4)           01/11/2007
c       3. Correct computation of jacobian xsj(4)           05/11/2007
c       4. Add shape for 16-node surface quadrilateral      07/02/2009
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Compute 2-d element shape functions

c      Inputs:
c        sg(3)   - Parent coordinates
c        xl(3,*) - Nodal coordinates
c        nef     - Number of face nodes

c      Outputs:
c        shps(*) - Surface shape functions
c        xsj(4)  - Surface normal and area
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'pconstant.h'

      integer   j,k,nef, xi1(16),xi2(16)
      real*8    xi1s2,xi2s2,xi1s9,xi2s9
      real*8    sg(3),xl(3,*),shps(*),xsj(4),x1(3),x2(3),xm(3)
      real*8    n1(7),n2(7), d1(8),d2(8), x7(3)

      save

      data       xi1/1,2,2,1,3,4,2,2,4,3,1,1,3,4,4,3/
      data       xi2/1,1,2,2,1,1,3,4,2,2,4,3,3,3,4,4/

c     Linear triangle

      if(nef.eq.3) then
        do j = 1,3
          x1(j)   = xl(j,1) - xl(j,3)
          x2(j)   = xl(j,2) - xl(j,3)
        end do ! j
        shps(1) = sg(1)
        shps(2) = sg(2)
        shps(3) = 1.d0 - sg(1) - sg(2)

c     Quadratic triangle

      elseif(nef.eq.6 .or. nef.eq.7) then

        d1(1) =  4.0d0 * sg(1)
        d1(2) =  4.0d0 * sg(2)
        d1(3) =  4.0d0 - d1(1) - d1(2)

        n1(1) =  1.0d0 - d1(1)
        n1(2) =  0.0d0
        n1(3) = -1.0d0 + d1(3)
        n1(4) = -d1(2)
        n1(5) =  d1(2)
        n1(6) =  d1(1) - d1(3)

        n2(1) =  0.0d0
        n2(2) =  1.0d0 - d1(2)
        n2(3) = -1.0d0 + d1(3)
        n2(4) = -d1(1)
        n2(5) =  d1(2) - d1(3)
        n2(6) =  d1(1)

c       Compute gradients

        if(nef.eq.6) then
          do j = 1,3
            x1(j) = n1(1)*xl(j,1) + n1(2)*xl(j,2) + n1(3)*xl(j,3)
     &            + n1(4)*xl(j,4) + n1(5)*xl(j,5) + n1(6)*xl(j,6)
            x2(j) = n2(1)*xl(j,1) + n2(2)*xl(j,2) + n2(3)*xl(j,3)
     &            + n2(4)*xl(j,4) + n2(5)*xl(j,5) + n2(6)*xl(j,6)
          end do ! j
          shps(7) = 0.0d0

c       Face node shape function and gradient

        else
          n1(7) = 1.6875d0*d1(2)*(d1(1) - d1(3))
          n2(7) = 1.6875d0*d1(1)*(d1(2) - d1(3))
          do j = 1,3
            x7(j) = xl(j,7) + one9 *(xl(j,1) + xl(j,2) + xl(j,3))
     &                      - four9*(xl(j,4) + xl(j,5) + xl(j,6))
          end do ! j
          do j = 1,3
            x1(j) = n1(1)*xl(j,1) + n1(2)*xl(j,2) + n1(3)*xl(j,3)
     &            + n1(4)*xl(j,4) + n1(5)*xl(j,5) + n1(6)*xl(j,6)
     &            + n1(7)*x7(j)
            x2(j) = n2(1)*xl(j,1) + n2(2)*xl(j,2) + n2(3)*xl(j,3)
     &            + n2(4)*xl(j,4) + n2(5)*xl(j,5) + n2(6)*xl(j,6)
     &            + n2(7)*x7(j)
          end do ! j
          shps(7) = 27.0d0*sg(1)*sg(2)*(1.d0 - sg(1) - sg(2))
        endif

        d1(1) =  sg(1)
        d1(2) =  sg(2)
        d1(3) =  1.d0 - d1(1) - d1(2)

        do j = 1,3
          k         = mod(j,3) + 1
          shps(j  ) = 2.0d0*d1(j)*d1(j) - d1(j) + one9 *shps(7)
          shps(j+3) = 4.0d0*d1(j)*d1(k)         - four9*shps(7)
        end do ! j

c     Linear quadrilateral

      elseif(nef.eq.4) then
        do j = 1,3
          xm(j) =   xl(j,1) - xl(j,2) + xl(j,3) - xl(j,4)
          x1(j) = (-xl(j,1) + xl(j,2) + xl(j,3) - xl(j,4)
     &          +   xm(j)*sg(2))*0.25d0
          x2(j) = (-xl(j,1) - xl(j,2) + xl(j,3) + xl(j,4)
     &          +   xm(j)*sg(1))*0.25d0
        end do ! j

        shps(1) = (0.5d0 - 0.5d0*sg(1))*(0.5d0 - 0.5d0*sg(2))
        shps(2) = (0.5d0 + 0.5d0*sg(1))*(0.5d0 - 0.5d0*sg(2))
        shps(3) = (0.5d0 + 0.5d0*sg(1))*(0.5d0 + 0.5d0*sg(2))
        shps(4) = (0.5d0 - 0.5d0*sg(1))*(0.5d0 + 0.5d0*sg(2))

c     Quadratic 8-node quadrilateral

      elseif(nef.eq.8) then

        shps(1) = (0.5d0 - 0.5d0*sg(1))*(0.5d0 - 0.5d0*sg(2))
     &          * (-sg(1) - sg(2) -1.d0)
        shps(2) = (0.5d0 + 0.5d0*sg(1))*(0.5d0 - 0.5d0*sg(2))
     &          * ( sg(1) - sg(2) -1.d0)
        shps(3) = (0.5d0 + 0.5d0*sg(1))*(0.5d0 + 0.5d0*sg(2))
     &          * ( sg(1) + sg(2) -1.d0)
        shps(4) = (0.5d0 - 0.5d0*sg(1))*(0.5d0 + 0.5d0*sg(2))
     &          * (-sg(1) + sg(2) -1.d0)
        shps(5) = (0.5d0 - 0.5d0*sg(1)*sg(1))*(1.d0 - sg(2))
        shps(6) = (0.5d0 - 0.5d0*sg(2)*sg(2))*(1.d0 + sg(1))
        shps(7) = (0.5d0 - 0.5d0*sg(1)*sg(1))*(1.d0 + sg(2))
        shps(8) = (0.5d0 - 0.5d0*sg(2)*sg(2))*(1.d0 - sg(1))

        d1(1)   = (sg(1) + 0.5d0*sg(2))*(0.5d0 - 0.5d0*sg(2))
        d1(2)   = (sg(1) - 0.5d0*sg(2))*(0.5d0 - 0.5d0*sg(2))
        d1(3)   = (sg(1) + 0.5d0*sg(2))*(0.5d0 + 0.5d0*sg(2))
        d1(4)   = (sg(1) - 0.5d0*sg(2))*(0.5d0 + 0.5d0*sg(2))
        d1(5)   = -sg(1)*(1.0d0 - sg(2))
        d1(6)   =  (0.5d0 - 0.5d0*sg(2)*sg(2))
        d1(7)   = -sg(1)*(1.0d0 + sg(2))
        d1(8)   = -(0.5d0 - 0.5d0*sg(2)*sg(2))

        d2(1)   = (sg(2) + 0.5d0*sg(1))*(0.5d0 - 0.5d0*sg(1))
        d2(2)   = (sg(2) - 0.5d0*sg(1))*(0.5d0 + 0.5d0*sg(1))
        d2(3)   = (sg(2) + 0.5d0*sg(1))*(0.5d0 + 0.5d0*sg(1))
        d2(4)   = (sg(2) - 0.5d0*sg(1))*(0.5d0 - 0.5d0*sg(1))
        d2(5)   = -(0.5d0 - 0.5d0*sg(1)*sg(1))
        d2(6)   = -sg(2)*(1.0d0 + sg(1))
        d2(7)   =  (0.5d0 - 0.5d0*sg(1)*sg(1))
        d2(8)   = -sg(2)*(1.0d0 - sg(1))

        do j = 1,3
          x1(j) = xl(j,1)*d1(1) + xl(j,2)*d1(2) + xl(j,3)*d1(3)
     &          + xl(j,4)*d1(4) + xl(j,5)*d1(5) + xl(j,6)*d1(6)
     &          + xl(j,7)*d1(7) + xl(j,8)*d1(8)
          x2(j) = xl(j,1)*d2(1) + xl(j,2)*d2(2) + xl(j,3)*d2(3)
     &          + xl(j,4)*d2(4) + xl(j,5)*d2(5) + xl(j,6)*d2(6)
     &          + xl(j,7)*d2(7) + xl(j,8)*d2(8)
        end do ! j

c     Quadratic 9-node quadrilateral

      elseif(nef.eq.9) then

        n1(1) =  0.5d0*sg(1)*(sg(1) - 1.d0)
        n1(2) =  0.5d0*sg(1)*(sg(1) + 1.d0)
        n1(3) =  1.0d0 - sg(1)*sg(1)

        n2(1) =  0.5d0*sg(2)*(sg(2) - 1.d0)
        n2(2) =  0.5d0*sg(2)*(sg(2) + 1.d0)
        n2(3) =  1.0d0 - sg(2)*sg(2)

        d1(1) =  sg(1) - 0.5d0
        d1(2) =  sg(1) + 0.5d0
        d1(3) = -sg(1) - sg(1)

        d2(1) =  sg(2) - 0.5d0
        d2(2) =  sg(2) + 0.5d0
        d2(3) = -sg(2) - sg(2)

        do j = 1,3
          x1(j) = (xl(j,1)*d1(1) + xl(j,2)*d1(2) + xl(j,5)*d1(3))*n2(1)
     &          + (xl(j,4)*d1(1) + xl(j,3)*d1(2) + xl(j,7)*d1(3))*n2(2)
     &          + (xl(j,8)*d1(1) + xl(j,6)*d1(2) + xl(j,9)*d1(3))*n2(3)

          x2(j) = (xl(j,1)*n1(1) + xl(j,2)*n1(2) + xl(j,5)*n1(3))*d2(1)
     &          + (xl(j,4)*n1(1) + xl(j,3)*n1(2) + xl(j,7)*n1(3))*d2(2)
     &          + (xl(j,8)*n1(1) + xl(j,6)*n1(2) + xl(j,9)*n1(3))*d2(3)
        end do ! j

        shps(1) = n1(1)*n2(1)
        shps(2) = n1(2)*n2(1)
        shps(3) = n1(2)*n2(2)
        shps(4) = n1(1)*n2(2)
        shps(5) = n1(3)*n2(1)
        shps(6) = n1(2)*n2(3)
        shps(7) = n1(3)*n2(2)
        shps(8) = n1(1)*n2(3)
        shps(9) = n1(3)*n2(3)

c     Cubic 16-node quadrilateral

      elseif(nef.eq.16) then

        xi1s9  = one9 - sg(1)*sg(1)
        xi2s9  = one9 - sg(2)*sg(2)
        xi1s2  = 1.d0 - sg(1)*sg(1)
        xi2s2  = 1.d0 - sg(2)*sg(2)

        n1(1)  = -9.d0*(1.d0 - sg(1))*xi1s9*0.0625d0
        n1(2)  = -9.d0*(1.d0 + sg(1))*xi1s9*0.0625d0
        n1(3)  = 27.d0*xi1s2*(one3 - sg(1))*0.0625d0
        n1(4)  = 27.d0*xi1s2*(one3 + sg(1))*0.0625d0

        n2(1)  = -9.d0*(1.d0 - sg(2))*xi2s9*0.0625d0
        n2(2)  = -9.d0*(1.d0 + sg(2))*xi2s9*0.0625d0
        n2(3)  = 27.d0*xi2s2*(one3 - sg(2))*0.0625d0
        n2(4)  = 27.d0*xi2s2*(one3 + sg(2))*0.0625d0

        do k = 1,16
          shps(k) = n1(xi1(k))*n2(xi2(k))
        end do ! k

        d1(1) = (  1.d0 + (18.d0 - 27.d0*sg(1))*sg(1))*0.0625d0
        d1(2) = ( -1.d0 + (18.d0 + 27.d0*sg(1))*sg(1))*0.0625d0
        d1(3) = (-27.d0 - (18.d0 - 81.d0*sg(1))*sg(1))*0.0625d0
        d1(4) = ( 27.d0 - (18.d0 + 81.d0*sg(1))*sg(1))*0.0625d0

        d2(1) = (  1.d0 + (18.d0 - 27.d0*sg(2))*sg(2))*0.0625d0
        d2(2) = ( -1.d0 + (18.d0 + 27.d0*sg(2))*sg(2))*0.0625d0
        d2(3) = (-27.d0 - (18.d0 - 81.d0*sg(2))*sg(2))*0.0625d0
        d2(4) = ( 27.d0 - (18.d0 + 81.d0*sg(2))*sg(2))*0.0625d0

        do j = 1,3
          x1(j) = 0.0d0
          x2(j) = 0.0d0
          do k = 1,16
            x1(j) = x1(j) + xl(j,k)*d1(xi1(k))*n2(xi2(k))
            x2(j) = x2(j) + xl(j,k)*n1(xi1(k))*d2(xi2(k))
          end do ! k
        end do ! j

      endif

      xsj(1)  = x1(2)*x2(3) - x1(3)*x2(2)
      xsj(2)  = x1(3)*x2(1) - x1(1)*x2(3)
      xsj(3)  = x1(1)*x2(2) - x1(2)*x2(1)

      xsj(4)  = sqrt(xsj(1)*xsj(1) + xsj(2)*xsj(2) + xsj(3)*xsj(3))

      end
