c$Id:$
      subroutine cnproj(ids, xs,xm,  xi,xp,t1,t2,gp, its)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Coded/modified by:                           Date:          rel.:
c               Robert L. Taylor             October 1998            1.0

c     Compute the projected point on quadrilateral

c     Inputs: ids(3)    - Boundary Condition indicators
c             xs (3)    - slave  postion at curent time, t_n+1
c             xm (3,4)  - master postion at curent time, t_n+1

c     Output: xi(2)     - Local projection points xi_1 and xi_2
c             xp(3)     - Master surface projection point at t_n+1
c             t1(3)     - Tangent vector 1
c             t2(3)     - Tangent vector 2
c             gp(3)     - Gap vector
c             its       - number of iterations to obtain solution
c-----[--.----+----.----+----.-----------------------------------------]
      implicit none

      include 'c_0.h'
      include 'c_contac.h'

      logical  noconv,   constr
      integer  i,j,      its,      nits,    ids(3)
      real*8   rdet,     tol,      lamtol,   nrmdlm
      real*8   xs(3),    xp(3),    xm(3,4),  shp(4),  dshp(3,4)
      real*8   t1(3),    t2(3),    gp(3),    dt(3),   xj(2),   xk(2)
      real*8   r1(2),    xi(2),    dxi(2),   lam(3),  a11(2,2)
      real*8   a12(2,3), a22(3,3), ai(2,3), bb(3,3), bi(3,3), dlam(3)
      real*8   r2(3)

      data     tol/0.5d-08/, nits/25/

      save

c     Initialization: Use mid point

      xi(1)  =  0.0d0
      xi(2)  =  0.0d0
      lam(1) =  0.0d0
      lam(2) =  0.0d0
      lam(3) =  0.0d0
      nrmdlm =  0.0d0
      lamtol =  1.0d0

c     Perform Newton iteration

      its    =  0
      noconv = .true.

      do while (noconv)

        its = its + 1

c       Surface shape functions and their derivatives
c       N.B. Range for natural coordinates:  -1 < xi < 1

        xj(1)     = 0.5d0 - 0.5d0*xi(1)
        xj(2)     = 0.5d0 - 0.5d0*xi(2)
        xk(1)     = 0.5d0 + 0.5d0*xi(1)
        xk(2)     = 0.5d0 + 0.5d0*xi(2)

        shp(1)    = xj(1)*xj(2)
        shp(2)    = xk(1)*xj(2)
        shp(3)    = xk(1)*xk(2)
        shp(4)    = xj(1)*xk(2)

        dshp(1,1) = -0.5d0*xj(2)
        dshp(1,2) = -dshp(1,1)
        dshp(1,3) =  0.5d0*xk(2)
        dshp(1,4) = -dshp(1,3)

        dshp(2,1) = -0.5d0*xj(1)
        dshp(2,2) = -0.5d0*xk(1)
        dshp(2,3) = -dshp(2,2)
        dshp(2,4) = -dshp(2,1)

c       Form basis vectors and their derivatives

        do i =  1,3
          t1(i) = dshp(1,1)*xm(i,1) + dshp(1,2)*xm(i,2)
     &          + dshp(1,3)*xm(i,3) + dshp(1,4)*xm(i,4)
          t2(i) = dshp(2,1)*xm(i,1) + dshp(2,2)*xm(i,2)
     &          + dshp(2,3)*xm(i,3) + dshp(2,4)*xm(i,4)
          xp(i) = shp(1)*xm(i,1) + shp(2)*xm(i,2)
     &          + shp(3)*xm(i,3) + shp(4)*xm(i,4)
          gp(i) = xs(i) - xp(i)
          dt(i) = 0.25d0*(xm(i,1) - xm(i,2) + xm(i,3) - xm(i,4))
        end do ! i

c       Tangent matrix

        a11(1,1) = t1(1)*t1(1) + t1(2)*t1(2) + t1(3)*t1(3)
        a11(1,2) = t1(1)*t2(1) + t1(2)*t2(2) + t1(3)*t2(3)
     &           - gp(1)*dt(1) - gp(2)*dt(2) - gp(3)*dt(3)
        a11(2,2) = t2(1)*t2(1) + t2(2)*t2(2) + t2(3)*t2(3)

c       Residuals

        r1(1) = t1(1)*gp(1) + t1(2)*gp(2) + t1(3)*gp(3)
        r1(2) = t2(1)*gp(1) + t2(2)*gp(2) + t2(3)*gp(3)

c       Lagrange multiplier modification to tangent and residual

        constr = .false.
        if(max(ids(1),ids(2),ids(3)).gt.0) then
          do i = 1,3
            if(ids(i).le.0) then
              constr = .true.
              exit
            endif
          end do ! i
        endif

        if(constr) then
          do i = 1,3
            do j = 1,3
              a22(j,i) = 0.0d0
            end do ! j
          end do ! i
          do i = 1,3
            if(ids(i).le.0) then
              a12(1,i) = -t1(i)
              a12(2,i) = -t2(i)
              a11(1,2) =  a11(1,2) - lam(i)*dt(i)
              r1(1)    =  r1(1)    + lam(i)*t1(i)
              r1(2)    =  r1(2)    + lam(i)*t2(i)
              r2(i)    =  gp(i)
            else
              a12(1,i) =  0.0d0
              a12(2,i) =  0.0d0
              a22(i,i) =  1.0d0
              r2(i)    =  0.0d0
            endif
          end do ! i
        endif ! constr

c       Compute incremental solution for xi

        a11(2,1) = a11(1,2)
        rdet     = 1.d0/(a11(1,1)*a11(2,2) - a11(1,2)*a11(2,1))

        dxi(1)   = ( a11(2,2)*r1(1) - a11(1,2)*r1(2))*rdet
        dxi(2)   = (-a11(2,1)*r1(1) + a11(1,1)*r1(2))*rdet

c       Lagrange multiplier modification to solution

        if(constr) then

          do i = 1,3
            ai(1,i) = ( a11(2,2)*a12(1,i) - a11(1,2)*a12(2,i))*rdet
            ai(2,i) = (-a11(2,1)*a12(1,i) + a11(1,1)*a12(2,i))*rdet
          end do

          do i = 1,3
            r2(i) = r2(i) + a12(1,i)*dxi(1) + a12(2,i)*dxi(2)
            do j = 1,3
              bb(j,i) = a22(j,i) + a12(1,j)*ai(1,i) + a12(2,j)*ai(2,i)
            end do ! j
          end do ! i

c         Adjoints for symmetric inverse of bb

          bi(1,1) = bb(2,2)*bb(3,3) - bb(2,3)*bb(3,2)
          bi(2,1) = bb(2,3)*bb(3,1) - bb(2,1)*bb(3,3)
          bi(3,1) = bb(2,1)*bb(3,2) - bb(2,2)*bb(3,1)

          bi(1,2) = bi(2,1)
          bi(2,2) = bb(3,3)*bb(1,1) - bb(3,1)*bb(1,3)
          bi(3,2) = bb(3,1)*bb(1,2) - bb(3,2)*bb(1,1)

          bi(1,3) = bi(3,1)
          bi(2,3) = bi(3,2)
          bi(3,3) = bb(1,1)*bb(2,2) - bb(1,2)*bb(2,1)

          rdet = 1.d0/(bi(1,1)*bb(1,1)+bi(1,2)*bb(2,1)+bi(1,3)*bb(3,1))

c         Compute Lagrange multiplier updates

          do i = 1,3
            dlam(i) = (bi(i,1)*r2(1)+bi(i,2)*r2(2)+bi(i,3)*r2(3))*rdet
            lam(i)  = lam(i) + dlam(i)
          end do ! i

          nrmdlm =      abs(dlam(1)) + abs(dlam(2)) + abs(dlam(3))
          lamtol = tol*(abs(lam(1))  + abs(lam(2))  + abs(lam(3)))

c         Update solution to xi from Lagrange multipliers

          dxi(1) = dxi(1) + ai(1,1)*dlam(1)
     &                    + ai(1,2)*dlam(2)
     &                    + ai(1,3)*dlam(3)
          dxi(2) = dxi(2) + ai(2,1)*dlam(1)
     &                    + ai(2,2)*dlam(2)
     &                    + ai(2,3)*dlam(3)

        endif

c       Update natural coordinates

        xi(1)  = xi(1) + dxi(1)
        xi(2)  = xi(2) + dxi(2)

c       Check if leaving facet

        if(max(abs(xi(1)),abs(xi(2))).gt.2.5d0) then
          its = nits
        endif

c       Check convergence

        if((abs(dxi(1))+abs(dxi(2)).lt.tol .and. nrmdlm.lt.lamtol)
     &                                     .or.its.ge.nits ) then
          noconv = .false.
        endif

      end do ! while

      end
