c$Id:$
      subroutine tmat7(rotv,rota, tl,h)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose:  Compute T and H matrices for rigid body tangent.

c      Inputs:
c         rotv(3)  - Rotational velocity
c         rota(3)  - Rotational acceleration

c      Outputs:
c         tl(3,3)  - T matrix
c         h(3,3)   - H matrix
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      integer   i, j
      real*8    norm, cons1, cons2, consm, t, tt, nn, h1,h2, t1,t2
      real*8    rotv(3), rota(3), h(3,3), tl(3,3), dyad(3,3)

      save

      data      h1,h2/ 1.6666666666667d-01, 8.3333333333333d-03/
      data      t1,t2/ 3.3333333333333d-01, 1.3333333333333d-01/

c     Compute T matrix

      call dyadic(rota,norm,dyad)

      if (norm.lt.1.d-15) then

        do j=1,3
          do i=1,3
            tl(i,j) = 0.0d0
          end do ! i
          tl(j,j) = 1.0d0
        end do ! j

      else

        t = norm*0.5d0
        if (norm.gt.1.d-03) then
          cons1 =t/tan(t)
        else
          tt    = t*t
          cons1 = 1.d0/(1.d0+tt*(t1 + tt*t2))
        endif

        consm = 1.d0 - cons1

        do j=1,3
          do i=1,3
            tl(i,j) = consm*dyad(i,j)
          end do ! i
        end do ! j

        tl(1,1) = tl(1,1) + cons1
        tl(1,2) = tl(1,2) + 0.5d0*rota(3)
        tl(1,3) = tl(1,3) - 0.5d0*rota(2)
        tl(2,1) = tl(2,1) - 0.5d0*rota(3)
        tl(2,2) = tl(2,2) + cons1
        tl(2,3) = tl(2,3) + 0.5d0*rota(1)
        tl(3,1) = tl(3,1) + 0.5d0*rota(2)
        tl(3,2) = tl(3,2) - 0.5d0*rota(1)
        tl(3,3) = tl(3,3) + cons1

      endif

c     Compute H matrix

      call dyadic(rotv,norm,dyad)

      t = norm*0.5d0

      if (norm.lt.0.001) then
        nn = norm**2
        tt = t**2
        cons1 = 1.d0 - nn*(h1 - nn*h2)
        cons2 = 1.d0 - tt*(h1 - tt*h2)
      else
        cons1 = sin(norm)/norm
        cons2 = 0.5d0*(sin(t)/t)**2
      endif

      consm = 1.d0 - cons1

      do j=1,3
        do i=1,3
          h(i,j) = consm*dyad(i,j)
        end do ! i
      end do ! j

      h(1,1) = h(1,1) + cons1
      h(1,2) = h(1,2) - cons2*rotv(3)
      h(1,3) = h(1,3) + cons2*rotv(2)
      h(2,1) = h(2,1) + cons2*rotv(3)
      h(2,2) = h(2,2) + cons1
      h(2,3) = h(2,3) - cons2*rotv(1)
      h(3,1) = h(3,1) - cons2*rotv(2)
      h(3,2) = h(3,2) + cons2*rotv(1)
      h(3,3) = h(3,3) + cons1

      end
