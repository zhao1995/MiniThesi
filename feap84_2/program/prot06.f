c$Id:$
      subroutine prot06 (tn,ta,t1,du,ua,vn,v1,an,a1,dm,isw)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Prototype rotational update routine

c      Inputs:
c            tn = lambda_n
c            vn = omega_n
c            an = dot-omega_n
c            dm =
c            isw = switch: 0=initialize; 1 = start step; 2 = update.

c      Outputs:
c            t1 = lambda_n+1
c            v1 = omega_n+1
c            a1 = dot-omega_n+1
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'tdata.h'
      include  'tdato.h'

      integer   isw, i,j
      real*8    th,th2,th4,th6, sinth,costh,fac3,fac30, dtbar
      real*8    tn(3,3) , ta(3,3) , t1(3,3) , dl(3,3), dm(3)
      real*8    du(3)   , ua(3)   , vn(3)   , v1(3)   , an(3) , a1(3)

      save

c     Initialize rotational matrix to identity at t_0

      if(isw.eq.0) then

        do j = 1,3
          do i = 1,3
            tn(i,j) = 0.d0
            t1(i,j) = 0.d0
            ta(i,j) = 0.d0
          end do ! i
          tn(j,j) = 1.d0
          t1(j,j) = 1.d0
          ta(j,j) = 1.d0
          du(j)   = 0.d0
          ua(j)   = 0.d0
          vn(j)   = 0.d0
          v1(j)   = 0.d0
          an(j)   = 0.d0
          a1(j)   = 0.d0
        end do ! j

c     At start of new time step

      elseif(isw.eq.1) then

c       Mean time step

        dtbar = 0.5d0*(dt + dtold)

c       Compute velocity at t_n+1/2 and solution state at t_n+1

        do j = 1,3
          v1(j) = vn(j) + dtbar*an(j)
          a1(j) = 0.0d0
          ua(j) = dt*v1(j)
        end do ! j

c       Compute incremental rotation matrix

        th2   = ua(1)**2 + ua(2)**2 + ua(3)**2
        th    = sqrt(th2)
        costh = cos(th)
        if(abs(th).gt.1.0d-02) then
          sinth = sin(th)/th
          fac30 = (1.d0 - costh)/(th*th)
        else
          th4   = th2*th2
          th6   = th2*th4
          sinth = 1.0d0 - th2/6.d0  + th4/120.d0 - th6/5040.d0
          fac30 = 0.5d0 - th2/24.d0 + th4/720.d0 + th6/40320.d0
        endif
        do j = 1,3
          fac3 = fac30*ua(j)
          do i = 1,3
            dl(i,j) = fac3*ua(i)
          end do ! i
          dl(j,j) = dl(j,j) + costh
        end do ! j
        dl(1,2) = dl(1,2) - sinth*ua(3)
        dl(1,3) = dl(1,3) + sinth*ua(2)
        dl(2,1) = dl(2,1) + sinth*ua(3)
        dl(2,3) = dl(2,3) - sinth*ua(1)
        dl(3,1) = dl(3,1) - sinth*ua(2)
        dl(3,2) = dl(3,2) + sinth*ua(1)

c       Multiply increment to form state at t_n+1

        do j = 1,3
          do i = 1,3
            t1(i,j) = dl(i,1)*tn(1,j)
     &              + dl(i,2)*tn(2,j)
     &              + dl(i,3)*tn(3,j)
          end do ! i
        end do ! j

c     Update after solution for acceleration

      elseif(isw.eq.2) then

c       Translational and rotational accelerations updates

        do i = 1,3
          a1(i) = a1(i) + dm(i)
        end do ! i

      endif

      end
