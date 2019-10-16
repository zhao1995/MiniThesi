c$Id:$
      subroutine jslide(lam,rlam,ebig,alp1,x1,x2h,x2,
     &                  r1h,r1,r2h,r2,n1,n2,pr,sr)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Form SLIDER joint arrays for connecting rigid bodies.

c      Inputs:
c         lam(4)    - Value of lagrange multiplier force on joint
c         rlam(3,*) - Rotation array for rigid bodies
c         ebig(3,3) - Orientation of slider axis
c         alp1      - alpha value for t_n+alpha
c         x1(3)     - Location of body 1 connection point at t_n+1
c         x2h(3)    - Location of body 2 connection point at t_n+alpha
c         x2(3)     - Location of body 2 connection point at t_n+1
c         r1h(3)    - Distance of body 1 mass center to x1 at t_n+alpha
c         r1(3)     - Distance of body 1 mass center to x1 at t_n+1
c         r2h(3)    - Distance of body 2 mass center to x2 at t_n+alpha
c         r2(3)     - Distance from body 2 mass center to x2 at t_n+1
c         n1        - Number of rigid body 1
c         n2        - Number of rigid body 2

c      Outputs:
c         pr(*)    - Residual array for joint
c         sr(*,*)  - Tangent  array for joint
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      integer   i,j,k,l,n1,n2
      real*8    alp1, alpn, lam(4), ebig(3,3), rlam(9,6,*)
      real*8    p(3,6), s(17,17), pr(3,7), sr(21,21)
      real*8    ean(3,3),ea1(3,3),eah(3,3),ebn(3,2),eb1(3,2),ebh(3,2)
      real*8    x1(3),x2h(3),x2(3),r1h(3),r1(3),r2h(3),r2(3)

      save

c     Form updated basis vectors at t_n and at t_n+1

      do i = 1,2
        call quavec(rlam(1,1,n2),ebig(1,i),ebn(1,i))
        call quavec(rlam(1,3,n2),ebig(1,i),eb1(1,i))
      end do ! i

      do i = 1,3
        call quavec(rlam(1,1,n1),ebig(1,i),ean(1,i))
        call quavec(rlam(1,3,n1),ebig(1,i),ea1(1,i))
      end do ! i

c     Interpolate the base vectors to t_n+a
      alpn = 1.d0 - alp1
      do i = 1,3
        do j = 1,2
          eah(i,j) = alpn*ean(i,j) + alp1*ea1(i,j)
          ebh(i,j) = alpn*ebn(i,j) + alp1*eb1(i,j)
        end do ! j
        eah(i,3) = alpn*ean(i,3) + alp1*ea1(i,3)
      end do ! i

      do i = 1,2

        call jtrots(ea1(1,3),eb1(1,i),eah(1,3),ebh(1,i),lam(i),p,s)
        do j = 1,3
          do k = 1,3
            sr(j+3,k+3) = sr(j+3,k+3) + s(j,k)
            sr(j+3,k+9) = sr(j+3,k+9) + s(j,k+3)
            sr(j+9,k+3) = sr(j+9,k+3) + s(j+3,k)
            sr(j+9,k+9) = sr(j+9,k+9) + s(j+3,k+3)
          end do ! k
          sr(j+3,12+i) = sr(j+3,12+i) + s(j,7)
          sr(j+9,12+i) = sr(j+9,12+i) + s(j+3,7)
          sr(12+i,j+3) = sr(12+i,j+3) + s(7,j)
          sr(12+i,j+9) = sr(12+i,j+9) + s(7,j+3)
        end do ! j
        do j = 1,3
          pr(j,2) = pr(j,2) + p(j,1)
          pr(j,4) = pr(j,4) + p(j,2)
        end do ! j
        pr(i,5) = pr(i,5) + p(1,3)

        l = i + 2
        call jtmix(ea1(1,i),eah(1,i),x1,x2,r1,r2,x2h,r1h,r2h,lam(l),p,s)
        do j = 1,12
          do k = 1,12
            sr(j,k) = sr(j,k) + s(j,k)
          end do ! k
          sr(j,12+l) = sr(j,12+l) + s(j,13)
          sr(12+l,j) = sr(12+l,j) + s(13,j)
        end do ! j
        do j = 1,3
          do k = 1,4
            pr(j,k) = pr(j,k) + p(j,k)
          end do ! k
        end do ! j
        if(i.eq.1) then
          pr(l,5) = pr(l,5) + p(1,5)
        elseif(i.eq.2) then
          pr(1,6) = pr(1,6) + p(1,5)
        endif

      end do ! i

      end
