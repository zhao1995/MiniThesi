c$Id:$
      subroutine inerbx(mass,rcg,rlam,jj,ctan,nst, sr,pr,iomeg,omega)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c     Purpose: Compute central difference explicit rigid body equation

c     Inputs:

c     Outputs:
c-----[--.----+----.----+----.-----------------------------------------]
      implicit none

      integer  i,j,nst
      real*8   mass, rcg(3,11),rlam(3,3,6),ctan(3)
      real*8   jj(3,3),sr(nst,*),pr(6), jrlam(3,3), omega(3), iomeg(3)

      save

c     Compute complete inertia mass/dyadic in current configuration

      do j = 1,3
        do i = 1,3
          jrlam(i,j) = jj(i,1)*rlam(j,1,3)
     &               + jj(i,2)*rlam(j,2,3)
     &               + jj(i,3)*rlam(j,3,3)
        end do ! i
      end do ! j

      do i = 1,6
        do j = 1,6
          sr(j,i) = 0.0d0
        end do ! j
      end do ! i

      do j = 1,3
        sr(j,j) = mass
        do i = 1,3
          sr(i+3,j+3) = rlam(i,1,3)*jrlam(1,j)
     &                + rlam(i,2,3)*jrlam(2,j)
     &                + rlam(i,3,3)*jrlam(3,j)
        end do ! i
      end do ! j

c     Compute residual

      do i = 1,3
        omega(i) = rlam(i,2,5)
      end do ! i

      do i = 1,3
        pr(i) = - mass*rcg(i,6)
      end do ! i

c     Subtract Coriolis term to residual

      do i = 1,3
        iomeg(i) = sr(i+3,4)*omega(1)
     &           + sr(i+3,5)*omega(2)
     &           + sr(i+3,6)*omega(3)
      end do ! i


      pr(4) = - omega(2)*iomeg(3) + omega(3)*iomeg(2)
     &        - sr(4,4)*rlam(1,3,5)
     &        - sr(4,5)*rlam(2,3,5)
     &        - sr(4,6)*rlam(3,3,5)
      pr(5) = - omega(3)*iomeg(1) + omega(1)*iomeg(3)
     &        - sr(5,4)*rlam(1,3,5)
     &        - sr(5,5)*rlam(2,3,5)
     &        - sr(5,6)*rlam(3,3,5)
      pr(6) = - omega(1)*iomeg(2) + omega(2)*iomeg(1)
     &        - sr(6,4)*rlam(1,3,5)
     &        - sr(6,5)*rlam(2,3,5)
     &        - sr(6,6)*rlam(3,3,5)

      do i = 1,6
        do j = 1,6
          sr(j,i) = sr(j,i)*ctan(3)
        end do ! j
      end do ! i

      end
