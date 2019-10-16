c$Id:$
      subroutine stcon3f(h,g21,dui,ndf,ndfi,nel,nst,nc,s,p)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: 3-D finite deformation enhanced strain element

c      Inputs:

c      Outputs:
c-----[--.----+----.----+----.-----------------------------------------]
      implicit none

      integer ndf,ndfi,nel,nst,nc,ns, i,j, i1,i2,j1,j2,  ii
      real*8  h(12,12), g21(12,24), s(nst,*), p(ndf,*)
      real*8  dui(12),tt(12,24), dot

c     Perform static condensation of internal modes

      ns = 8*ndfi

c     Reduce load vector: p - K12 * K22-inverse * b
c               (K12 = K21^T)
      i2 = 0
      do i = 1,nel
        do j = 1,ndfi
          p(j,i) = p(j,i) - dot(g21(1,i2+j),dui,nc)
        end do ! j
        i2 = i2 + ndfi
      end do ! i

c     Compute K12 * K22-inverse (store transpose)

      do i = 1,ns
        do j = 1,nc
          tt(j,i) = dot(g21(1,i),h(1,j),nc)
        end do ! j
      end do ! i

c     Compute K12 * K22-inverse * K21

      i1 = 0
      i2 = 0
      do i = 1,nel
        j1 = 0
        j2 = 0
        do j = 1,i
          do ii = 1,ndfi
            s(i1+ii,j1+1) = s(i1+ii,j1+1)
     &                    - dot(tt(1,i2+ii),g21(1,j2+1),nc)
            s(i1+ii,j1+2) = s(i1+ii,j1+2)
     &                    - dot(tt(1,i2+ii),g21(1,j2+2),nc)
            s(i1+ii,j1+3) = s(i1+ii,j1+3)
     &                    - dot(tt(1,i2+ii),g21(1,j2+3),nc)
          end do ! ii
          j1 = j1 + ndf
          j2 = j2 + ndfi
        end do ! j
        i1 = i1 + ndf
        i2 = i2 + ndfi
      end do ! i

      end
