c$Id:$
      subroutine pfolel(xl,ul,fld,p,s,prop,ndm,ndf,nst)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c     Purpose:

c     Inputs:

c     Outputs:
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      integer   ndm,ndf,nst, i,j,j1
      real*8    xnorm,force,fld,prop, xg(3,2),nn(3)
      real*8    xl(ndm,2),ul(ndf,2),p(*),s(nst,*)

      save

      do i = 1,2
        do j = 1,ndm
          xg(j,i) = xl(j,i) + ul(j,i)
        end do ! j
      end do ! i

c     Compute unit normal vector

      xnorm = 0.0d0
      do j = 1,ndm
        nn(j) = xg(j,1) - xg(j,2)
        xnorm = xnorm + nn(j)**2
      end do ! j
      xnorm = 1.d0/sqrt(xnorm)
      do j = 1,ndm
        nn(j) = nn(j)*xnorm
      end do ! j

      force = fld*prop

c     Compute residual/tangent

      xnorm = xnorm*force

      do j = 1,nst
        do i = 1,nst
          s(i,j) = 0.0d0
        end do ! i
        p(j) = 0.0d0
      end do ! j

      do j = 1,ndm
        p(j) = - force*nn(j)
        j1 = j + ndf
        do i = 1,ndm
          s(i,j ) = -nn(i)*nn(j)*xnorm
          s(i,j1) = -s(i,j)
        end do ! i
        s(j,j)  = s(j,j)  + xnorm
        s(j,j1) = s(j,j1) - xnorm
      end do ! j

      end
