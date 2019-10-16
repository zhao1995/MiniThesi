c$Id:$
      subroutine bernshp(ss,nn,xl,ndm, xx)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Generate nodes using Bernstein polynomials

c      Inputs:
c         ss        - interpolation point
c         nn        - Number of control points
c         xl(ndm,*) - Control point coordinates
c         ndm       - Mesh dimension

c      Outputs:
c         xx(ndm)   - Node coordinates
c-----[--.----+----.----+----.-----------------------------------------]
      implicit   none

      integer    nn,ndm, no,n, i
      real*8     ss, xl(ndm,*), xx(*)
      real*8     fm, fp, fact(11), bi(11), bern(11)

      call binomial(nn, fact, bi)

c     Generate the Bernstein functions

      fm = 0.5d0 - 0.5d0*ss
      fp = 0.5d0 + 0.5d0*ss
      no = nn - 1
      bern( 1) = fm**no
      bern(nn) = fp**no
      do n = 1,no-1
        bern(n+1) = bi(n+1)*fm**(no-n)*fp**n
      end do ! n

      do i = 1,ndm
        xx(i) = 0.0d0
        do n = 1,nn
          xx(i) = xx(i) + bern(n)*xl(i,n)
        end do ! n
      end do ! i

      end
