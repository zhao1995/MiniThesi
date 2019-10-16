c$Id:$
      subroutine seigen(s,nn,nst, dr, smax)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    17/01/2008
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose:  Compute the maximum eigenvalue of symmetric s(nst,nst)

c      Inputs:
c         s(nst,*) - Element array
c         nn       - Size of upper block
c         nst      - Dimension of element array

c      Working:
c         dr(*)    - Temporary work space

c      Outputs:
c         smax     - Maximum eigenvalue of s
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      integer   nst, i, j, j1, nn, n1,n2,n3
      real*8    s(nst,nst), dr(*), smax

      save

c     Set pointers                  ! Eigenvalues in dr(1:nst)

      n1 =  1 +  nn               ! Triangular matrix of s
      n2 = n1 + (nn*(nn+1))/2     ! Eigenvalues of s
      n3 = n2 +  nn*nn            ! Working space

c     Load stiffness terms into triangular matrix

      j1 = -1
      do j = 1,nn
        do i = 1,j
          dr(n1+j1+i) = s(i,j)
        end do ! i
        j1 = j1 + j
      end do ! j

c     Compute eigenpairs for element

      call eisql(dr(n1),dr(1),dr(n3),dr(n2),nn,j1)

c     Compute maximum eigenvalue

      smax = 0.0d0
      do i = 1,nn
        smax = max(smax,abs(dr(i)))
      end do ! i

      end
