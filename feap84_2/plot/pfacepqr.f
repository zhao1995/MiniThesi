c$Id:$
        subroutine pfacepqr(np,iu,nfac)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    06/02/2009
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Set plot data for brick elements

c      Inputs:
c         np      - Order of element (2 for 27 node; 3 for 64 node; etc)

c      Output:
c         iu(4,*) - 4-node quadrilateral face
c         nfac    - Number of faces
c-----[--.----+----.----+----.-----------------------------------------]
      implicit   none

      include   'ublk1.h'

      integer    np, nfac
      integer    iu(4,*)
      integer    ne(3),ns(3), i,j, ii

      do i = 1,3
        ne(i) = np
        ns(i) = ne(i) + 1
      end do ! i

c     Negative 1-face

      nfac = 0
      do j = 1,ne(3)
        ii = ns(1)*ns(2)*(j-1)
        do i = 1,ne(2)
          nfac       = nfac + 1
          iu(1,nfac) = ns(1)*(i-1) + ii + 1
          iu(2,nfac) = iu(1,nfac) + ns(1)*ns(2)
          iu(3,nfac) = iu(2,nfac) + ns(1)
          iu(4,nfac) = iu(1,nfac) + ns(1)
        end do ! i
      end do ! j

c     Positive 1-face

      do j = 1,ne(3)
        ii = ns(1)*ns(2)*(j-1) + ne(1)
        do i = 1,ne(2)
          nfac       = nfac + 1
          iu(1,nfac) = ns(1)*(i-1) + ii + 1
          iu(2,nfac) = iu(1,nfac) + ns(1)
          iu(4,nfac) = iu(1,nfac) + ns(1)*ns(2)
          iu(3,nfac) = iu(4,nfac) + ns(1)
        end do ! i
      end do ! j

c     Negative 2-face

      do j = 1,ne(3)
        ii = ns(1)*ns(2)*(j-1)
        do i = 1,ne(1)
          nfac       = nfac + 1
          iu(1,nfac) = ii + i
          iu(2,nfac) = iu(1,nfac) + 1
          iu(4,nfac) = iu(1,nfac) + ns(1)*ns(2)
          iu(3,nfac) = iu(4,nfac) + 1
        end do ! i
      end do ! j

c     Positive 2-face

      do j = 1,ne(3)
        ii = ns(1)*ns(2)*j - ns(1)
        do i = 1,ne(1)
          nfac       = nfac + 1
          iu(1,nfac) = ii + i
          iu(2,nfac) = iu(1,nfac) + ns(1)*ns(2)
          iu(3,nfac) = iu(2,nfac) + 1
          iu(4,nfac) = iu(1,nfac) + 1
        end do ! i
      end do ! j

c     Negative 3-face

      do j = 1,ne(2)
        ii = ns(1)*(j-1)
        do i = 1,ne(1)
          nfac       = nfac + 1
          iu(1,nfac) = ii + i
          iu(2,nfac) = iu(1,nfac) + ns(1)
          iu(3,nfac) = iu(2,nfac) + 1
          iu(4,nfac) = iu(1,nfac) + 1
        end do ! i
      end do ! j

c     Positive 3-face

      do j = 1,ne(2)
        ii = ns(1)*(j-1) + ns(1)*ns(2)*(ns(3)-1)
        do i = 1,ne(1)
          nfac       = nfac + 1
          iu(1,nfac) = ii + i
          iu(2,nfac) = iu(1,nfac) + 1
          iu(4,nfac) = iu(1,nfac) + ns(1)
          iu(3,nfac) = iu(4,nfac) + 1
        end do ! i
      end do ! j

      end
