c$Id:$
      subroutine setclp(x,ndm,numnp)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Compute clip range from mesh data

c      Inputs:
c         x(ndm,*)  - Nodal coordinates of mesh
c         ndm       - Spatial dimension of region
c         numnp     - Number of nodes in mesh

c      Outputs:
c         none      - Output through common /plclip/
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'plclip.h'
      include  'pointer.h'
      include  'comblk.h'

      include  'p_int.h'

      logical   flag
      integer   i,n,ndm,numnp
      real*8    x(ndm,*)

      save

      flag  = .true.
      fp(1) =  npty - 1
      do n = 1,numnp
        if(mr(fp(1)+n).ge.0) then

c         Set defaults to an existing coordinate

          if( flag ) then
            do i = 1,ndm
              cmin(i) = x(i,n)
              cmax(i) = x(i,n)
            end do ! i
            flag = .false.
          end if

c         Set min / max for other coordinates

          do i = 1,ndm
            cmin(i) = min( cmin(i), x(i,n))
            cmax(i) = max( cmax(i), x(i,n))
          end do ! i
        end if
      end do ! n

      end
