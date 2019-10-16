c$Id:$
      integer function nbuck(x,xmd,xmn,ndm,nsize)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Determine bucket number for bucket search routine

c      Inputs:
c         x(3)     - Coordinate to determine bucket number
c         xmd(3)   - Size of bucket
c         xmn      - Minimum coordinate of bucket
c         ndm      - Spatial dimension of mesh
c         nsize(4) - Number of buckets/direction: (4) indicates
c                    if all points to be placed in one bucket.

c      Outputs:
c         nbuck    - Bucket number containing data point x
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      integer   ndm,nn,ni,i, nsize(4)
      real*8    x(3),xmd(3),xmn(3)

      save

      if(nsize(4).gt.1) then
        nn = nint((x(1) - xmn(1))/xmd(1))
        nn = max(0,min(nsize(1)-1,nn))
        do i = 2,ndm
          ni = nint((x(i) - xmn(i))/xmd(i))
          nn = nsize(i)*nn + max(0,min(nsize(i)-1,ni))
        end do ! i
        nbuck = nn + 1
      else
        nbuck = 1
      endif

      end
