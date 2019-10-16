c$Id:$
      subroutine uptriad(ia,t,ul,ndf,isw)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Update for triad  angle 3-d solutions

c      Inputs:

c      Outputs:
c-----[--.----+----.----+----.-----------------------------------------]
      implicit   none

      integer    i,j, ia(3), ndf,isw
      real*8     t(3,3),ul(ndf,*), cs(3)

c     Loop over nodes

      if(ia(1).ne.0 .and. ia(2).ne.0 .and. ia(3).ne.0) then

c       Rotate displacements to cartesian coordinates

        if(isw.eq.1) then
          do j = 1,3
            cs(j) = t(1,j)*ul(ia(1),1)
     &            + t(2,j)*ul(ia(2),1)
     &            + t(3,j)*ul(ia(3),1)
          end do ! j
          do j = 1,3
            ul(ia(j),1) = cs(j)
          end do !j
        else
          do i = 1,3
            do j = 1,3
              cs(j) = t(j,1)*ul(ia(1),i)
     &              + t(j,2)*ul(ia(2),i)
     &              + t(j,3)*ul(ia(3),i)
            end do ! j
            do j = 1,3
              ul(ia(j),i) = cs(j)
            end do !j
          end do ! i
        endif

      endif

      end
