c$Id:$
      subroutine upeul(ia,eang,ul,ndf,isw)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Update for euler angle 3-d solutions

c      Inputs:

c      Outputs:
c-----[--.----+----.----+----.-----------------------------------------]
      implicit   none

      logical    etrans
      integer    i,j, ia(3), ndf,isw
      real*8     eang(3),ul(ndf,*), sn(3),cs(3), t(3,3)

c     Loop over nodes

      if(ia(1).ne.0 .and. ia(2).ne.0 .and. ia(3).ne.0) then
        etrans = .false.
        do j = 1,3
          if(eang(j).ne.0.0d0) then
            call pdegree(eang(j), sn(j),cs(j))
            etrans = .true.
          else
            sn(j) = 0.0d0
            cs(j) = 1.0d0
          endif
        end do ! j

c       Perform transformations for non-zero rotations

        if(etrans) then
          t(1,1) =  cs(1)*cs(3) - sn(1)*sn(2)*sn(3)
          t(1,2) = -sn(1)*cs(2)
          t(1,3) =  cs(1)*sn(3) + sn(1)*sn(2)*cs(3)
          t(2,1) =  sn(1)*cs(3) + cs(1)*sn(2)*sn(3)
          t(2,2) =  cs(1)*cs(2)
          t(2,3) =  sn(1)*sn(3) - cs(1)*sn(2)*cs(3)
          t(3,1) = -cs(2)*sn(3)
          t(3,2) =  sn(2)
          t(3,3) =  cs(2)*cs(3)

c         Rotate displacements to cartesian coordinates

          if(isw.eq.1) then
            do j = 1,3
              cs(j) = t(1,j)*ul(ia(1),1)
     &              + t(2,j)*ul(ia(2),1)
     &              + t(3,j)*ul(ia(3),1)
            end do ! j
            do j = 1,3
              ul(ia(j),1) = cs(j)
            end do !j
          else
            do i = 1,3
              do j = 1,3
                cs(j) = t(j,1)*ul(ia(1),i)
     &                + t(j,2)*ul(ia(2),i)
     &                + t(j,3)*ul(ia(3),i)
              end do ! j
              do j = 1,3
                ul(ia(j),i) = cs(j)
              end do !j
            end do ! i
          endif
        endif
      endif

      end
