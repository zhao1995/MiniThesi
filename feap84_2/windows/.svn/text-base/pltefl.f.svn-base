c$Id:$
      subroutine pltefl (nel,ic,x,v,vc,nc)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Dimension v(*)                                   01/05/2012
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Element fill routine for contour panel plots

c      Inputs:
c         nel       - Number of nodes/element
c         ic(*)     - Contour value at nodes
c         x(3,*)    - Nodal coordinates of element
c         v(*)      - Nodal values of contour
c         nc        - Number of values for contours

c      Outputs:
c         none      - Plot outputs to screen/file
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'prmptd.h'

      integer   i, j, k, icol, icx, icn, nel, nc, ic(*)
      real*8    vci, s, xp(10),yp(10),zp(10),x(3,*),v(*),vc(*)

      save

c     color pallet

      integer   ipal(12)
      data      ipal  / 13, 4,12, 6,14, 3,10, 5,11, 8, 9, 2/
      integer   rpal(12)
      data      rpal  /  2, 9, 8,11, 5,10, 3,14, 6,12, 4,13/

      vc(nc+1) = max ( vc(1),vc(nc) )
      vci = 0.0d0
      icx = ic(1)
      icn = ic(1)
      do i = 1,nel
        vci      = max(vci,abs(v(i)))
        vc(nc+1) = max(vc(nc+1),v(i))
        icx      = max(ic(i),icx)
        icn      = min(ic(i),icn)
      end do ! i
      vc(nc+1)   = vc(nc+1)*1.001 + vci + 1.0d-8
      do icol = icn,icx
        if(psrevs) then
          i = rpal(icol)
        else
          i = ipal(icol)
        endif
        call pppcol(i,2)
        k = 0
        i = nel
        do j = 1,nel
          if((ic(j).ge.icol .and. ic(i).le.icol)
     &                      .and. (ic(i).ne.ic(j))) then
            if(icol-1.ge.ic(i)) then
              s = (vc(icol-1)-v(i))/(v(j)-v(i))
              k = k + 1
              xp(k) = x(1,i) + (x(1,j)-x(1,i))*s
              yp(k) = x(2,i) + (x(2,j)-x(2,i))*s
              zp(k) = x(3,i) + (x(3,j)-x(3,i))*s
            endif
            s = (vc(icol)-v(i))/(v(j)-v(i))
            if(s.le.1.0d0) then
              k = k + 1
              xp(k) = x(1,i) + (x(1,j)-x(1,i))*s
              yp(k) = x(2,i) + (x(2,j)-x(2,i))*s
              zp(k) = x(3,i) + (x(3,j)-x(3,i))*s
            endif
          elseif((ic(i).ge.icol .and. ic(j).le.icol)
     &                          .and. (ic(i).ne.ic(j))) then
            s = (vc(icol)-v(i))/(v(j)-v(i))
            if(s.ge.0.0d0) then
              k = k + 1
              xp(k) = x(1,i) + (x(1,j)-x(1,i))*s
              yp(k) = x(2,i) + (x(2,j)-x(2,i))*s
              zp(k) = x(3,i) + (x(3,j)-x(3,i))*s
            endif
            if(icol-1.ge.ic(j)) then
              s = (vc(icol-1)-v(i))/(v(j)-v(i))
              k = k + 1
              xp(k) = x(1,i) + (x(1,j)-x(1,i))*s
              yp(k) = x(2,i) + (x(2,j)-x(2,i))*s
              zp(k) = x(3,i) + (x(3,j)-x(3,i))*s
            endif
          endif
          if(ic(j).eq.icol) then
            k = k + 1
            xp(k) = x(1,j)
            yp(k) = x(2,j)
            zp(k) = x(3,j)
          endif

c         Check for duplicate points

          if(k.gt.1) then
            if(xp(k).eq.xp(k-1) .and.
     &         yp(k).eq.yp(k-1) .and.
     &         zp(k).eq.zp(k-1) ) then
              k = k - 1
            endif
          endif
          i = j
        end do ! j

c       Plot panel of this color

        call plotl(xp(1),yp(1),zp(1),1)
        do j = 2,k
          call plotl(xp(j),yp(j),zp(j),2)
        end do ! j
        call clpan
      end do ! icol

      end
