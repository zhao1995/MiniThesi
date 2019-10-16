c$Id:$
      subroutine pltlfl(nel,xl,v,vc,nc)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Plot of mesh contours: With inter-element smoothing

c      Inputs:
c         nel      - Number of nodes on segment
c         xl(3,*)  - Nodal coordinates of line segment
c         v(*)     - Nodal solution values
c         vc(*)    - Contour values
c         nc       - Number of contour intervals

c      Outputs:
c         none      - Plot outputs to screen/file
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'prmptd.h'

      integer   nel,nc,icol,icx,icn,i,j,k, ilc(9)
      real*8    s,vci, xl(3,*),v(*),vc(*) ,xp(9),yp(9),zp(9)

      save

c     Color pallet

      integer   ipal(12)
      data      ipal  / 13, 4,12, 6,14, 3,10, 5,11, 8, 9, 2/
      integer   rpal(12)
      data      rpal  /  2, 9, 8,11, 5,10, 3,14, 6,12, 4,13/

      call pltcor(nel,ilc,v,vc,nc)

      vc(nc+1) = max ( vc(1),vc(nc) )
      vci = 0.0d0
      icx = ilc(1)
      icn = ilc(1)
      do i = 1,nel
        vci      = max(vci,abs(v(i)))
        vc(nc+1) = max(vc(nc+1),v(i))
        icx      = max(ilc(i),icx)
        icn      = min(ilc(i),icn)
      end do ! i
      vc(nc+1)   = vc(nc+1)*1.001 + vci + 1.0d-8
      do icol = icn,icx
        if(psrevs) then
          i = rpal(icol)
        else
          i = ipal(icol)
        endif
        call pppcol(i,0)
        k = 0
        i = nel
        do j = 1,nel
          if((ilc(j).ge.icol .and. ilc(i).le.icol)
     &                      .and. (ilc(i).ne.ilc(j))) then
            if(icol-1.ge.ilc(i)) then
              s = (vc(icol-1)-v(i))/(v(j)-v(i))
              k = k + 1
              xp(k) = xl(1,i) + (xl(1,j)-xl(1,i))*s
              yp(k) = xl(2,i) + (xl(2,j)-xl(2,i))*s
              zp(k) = xl(3,i) + (xl(3,j)-xl(3,i))*s
            endif
            s = (vc(icol)-v(i))/(v(j)-v(i))
            if(s.le.1.0d0) then
              k = k + 1
              xp(k) = xl(1,i) + (xl(1,j)-xl(1,i))*s
              yp(k) = xl(2,i) + (xl(2,j)-xl(2,i))*s
              zp(k) = xl(3,i) + (xl(3,j)-xl(3,i))*s
            endif
          elseif((ilc(i).ge.icol .and. ilc(j).le.icol)
     &                          .and. (ilc(i).ne.ilc(j))) then
            s = (vc(icol)-v(i))/(v(j)-v(i))
            if(s.ge.0.0d0) then
              k = k + 1
              xp(k) = xl(1,i) + (xl(1,j)-xl(1,i))*s
              yp(k) = xl(2,i) + (xl(2,j)-xl(2,i))*s
              zp(k) = xl(3,i) + (xl(3,j)-xl(3,i))*s
            endif
            if(icol-1.ge.ilc(j)) then
              s = (vc(icol-1)-v(i))/(v(j)-v(i))
              k = k + 1
              xp(k) = xl(1,i) + (xl(1,j)-xl(1,i))*s
              yp(k) = xl(2,i) + (xl(2,j)-xl(2,i))*s
              zp(k) = xl(3,i) + (xl(3,j)-xl(3,i))*s
            endif
          endif
          if(ilc(j).eq.icol) then
            k = k + 1
            xp(k) = xl(1,j)
            yp(k) = xl(2,j)
            zp(k) = xl(3,j)
          endif

c         Check for duplicate points

          if(k.gt.1) then
            if(xp(k).eq.xp(k-1) .and.
     &         yp(k).eq.yp(k-1) .and.
     &         zp(k).eq.zp(k-1)) then
              k = k - 1
            endif
          endif
          i = j
        end do ! j

c       Plot line of this color

        call plotl(xp(1),yp(1),zp(1),3)
        do j = 2,k
          call plotl(xp(j),yp(j),zp(j),2)
        end do ! j
        call pppcol(0,0)
      end do ! icol

      end
