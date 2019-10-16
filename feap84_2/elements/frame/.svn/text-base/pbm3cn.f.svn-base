c$Id:$
      subroutine pbm3cn(mx,xbm,sv,wtb, yy,zz,sig,jj,nfc,nqy,nqz)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose:

c      Inputs:

c      Outputs:
c-----[--.----+----.----+----.-----------------------------------------]
      implicit none

      include 'cdata.h'
      include 'eldata.h'
      include 'pview.h'

      real*8         xll     ,rot0     ,rot1
      integer                                     mel,mxl
      common /bm3cn/ xll(3,3),rot0(3,3),rot1(3,3),mel,mxl(4)

      integer  i,m,m1,m2,nf,jj,nfc,nqy,nqz, mx(5,nfc,numel)
      real*8   dx(3), yy,zz, sig
      real*8   xbm(3,nfc,numnp),sv(nfc,numnp),wtb(numnp)

      save

c     If deformed update the rotations

      if(cs.gt.0.0d0) then
        do i = 1,3
          dx(i) =  rot1(i,2)*yy + rot1(i,3)*zz
        end do ! i
      else
        do i = 1,3
          dx(i) =  rot0(i,2)*yy + rot0(i,3)*zz
        end do ! i
      endif

c     Set nodal quantities

      do m = 1,mel
        if(jj.eq.1) then
          wtb(mxl(m)) = wtb(mxl(m)) + 1.d0
        endif
        do i = 1,3
          xbm(i,jj,mxl(m)) = xbm(i,jj,mxl(m)) + xll(i,m) + dx(i)
        end do ! i
        sv(jj,mxl(m)) = sv(jj,mxl(m)) + sig
      end do ! m

c     Set element quantities

      if (jj.eq.1) then
        nf = 0
        m1 = nfc*mxl(1) - nfc
        m2 = nfc*mxl(2) - nfc
        do m = 1,nqy-1
          nf = nf+1
          mx(1,nf,n) = m1 + m
          mx(2,nf,n) = m1 + m + 1
          mx(3,nf,n) = m2 + m + 1
          mx(4,nf,n) = m2 + m
          mx(5,nf,n) = mxl(4)
        end do ! m

        nf = nf + 1
        mx(1,nf,n) = m1 + nqy + 1
        mx(2,nf,n) = m1 + 1
        mx(3,nf,n) = m2 + 1
        mx(4,nf,n) = m2 + nqy + 1
        mx(5,nf,n) = mxl(4)

        if (nqz.gt.2) then
        nf = nf + 1
        mx(1,nf,n) = m1 + nqy
        mx(2,nf,n) = m1 + nqy + 2
        mx(3,nf,n) = m2 + nqy + 2
        mx(4,nf,n) = m2 + nqy
        mx(5,nf,n) = mxl(4)
        end if

        do m = 2,nqz-2
          nf = nf+1
          mx(1,nf,n) = m1 + nqy + 2*m - 1
          mx(2,nf,n) = m1 + nqy + 2*m - 3
          mx(3,nf,n) = m2 + nqy + 2*m - 3
          mx(4,nf,n) = m2 + nqy + 2*m - 1
          mx(5,nf,n) = mxl(4)
          nf = nf+1
          mx(1,nf,n) = m1 + nqy + 2*m - 2
          mx(2,nf,n) = m1 + nqy + 2*m
          mx(3,nf,n) = m2 + nqy + 2*m
          mx(4,nf,n) = m2 + nqy + 2*m - 2
          mx(5,nf,n) = mxl(4)
        end do ! m

        if (nqz.gt.2) then
          nf = nf+1
          mx(1,nf,n) = m1 + nqy + 2*nqz - 3
          mx(2,nf,n) = m1 + nqy + 2*nqz - 5
          mx(3,nf,n) = m2 + nqy + 2*nqz - 5
          mx(4,nf,n) = m2 + nqy + 2*nqz - 3
          mx(5,nf,n) = mxl(4)
        end if

          nf = nf+1
          mx(1,nf,n) = m1 + nqy + 2*nqz - 4
          mx(2,nf,n) = m1 + 2*(nqy + nqz) - 4
          mx(3,nf,n) = m2 + 2*(nqy + nqz) - 4
          mx(4,nf,n) = m2 + nqy + 2*nqz - 4
          mx(5,nf,n) = mxl(4)

        do m = 1,nqy-1
          nf = nf+1
          mx(1,nf,n) = m1 + nqy + 2*nqz + m - 3
          mx(2,nf,n) = m1 + nqy + 2*nqz + m - 4
          mx(3,nf,n) = m2 + nqy + 2*nqz + m - 4
          mx(4,nf,n) = m2 + nqy + 2*nqz + m - 3
          mx(5,nf,n) = mxl(4)
        end do ! m
      endif ! jj test

      end
