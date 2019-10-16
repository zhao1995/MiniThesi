c$Id:$
      subroutine pbm3sc(mx,xbm,sv,wtb, yy,zz,sig,ii,nqudr)

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

      include 'bm3cn.h'
      include 'bm3fd.h'
      include 'cdata.h'
      include 'eldata.h'
      include 'pview.h'

      integer  i,m,m1,m2,nf,ii,nqudr, mx(5,nqudr,numel)
      real*8   xbm(3,nqudr,numnp),sv(nqudr,numnp),wtb(numnp)
      real*8   dx(3), yy,zz, sig

      save

c     Deformed update with current rotations

      if(cs.gt.0.0d0) then
        do i = 1,3
          dx(i) = rot1(i,2)*yy + rot1(i,3)*zz
        end do ! i

c     Undeformed update with inital rotations

      else
        do i = 1,3
          dx(i) = rot0(i,2)*yy + rot0(i,3)*zz
        end do ! i
      endif

c     Set nodal quantities

      do m = 1,mel
        if(ii.eq.1) then
          wtb(mxl(m)) = wtb(mxl(m)) + 1.d0
        endif
        do i = 1,3
          xbm(i,ii,mxl(m)) = xbm(i,ii,mxl(m)) + xll(i,m) + dx(i)
        end do ! i
        sv(ii,mxl(m)) = sv(ii,mxl(m)) + sig
      end do ! m

c     Set element quantities

      if (ii.eq.1) then
        m1 = nqudr*mxl(1) - nqudr
        m2 = nqudr*mxl(2) - nqudr
        nf = 1
        mx(1,nf,n) = m1 + 1
        mx(2,nf,n) = m2 + 1
        mx(3,nf,n) = m2 + nqudr
        mx(4,nf,n) = m1 + nqudr
        mx(5,nf,n) = mxl(4)
        do nf = 2,nqudr
          mx(1,nf,n) = m1 + nf
          mx(2,nf,n) = m2 + nf
          mx(3,nf,n) = m2 + nf - 1
          mx(4,nf,n) = m1 + nf - 1
          mx(5,nf,n) = mxl(4)
        end do ! nf
      endif ! ii test

      end
