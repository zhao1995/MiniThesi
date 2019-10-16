c$Id:$
      subroutine plotelm(ie,ix, x, num, axs)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Rewrite for current use of 'pstyp' specification 14/02/2007
c-----[--+---------+---------+---------+---------+---------+---------+-]
c     Purpose: Plot curved edges on individual elements

c     Inputs:
c       ie(nie,*) - Element type descriptions
c       ix(nen)   - List of nodes on element
c       x(ndm,*)  - Nodal coordinates of element
c       num       - Plot if true
c       axs       - Axes at node axs


c     Output:
c       Graphical plot of element edges
c-----[--+---------+---------+---------+---------+---------+---------+-]
      implicit   none

      include   'cdata.h'
      include   'cdat1.h'
      include   'pdata1.h'
      include   'pdata4.h'
      include   'ppers.h'
      include   'sdata.h'

      logical    num
      integer    ma,pstyp,iel, axs, nel, n,ns,i, iu,iutot
      integer    ie(nie,*), ix(nen), iplt(50), ixl(50)
      real*8     x(ndm,*), xl(3,50), xq(3), xx0(3)
      real*8     os0(2),osx(2),odx(2), oscale,oscaleg,ofact, sz
      real*8     oxmin(3),oxmax(3)

      save

c     Determine material number of element and plot type

      ma    = ix(nen1)
      pstyp = ie(1    ,ma)
      iel   = ie(nie-1,ma)

c     Determine number of nodes on element

      do n = 1,nen
        if(ix(n).ne.0) then
          nel = n
        endif
      end do ! n

      call plftyp(pstyp,nel,iel)
      call pltord(ix,iel, iu,iplt)

c     Set plot coordinates

      iutot = 0
      do n = 1,iu
        ns = iplt(n)
        if(ns.le.nel) then
          ixl(n) = ix(ns)
          if(ixl(n).gt.0) then
            iutot    = iutot + 1
            do i = 1,ndm
              xl(i,iutot) = x(i,ixl(n))
            end do ! i
            do i = ndm+1,3
              xl(i,iutot) = 0.0d0
            end do ! i
          endif
        endif
      end do ! n

c     Scale element

      do i = 1,ndm
        xx0(i) = xl(i,1)
      end do ! i
      sz = 0.0d0
      do n = 1,iutot
        do i = 1,ndm
          xl(i,n) = xl(i,n) - xx0(i)
          sz      = max(sz,abs(xl(i,n)))
        end do ! i
      end do ! n

c     Save old scaling

      oscale  = scale
      oscaleg = scaleg
      ofact   = fact
      do i = 1,2
        os0(i) = s0(i)
        odx(i) = dx(i)
        osx(i) = sx(i)
      end do ! i
      do i = 1,3
        oxmin(i) = xmin(i)
        oxmax(i) = xmax(i)
      end do ! i

c     Compute new scaling

      if(kpers.ne.0) then
        call frame(xl,3,iutot,-1)
      else
        call frame(xl,3,iutot, 1)
      endif

c     Line elements for edges

      xq(1) = xl(1,1)
      xq(2) = xl(2,1)
      xq(3) = xl(3,1)
      do i = 1,iutot-1

        call plotl(xq(1),xq(2),xq(3),3)

        xq(1) = xl(1,i+1)
        xq(2) = xl(2,i+1)
        xq(3) = xl(3,i+1)
        call plotl(xq(1),xq(2),xq(3),2)
      end do ! i

c     Add axes

      call pppcol (3,0)
      if(axs.gt.0) then
        call pltaxs(xl(1,axs),ndm,0.1d0*sz)
      endif

c     Add nodes and numbers

      call pppcol (5,0)
      call pltelnd(xl,ixl,3,iutot-1, num)

c     Check for any added nodes

      if(pstyp.eq.2) then
        if(nel.eq.9) then
          do i = 1,ndm
            xl(i,1) = x(i,ix(9)) - xx0(i)
          end do ! i
          do i = ndm+1,3
            xl(i,1) = 0.0d0
          end do ! i
          call pltelnd(xl,ix(9),3,1, num)
        elseif(nel.eq.16) then
          do n = 13,16
            ns = ix(n)
            if(ns.gt.0) then
              do i = 1,ndm
                xl(i,n) = x(i,ns) - xx0(i)
              end do ! i
              do i = ndm+1,3
                xl(i,n) = 0.0d0
              end do ! i
            endif
          end do ! n
          call pltelnd(xl(1,13),ix(13),3,4, num)
        endif
      elseif(pstyp.eq.3) then
        if(nel.eq.27) then
          do n = 21,27
            ns = ix(n)
            if(ns.gt.0) then
              do i = 1,ndm
                xl(i,n) = x(i,ns) - xx0(i)
              end do ! i
              do i = ndm+1,3
                xl(i,n) = 0.0d0
              end do ! i
            endif
          end do ! n
          call pltelnd(xl(1,21),ix(21),3,7, num)
        endif
      endif

c     Restore old scaling

      scale  = oscale
      scaleg = oscaleg
      fact   = ofact
      do i = 1,2
        s0(i) = os0(i)
        dx(i) = odx(i)
        sx(i) = osx(i)
      end do ! i
      do i = 1,3
        xmin(i) = oxmin(i)
        xmax(i) = oxmax(i)
      end do ! i

      end
