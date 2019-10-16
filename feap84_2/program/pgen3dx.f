c$Id:$
      subroutine pgen3dx(x,ndtyp, x2,x0,xp, numn2,nsty,nseg,ndir,nn,prt)

      implicit   none

      include   'bdata.h'
      include   'iofile.h'
      include   'sdata.h'
      include   'pconstant.h'
      include   'trdata.h'

      logical    prt
      integer    i, numn2,nsty,nseg, ndir, nn, n,ns, ndtyp(*)

      real*8     x2(2,numn2),x0(3),xp(3,*), xl(3), x(ndm,*)
      real*8     dpi

      if(prt) then
        write(iow,2000) head
      endif
      dpi = pi/180.0d0
      nn  = 0
      do ns = 1,nseg+1
        do n = 1,numn2
          if(nsty.eq.1) then
            xl(1) = x0(1) + xp(1,ns) + x2(1,n)
            xl(2) = x0(2) + xp(2,ns) + x2(2,n)
            xl(3) = x0(3) + xp(3,ns)
          else
            if(ndir.eq.1) then
              xl(1) = x0(1) +  xp(3,ns) + x2(1,n)
              xl(2) = x0(2) + (xp(1,ns) + x2(2,n))*cos(dpi*xp(2,ns))
              xl(3) = x0(3) + (xp(1,ns) + x2(2,n))*sin(dpi*xp(2,ns))
            else
              xl(1) = x0(1) + (xp(1,ns) + x2(1,n))*cos(dpi*xp(2,ns))
              xl(2) = x0(2) - (xp(1,ns) + x2(1,n))*sin(dpi*xp(2,ns))
              xl(3) = x0(3) +  xp(3,ns) + x2(2,n)
            endif
          endif

c         Transform to global coordinates

          nn = nn + 1
          do i = 1,3
            x(i,nn) = xr(i) + tr(i,1)*xl(1)
     &                      + tr(i,2)*xl(2)
     &                      + tr(i,3)*xl(3)
          end do ! i

c         Set node to active

          ndtyp(nn) = 0

c         Output to file

          if(prt) then
            write(iow,2001) nn,(x(i,nn),i=1,3)
          endif
        end do ! n
      end do ! ns

c     Formats

2000  format(/1x,20a4//'   N o d a l   C o o r d i n a t e s'//
     &      '      Node     1-Coord     2-Coord     3-Coord')

2001  format(i10,1p,3e12.4)

      end
