      subroutine  dswap (n,dx,incx,dy,incy)

c     interchanges two vectors.
c     uses unrolled loops for increments equal one.
c     jack dongarra, linpack, 3/11/78.
      implicit   none

      double precision dx(1),dy(1),dtemp
      integer i,incx,incy,ix,iy,m,mp1,n

      if(n.le.0)return

c     code for both increments equal to 1

      if(incx.eq.1.and.incy.eq.1) then

c       clean-up loop

        m = mod(n,3)
        if( m .gt. 0 ) then
          do i = 1,m
            dtemp = dx(i)
            dx(i) = dy(i)
            dy(i) = dtemp
          end do ! i
          if( n .lt. 3 ) return
        endif
        mp1 = m + 1
        do i = mp1,n,3
          dtemp     = dx(i)
          dx(i    ) = dy(i)
          dy(i    ) = dtemp
          dtemp     = dx(i + 1)
          dx(i + 1) = dy(i + 1)
          dy(i + 1) = dtemp
          dtemp     = dx(i + 2)
          dx(i + 2) = dy(i + 2)
          dy(i + 2) = dtemp
        end do ! i

c     code for unequal increments or equal increments not equal to 1

      else
        ix = 1
        iy = 1
        if(incx.lt.0) ix = (-n+1)*incx + 1
        if(incy.lt.0) iy = (-n+1)*incy + 1
        do i = 1,n
          dtemp  = dx(ix)
          dx(ix) = dy(iy)
          dy(iy) = dtemp
          ix     = ix + incx
          iy     = iy + incy
        end do ! i
      endif

      end
