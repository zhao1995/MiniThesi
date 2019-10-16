      subroutine  dcopy(n,dx,incx,dy,incy)

c     copies a vector, x, to a vector, y.
c     uses unrolled loops for increments equal to one.
c     jack dongarra, linpack, 3/11/78.
      implicit   none

      double precision dx(1),dy(1)
      integer i,incx,incy,ix,iy,m,mp1,n

      if(n.le.0)return

c     code for both increments equal to 1

      if(incx.eq.1.and.incy.eq.1) then

c       clean-up loop

        m = mod(n,7)
        if( m .gt. 0 ) then
          do i = 1,m
            dy(i) = dx(i)
          end do ! i
          if( n .lt. 7 ) return
        endif
        mp1 = m + 1
        do i = mp1,n,7
          dy(i) = dx(i)
          dy(i + 1) = dx(i + 1)
          dy(i + 2) = dx(i + 2)
          dy(i + 3) = dx(i + 3)
          dy(i + 4) = dx(i + 4)
          dy(i + 5) = dx(i + 5)
          dy(i + 6) = dx(i + 6)
        end do ! i

c        code for unequal increments or equal increments
c          not equal to 1

      else
        ix = 1
        iy = 1
        if(incx.lt.0)ix = (-n+1)*incx + 1
        if(incy.lt.0)iy = (-n+1)*incy + 1
        do i = 1,n
          dy(iy) = dx(ix)
          ix = ix + incx
          iy = iy + incy
        end do ! i
      endif

      end
