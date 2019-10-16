      subroutine daxpy(n,da,dx,incx,dy,incy)

c     constant times a vector plus a vector.
c     uses unrolled loops for increments equal to one.
c     jack dongarra, linpack, 3/11/78.
      implicit   none

      double precision dx(1),dy(1),da
      integer i,incx,incy,ix,iy,m,mp1,n

      if(n.le.0)return
      if (da .eq. 0.0d0) return

c      code for both increments equal to 1

      if(incx.eq.1.and.incy.eq.1) then

c       clean-up loop

        m = mod(n,4)
        if( m .gt. 0 ) then
          do i = 1,m
            dy(i) = dy(i) + da*dx(i)
          end do ! i
          if( n .lt. 4 ) return
        endif
        mp1 = m + 1
        do i = mp1,n,4
          dy(i) = dy(i) + da*dx(i)
          dy(i + 1) = dy(i + 1) + da*dx(i + 1)
          dy(i + 2) = dy(i + 2) + da*dx(i + 2)
          dy(i + 3) = dy(i + 3) + da*dx(i + 3)
        end do ! i

c     code for unequal increments or equal increments
c     not equal to 1

      else

        ix = 1
        iy = 1
        if(incx.lt.0)ix = (-n+1)*incx + 1
        if(incy.lt.0)iy = (-n+1)*incy + 1
        do i = 1,n
          dy(iy) = dy(iy) + da*dx(ix)
          ix = ix + incx
          iy = iy + incy
        end do ! i
      endif

      end
