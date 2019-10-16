      subroutine  dscal(n,da,dx,incx)

c     scales a vector by a constant.
c     uses unrolled loops for increment equal to one.
c     jack dongarra, linpack, 3/11/78.
c     modified to correct problem with negative increment, 8/21/90.
      implicit   none

      double precision da,dx(1)
      integer i,incx,ix,m,mp1,n

      if(n.le.0)return
c     code for increment equal to 1
      if(incx.eq.1) then

c        clean-up loop

         m = mod(n,5)
         if( m .gt. 0 ) then
           do i = 1,m
             dx(i) = da*dx(i)
           end do ! i
           if( n .lt. 5 ) return
         endif
         mp1 = m + 1
         do i = mp1,n,5
           dx(i    ) = da*dx(i)
           dx(i + 1) = da*dx(i + 1)
           dx(i + 2) = da*dx(i + 2)
           dx(i + 3) = da*dx(i + 3)
           dx(i + 4) = da*dx(i + 4)
         end do ! i

c     code for increment not equal to 1

      else
        ix = 1
        if(incx.lt.0)ix = (-n+1)*incx + 1
        do i = 1,n
          dx(ix) = da*dx(ix)
          ix     = ix + incx
        end do ! i
      endif

      end
