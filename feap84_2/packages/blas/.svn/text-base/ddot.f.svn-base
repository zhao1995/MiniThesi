      double precision function ddot(n,dx,incx,dy,incy)

c     forms the dot product of two vectors.
c     uses unrolled loops for increments equal to one.
c     jack dongarra, linpack, 3/11/78.
      implicit   none

      double precision dx(1),dy(1),dtemp
      integer i,incx,incy,ix,iy,m,mp1,n

      ddot  = 0.0d0
      if(n.le.0)return
      dtemp = 0.0d0

c     code for both increments equal to 1

      if(incx.eq.1.and.incy.eq.1) then

c       clean-up loop

        m = mod(n,5)
        if( m .gt. 0 ) then
          do i = 1,m
            dtemp = dtemp + dx(i)*dy(i)
          end do ! i
          if( n .lt. 5 ) go to 100
        endif
        mp1 = m + 1
        do i = mp1,n,5
          dtemp = dtemp + dx(i  )*dy(i  ) + dx(i+1)*dy(i+1)
     *                  + dx(i+2)*dy(i+2) + dx(i+3)*dy(i+3)
     &                  + dx(i+4)*dy(i+4)
        end do ! i

c     code for unequal increments or equal increments not equal to 1

      else
        ix = 1
        iy = 1
        if(incx.lt.0)ix = (-n+1)*incx + 1
        if(incy.lt.0)iy = (-n+1)*incy + 1
        do i = 1,n
          dtemp = dtemp + dx(ix)*dy(iy)
          ix    = ix + incx
          iy    = iy + incy
        end do ! i
      endif

  100 ddot = dtemp

      end
