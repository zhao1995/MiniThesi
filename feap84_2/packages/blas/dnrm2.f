      double precision function dnrm2 ( n, dx, incx)
      implicit   none
      integer i, incx, ix, j, n, next
      double precision   dx(1), cutlo, cuthi, hitest, sum, xmax,zero,one
      data   zero, one /0.0d0, 1.0d0/

c     euclidean norm of the n-vector stored in dx() with storage
c     increment incx .
c     if    n .le. 0 return with result = 0.
c     if n .ge. 1 then incx must be .ge. 1

c           c.l.lawson, 1978 jan 08
c     modified to correct problem with negative increment, 8/21/90.

c     four phase method     using two built-in constants that are
c     hopefully applicable to all machines.
c         cutlo = maximum of  dsqrt(u/eps)  over all known machines.
c         cuthi = minimum of  dsqrt(v)      over all known machines.
c     where
c         eps = smallest no. such that eps + 1. .gt. 1.
c         u   = smallest positive no.   (underflow limit)
c         v   = largest  no.            (overflow  limit)

c     brief outline of algorithm..

c     phase 1    scans zero components.
c     move to phase 2 when a component is nonzero and .le. cutlo
c     move to phase 3 when a component is .gt. cutlo
c     move to phase 4 when a component is .ge. cuthi/m
c     where m = n for x() real and m = 2*n for complex.

c     values for cutlo and cuthi..
c     from the environmental parameters listed in the imsl converter
c     document the limiting values are as follows..
c     cutlo, s.p.   u/eps = 2**(-102) for  honeywell.  close seconds are
c                   univac and dec at 2**(-103)
c                   thus cutlo = 2**(-51) = 4.44089e-16
c     cuthi, s.p.   v = 2**127 for univac, honeywell, and dec.
c                   thus cuthi = 2**(63.5) = 1.30438e19
c     cutlo, d.p.   u/eps = 2**(-67) for honeywell and dec.
c                   thus cutlo = 2**(-33.5) = 8.23181d-11
c     cuthi, d.p.   same as s.p.  cuthi = 1.30438d19
c     data cutlo, cuthi / 8.232d-11,  1.304d19 /
c     data cutlo, cuthi / 4.441e-16,  1.304e19 /
c     data cutlo, cuthi / 8.232d-11,  1.304d19 /

c     From Ed Anderson

      data cutlo /   0.1415686533102923D-145 /
      data cuthi /   0.1340780792994260D+155 /

      if(n .le. 0)  then
         dnrm2  = zero
      else
        next = 1
        sum  = zero
        xmax = zero
        i = 1
        if( incx .lt. 0 )i = (-n+1)*incx + 1
        ix = 1

c       Begin main loop

   10   go to (20, 30, 60, 70), next

   20   if( dabs(dx(i)) .gt. cutlo) go to 100
        next = 2
        xmax = zero

c       Phase 1.  sum is zero

   30   if( dx(i) .eq. zero) go to 110
        if( dabs(dx(i)) .gt. cutlo) go to 100

c       Prepare for phase 2.
        next = 3
        go to 50

c       Prepare for phase 4.

   40   continue
        next = 4
        sum  = (sum / dx(i)) / dx(i)
   50   xmax = dabs(dx(i))
        go to 80

c       Phase 2.  Sum is small.
c                 Scale to avoid destructive underflow.

   60   if( dabs(dx(i)) .gt. cutlo ) go to 90

c       Common code for phases 2 and 4.
c       In phase 4 sum is large.  Scale to avoid overflow.

   70   if( dabs(dx(i)) .le. xmax ) go to 80
           sum = one + sum * (xmax / dx(i))**2
           xmax = dabs(dx(i))
           go to 110

   80   sum = sum + (dx(i)/xmax)**2
        go to 110

c       Prepare for phase 3.

   90   sum = (sum * xmax) * xmax

c       For real or d.p. set hitest = cuthi/n
c       For complex      set hitest = cuthi/(2*n)

  100   hitest = cuthi/float( n )

c       Phase 3.  sum is mid-range.  no scaling.

        do j = ix,n
          if(dabs(dx(i)) .ge. hitest) go to 40
           sum = sum + dx(i)**2
           i = i + incx
        end do ! j
        dnrm2 = dsqrt( sum )
        go to 120

  110   continue
        ix = ix + 1
        i  = i + incx
        if( ix .le. n ) go to 10

c       End of main loop.

c       Compute square root and adjust for scaling.

        dnrm2 = xmax * dsqrt(sum)
      endif

  120 continue

      end
