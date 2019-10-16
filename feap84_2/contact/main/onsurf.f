c$Id:$
      integer function onsurf(zp,zl,dz,gap,dnorm,nel)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Define surface nodes for 3-D contact slidelines

c      Inputs:
c         zp(3)       - Transformed coordinates for node point
c         zl(3,*)     - Surface patch coordinates
c         dz(3)       - Size of bounding box
c         gap         - Size of gap to use in search
c         nel         - Number of nodes on patch

c      Scratch:
c         onsurf      - 1 = on surface, 0 = not on surface
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'iofile.h'

      logical   noconv
      integer   nel, i,j, count,maxc
      real*8    gap, len2, f1,f2, a11,a12,a22, det, tol, tolc
      real*8    zp(3), zl(3,*), dz(3), zs(3), z1(2,3), z2(3,3)
      real*8    xi(2), dxi(2), dnorm(3), shp(9), dshp(2,9), d2shp(3,9)

      save

      data      maxc /20/, tol /1.d-20/, tolc /1.d-03/

c     Compute initial values of surface coordinates

      xi(1) = zp(1)/dz(1)
      xi(2) = zp(2)/dz(2)

      xi(1) = 0.0d0
      xi(2) = 0.0d0
      noconv = .true.
      count  =  0

      do while (noconv .and. count.le.maxc)

c       Compute the shape functions and derivatives

        if(nel.eq.4) then
          call surf4(xi, shp, dshp, d2shp )
        else
          call surf9(xi, shp, dshp, d2shp )
        end if

c       Compute points and derivatives

        do i = 1,3
          zs(i)   = 0.0d0
          z1(1,i) = 0.0d0
          z1(2,i) = 0.0d0
          z2(3,i) = 0.0d0
          do j = 1,nel
            zs(i)   = zs(i)   + shp(j)*zl(i,j)
            z1(1,i) = z1(1,i) + dshp(1,j)*zl(i,j)
            z1(2,i) = z1(2,i) + dshp(2,j)*zl(i,j)
            z2(3,i) = z2(3,i) + d2shp(3,j)*zl(i,j)
          end do

          if(nel.eq.9) then
            z2(1,i) = d2shp(1,5)*zl(i,5) + d2shp(1,7)*zl(i,7)
     &              + d2shp(1,9)*zl(i,9)
            z2(2,i) = d2shp(2,6)*zl(i,6) + d2shp(2,8)*zl(i,8)
     &              + d2shp(2,9)*zl(i,9)
          endif
        end do

c       Compute residual

        f1 = (zs(1)-zp(1))*z1(1,1) + (zs(2)-zp(2))*z1(1,2)
     &     + (zs(3)-zp(3))*z1(1,3)
        f2 = (zs(1)-zp(1))*z1(2,1) + (zs(2)-zp(2))*z1(2,2)
     &     + (zs(3)-zp(3))*z1(2,3)

c       Compute tangent

        a11 = z1(1,1)*z1(1,1) + z1(1,2)*z1(1,2) + z1(1,3)*z1(1,3)
        a12 = z1(1,1)*z1(2,1) + z1(1,2)*z1(2,2) + z1(1,3)*z1(2,3)
        a22 = z1(2,1)*z1(2,1) + z1(2,2)*z1(2,2) + z1(2,3)*z1(2,3)

        if(count.gt.3) then
        a11 = a11
     &      + (zs(1)-zp(1))*z2(1,1) + (zs(2)-zp(2))*z2(1,2)
     &      + (zs(3)-zp(3))*z2(1,3)

        a12 = a12
     &      + (zs(1)-zp(1))*z2(3,1) + (zs(2)-zp(2))*z2(3,2)
     &      + (zs(3)-zp(3))*z2(3,3)

        a22 = a22
     &      + (zs(1)-zp(1))*z2(2,1) + (zs(2)-zp(2))*z2(2,2)
     &      + (zs(3)-zp(3))*z2(2,3)
        endif

c       Solve Newton step

        det = a11*a22 - a12*a12
        if(det.ne.0.0d0) then
          det = 1.d0/det
          dxi(1) = (-a22*f1 + a12*f2)*det
          dxi(2) = ( a12*f1 - a11*f2)*det

          xi(1)  = xi(1) + dxi(1)
          xi(2)  = xi(2) + dxi(2)

c         Check convergence

          count = count + 1

          len2  = xi(1)**2 + xi(2)**2

          if(dxi(1)**2 + dxi(2)**2 .lt. tol) then
            noconv = .false.
          endif
        else
          write(  *,*) 'XI',xi(1),xi(2)
          write(  *,*) 'fi',f1,f2
          write(  *,*) 'ai',a11,a12,a22
          write(  *,3000)
          write(ilg,*) 'XI',xi(1),xi(2)
          write(ilg,*) 'fi',f1,f2
          write(ilg,*) 'ai',a11,a12,a22
          write(ilg,3000)
          call plstop()
        endif

      end do

      len2 = (zp(1)-zs(1))**2 + (zp(2)-zs(2))**2 + (zp(3)-zs(3))**2

c     Set value of function for return

      if(len2 .le. gap*gap .and. abs(xi(1)).le.1.0d0 + tolc
     &                     .and. abs(xi(2)).le.1.0d0 + tolc) then
        onsurf = 1
        dnorm(1) = z1(1,2)*z1(2,3) - z1(1,3)*z1(2,2)
        dnorm(2) = z1(1,3)*z1(2,1) - z1(1,1)*z1(2,3)
        dnorm(3) = z1(1,1)*z1(2,2) - z1(1,2)*z1(2,1)

        len2     = 1.d0/sqrt(dnorm(1)**2 + dnorm(2)**2 + dnorm(3)**2)

        dnorm(1) = dnorm(1)*len2
        dnorm(2) = dnorm(2)*len2
        dnorm(3) = dnorm(3)*len2
      else
        onsurf = 0

        dnorm(1) = 0.0d0
        dnorm(2) = 0.0d0
        dnorm(3) = 0.0d0

      endif

c     Format

3000  format(' *ERROR* ONSURF: Zero Determinant in ONSURF')

      end
