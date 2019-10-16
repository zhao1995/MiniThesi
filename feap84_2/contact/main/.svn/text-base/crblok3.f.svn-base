c$Id:$
      subroutine crblok3(nope,dnope,ics,neps,nsopt,scdat,polfl,
     &           ix,ip,x,norm,ndm,nen,nen1,numnp,numel,prt,prth)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c     1. Use 'pi' from 'pconstant.h'                        14/11/2006
c     2. Add ctface to permit triangular facets             30/03/2012
c     3. Correct set of angle for polar faces               26/02/2013
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Define surface nodes for 3-D contact slidelines

c      Inputs:
c         nope        - Nodes per element
c         dnope       - Dimension Nodes per element
c         nsopt       - Subcommand # 1=gap,2=segm,3=pola,4=cart,5=regi
c         scdat(*)    - Subcommand data
c         ix(nen1,*)  - Element nodal connection data
c         x(ndm,*)    - Nodal coordinates of mesh
c         norm(3,*)   - Nodal normals of mesh
c         ndm         - Spatial dimension of mesh
c         nen         - Maximum number of nodes/element
c         nen1        - Dimension of ix array
c         numnp       - Number of nodes in mesh
c         numel       - Number of elements in mesh
c         prt         - Output results if true
c         prth        - Output title/header data if true

c      Scratch:
c         ip(*)       - Nodal integer list storage

c      Output:
c         ics(*)      - Contact facets
c         neps        - Last facet number
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'iofile.h'
      include  'pconstant.h'

      logical   prt,prth,errck,gettd, polfl, noconv
      integer   nsopt,nope,dnope,neps,ndm,nen,nen1,numnp,numel
      integer   i,j,k,n,nits,numprt, count, onsurf, regn
      integer   ics(dnope,*),ix(nen1,*),ip(*),ld(9)
      real*8    gap,gap0, r1,r2, th, dx,dy,cn,sn
      real*8    c1,c2,c3,c4,c5,c6,c7,c8,c9, a11,a12,a22, det, z3
      real*8    xi1,xi2, di1,di2, x12,x1s,x2s, tol, tol0
      real*8    scdat(14),norm(3,*), x(ndm,*),xl(3,9),x0(3), dnorm(3)
      real*8    v1(3),v2(3), e1(3),e2(3),e3(3), z0(3), zl(3,9)
      real*8    zp(3),dz(3), zmax(3),zmin(3), td(4), xx(3)

      save

      data      nits / 50 /, tol /1.d-12/, tol0 /1.d-5/
      data      gap0  / 0.0d0 /
      data      regn  / 0 /

      call cdebug0 ('        crblok3',-1)

c     Coordinate contact surface descriptions

      numprt = 0


      if( nsopt.eq.1 ) then
        gap0 = scdat(1)
        write(iow,*) '      GAP =',gap0
        if(ior.lt.0) then
          write(*,*) '      GAP =',gap0
        end if

      elseif( nsopt.eq.2 ) then

c       Input coordinates for surface search

        do j = 1,9
          ld(j) = 0
        end do
        do j = 1,9
          errck = gettd(td,4,'skip')
          if(errck) then
            backspace ior
            go to 10
          endif
          n     = int(td(1))
          if(n.eq.0) go to 10
          ld(n) = 1
          do i = 1,3
            xl(i,n) = td(i+1)
          end do
        end do

c       Fill missing mid-side nodes

10      do i = 1,4
          if(ld(i+4).eq.0) then
            j = mod(i,4) + 1
            xl(1,i+4) = 0.5d0*(xl(1,i) + xl(1,j))
            xl(2,i+4) = 0.5d0*(xl(2,i) + xl(2,j))
            xl(3,i+4) = 0.5d0*(xl(3,i) + xl(3,j))
          end if
        end do

c       Fill central node

        if(ld(9).eq.0) then
           xl(1,9) = 0.50d0*(xl(1,5) + xl(1,6) + xl(1,7) + xl(1,8))
     &             - 0.25d0*(xl(1,1) + xl(1,2) + xl(1,3) + xl(1,4))
           xl(2,9) = 0.50d0*(xl(2,5) + xl(2,6) + xl(2,7) + xl(2,8))
     &             - 0.25d0*(xl(2,1) + xl(2,2) + xl(2,3) + xl(2,4))
           xl(3,9) = 0.50d0*(xl(3,5) + xl(3,6) + xl(3,7) + xl(3,8))
     &             - 0.25d0*(xl(3,1) + xl(3,2) + xl(3,3) + xl(3,4))
        end if

c       Output data

        if(prt) then
          call prtitl(prth)
          write(iow,2000) (i,xl(1,i),xl(2,i),xl(3,i),i=1,9)
          if(ior.lt.0) then
            write(*,2000) (i,xl(1,i),xl(2,i),xl(3,i),i=1,9)
          end if
        endif

c       Compute bounding box

        do i = 1,3
          v1(i) = xl(i,3) - xl(i,1)
          v2(i) = xl(i,4) - xl(i,2)
        end do

        r1 = 1.d0/sqrt(v1(1)**2 + v1(2)**2 + v1(3)**2)
        r2 = 1.d0/sqrt(v2(1)**2 + v2(2)**2 + v2(3)**2)

        do i = 1,3
          v1(i) = v1(i)*r1
          v2(i) = v2(i)*r2
        end do

        do i = 1,3
          e1(i) = v1(i) - v2(i)
          e2(i) = v1(i) + v2(i)
        end do

        r1 = 1.d0/sqrt(e1(1)**2 + e1(2)**2 + e1(3)**2)
        r2 = 1.d0/sqrt(e2(1)**2 + e2(2)**2 + e2(3)**2)

        do i = 1,3
          e1(i) = e1(i)*r1
          e2(i) = e2(i)*r2
        end do

        e3(1) = e1(2)*e2(3) - e1(3)*e2(2)
        e3(2) = e1(3)*e2(1) - e1(1)*e2(3)
        e3(3) = e1(1)*e2(2) - e1(2)*e2(1)

        do j = 1,9
          zl(1,j) = e1(1)*xl(1,j) + e1(2)*xl(2,j) + e1(3)*xl(3,j)
          zl(2,j) = e2(1)*xl(1,j) + e2(2)*xl(2,j) + e2(3)*xl(3,j)
          zl(3,j) = e3(1)*xl(1,j) + e3(2)*xl(2,j) + e3(3)*xl(3,j)
        end do

c       Set max/min based on nodes of patch

        do i = 1,3
          zmin(i) = zl(i,1)
          zmax(i) = zl(i,1)
          do j = 1,9
            zmin(i) = min(zmin(i),zl(i,j))
            zmax(i) = max(zmax(i),zl(i,j))
          end do
        end do

c       Check for max/min on edges

        do i = 1,3
          do j = 1,4
            k = mod(j,4) + 1
            r1 = 0.5d0*(zl(i,k) - zl(i,j))
            r2 = 2.d0*zl(i,j+4) - zl(i,j) - zl(i,k)
            if(abs(r2).gt.abs(r1)) then
              zmax(i) = max(zmax(i),zl(i,j+4) + 0.5d0*r1*r1/r2)
              zmin(i) = min(zmin(i),zl(i,j+4) + 0.5d0*r1*r1/r2)
            endif
          end do
        end do

c       Convert coordinates to hierarchical form

        do i = 1,3
          zl(i,9) = zl(i,9) - 0.50d0*(zl(i,5)+zl(i,6)+zl(i,7)+zl(i,8))
     &                      + 0.25d0*(zl(i,1)+zl(i,2)+zl(i,3)+zl(i,4))
          do n = 1,4
            j = mod(n,4) + 1
            zl(i,n+4) = zl(i,n+4) - 0.5d0*(zl(i,n) + zl(i,j))
          end do
        end do

c       Check interior for maximum/minimums of z_3 coordinate.

        c1 = 0.25d0*( zl(3,1)+zl(3,2)+zl(3,3)+zl(3,4))
        c2 = 0.25d0*(-zl(3,1)+zl(3,2)+zl(3,3)-zl(3,4))
        c3 = 0.25d0*(-zl(3,1)-zl(3,2)+zl(3,3)+zl(3,4))
        c4 = 0.25d0*( zl(3,1)-zl(3,2)+zl(3,3)-zl(3,4))

        c5 =         zl(3,7) + zl(3,5)
        c6 = 0.50d0*(zl(3,7) - zl(3,5))

        c7 =         zl(3,6) + zl(3,8)
        c8 = 0.50d0*(zl(3,6) - zl(3,8))

        c9 = 2.00d0* zl(3,9)

        xi1 = 0.0d0
        xi2 = 0.0d0

        noconv = .true.
        count  = 0

        do while (noconv .and. count.lt.nits)

          count = count + 1

          x1s = 1.d0 - xi1*xi1
          x2s = 1.d0 - xi2*xi2
          x12 = 2.d0*xi1*xi2

          r1 = c1 - c5*xi1 + c4*xi2  + x2s*(c8 - xi1*c9) - x12*c6
          r2 = c2 + c4*xi1 - c7*xi2  + x1s*(c6 - xi2*c9) - x12*c8

          a11 = -c5 -2.d0*xi2*c6 - (1.d0-xi2*xi2)*c9
          a22 = -c7 -2.d0*xi1*c8 - (1.d0-xi1*xi1)*c9
          a12 =  c4 -2.d0*(xi1*c6 + xi2*(c8 - xi1*c9))

          det = a11*a22 - a12*a12

          if(det.gt.tol) then
            di1 = (-a22*r1 + a12*r2)/det
            di2 = ( a12*r1 - a11*r2)/det

            xi1 = xi1 + di1
            xi2 = xi2 + di2

            if(di1**2+di2**2 .le. tol*(xi1**2+xi2**2)) then
              noconv = .false.
              if(abs(xi1).lt.1.d0 .and. abs(xi2).lt.1.d0) then
               z3 = c1 + c2*xi1 + c3*xi2 + c4*x12
     &            + x1s*(0.5d0*c5 + xi2*c6) + x2s*(0.5d0*c7 + xi1*c8)
     &            + x1s*x2s*zl(3,9)
               zmax(3) = max(zmax(3),z3)
               zmin(3) = min(zmin(3),z3)
              endif
            endif
          else
            noconv = .false.
          end if
        end do

c       Set center of bounding box and size of half lengths

        do i = 1,3
          z0(i) = 0.5d0*(zmax(i) + zmin(i))
          dz(i) = 0.5d0*(zmax(i) - zmin(i))
        end do

c       Set gap values and increase bounding box to capture nodes

        gap = max(gap0,tol0*max(dz(1),dz(2),dz(3)))

        do i = 1,3
          dz(i) = dz(i) + 1.5*gap
        end do

c       Shift vertex nodes of patch to relative coordinates to z0.

        do j = 1,4
          do i = 1,3
            zl(i,j) = zl(i,j) - z0(i)
          end do
        end do

c       Find nodes within bounding box and check if in surface

        if(polfl) then
          cn = cos(z0(2)/th)
          sn = sin(z0(2)/th)
        endif
        do n = 1,numnp
          ip(n) = 0
          norm(1,n) = 0.0d0
          norm(2,n) = 0.0d0
          norm(3,n) = 0.0d0
          xx(1) = x(1,n)
          xx(2) = x(2,n)
          xx(3) = x(3,n)

c         Polar adjustment

          if(polfl) then
            r1    =  sqrt ((xx(1)-x0(1))**2 + (xx(2)-x0(2))**2)
            dx    =  cn*(xx(1) - x0(1)) + sn*(xx(2) - x0(2))
            dy    = -sn*(xx(1) - x0(1)) + cn*(xx(2) - x0(2))
            r2    = atan2(dy,dx)*th + z0(2)
            xx(1) = r1
            xx(2) = r2
            xx(3) = xx(3) - x0(3)
            zp(1) = e1(1)*xx(1) + e1(2)*xx(2) + e1(3)*xx(3)

            if(zmin(1).lt.180.0d0 .and. zmax(1).gt.180 .and.
     &                                  zp(1).lt.0.0d0) then
              zp(1) = zp(1) + z0(1)
            else
              zp(1) = zp(1) - z0(1)
            endif

c         Cartesian adjustment

          else
            zp(1) = e1(1)*xx(1) + e1(2)*xx(2) + e1(3)*xx(3) - z0(1)
          endif

          if(abs(zp(1)).le.dz(1)) then
            zp(2) = e2(1)*xx(1) + e2(2)*xx(2) + e2(3)*xx(3) - z0(2)
            if(abs(zp(2)).le.dz(2)) then
              zp(3) = e3(1)*xx(1) + e3(2)*xx(2) + e3(3)*xx(3) -z0(3)
              if(abs(zp(3)).le.dz(3)) then
                ip(n)     = onsurf(zp,zl,dz,gap,dnorm,9)
                norm(1,n) = e1(1)*dnorm(1)
     &                    + e2(1)*dnorm(2)
     &                    + e3(1)*dnorm(3)
                norm(2,n) = e1(2)*dnorm(1)
     &                    + e2(2)*dnorm(2)
     &                    + e3(2)*dnorm(3)
                norm(3,n) = e1(3)*dnorm(1)
     &                    + e2(3)*dnorm(2)
     &                    + e3(3)*dnorm(3)
              end if
            end if
          end if
        end do ! n

c       Convert polar normals to cartesian form

        if(polfl) then
          do n = 1,numnp
            if(ip(n).ne.0) then
              r2        = atan2(x(2,n)-x0(2),x(1,n)-x0(1))
              r1        = cos(r2)
              r2        = sin(r2)
              c1        = r1*norm(1,n) - r2*norm(2,n)
              norm(2,n) = r2*norm(1,n) + r1*norm(2,n)
              norm(1,n) = c1
            endif
          end do ! n
        endif

c       Compute face lists

        if(nope.eq.1) then
          do n = 1,numnp
            if(ip(n).gt.0) then
              neps = neps + 1
              ics(1,neps) = n
            endif
          end do
        else
          do n = 1,numel
            if(ix(nen1-1,n).eq.regn .or. regn.eq.0) then
              if(nope.eq.3) then
                call ctface3(dnope,ics,ix(1,n),ip,x,norm,
     &                       nen,ndm,neps)
              else
                call crface3(nope,dnope,ics,ix(1,n),ip,x,norm,
     &                       nen,ndm,neps)
              endif
            endif
          end do
        endif

c     Segment faces to have coordinates in polar form

      elseif( nsopt.eq.3 ) then

        polfl = .true.
        do i = 1,3
          x0(i) = scdat(i)
        end do
        th = 180.0d0/pi
        write(iow,2001) 'POLAR',x0

c     Segment faces to have coordinates in Cartesian form

      elseif( nsopt.eq.4 ) then

        polfl = .false.
        write(iow,2001) 'CARTESIAN',x0

c     Set region number

      elseif( nsopt.eq.5 ) then

        regn = nint(scdat(1))

      endif

c     Formats

2000  format('   C o n t a c t   S u r f a c e   P a t c h',
     &       '   C o o r d i n a t e s'//
     &      ('         Node    X-coord.     Y-coord.     Z-coord.')//
     &        (i12,1p,3e13.4))

2001  format(5x,'Block coordinates in ',a,' form'/
     &       7x,'Reference location:'/
     &       9x,'x_1 =',1p,1e12.4,' x_2 =',1p,1e12.4,' x_3 =',1p,1e12.4)

      end
