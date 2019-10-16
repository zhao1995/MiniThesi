c$Id:$
      subroutine crblen3(nope,dnope,ics,neps,nsopt,scdat,xs,ix,ip,x,
     &                   norm,ndm,nen,nen1,numnp,numel,prt,prth)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Add ctface to permit triangular facets           30/03/2012
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Define blended surface nodes for 3-D contact slideline

c      Inputs:
c         nope        - Nodes per element
c         dnope       - Dimension Nodes per element
c         nsopt       - Subcommand # 1=gap,2=segm,3=pola,4=cart
c         scdat(*)    - Subcommand data
c         xs(3,*)     - Supernode coordinates
c         ix(nen1,*)  - Element nodal connection data
c         x(ndm,*)    - Nodal coordinates of mesh
c         norm(3,*)   - Normal to nodes on surface
c         ndm         - Spatial dimension of mesh
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
      include  'pointer.h'
      include  'comblk.h'

      include  'p_int.h'

      logical   prt,prth,errck,gettd, wdflg, extfl
      integer   nsopt,nope,dnope,neps,ndm,nen,nen1,numnp,numel
      integer   i,j,n,numprt, regn, ii,i1,jj,j1,kk, onsurf
      integer   ics(dnope,*),ix(nen1,*),ip(*)
      integer   tblend(20),is(16),iside(4),isd
      real*8    gap,gap0,tol0, dxi, xi1,xi2,r1,r2
      real*8    scdat(14), x(ndm,*),xl(3,9),norm(3,*)
      real*8    td(4), zp(3),zl(3,4),zmin(3),zmax(3),z0(3),dz(3),xx(3)
      real*8    v1(3),v2(3),e1(3),e2(3),e3(3),dnorm(3)
      real*8    shp(181),xr(3,181,4),rt(3,13), xs(3,*),tr(3,4)

      save

      data      tr    /1.d0,3*0.d0,1.d0,3*0.d0,1.d0,3*0.d0 /
      data      gap0  / 0.0d0   /
      data      tol0  / 1.0d-05 /
      data      regn  / 0       /
      data      isd   / 16      /

      data      extfl / .false. /

      call cdebug0 ('        crblen3',-1)

c     Coordinate contact surface descriptions

      numprt = 0

      if( nsopt.eq.1 ) then
        gap0 = scdat(1)
        write(iow,*) '      GAP =',gap0
        if(ior.lt.0) then
          write(*,*) '      GAP =',gap0
        end if

      elseif( nsopt.eq.2 ) then

        do n = 1,numnp
          ip(n)     = 0
          norm(1,n) = 0.0d0
          norm(2,n) = 0.0d0
          norm(3,n) = 0.0d0
        end do ! n

c       Define segments with supernodes for surface search

        if(prth) then
          write(iow,2000)
        endif
        wdflg = .true.
        do while (wdflg)

c         Define list of segment nodes in IS array

          errck = gettd(td,15,'skip')
          if(.not.errck) then
            j = 0
            do i = 1,15
              is(i) = nint(td(i))
              if(is(i).gt.0) j = i
            end do
            if(j.eq.0) errck = .true.
          endif

c         Check for error

          if(errck) then
            backspace ior
            wdflg = .false.

c         Do search to find nodal points within gap of blend

          else
            if(prt) then
              write(iow,2001) (is(i),i=1,j)
            endif

            do i = 1,20
              tblend(i) = 0
            end do ! i
            do i = 1,j
              tblend(i+10) = is(i)
            end do ! i

c           Set up data for the 2-d surface blending

            call pblend2a(tblend,iside,isd)

c           Establish side coordinates for grid on blend surface

            do i = 1,4
              fp(1) = np(162) + isd*(iside(i)-1)
              do ii = 1,isd-1
                if(mr(fp(1)+ii).eq.0) go to 100
              end do ! ii
100           jj = 181
              j1 = jj + 1
              call pside1(jj-1,xs,tr,1,mr(fp(1)+1),ii-1,ndm,shp,
     &                    rt,xr(1,1,i),mr(fp(1)))
            end do ! i

c           Check that each side is properly directed

            fp(1) = np(162) + isd*(iside(1)-1) + 1
            if(is(1).ne.mr(fp(1))) then
              do ii = 1,jj/2
                do i = 1,3
                  xi1           = xr(i,   ii,1)
                  xr(i,   ii,1) = xr(i,j1-ii,1)
                  xr(i,j1-ii,1) = xi1
                end do ! i
              end do ! ii
            endif

            fp(1) = np(162) + isd*(iside(4)-1) + 1
            if(is(1).ne.mr(fp(1))) then
              do ii = 1,jj/2
                do i = 1,3
                  xi1           = xr(i,   ii,4)
                  xr(i,   ii,4) = xr(i,j1-ii,4)
                  xr(i,j1-ii,4) = xi1
                end do ! i
              end do ! ii
            endif

            fp(1) = np(162) + isd*(iside(2)-1) + 1
            if(is(3).eq.mr(fp(1))) then
              do ii = 1,jj/2
                do i = 1,3
                  xi1           = xr(i,   ii,2)
                  xr(i,   ii,2) = xr(i,j1-ii,2)
                  xr(i,j1-ii,2) = xi1
                end do ! i
              end do ! ii
            endif

            fp(1) = np(162) + isd*(iside(3)-1) + 1
            if(is(3).eq.mr(fp(1))) then
              do ii = 1,jj/2
                do i = 1,3
                  xi1           = xr(i,   ii,3)
                  xr(i,   ii,3) = xr(i,j1-ii,3)
                  xr(i,j1-ii,3) = xi1
                end do ! i
              end do ! ii
            endif

c           Get coordinates for corners of search patch

            dxi =  2.d0/dble(jj-1)
            xi2 = -1.d0
            do j1 = 1,jj-1
              xi1 = -1.d0
              do i1 = 1,jj-1
                call phibl(xr,xi1,xi2,dxi, xl, i1,j1,jj)

c               Compute bounding box

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

                do j = 1,4
                  zl(1,j) = e1(1)*xl(1,j)
     &                    + e1(2)*xl(2,j)
     &                    + e1(3)*xl(3,j)
                  zl(2,j) = e2(1)*xl(1,j)
     &                    + e2(2)*xl(2,j)
     &                    + e2(3)*xl(3,j)
                  zl(3,j) = e3(1)*xl(1,j)
     &                    + e3(2)*xl(2,j)
     &                    + e3(3)*xl(3,j)
                end do

c               Set max/min based on nodes of patch

                do i = 1,3
                  zmin(i) = zl(i,1)
                  zmax(i) = zl(i,1)
                  do j = 1,4
                    zmin(i) = min(zmin(i),zl(i,j))
                    zmax(i) = max(zmax(i),zl(i,j))
                  end do
                end do

c               Set center of bounding box and size of half lengths

                do i = 1,3
                  z0(i) = 0.5d0*(zmax(i) + zmin(i))
                  dz(i) = 0.5d0*(zmax(i) - zmin(i))
                end do

c               Set gap values & increase bounding box to capture nodes

                gap = max(gap0,tol0*max(dz(1),dz(2),dz(3)))

                do i = 1,3
                  dz(i) = dz(i) + 1.5*gap
                end do

c               Shift vertex nodes of patch to relative coordinates z0.

                do j = 1,4
                  do i = 1,3
                    zl(i,j) = zl(i,j) - z0(i)
                  end do
                end do

c               Find nodes within bounding box and check if on surface

                do n = 1,numnp
                  if(ip(n).eq.0) then
                    xx(1) = x(1,n)
                    xx(2) = x(2,n)
                    xx(3) = x(3,n)
                    zp(1) = e1(1)*xx(1)
     &                    + e1(2)*xx(2)
     &                    + e1(3)*xx(3) - z0(1)
                    if(abs(zp(1)).le.dz(1)) then
                      zp(2) = e2(1)*xx(1)
     &                      + e2(2)*xx(2)
     &                      + e2(3)*xx(3) - z0(2)
                      if(abs(zp(2)).le.dz(2)) then
                        zp(3) = e3(1)*xx(1)
     &                        + e3(2)*xx(2)
     &                        + e3(3)*xx(3) - z0(3)
                        if(abs(zp(3)).le.dz(3)) then
                          kk        = onsurf(zp,zl,dz,gap,dnorm,4)
                          if(kk.eq.1) then
                            ip(n) = kk
                            norm(1,n) = e1(1)*dnorm(1)
     &                                + e2(1)*dnorm(2)
     &                                + e3(1)*dnorm(3)
                            norm(2,n) = e1(2)*dnorm(1)
     &                                + e2(2)*dnorm(2)
     &                                + e3(2)*dnorm(3)
                            norm(3,n) = e1(3)*dnorm(1)
     &                                + e2(3)*dnorm(2)
     &                                + e3(3)*dnorm(3)
                          end if
                        end if
                      end if
                    end if
                  endif ! ip(n).eq.0
                end do ! n

                xi1 = xi1 + dxi
              end do ! i
              xi2 = xi2 + dxi
            end do ! j

          endif
        end do ! while

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

c     Set region number for patch

      elseif( nsopt.eq.5 ) then

        regn = nint(scdat(1))

c     Set external node indicator

      elseif( nsopt.eq.3 ) then

        if(nint(scdat(1)).eq.0) then
          if(prt) write(iow,2002)
          extfl = .true.
        else
          if(prt) write(iow,2003)
          extfl = .false.
        endif

      endif

c     Formats

2000  format(
     &/5x,'C o n t a c t   S u r f a c e   S u p e r   N o d e s',//
     & '       1-SNode 2-SNode 3-SNode 4-SNode 5-SNode 6-SNode')

2001  format(i10,6i8/(10x,6i8))

2002  format(5x,'Search restricted to external surface nodes only.')

2003  format(5x,'Search over all nodes.')

      end
