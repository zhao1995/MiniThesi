c$Id:$
      subroutine crblok(nope,dnope,ics,neps,nsopt,scdat,polfl,
     &                  ix,ip,ep,xin,x,ndm,nen,nen1,numnp,numel)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Coded/modified by:                           Date:          rel.:
c               Robert L. Taylor           April 10, 1996            1.0

c      Acronym: Contact Read BLOcK

c      Purpose: Automatic generation of contact surfaces

c      Inputs:
c         nope    - # of NOde Per Element
c         dnope   - Dimension of NOde Per Element
c         nsopt   - Subcommand #: 1=gap,2=segm,3=pola,4=cart
c         scdat(*)- Subcommand data
c         ix(*)   - Element nodal connection list
c         x(*)    - nodal coordinates
c         ndm     - Space dimension of mesh
c         nen     - Maximum number of nodes/element
c         nen1    - Dimension of ix array
c         numnp   - Number of nodes in mesh
c         numel   - Number of elements in mesh

c      Scratch:
c         ip(*)   - Nodal integer list storage
c         ep(*)   - Element list storage
c         xin(*)  - Nodal real list storage

c      Outputs:
c         ics(*)  - Contact element nodal connection array
c         neps    - # of Elements Per Surface
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'iofile.h'

      logical   errck,gettd, polfl, wdflg
      integer   nope,dnope,neps,nsopt, i,i1,i2,j,m,n,nlast,regn
      integer   ndm,nen,nen1,numnp,numel,jsw
      integer   ics(dnope,*),ix(nen1,*),ip(*),ep(numel,*),nend(2,2)
      real*8    cn,sn, tol,tol0,tolxi,d,gap0(2),ximin,ximax
      real*8    scdat(14),xin(*),x(ndm,*)
      real*8    xl(3,3), x0(3),xld(2),xqd(2), td(30)

      save

      data      tol0 / 1.d-5/, tolxi / 1.0d-8 /
      data      gap0 /0.0d0,1.d0/
      data      regn /0/

c     Set gap value for search

      if( nsopt.eq.1) then
        gap0(1) = scdat(1)
        gap0(2) = scdat(2)
        if(gap0(2).eq.0.0d0) then
          gap0(2) = 1.0d0
        endif
        write(iow,*) '  GAP =',gap0
        if(ior.lt.0) then
          write(*,*) '  GAP =',gap0
        end if

c     Set segment end and center coordinate values

      elseif( nsopt.eq.2 ) then

c       Set coordinates for surface search

        write(iow,2000)
        jsw   = 1
        wdflg = .true.
        do while (wdflg)
          errck = gettd(td,3,'skip')
          if(errck) then
            backspace ior
            wdflg = .false.
          else
            i = nint(td(1))
            if(i.ge.1 .and. i.le.3) then
              xl(1,i) = td(2)
              xl(2,i) = td(3)
              xl(3,i) = 0.0d0
              if(i.eq.3) jsw = 2
              write(iow,2001) i,xl(1,i),xl(2,i)
            else
              wdflg = .false.
            endif
          endif
        end do

c       Set segment data for cartesian linear edge

        if(jsw.eq.1) then
          xl(1,3) = 0.5d0*(xl(1,1) + xl(1,2))
          xl(2,3) = 0.5d0*(xl(2,1) + xl(2,2))
          xl(3,3) = 0.0d0
        endif

c       Projections nodes to points

        errck = .false.
        tol   =  tol0

100     call prj2dl(gap0,tol0,tol,x0,ip,xin,x,xl,ndm,numnp,polfl)

c       Check if end points found

        ximin = +1.d0
        ximax = -1.d0
        do n = 1,numnp
          if(xin(n).ne.0.0d0) then
            ximin = min(ximin,xin(n))
            ximax = max(ximax,xin(n))
          end if
        end do

c       For cylindrical coordinates check number of 1 and -1 values

        if(polfl) then
          i1 = 0
          i2 = 0
          do n = 1,numnp
            if(xin(n).ge.  1.d0 - tol0) then
              i1 = i1 + 1
            endif
            if(xin(n).le. -1.d0 + tol0) then
              i2 = i2 + 1
            endif
          end do ! n
          if(i1.ge.2) then
            ximin = -1.d0
          endif
          if(i2.ge.2) then
            ximax =  1.d0
          endif
        endif

        tol = tol0 + max(1.d0+ximin,1.d0-ximax,0.0d0)

        errck = .not.errck
        if(errck .and. tol.gt.tol0+tolxi) go to 100

c       Loop through elements to find adjacent nodes

        do n = 1,numel
          if(ix(nen1-1,n).eq.regn .or. regn.eq.0) then
            ep(n,1) = 0
            ep(n,2) = 0
            ep(n,3) = 0
            do j = 1,min(nen,4)
              i1 = ix(j,n)
              if(i1.ne.0) then
                if(ip(i1).ne.0) then
                  if(j.lt.min(nen,4)) then
                    i2 = ix(j+1,n)
                    if(i2.eq.0) then
                      i2 = ix(1,n)
                    endif
                  else
                    i2 = ix(1,n)
                  endif
                  if(ip(i2).gt.0 .and. i1.ne.i2) then

c                   Check if already assigned on 2-node element type.

                    if(ep(n,1).le.0 .and. ep(n,2).le.0) then
                      ep(n,1) = i1
                      ep(n,2) = i2
                      if(j+4.le.nen) then
                        if(ip(ix(j+4,n)).gt.0) ep(n,3) = ix(j+4,n)
                      endif
                    endif
                  endif
                endif
              endif
            end do
          endif
        end do

c       Remove duplicates

        do n = 1,numel
          if(ix(nen1-1,n).eq.regn .or. regn.eq.0) then
            if(ep(n,1).gt.0) then
              do m = n+1,numel
                if(ix(nen1-1,m).eq.regn .or. regn.eq.0) then
                  if(ep(m,1).eq.ep(n,1) .and. ep(m,2).eq.ep(n,2) .or.
     &               ep(m,1).eq.ep(n,2) .and. ep(m,2).eq.ep(n,1)) then
                    ep(m,1) = 0
                    ep(m,2) = 0
                    ep(m,3) = 0
                  end if
                end if
              end do
            end if
          end if
        end do

c       Compute surface segments

c       Find end points

        do n = 1,numnp
          ip(n) = 0
        end do

        xld(1) = 0.5d0*(xl(1,2) - xl(1,1))
        xld(2) = 0.5d0*(xl(2,2) - xl(2,1))
        xqd(1) = xl(1,1) + xl(1,2) - 2.d0*xl(1,3)
        xqd(2) = xl(2,1) + xl(2,2) - 2.d0*xl(2,3)
        do n = 1,numel
          if(ix(nen1-1,n).eq.regn .or. regn.eq.0) then
            if(ep(n,1).gt.0) then

c             Order segments along master direction

              if(polfl) then
                d  = xl(2,3) + xin(ep(n,1))*(xld(2)
     &                       + 0.5d0*xin(ep(n,1))*xqd(2))
                call pdegree(d, cn,sn)
                cn = -cn
              else
                cn = xld(1) + xin(ep(n,1))*xqd(1)
                sn = xld(2) + xin(ep(n,1))*xqd(2)
c               d  = 1.d0/sqrt(cn*cn + sn*sn)
c               sn = sn*d
c               cn = cn*d
              endif

              d = (x(1,ep(n,2)) - x(1,ep(n,1)))*cn
     &          + (x(2,ep(n,2)) - x(2,ep(n,1)))*sn

              if(polfl.and.xld(2).lt.0.0d0) then
                d = -d
              endif

c             Delete edges going wrong direction

              if(d.lt.0.0d0) then
                ep(n,1) = 0
                ep(n,2) = 0
                ep(n,3) = 0

c             Count occurances of node

              else
                ip(ep(n,1)) = ip(ep(n,1)) + 1
                ip(ep(n,2)) = ip(ep(n,2)) + 1
              end if
            end if
          end if
        end do

c       Set end point array

        j         = 0
        nend(1,1) = 0
        nend(2,1) = 0
        do n = 1,numnp
          if(ip(n).eq.1) then
            j         = j + 1
            nend(j,1) = n
          end if
        end do

c       Check for error

        if(nend(1,1).eq.0 .or. nend(2,1).eq.0) then
          write(ilg,3000)
          write(iow,3000)
          if(ior.lt.0) then
            write(*,3000)
          endif
        endif

c       Output lists

        nlast = 0
300     n     = 0
        ximin = 2.0d0

        do i = 1,numel
          i1 = ep(i,1)
          if(i1.gt.0) then
            if(xin(i1).lt.ximin) then
              ximin = xin(i1)
              n     = i
            end if
          end if
        end do

        if(n.gt.0) then

          neps = neps + 1
          ics(1,neps) = ep(n,1)
          if(nope.eq.1) nlast       = n
          if(nope.gt.1) ics(2,neps) = ep(n,2)
          if(nope.gt.2) ics(3,neps) = ep(n,3)

          xin(ep(n,1)) = 3.0d0

          go to 300

c       Set last point for 'point' segments

        elseif(nlast.gt.0) then

          neps        = neps + 1
          ics(1,neps) = ep(nlast,2)

        end if

c     Specify that segment coordinates are input in polar form

      elseif( nsopt.eq.3 ) then
        polfl = .true.
        do i = 1,2
          x0(i) = scdat(i)
        end do
        x0(3) = 0.0d0

c     Specify that segment coordinates are input in cartesian form

      elseif( nsopt.eq.4 ) then

        polfl = .false.

c     Specify region number to restrict search

      elseif( nsopt.eq.5 ) then

        regn = nint(scdat(1))

      endif

c     Formats

2000  format(/5x,'C o n t a c t    S u r f a c e    ',
     & 'C o o r d i n a t e s'//'      Node    X-Coord.    Y-Coord.')

2001  format(i10,1p,2e13.4)

3000  format(' *ERROR* CRBLOK: No surface located')

      end
