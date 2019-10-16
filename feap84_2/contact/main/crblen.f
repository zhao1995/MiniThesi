c$Id:$
      subroutine crblen(nope,dnope,ics,neps,nsopt,scdat,xs,
     &                  ix,ip,x)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Add 'oname' and choice by word type 'cart', etc. 11/01/2007
c-----[--.----+----.----+----.-----------------------------------------]
c      Coded/modified by:                           Date:          rel.:
c               Robert L. Taylor           March 25, 1997            1.0

c      Acronym: Contact Read BLOcK

c      Purpose: Automatic blending generation of contact surfaces

c      Inputs:
c         nope    - # of NOde Per Element
c         dnope   - Dimension of NOde Per Element
c         nsopt   - Subcommand #: 1=gap,2=segm,3=pola,4=cart
c         scdat(*)- Subcommand data
c         xs(3,*) - Super node coordinates
c         ix(*)   - Element nodal connection list
c         x(*)    - nodal coordinates

c      Scratch:
c         ip(*)   - Nodal integer list storage

c      Outputs:
c         ics(*)  - Contact element nodal connection array
c         neps    - # of Elements Per Surface
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'cblend.h'
      include  'cdata.h'
      include  'iofile.h'
      include  'sdata.h'

      character oname(2)*15
      logical   gettxtd, wdflg,gapfl, errck, pcomp
      integer   nope,dnope,neps,nsopt, i,i1,i2,i2last,j,n,jsw, styp
      integer   ics(dnope,*),ix(nen1,*),ip(*), is(15)
      real*8    y1a,y2a,y1b,y2b, tolgp,d,gap,gap0,s
      real*8    scdat(14),x(ndm,*),xs(3,*),tr(3,4), xx(2),dx(2),td(30)
      real*8    shp(181),xr(2,181),rt(3,13)

      save

      data      tolgp / 1.0d-3 /
      data      gap0 /0.0d0/
      data      tr   /1.d0,3*0.d0,1.d0,3*0.d0,1.d0,3*0.d0/

c     Coordinate contact surface descriptions

c     Set gap value for search

      if( nsopt.eq.1) then
        gap0 = scdat(1)
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
          errck = gettxtd(oname,1,td,15,'skip')
          if(pcomp(oname(1),'cart',4)) then
            styp = 0
          elseif(pcomp(oname(1),'pola',4)) then
            styp = 1
          elseif(pcomp(oname(1),'segm',4)) then
            styp = 2
          elseif(pcomp(oname(1),'elip',4)) then
            styp = 3
          elseif(.not.pcomp(oname(1),'    ',4)) then
            write(iow,3000) oname(1)
            call plstop()
          endif

          if(.not.errck) then
            i2 = 0
            do i = 1,15
              is(i) = int(td(i))
              if(is(i).gt.0) i2 = i
            end do
            if(i2.eq.0) errck = .true.
          endif
          if(errck) then
            backspace ior
            wdflg = .false.
          else
            write(iow,2001) oname(1)(1:5),(is(i),i=1,i2)
            j = 180
            call pside1(j,xs,tr,1,is,i2,ndm,shp,rt,xr, styp)

c           Compute bounding box

            y1a = xr(1,1)
            y1b = xr(1,1)
            y2a = xr(2,1)
            y2b = xr(2,1)
            do i = 2,j+1
              y1a = min(y1a,xr(1,i))
              y2a = min(y2a,xr(2,i))
              y1b = max(y1b,xr(1,i))
              y2b = max(y2b,xr(2,i))
            end do
            y1b = 0.5d0*sqrt((y1b - y1a)**2 + (y2b - y2a)**2)

c           Do a double pass if necessary to get end points

            gap = max(tolgp*y1b,gap0)

c           Find nodes on surface (within gap)

            do n = 1,numnp
              ip(n) = 0
            end do

            do n = 1,numnp

c             Compute location of closest point on surface

              gapfl = .true.
              i     = 1
              do while(gapfl .and. i.le.j)

                xx(1) = xr(1,i  )
                xx(2) = xr(2,i  )
                dx(1) = xr(1,i+1) - xx(1)
                dx(2) = xr(2,i+1) - xx(2)
                d =  dx(1)**2 + dx(2)**2
                s = ((x(1,n)-xx(1))*dx(1) + (x(2,n)-xx(2))*dx(2))/d
                if(abs(s-0.5d0).le.0.5d0+1.d-5) then
                  xx(1) = xx(1) + s*dx(1)
                  xx(2) = xx(2) + s*dx(2)
                  d     = sqrt((x(1,n)-xx(1))**2 + (x(2,n)-xx(2))**2)
                  gapfl = d.gt.gap
                endif

                i = i + 1
              end do
              if(gapfl) then
                ip(n) = 0
              else
                ip(n) = i
              endif
            end do
          endif
        end do

c       Set up contact segments list

        i2last = 0
        do n = 1,numel
          do j = 1,min(nen,4)
            i1 = ix(j,n)
            if(i1.ne.0) then
              if(ip(i1).gt.0) then
                if(j.lt.min(nen,4)) then
                  i2 = ix(j+1,n)
                  if(i2.eq.0) then
                    i2 = ix(1,n)
                  endif
                else
                  i2 = ix(1,n)
                endif
                if(ip(i2).gt.0 .and. i1.ne.i2) then

c                 Check direction cosine for tangents

                  d = (x(1,i2)-x(1,i1))*(xr(1,ip(i2))-xr(1,ip(i1)))
     &              + (x(2,i2)-x(2,i1))*(xr(2,ip(i2))-xr(2,ip(i1)))

c                 Insert nodes into list

                  if(d.gt.0.0d0) then
                    neps = neps + 1
                    ics(1,neps) = i1
                    if(nope.eq.1)  i2last      = i2
                    if(nope.gt.1 ) ics(2,neps) = i2
                    if(nope.gt.2 ) ics(3,neps) = ix(j+4,n)
                  endif
                endif
              endif
            endif
          end do
        end do

c       Set last value for 'point' surfaces

        if(i2last.gt.0) then
          neps        = neps + 1
          ics(1,neps) = i2last
        endif

      endif

c     Formats

2000  format(/5x,'C o n t a c t    S u r f a c e    S i d e s',//
     & '    Type    1-SNode 2-SNode 3-SNode 4-SNode 5-SNode 6-SNode')

2001  format(5x,a5,6i8/(10x,6i8))

3000  format(' ->ERROR: Incorrect contact blend option ',a/
     &       '     USE: CARTesian, POLAr, SEGMent, ELLIpse')

      end
