c$Id:$
      subroutine polar(nty,x,ndm,prt,prth)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c     1. Use 'pi' from 'pconstant.h'                        14/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Convert polar/cylindrical to cartesian coordinates

c      Inputs:
c         nty(*)    - Node type
c         x(ndm,*)  - Polar coordinate values
c         ndm       - Spatial dimension of mesh
c         prt       - Output generated data if true
c         prth      - Output title/header data if true

c      Outputs:
c         x(ndm,*)  - Cartesian coordinate values
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'cdata.h'
      include  'cdat2.h'
      include  'crotas.h'
      include  'dstars.h'
      include  'iofile.h'
      include  'pconstant.h'
      include  'trdata.h'
      include  'pointer.h'
      include  'comblk.h'

      logical   prt,prth, errck, tinput, pcomp
      character type*15
      integer   i,n, inc,mct,ne,ni,ndm, nty(*)
      real*8    r,th, x(ndm,*),td(6), xx(3),xl(3)

      save

      if(ndm.eq.1) return
      mct = 0
      th = pi/180.0d0
100   if(ior.lt.0) then
        write(*,4000)
        call pprint('   >')
      endif
      errck = tinput(type,1,td,6)
      if(errck) go to 100
      if(pcomp(type,'node',4)) then
        if(nint(td(1)).le.0) return
        ni  = nint(td(1)) + starnd
        ne  = nint(td(2)) + starnd
        if(ni.gt.numnp.or.ne.gt.numnp) go to 300
        if(ne.eq.0) ne = ni
        inc = nint(td(3))
        inc = sign(max(abs(inc),1),ne-ni)
      elseif(pcomp(type,'all',3)) then
        ni  = 1
        ne  = numnp
        inc = 1
      else
        return
      endif
      n = ni
200   if(nty(n).eq.0) then
        nty(n) = 1
        xl(1) = x(1,n)
        xl(2) = x(2,n)
        if(ndm.eq.3) then
          xl(3) = x(3,n)
        else
          xl(3) = 0.0d0
        endif
        do i = 1,3
          xx(i)  = xr(i) + tr(i,1)*xl(1) + tr(i,2)*xl(2) + tr(i,3)*xl(3)
        end do ! i
        r      = xx(1)
        x(1,n) = x0(1) + r*cos(xx(2)*th) + td(4)
        x(2,n) = x0(2) + r*sin(xx(2)*th) + td(5)
        if(ndm.eq.3) x(3,n) = x0(3) + xx(3) + td(6)

c       Do shell surface modifications

        if(frotas) call poldir(hr(np(82)),n,x0(1),x0(2),th)

        if(mct.le.0) then
          if(prt) then
            call prtitl(prth)
            write(iow,2000) x0,td(4),td(5),td(6),(i,i=1,ndm)
            if(ior.lt.0) then
              write(*,2000) x0,td(4),td(5),td(6),(i,i=1,ndm)
            endif
          endif
          mct = 50
        endif
        if(prt) then
          write(iow,2001) n,(x(i,n),i=1,ndm)
          if(ior.lt.0) write(*,2001) n,(x(i,n),i=1,ndm)
        endif
        mct = mct - 1
      elseif(nty(n).gt.0) then
        write(ilg,3001) n
        write(iow,3001) n
        if(ior.lt.0) then
          write(iow,3001) n
        endif
      else
        write(ilg,3002) n
        write(iow,3002) n
        if(ior.lt.0) then
          write(iow,3002) n
        endif
      endif
      n = n + inc
      if((ne-n)*inc.ge.0) go to 200
      if(mod(ne-ni,inc).eq.0) go to 100
      ni = ne
      n = ne
      go to 200

c     Error

300   write(ilg,3000) ni,ne
      write(iow,3000) ni,ne
      if(ior.lt.0) write(*,3000) ni,ne
      call plstop()

c     Formats

2000  format('  Cartesian coordinates computed from polar inputs.'//
     & '  Global: x0 =',1p,1e12.4,' y0 =',1p,1e12.4,' z0 =',1p,1e12.4/
     & '  Local:  x0 =',1p,1e12.4,' y0 =',1p,1e12.4,' z0 =',1p,1e12.4/
     &   /4x,'Node',6(i6,'-Coord')/(8x,6(i6,'-Coord')))

2001  format(i8,6f12.4)

3000  format(' *ERROR* POLAR: Attempt to convert nodes ni =',i8,
     &       ' to ne =',i8)

3001  format(' *ERROR* POLAR: Attempt to convert cartesian node',i9)

3002  format(' *ERROR* POLAR: Attempt to convert undefined node',i9)

4000  format(' Input: node-1,node-2,inc, x0, y0, z0')

      end
