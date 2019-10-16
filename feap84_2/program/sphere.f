c$Id:$
      subroutine sphere(nty,x,ndm,prt,prth)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Converts spherical to cartesian coordinates

c      Inputs:
c         nty(*)    - Node type
c         x(ndm,*)  - Spherical coordinates of point
c         ndm       - Spatial dimension of mesh
c         prt       - Flag, output results if true
c         prth      - Flag, output title header if true

c      Outputs:
c         x(ndm,*)  - Cartesian coordinates of point
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'cdata.h'
      include  'cdat2.h'
      include  'crotas.h'
      include  'iofile.h'
      include  'pointer.h'
      include  'comblk.h'

      logical   prt,prth, errck, pinput
      integer   i,inc, mct, n,ne,ni,ndm,nty(*)
      real*8    r,sn2,cn2,sn3,cn3, x(ndm,*),td(6)

      save

      if(ndm.lt.3) then
        write(iow,4000)
        if(ior.lt.0) write(*,4000)
        return
      else
        mct = 0
100     if(ior.lt.0) then
          write(*,4001)
          call pprint('   >')
        endif
        errck = pinput(td,6)
        if(errck) go to 100
        ni  = nint(td(1))
        ne  = nint(td(2))
        inc = nint(td(3))
        if(ni.le.0) return
        if(ni.gt.numnp.or.ne.gt.numnp) go to 300
        inc = sign(max(abs(inc),1),ne-ni)
        if(ne.eq.0) ne = ni
        n = ni
200     if(nty(n).eq.0) then
          call pdegree(x(2,n), sn2,cn2)
          call pdegree(x(3,n), sn3,cn3)
          nty(n) = 2
          r      = x(1,n)
          x(1,n) = x0(1) + r*cn2*sn3 + td(4)
          x(2,n) = x0(2) + r*sn2*sn3 + td(5)
          x(3,n) = x0(3) + r*cn3     + td(6)

c         Set shell surface directors

          if(frotas) call sphdir(hr(np(82)),n,x0(1),x0(2),x0(3),sn3)

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
          write(iow,3001) n
          if(ior.lt.0) then
            write(iow,3001) n
          endif
        else
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
      endif

c     Error

300   if(ior.gt.0) then
        write(iow,3000) ni,ne
        call plstop()
      endif
      write(*,3000) ni,ne

c     Formats

2000  format('  Cartesian coordinates computed from spherical input.',
     & '  Global: x0 =',1p,1e12.4,' y0 =',1p,1e12.4,' z0 =',1p,1e12.4/
     & '  Local:  x0 =',1p,1e12.4,' y0 =',1p,1e12.4,' z0 =',1p,1e12.4/
     &   /4x,'Node',6(i6,'-Coord')/(8x,6(i6,'-Coord')))

2001  format(i8,6f12.4/(8x,6f12.4))

3000  format(' *ERROR* Attempt to convert nodes ni = ',i6,
     & ' to ne = ',i6)

3001  format(' *ERROR* Attempt to convert cartesian node',i9)

3002  format(' *ERROR* Attempt to convert undefined node',i9)

4001  format(' Input: node-1,node-2,inc, x0, y0, z0')

4000  format(' *ERROR* Attempt to convert spherical coordinates'/
     &       '         for a problem with less than 3-dimensions.')

      end
