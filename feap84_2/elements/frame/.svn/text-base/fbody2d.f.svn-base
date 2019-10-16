c$Id:$
      subroutine fbody2d(d,xl,ul, r,s, ndm,ndf,nst, isw)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Add fixed end moments for cubic form;            10/03/2007
c          Follower and normal loading for 3-node element.
c       2. Add fixed end moments for gravity load term      05/09/2008
c-----[--.----+----.----+----.-----------------------------------------]
c     Purpose: Compute body force loads

c     Inputs:
c       d(*)     - Loading parameters
c       xl(ndm,*) - Reference element coordinates
c       ul(ndf,*) - Element displacements
c       ndm       - Mesh spatial dimension
c       ndf       - Degree of freedoms/node
c       nst       - Element matrix dimension
c       isw       - Switch option

c     Outputs:
c       r(ndf,*)  - Nodal loading vector
c       s(nst,*)  - Loading tangent matrix (follower load only)
c-----[--.----+----.----+----.-----------------------------------------]
      implicit   none

      include   'elbody.h'
      include   'eldata.h'
      include   'pconstant.h'

      integer    ndm,ndf,nst,isw, ii,ll,lint, nd2
      real*8     d(*),xl(ndm,*),ul(ndf,*),r(ndf,*),s(nst,nst)
      real*8     shp(2,3), sg(2,3), body(3)
      real*8     le,dx,ddx,dy,ddy,pn

c     Set body force loading factors

      if(isw.eq.15) then
        body(1) = bodyf(1)
        body(2) = bodyf(2)
        pn     = 0.5d0*d(10)
      else
        call sbodyf(d, body)
        pn     = 0.5d0*d(10)*dm
      endif

c     2-node element

      if(nel.eq.2) then

        le     = sqrt((xl(1,2)-xl(1,1))**2
     &              + (xl(2,2)-xl(2,1))**2)*0.5d0
        dx     = (xl(1,2) - xl(1,1))
     &         + (ul(1,2) - ul(1,1))*d(68)
        dy     = (xl(2,2) - xl(2,1))
     &         + (ul(2,2) - ul(2,1))*d(68)
        r(1,1) = r(1,1)   + body(1)*le - pn*dy
        r(2,1) = r(2,1)   + body(2)*le + pn*dx
        r(1,2) = r(1,2)   + body(1)*le - pn*dy
        r(2,2) = r(2,2)   + body(2)*le + pn*dx

c       Compute fixed end moments for cubic displacements

        if(d(79).gt.0.0d0) then

c         x and y body loading

          dx     = (xl(1,2) - xl(1,1))*abs(xl(1,2) - xl(1,1))*one12
          dy     = (xl(2,2) - xl(2,1))*abs(xl(2,2) - xl(2,1))*one12
          r(3,1) = r(3,1)   + body(2)*dx - body(1)*dy + pn*le*le*two3
          r(3,2) = r(3,2)   - body(2)*dx + body(1)*dy - pn*le*le*two3

        endif ! End fixed end moment

c       Tangent for follower loading (unsymmetric)

        if(d(68).gt.0.0d0) then
          s(1    ,2    ) = -pn
          s(1    ,2+ndf) =  pn
          s(2    ,1    ) =  pn
          s(2    ,1+ndf) = -pn
          s(1+ndf,2    ) = -pn
          s(1+ndf,2+ndf) =  pn
          s(2+ndf,1    ) =  pn
          s(2+ndf,1+ndf) = -pn
        endif

c     3-node element

      else
        pn   = pn*2.d0
        lint = nel
        if(nint(d(182)).gt.0) then
          call int1dn(lint, sg)
        else
          call int1d(lint, sg)
        endif
        do ll = 1,lint
          call shp1d(sg(1,ll),xl,shp,ndm,nel,dx)
          dx = sg(2,ll)*dx
          do ii = 1,nel
            r(1,ii)   = r(1,ii) + body(1)*shp(2,ii)*dx
            r(2,ii)   = r(2,ii) + body(2)*shp(2,ii)*dx
          end do ! ii
        end do ! ll

c       Follower loading

        if(d(68).ne.0.0d0) then

          dx     = (xl(1,2) - xl(1,1))
     &           + (ul(1,2) - ul(1,1))
          ddx    = (xl(1,1) + xl(1,2) - 2.d0*xl(1,3))
     &           + (ul(1,1) + ul(1,2) - 2.d0*ul(1,3))
          dy     = (xl(2,2) - xl(2,1))
     &           + (ul(2,2) - ul(2,1))
          ddy    = (xl(2,1) + xl(2,2) - 2.d0*xl(2,3))
     &           + (ul(2,1) + ul(2,2) - 2.d0*ul(2,3))

          nd2            =  ndf + ndf
          s(1    ,2    ) = -one2*pn
          s(1    ,2+ndf) = -one6*pn
          s(1    ,2+nd2) =  two3*pn
          s(2    ,1    ) =  one2*pn
          s(2    ,1+ndf) =  one6*pn
          s(2    ,1+nd2) = -two3*pn

          s(1+ndf,2    ) =  one6*pn
          s(1+ndf,2+ndf) =  one2*pn
          s(1+ndf,2+nd2) = -two3*pn
          s(2+ndf,1    ) = -one6*pn
          s(2+ndf,1+ndf) = -one2*pn
          s(2+ndf,1+nd2) =  two3*pn

          s(1+nd2,2    ) = -two3*pn
          s(1+nd2,2+ndf) =  two3*pn
          s(2+nd2,1    ) =  two3*pn
          s(2+nd2,1+ndf) = -two3*pn
        else
          dx     = (xl(1,2) - xl(1,1))
          ddx    = (xl(1,1) + xl(1,2) - 2.d0*xl(1,3))
          dy     = (xl(2,2) - xl(2,1))
          ddy    = (xl(2,1) + xl(2,2) - 2.d0*xl(2,3))
        endif

c       Normal load

        r(1,1) = r(1,1) - pn*(one6*dy - one3*ddy)
        r(2,1) = r(2,1) + pn*(one6*dx + one3*ddx)
        r(1,2) = r(1,2) - pn*(one6*dy + one3*ddy)
        r(2,2) = r(2,2) + pn*(one6*dx - one3*ddx)
        r(2,3) = r(2,2) - pn*(two3*dx)
        r(1,3) = r(1,2) + pn*(two3*dy)

      endif

      end
