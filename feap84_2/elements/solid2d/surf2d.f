c$Id:$
      subroutine surf2d(d,ul,xl,ma,ndf,ndm,nel,mct,nst, p,s)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose:  Surface Loading traction load computations

c     Input:
c        d(mct)       - Pressure parameters
c        ul(ndf,nel)  - Nodal displacements
c        xl(ndf,nel)  - Reference Nodal coordinates
c        ma           - Loading Type
c                       (ma.eq.1) Constant Pressure
c        ndf          - Degrees of freedom per node
c        ndm          - Spatial dimensions of mesh
c        nel          - Number of nodes on element surface
c        mct          - Number of parameter to define loading
c        nst          - nel*ndf

c     Output:
c        p(nst)       - Residual
c        s(nst,nst)   - Tangent
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'iofile.h'
      include  'pmod2d.h'

      integer   ma,ndf,ndm,nel,mct,nst
      real*8    dx,dy,dpn,dpt,dx31,dx32,dn1,dn2,dt1,dt2
      real*8    d(6),ul(ndf,*),xl(ndm,*),p(*),s(nst,*)

      save

c     Small Deformation load

      if(dtype.gt.0) then
        if(ma.eq.1) then

          dx = xl(1,1) - xl(1,2)
          dy = xl(2,2) - xl(2,1)
          if(mct.eq.1) then
            p(1)     = d(1)*dy/2.
            p(2)     = d(1)*dx/2.
            p(ndf+1) = p(1)
            p(ndf+2) = p(2)
          elseif(mct.eq.2) then
            p(1)     = dy*(d(1)/3. + d(2)/6.)
            p(2)     = dx*(d(1)/3. + d(2)/6.)
            p(ndf+1) = dy*(d(1)/6. + d(2)/3.)
            p(ndf+2) = dx*(d(1)/6. + d(2)/3.)
          elseif(mct.ge.4) then
            dpn      = d(1)/3. + d(2)/6.
            dpt      = d(3)/3. + d(4)/6.
            p(1)     = dy*dpn - dx*dpt
            p(2)     = dx*dpn + dy*dpt
            dpn      = d(2)/3. + d(1)/6.
            dpt      = d(4)/3. + d(3)/6.
            p(ndf+1) = dy*dpn - dx*dpt
            p(ndf+2) = dx*dpn + dy*dpt
          endif
          if(nel .gt. 2) then
            p(2*ndf+1) = 2.*(p(1) + p(ndf+1))/3.
            p(2*ndf+2) = 2.*(p(2) + p(ndf+2))/3.
            p(1)       = p(1) - 0.5*p(2*ndf+1)
            p(2)       = p(2) - 0.5*p(2*ndf+2)
            p(ndf+1)   = p(ndf+1) - 0.5*p(2*ndf+1)
            p(ndf+2)   = p(ndf+2) - 0.5*p(2*ndf+2)
         endif

        elseif(ma.eq.2) then

          dx         = (xl(1,1)-xl(1,2))/6.
          dy         = (xl(2,2)-xl(2,1))/6.
          d(5)       = d(5) -0.5*(d(1) + d(2))
          d(6)       = d(6) - 0.5*(d(3) + d(4))
          dx31       = -2.*(xl(2,3)-0.5*(xl(2,1)+xl(2,2)))/3.0
          dx32       =  2.*(xl(1,3)-0.5*(xl(1,1)+xl(1,2)))/3.0
          dn1        = d(1) + 0.4*d(5)
          dt1        = d(3) + 0.4*d(6)
          dn2        = d(5) + 2*d(1) + 0.5*d(2)
          dt2        = d(6) + 2*d(3) + 0.5*d(4)
          p(1)       = dn1*dy - dt1*dx - 0.4*(dn2*dx31 - dt2*dx32)
          p(2)       = dn1*dx + dt1*dy - 0.4*(dt2*dx31 + dn2*dx32)
          dn1        = d(2) + 0.4*d(5)
          dt1        = d(4) + 0.4*d(6)
          dn2        = d(5) + 2*d(2) + 0.5*d(1)
          dt2        = d(6) + 2*d(4) + 0.5*d(3)

          p(ndf+1)   = dn1*dy - dt1*dx + 0.4*(dn2*dx31 - dt2*dx32)
          p(ndf+2)   = dn1*dx + dt1*dy + 0.4*(dt2*dx31 + dn2*dx32)

          dn1        = d(1) + d(2) + 1.6*d(5)
          dt1        = d(3) + d(4) + 1.6*d(6)

          p(2*ndf+1) = 2.*(dn1*dy - dt1*dx)
     &               + 0.4*((d(2)-d(1))*dx31 - (d(4)-d(3))*dx32)
          p(2*ndf+2) = 2.*(dn1*dx + dt1*dy)
     &               + 0.4*((d(4)-d(3))*dx31 + (d(2)-d(1))*dx32)

        endif

c     Finite Deformation follower load

      elseif(dtype.lt.0) then

        if(ma.ne.1) then

          if(ior.lt.0) then
            write(*,*) 'ERROR: Surface load type',ma,' unavailable.'
          endif
          write(iow,*) 'ERROR: Surface load type',ma,' unavailable.'
          return

        else

          if(mct.ne.1) then
            if(ior.lt.0) then
              write(*,*) 'ERROR:',mct,' surface parameter not allowed.'
            endif
            write(iow,*) 'ERROR:',mct,' surface parameter not allowed.'
            return
          endif

          dn1      = d(1)*0.5d0

          p(1)     =  dn1*(xl(2,2)+ul(2,2) - xl(2,1) - ul(2,1))
          p(2)     =  dn1*(xl(1,1)+ul(1,1) - xl(1,2) - ul(1,2))
          p(ndf+1) =  p(1)
          p(ndf+2) =  p(2)

          s(1,2)   =  dn1
          s(1,4)   = -dn1
          s(2,1)   = -dn1
          s(2,3)   =  dn1
          s(3,2)   =  dn1
          s(3,4)   = -dn1
          s(4,1)   = -dn1
          s(4,3)   =  dn1
        endif

      endif

      end
