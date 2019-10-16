c$Id:$
      subroutine ipenck1 (ndm,x,ns,ix2,cxs,ch2)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Correct computation for 'csi' and 'gt'           08/03/2013
c-----[--.----+----.----+----.-----------------------------------------]
c      Coded/modified by:                           Date:          rel.:
c               Giorgio Zavarise           April 10, 1996            1.0

c      Acronym: Initial PENetration ChecK

c      Purpose: Verify the contact status of the pair

c      Inputs :
c         ndm     - Space dimension of mesh
c         x(*)    - nodal coordinates
c         ns      - Slavenode global number
c         ix2(*)  - Element nodal connection list for surface 2
c         cxs(*)  - Slave node coordinates
c         ch2(*)  - Contact history variables (current)

c      Outputs:
c         x(*)    - nodal coordinates
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'c_0.h'
      include  'c_geom.h'
      include  'c_keyh.h'
      include  'c_tole.h'
      include  'iofile.h'
      include  'print.h'

      integer   ndm,ns,ix2(dnope2,*), masts,n1,n2
      real*8    d21,s21,c21,gt,csi,gn
      real*8    x(ndm,*),cxs(*),ch2(*), cx1(3),cx2(3)

      save

      call cdebug0 ('      ipenck1',-1)

c     Set main geom variables

      masts = nint(ch2(p1(1)))
      if(masts.gt.0) then
        n1     = ix2(1,masts)
        n2     = ix2(2,masts)
        cx1(1) = x(1,n1)
        cx1(2) = x(2,n1)
        cx2(1) = x(1,n2)
        cx2(2) = x(2,n2)

c       Length of master segment

        d21 = sqrt((cx2(1) - cx1(1))*(cx2(1) - cx1(1)) +
     &             (cx2(2) - cx1(2))*(cx2(2) - cx1(2)))

c       Tangent vector
c       s21 and c21 normalized later to minimize numerical errors

        s21 = (cx2(2) - cx1(2))
        c21 = (cx2(1) - cx1(1))

c       Distance S - 21  (positive when gap is open)

        gn = (cxs(1) - cx1(1)) * s21 + (cxs(2) - cx1(2)) * (-c21)
        gn = gn / d21

c       Check initial distance and move slave node if necessary

        if (gn.lt.(-tlipen)) then

c         Compute closest point location

          csi = ((cxs(1)-cx1(1))*c21 + (cxs(2)-cx1(2))*s21)/(d21*d21)
          gt  = csi

          if ((csi.gt.0.d0) .and. (csi.lt.1.d0)) then
            cxs(1) = cx1(1) + csi*(cx2(1)-cx1(1))
            cxs(2) = cx1(2) + csi*(cx2(2)-cx1(2))

c           Outputs

            if(ior.lt.0) then
              write (*,3000) ns
              write (*,3001) cxs(1),cxs(2),x(1,ns),x(2,ns)
            endif
            write (iow,3000) ns
            if (prt) then
              write (iow,3001) cxs(1),cxs(2),x(1,ns),x(2,ns)
            endif

c           Update

            x(1,ns) = cxs(1)
            x(2,ns) = cxs(2)
          endif
        endif
      endif

c     Formats

3000  format (/,'WARNING - Initial coord. changed for node ',i5)
3001  format (1x,'New coordinates',1p,2e20.10/
     +        1x,'Old coordinates',1p,2e20.10)

      end
