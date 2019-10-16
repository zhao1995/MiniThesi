c$Id:$
      subroutine gloscln (ndm,ndf,x,u,ix2,ns,cxs,ch2)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Coded/modified by:                           Date:          rel.:
c               Giorgio Zavarise           April 10, 1996           1.00
c               Robert L. Taylor         October 12, 1996           1.01

c      Acronym: GLObal Search of CLosest Node

c      Purpose: Search closest node on whole master surface

c      Inputs :
c         ndm     - Space dimension of mesh
c         ndf     - Number dof/node
c         x(*)    - nodal coordinates
c         u(*)    - Current nodal solution vectors
c         ix2(*)  - Element nodal connection list for surface 2
c         ns      - Node number for slave node
c         cxs(*)  - Coordinates of slavenode S

c      Outputs:
c         ch2(*)  - Contact history variables (current)
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'c_0.h'
      include  'c_geom.h'
      include  'c_keyh.h'
      include  'c_tole.h'

      logical   firsf
      integer   ndm,ndf,ix2(dnope2,*), ns,nm,nel2,nod2,ec,lnc
      real*8    x(ndm,*),u(ndf,*),cxs(*),ch2(*),cxm(3),dms,dcsmin

      save

      call cdebug0 ('      gloscln',-1)

      firsf = .true.

c    Loop over  surface 2 elements starting from last one

      do nel2 = 1,neps2
        do nod2 = 1,nope2

          if((nod2.eq.1) .or. (ix2(dnope2,nel2).eq.0)) then

c           Set master node number and coordinates

            nm     = ix2(nod2,nel2)

c           If not same as slave number check

            if(nm.ne.ns) then
              cxm(1) = x(1,nm) + u(1,nm)
              cxm(2) = x(2,nm) + u(2,nm)

c             Compute distance

              dms = sqrt((cxs(1) - cxm(1))*(cxs(1) - cxm(1)) +
     &                   (cxs(2) - cxm(2))*(cxs(2) - cxm(2)))

c             Set first minimal distance

              if (firsf) then
                firsf  = .false.
                dcsmin = dms
                ec     = nel2
                lnc    = nod2
              endif

c             Update minimal distance

              if (dms.lt.dcsmin) then
                dcsmin = dms
                ec     = nel2
                lnc    = nod2
              endif
            endif ! ns .ne. nm
          endif ! one node/segment except last
        end do ! nope2
      end do ! nel2

c     Store values

      ch2(p1(1)) = ec
      ch2(p1(2)) = lnc

      end
