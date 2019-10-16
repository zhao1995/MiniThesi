c$Id:$
      subroutine mastseg (ndm,ndf,x,u,ix2,ns,cxs,ch2)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Coded/modified by:                           Date:          rel.:
c               Giorgio Zavarise           April 10, 1996            1.0

c      Acronym: MASTer SEGment

c      Purpose: Find master segment per each contact node

c      Inputs :
c         ndf     - Number dof/node
c         ndm     - Space dimension of mesh
c         x(*)    - nodal coordinates
c         u(*)    - Current nodal solution vectors
c         ix2(*)  - Element nodal connection list for surface 2
c         ns      - Slave node number
c         cxs(*)  - Coordinates of slavenode S

c      Outputs:
c         ch2(*)  - Contact history variables (current)
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'c_0.h'
      include  'c_geom.h'
      include  'c_keyh.h'

      integer   ndm,ndf,ns,ec,lnc,ix2(dnope2,*)
      integer   masts,nc,ec1,ec2,n1,n2,istgt
      real*8    x(ndm,*),u(ndf,*),cxs(*),ch2(*)
      real*8    dot1c,dot2c,cx1(3),cx2(2),cxc(3)

      save

      call cdebug0 ('      mastseg',-1)

c     Get data

      ec  = nint(ch2(p1(1)))
      lnc = nint(ch2(p1(2)))

c     Check if only one segment

      if(ix2(dnope2-1,1).eq.0 .and. ix2(dnope2,1).eq.0) then

c       Global nodes

        n1  = ix2(1,1)
        n2  = ix2(2,1)

c       Coordinates

        cx1(1) = x(1,n1) + u(1,n1)
        cx1(2) = x(2,n1) + u(2,n1)
        cx2(1) = x(1,n2) + u(1,n2)
        cx2(2) = x(2,n2) + u(2,n2)

c       Scalar product to find the right segment

        dot1c = (cx2(1) - cx1(1)) * (cxs(1) - cx1(1)) +
     &          (cx2(2) - cx1(2)) * (cxs(2) - cx1(2))

        if(dot1c.gt.0.0d0) then
          istgt = 1
        endif
        ch2(p1(1)) = 1
        ch2(p1(2)) = lnc
        ch2(p1(3)) = istgt
        return
      endif

c     Find two adjacent segments

      if(lnc.eq.1) then
        if(ix2(dnope2-1,ec).eq.0) then
          ec1 = ec
        else
          ec1 = ix2(dnope2-1,ec)
        endif
        ec2 = ec
      else
        ec1 = ec
        if(ix2(dnope2  ,ec).eq.0) then
          ec2 = ec
        else
          ec2 = ix2(dnope2  ,ec)
        endif
      endif

c     Determine global nodes

      nc  = ix2(lnc,ec)
      n1  = ix2(1,ec1)
      n2  = ix2(2,ec2)

c     No contact for single surface state

      if(ns.eq.nc .or. ns.eq.n1 .or. ns.eq.n2) then

        masts = 0
        lnc   = 0
        istgt = 0

      else

c       Coordinates

        cx1(1) = x(1,n1) + u(1,n1)
        cx1(2) = x(2,n1) + u(2,n1)
        cx2(1) = x(1,n2) + u(1,n2)
        cx2(2) = x(2,n2) + u(2,n2)
        cxc(1) = x(1,nc) + u(1,nc)
        cxc(2) = x(2,nc) + u(2,nc)

c       Scalar product to find the right segment

        dot1c = (cx1(1) - cxc(1)) * (cxs(1) - cxc(1)) +
     &          (cx1(2) - cxc(2)) * (cxs(2) - cxc(2))
        dot2c = (cx2(1) - cxc(1)) * (cxs(1) - cxc(1)) +
     &          (cx2(2) - cxc(2)) * (cxs(2) - cxc(2))

c       Limit case mastersegment is # left

        if ((ec.eq.ec1) .and. (lnc.eq.1)) then
          masts = ec
          if (dot2c.ge.0.d0) then
            istgt = 1
          else
            istgt = 4
          endif

c       Limit case mastersegment is # right

        elseif ((ec.eq.ec2) .and. (lnc.eq.2)) then
          masts = ec
          if (dot1c.ge.0.d0) then
            istgt = 1
          else
            istgt = 5
          endif

c       Normal case - find segments related to the closest node
c       Master segment is 1-C

        elseif ((dot1c.ge.0.d0) .and. (dot2c.le.0.d0)) then
          masts = ec1
          istgt = 1

c       Master segment is C-2

        elseif ((dot1c.le.0.d0) .and. (dot2c.ge.0.d0)) then
          masts = ec2
          istgt = 1

c       Corner position - no master segment exists - default choice ec1
c       Out of both

        elseif ((dot1c.le.0.d0) .and. (dot2c.le.0.d0)) then
          masts = ec1
          istgt = 2

c       In of both
c       1C dominant

        elseif (dot1c.gt.dot2c) then
          masts = ec1
          istgt = 3

c       2C dominant

        else
          masts = ec2
          istgt = 3
        endif

c       Update local definition of closest node (global unchanged)

        if (masts.ne.ec) then
          if (lnc.eq.1) then
            lnc = 2
          else
            lnc = 1
          endif
        endif

      endif

c     Store

      ch2(p1(1)) = masts
      ch2(p1(2)) = lnc
      ch2(p1(3)) = istgt

      end
