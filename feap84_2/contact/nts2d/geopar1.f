c$Id:$
      subroutine geopar1 (ndm,ndf,x,u,nel1,nod1,ix1,ix2,cxs,nc,
     &                    id,ch1,ch2,ch3)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Remove 'cp0' from routine -- unused              21/04/2007
c       2. Correct computation for 'csi' and 'gt'           08/03/2013
c-----[--.----+----.----+----.-----------------------------------------]
c      Coded/modified by:                           Date:          rel.:
c               Giorgio Zavarise           April 10, 1996            1.0

c      Acronym: GEOmetrical PARameters subroutine

c      Purpose: Compute and store all geometrical data

c      Inputs :
c         ndm     - Space dimension of mesh
c         ndf     - Number dof/node
c         x(*)    - nodal coordinates
c         u(*)    - Current nodal solution vectors
c         nel1    - Current contact element of surf. 1
c         nod1    - Current node of contact element nel1
c         ix1(*)  - Element nodal connection list for surface 1
c         ix2(*)  - Element nodal connection list for surface 2
c         cxs(*)  - Coordinates of slavenode S
c         nc      - Slave node number
c         ch1(*)  - Contact history variables (old)
c         ch2(*)  - Contact history variables (current)

c      Outputs:
c         ch2(*)  - Contact history variables (current)
c         ch3(*)  - Contact history variables (static)
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'c_0.h'
      include  'c_comnd.h'
      include  'c_contac.h'
      include  'c_geom.h'
      include  'c_keyh.h'
      include  'c_pair.h'
      include  'c_tole.h'
      include  'counts.h'

      integer   ndm,ndf,nel1,nod1,ix1(dnope1,*),ix2(dnope2,*),id(ndf,*)
      integer   nc,masts,lnc,istgt,n1,n2,ns,nb,na,n1s,n2s,ifal
      integer   istgn,istgno,istgni,ke,mastso
      real*8    x(ndm,*),u(ndf,*),cxs(*), augfn
      real*8    ch1(*),ch2(*),ch3(*),cx1(2),cx2(2),cxa(2),cxb(2),cxc(2)
      real*8    d21,s21,c21,gt,csi,gn,area,gnbar,s21c,c21c,csic
      real*8    dtd,d21o,d21oi,csio,gto,cx1s(2),cx2s(2),tdo

      save

      call cdebug0 ('      geopar1',-1)

c     Get main geometrical variables

      masts = nint(ch2(p1(1)))

      if(masts.eq.0) then
        istgn =-1
        s21   = 0.0d0
        c21   = 1.0d0
        d21   = 1.0d0
        csi   = 0.0d0
        gn    = 0.0d0
        gt    = 0.0d0
        area  = 1.0d0
        go to 100
      endif

c     Find geometrical basic parameters

      lnc    = nint(ch2(p1(2)))
      istgt  = nint(ch2(p1(3)))
      n1     = ix2(1,masts)
      n2     = ix2(2,masts)
      cx1(1) = x(1,n1) + u(1,n1)
      cx2(1) = x(1,n2) + u(1,n2)
      cx1(2) = x(2,n1) + u(2,n1)
      cx2(2) = x(2,n2) + u(2,n2)

c     Tangent vector t
c     REMARK - s21 & c21 normalized at end to minimize numerical error

      s21 = (cx2(2) - cx1(2))
      c21 = (cx2(1) - cx1(1))

c     Compute length of master segment & Normal projection onto 21
c     REMARK - gn is positive when gap open

      if(lnc.eq.1 .and. id(1,n1).ne.0 .and. id(1,nc).ne.0) then
        s21 = 0.0d0
        d21 = abs(c21)
        gn  = (cxs(2) - cx1(2)) * (-c21)
      elseif(lnc.eq.2 .and. id(1,n2).ne.0 .and. id(1,nc).ne.0) then
        s21 = 0.0d0
        d21 = abs(c21)
        gn  = (cxs(2) - cx2(2)) * (-c21)
      elseif(lnc.eq.1 .and. id(2,n1).ne.0 .and. id(2,nc).ne.0) then
        c21 = 0.0d0
        d21 = abs(s21)
        gn  = (cxs(1) - cx1(1)) * s21
      elseif(lnc.eq.2 .and. id(2,n2).ne.0 .and. id(2,nc).ne.0) then
        c21 = 0.0d0
        d21 = abs(s21)
        gn  = (cxs(1) - cx2(1)) * s21
      else
        d21 = sqrt(c21*c21 + s21*s21)
        gn  = (cxs(1) - cx1(1)) * s21 + (cxs(2) - cx1(2)) * (-c21)
      endif

      if(d21 .eq. 0.0d0) then
        istgn = -1
        go to 100
      endif

c     Compute closest point location

      csi = ((cxs(1)-cx1(1))*c21 + (cxs(2)-cx1(2))*s21)/(d21*d21)

c     Normal & tangential projection of S onto 21

      gn = gn / d21
      gt = csi

c     Now normalize s21 and c21

      s21 = s21 / d21
      c21 = c21 / d21

c     Slave node in special positions -> update parameters

      if (istgt.gt.1) then

c       Set csi to the closest node

        if (lnc.eq.1) then
          csic = 0.d0
          cxc(1) = cx1(1)
          cxc(2) = cx1(2)
        else
          csic = 1.d0
          cxc(1) = cx2(1)
          cxc(2) = cx2(2)
        endif

c       Compute new parameters

        gnbar = sqrt((cxs(1) - cxc(1))*(cxs(1) - cxc(1)) +
     &               (cxs(2) - cxc(2))*(cxs(2) - cxc(2)))

c       WARNING - this choice valid if the angle between two segment
c                 measured on the material side is > 90'
c                 (no realistic possibilities that angle < 90')

        gn = gnbar * sign(1.d0,gn)

c       If penetration and case 4 or 5 check if inside tolerance

        if ((istgt.eq.4) .and. (gn.lt.0.d0)) then

c         If outside tolerance change sign to consider it open

          if (abs(csi*d21).gt.tlopen) then
            gn = abs(gn)
          endif

        else if ((istgt.eq.5) .and. (gn.lt.0.d0)) then

c         If outside tolerance change sign to consider it open

          if (abs((csi-1.d0)*d21).gt.tlopen) then
            gn = abs(gn)
          endif
        endif

c       New tangent vector
c       REMARK - if gn=0 we are never here !
c       REMARK - signs of s21c and c21c always ok  because of gn
c       REMARK - DO NOT TOUCH ORIGINAL s21 and c21, needed for friction

        s21c =  (cxs(1) - cxc(1)) / gn
        c21c = -(cxs(2) - cxc(2)) / gn

c       Store private variables

        ch2(p1(13)) = s21c
        ch2(p1(14)) = c21c
        ch2(p1(15)) = csic

      endif

c     Correct distance for small opening
c     REMARK - only for one iteration

      if ((gn.gt.0.d0) .and. (gn.le.tlopen)) then
        istgni = nint(ch2(p1(9)))
        if ((istgni.ne.0) .or. (niter.eq.0)) then
          if(ifdb) then
            write (*,*) 'WARNING, gn set to 0'
          endif
          gn = 0.d0
        endif
      endif

c     Set gap-status flag
c     REMARK - if gn > 0.d0 gap considered open
c     Closed with contact force

c     if (gn.lt.0.d0 .or. ifadhe.eq.1 ) then
      if (gn.lt.0.d0) then
        istgn = 1

c     Closed but no contact force

      elseif (gn.eq.0.d0) then
        istgn = 0

c       Also check if Lagrange multiplier is active

        if (ifsolm.eq.2 .and. ch2(p1(21)).lt.0.0d0) then
          istgn = 1
        endif

c     Open

      elseif (gn.gt.0.d0) then
        istgn = -1

c       Correct if augmented forces are active

        if (ifaugm.ne.1) then
          augfn = ch2(p1(151))
          if (augfn.ne.0.d0) then
            istgn = 2
          endif
        endif

c       Also check if Lagrange multiplier is active

        if (ifsolm.eq.2 .and. ch2(p1(21)).lt.0.0d0
     &                  .and. abs(csi).le.1.d0+tlouts) then
          istgn = 1
        endif

      endif

c     Check for adhesion

      if(ifadhe.eq.1) then
        istgn = 1
      elseif(ifadhe.eq.2) then
        istgn = 1  ! N.B. Need to check stress for this case
      endif

c     Finally: If point is outside segment set open

      if (abs(csi).ge.1.d0+tlouts) then
        istgn = -1
      endif

c     Set tangential displacement if requested

      if (iffric.eq.1) then
        istgno = nint(ch1(p1(4)))
        mastso = nint(ch1(p1(1)))
        d21o   = ch1(p1(7))
        csio   = ch1(p1(8))
        gto    = ch1(p1(10))
        tdo    = ch1(p1(16))

c       Check tangential solution mode

c       Master segment not changed from the previous step

        if (masts.eq.mastso) then
          dtd = (csi-csio)*d21o
          d21oi = d21o

c       Master changed -> compute total sliding,
c       Set the old segm. length = to the new one for the iteration

        elseif (masts.gt.mastso) then
          if(ifdb) then
            write (*,*) 'WARNING 3 MASTS CHANGED',
     &                            masts,mastso,nstep,niter,istgn,istgt
          endif
          dtd = (1.d0 - csio)*d21o + csi*d21
          do ke = mastso+1,masts
            n1s = ix2(1,ke)
            n2s = ix2(2,ke)
            cx1s(1) = x(1,n1s) + u(1,n1s)
            cx1s(2) = x(2,n1s) + u(2,n1s)
            cx2s(1) = x(1,n2s) + u(1,n2s)
            cx2s(2) = x(2,n2s) + u(2,n2s)
            dtd = dtd + sqrt((cx2s(1)-cx1s(1))*(cx2s(1)-cx1s(1)) +
     &                       (cx2s(2)-cx1s(2))*(cx2s(2)-cx1s(2)))
            d21oi = d21
          end do
        else
            if(ifdb) then
            write (*,*) 'WARNING 4 MASTS CHANGED',
     &                            masts,mastso,nstep,niter,istgn,istgt
          endif
          dtd = -(csio*d21o + (1.d0 - csi)*d21)
          do ke = masts+1,mastso-1
            n1s = ix2(1,ke)
            n2s = ix2(2,ke)
            cx1s(1) = x(1,n1s) + u(1,n1s)
            cx1s(2) = x(2,n1s) + u(2,n1s)
            cx2s(1) = x(1,n2s) + u(1,n2s)
            cx2s(2) = x(2,n2s) + u(2,n2s)
            dtd = dtd - sqrt((cx2s(1)-cx1s(1))*(cx2s(1)-cx1s(1)) +
     &                       (cx2s(2)-cx1s(2))*(cx2s(2)-cx1s(2)))
            d21oi = d21
          end do
        endif

c       Store private parameters

        ch2(p1(12)) = dtd
        ch3(p3( 1)) = d21oi

c       Tangential initialization

      elseif (iffric.eq.2) then
        ch2(p1(17)) = csi
      endif

c     Set contact area
c     Set left and right nodes with respect to the slave one

      if(dnope1.gt.1) then
        if (nod1.eq.1) then
          na = ix1(dnope1-1,nel1)
          if (na.eq.0) then
            na = nel1
          endif
          na = ix1(1,na)
          nb = ix1(2,nel1)
        else
          na = ix1(1,nel1)
          nb = ix1(dnope1  ,nel1)
          if (nb.eq.0) then
            nb = nel1
          endif
          nb = ix1(2,nb)
        endif

c       If contact area is linearized pick current coord
c       REMARK - still to do --> constant contact area

        ifal = 0

c       Constant contact area (original)

        if ((ifal.eq.0) .or. (ndf.eq.0)) then
          cxa(1) = x(1,na)
          cxa(2) = x(2,na)
          cxb(1) = x(1,nb)
          cxb(2) = x(2,nb)

c         Correct node slave coordinates

          ns = ix1(nod1,nel1)
          cxs(1) = x(1,ns)
          cxs(2) = x(2,ns)

c         Variable contact area (current)

        else
          cxa(1) = x(1,na) + u(1,na)
          cxa(2) = x(2,na) + u(2,na)
          cxb(1) = x(1,nb) + u(1,nb)
          cxb(2) = x(2,nb) + u(2,nb)
        endif

c       Compute area

        area = (sqrt((cxa(1)-cxs(1))*(cxa(1)-cxs(1))   +
     &               (cxa(2)-cxs(2))*(cxa(2)-cxs(2)))  +
     &          sqrt((cxb(1)-cxs(1))*(cxb(1)-cxs(1))   +
     &               (cxb(2)-cxs(2))*(cxb(2)-cxs(2)))) * 0.5d0

      else
        area = 1.d0
      endif

c     Store variables

100   ch2(p1(4))  = istgn
      ch2(p1(5))  = s21
      ch2(p1(6))  = c21
      ch2(p1(7))  = d21
      ch2(p1(8))  = csi
c     ch2(p1(9))  = min(1.d-6,gn)
      ch2(p1(9))  = gn
      ch2(p1(10)) = gt
      ch2(p1(11)) = area

      end
