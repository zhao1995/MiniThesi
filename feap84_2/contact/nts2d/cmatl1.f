c$Id:$
      subroutine cmatl1 (cp0,cm,ch1,ch2,ch3)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Coded/modified by:                           Date:          rel.:
c               Giorgio Zavarise           April 10, 1996            1.0

c      Acronym: Contact MATerial Law 1
c      Purpose: compute parameters related to contact material model

c      Inputs :
c         cp0(*)  - Contactpair control data
c         cm(*)   - Contact materials data storage
c         ch1(*)  - Contact history variables (old)
c         ch2(*)  - Contact history variables (current)

c      Outputs:
c         ch2(*)  - Contact history variables (current)
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'c_0.h'
      include  'c_comnd.h'
      include  'c_contac.h'
      include  'c_keyh.h'
      include  'c_mate.h'
      include  'c_pair.h'
      include  'counts.h'
      include  'iofile.h'

      integer   istgn,istgno
      real*8    cp0(nr0,n0c3:*),cm(*),ch1(*),ch2(*),ch3(*)
      real*8    cn,gn,area,fn,dfngn,ftt,ft,dtdft,dgnft,tde,tdet
      real*8    tdp,tolslip,tdeo,dtd,ct,ftmax,mu,istfr,istfro
      real*8    csibb,csi,d21,fnmax,fuoi,fnr,fnmaxnew,fnresid,fnaug
      real*8    fnstif,fninf,gninf,frlim, augfn,augft,gnlim

      save

      data      tolslip /1.e-8/

      call cdebug0 ('      cmatl1',-1)

c     Get penalty for normal stiffness

      cn = cp0(3,2)

c     Get geometrical variables

      gn   = ch2(p1(9))
      area = ch2(p1(11))

c     Compute normal force and derivative

      dfngn = cn * area
      fn    = dfngn * gn
      fnmax = ch3(p3(4))

c     Add augmented force

      if (ifaugm.ne.1) then
        istgn   = nint(ch2(p1(4)))
        augfn = ch2(p1(151))
        if (niter.ne.0) then
          fn = fn + augfn
        else
          fn = fn + augfn
        endif
      endif

c     Add Lagrange multiplier force

      if (ifsolm.eq.2) then
        fn = fn + ch2(p1(21))
      endif

c     This has to go in init subroutine

      if (fnmax.eq.0.d0) then
        fnmax = -cp0(5,2)
        ch3(p3(4)) = fnmax
      endif

c     Check if algorithm is activated

      if (fnmax.ne.0.d0) then
        fuoi  = ch3(p3(3))
        fnaug = ch3(p3(6))

c       Limit update: perform only if current force is over limit

        frlim = cp0(8,2)
        if (frlim.eq.0) then
          frlim = 0.8
        endif

        if (fn.lt.fnmax) then

          if (fuoi.ne.0.d0) then
            fnr = fn/fuoi
            if (fnr.gt.1.d0) fnr = 1.d0/fnr
          else
            fnr = -1.d0
          endif

          if ((fnr.gt.0.d0) .and.  (fnr.gt.frlim)) then
            fnmaxnew = max(fnmax+0.5d0*(-cp0(5,2)),
     &                    (fnmax + (fn-fnmax)*fnr)*1.05)

            fnmaxnew = min(fnmaxnew,fnmax+0.5d0*(-cp0(5,2)))
            fnmax    = min(fnmax,fnmaxnew)
            if (ifdb .and. fnmax.eq.fnmaxnew) then
              write (iow,4001) fn,fnmaxnew,fnr
              write (  *,4001) fn,fnmaxnew,fnr
4001        format ('FMAX INCREASED  fn fnmax fnr', 3f15.8)
            endif
          endif

c         Variable fn with respect to penetration
c         store here because then changed

          ch3(p3(4)) = fnmax
          gnlim      = fnmax/(cn*area)
          if (gn.lt.gnlim) then
            fninf = -cp0(6,2)
            gninf = -cp0(7,2)
            fnmax = fnmax + (fninf-fnmax)*(gn-gnlim)/(gninf-gnlim)
          endif
        endif

c       Closed

        if (fn.lt.0.d0) then

c         Limited

          if (fn.lt.fnmax) then
            if(ifdb) then
              write (iow,*) 'FN DGNFN',fn,dfngn
              write (  *,*) 'FN DGNFN',fn,dfngn
            endif
            fnresid = fnmax
            fnaug   = fnmax

c           This choice should be modified with a smooth translation

            fnstif  = 0.d0
            dfngn   = fnmax/gn

            if(ifdb) then
              write (iow,*) 'LIMITED ',fnmax,dfngn
              write (  *,*) 'LIMITED ',fnmax,dfngn
            endif

c         Augmented (skipped now)

          elseif ((gn.lt.-1.d1) .and. (fnr.gt.0.7d0)) then
            if(ifdb) then
              write (iow,*) 'AUGMENTED',fn,fn+fnaug,gn
              write (  *,*) 'AUGMENTED',fn,fn+fnaug,gn
            endif
            if (fnaug.eq.0.d0) then
              fnresid    = fn
              ch3(p3(6)) = fn
            else
              fnresid    = max((fn + fnaug),fnaug*1.3d0)
              fn         = max((fn + fnaug),fnaug*1.3d0)
              ch3(p3(6)) = fn
            endif

c           Normal case

          else
            if(ifdb) then
              write (iow,*) 'CONSTANT',fn,dfngn
              write (  *,*) 'CONSTANT',fn,dfngn
            endif

            fnresid = fn
            fnaug   = fn
            fnstif  = fn
          endif

c       Open

        else
          gnlim = -fnmax/(cn*area)
          if (gn.lt.gnlim) then

            fnresid = fnaug + fn
            if(ifdb) then
              write (iow,*) 'FNRESID',fnresid,fn,fnaug
              write (  *,*) 'FNRESID',fnresid,fn,fnaug
            endif
            if (fnresid.lt.0.d0) then
              fnaug  = fnresid

              fnstif = 0.d0
              dfngn = 0.d0

              if(ifdb) then
                write (iow,*) 'WORKING',fnresid,dfngn,gn
                write (  *,*) 'WORKING',fnresid,dfngn,gn
              endif

              ch3(p3(4)) = -cp0(5,2)

            else
              fnstif  = 0.d0
              fnaug   = 0.d0
              fnresid = 0.d0
              fn      = 0.d0
              dfngn   = 0.d0
              if(ifdb) then
                write (iow,*) 'OPEN 1'
                write (  *,*) 'OPEN 1'
              endif
              ch3(p3(4)) = -cp0(5,2)
            endif
          else

            fn = ch3(p3(3))

            fnresid = ch3(p3(3))
            fnresid = fnstif
            fnaug   = 0.d0

            fnresid = 0.d0
            fnstif  = 0.d0
            dfngn = 0.d0

          endif
        endif

c       Store

        ch3(p3(3))  = fn
        ch3(p3(5))  = fnresid
        ch3(p3(6))  = fnaug
        ch2(p1(51)) = fnstif
        ch2(p1(52)) = dfngn

c     Algorithm non active

      else

c       Store

        ch3(p3(5))  = fn
        ch2(p1(51)) = fn
        ch2(p1(52)) = dfngn
      endif

c     Coulomb friction

      if (iffric.eq.1) then

c       Get penalty for tangential stiffness

        ct = cp0(4,2)

c       Get friction coefficient

        mu = cm(1)

c       Get geometrical variables

        dtd    = ch2(p1(12))
        istgno = nint(ch1(p1(4)))
        tdeo   = ch1(p1(56))
        istfro = ch1(p1(58))

c       Previous step open -> no force, no stiffness
c       REMARK - this is a personal choice
c       REMARK - you are here only if gap is closed (istgn = 0 or  1)
c       REMARK - interesting to find contact time and start there

        if (istgno.lt.0) then
          ft    = 0.d0
          dtdft = 0.d0
          dgnft = 0.d0
          tde   = 0.d0
          tdp   = 0.d0
          istfr = -1

c       Generic case

        else

c         Trial stick state

          tdet = tdeo + dtd

c         Trial and max possible force

          ftt   = ct * tdet * area
          ftmax = mu * abs(fn)

c         Add augmented force

          if (ifaugm.ne.1) then
            augft = ch2(p1(153))
            ftt = ftt + augft
          endif

c         Check status using tolerance of slip function

c         Stick with no force

          if (abs(ftt).le.tolslip) then
            ft     = 0.d0
            dtdft  = ct * area
            dgnft  = 0.d0
            tde    = tdet
            tdp    = 0.d0
            istfr  = 0

c         Stick state

          elseif (abs(ftt).lt.(ftmax+tolslip)) then
            ft     = ftt
            dtdft  = ct * area
            dgnft  = 0.d0
            tde    = tdet
            tdp    = 0.d0
            istfr  = 1

c         Slip

          else
            ft     = sign(1.d0,ftt) * ftmax
            dtdft  = 0.d0
            dgnft  = sign(1.d0,ftt) * mu * sign(1.d0,fn) * dfngn
            tde    = ft / (ct*area)
            tdp    = dtd - (tde-tdeo)
            istfr  = 2

c           If sliding and segment not changed reset csibb

            csibb = ch2(p1(17))
            if((csibb.ne.0.d0) .and. (csibb.ne.1.d0)) then
              csi = ch2(p1(8))
              d21 = ch2(p1(7))
              csibb = csi - ft/(ct*area * d21)
              if(ifdb) then
                write (iow,*) 'csibb old,new',ch2(p1(17)),csibb
                write (  *,*) 'csibb old,new',ch2(p1(17)),csibb
              endif
              ch2(p1(17)) = csibb
            endif
          endif
        endif

c       Store

        ch2(p1(53)) = ft
        ch2(p1(54)) = dtdft
        ch2(p1(55)) = dgnft
        ch2(p1(56)) = tde
        ch2(p1(57)) = tdp
        ch2(p1(58)) = istfr
      endif

      end
