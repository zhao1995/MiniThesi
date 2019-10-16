c$Id:$
      subroutine dynlib(u,urate, isw)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Change '3' to 'ip' for udynam call               06/07/2007
c       2. Add dyna10 for Euler implicit algorithm          23/07/2008
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Call dynamic integrator for update steps

c      Inputs:
c        u(ndf,numnp,*)      - Solution states
c        isw                 - Switch state: 1 = initialize step
c                                            2 = update within step
c                                            3 = back up to start

c      Outputs:
c        urate(ndf,numnp,*)  - Solution rate states
c-----[--+---------+---------+---------+---------+---------+---------+-]
      implicit   none

      include   'cdata.h'
      include   'ddata.h'
      include   'part0.h'
      include   'sdata.h'

      integer    isw, ip,ipl(3)
      real*8     u(ndf,numnp,*),urate(*)

      save

      data       ipl / 1, 3, 2 /

c     Set update parameter

      ip = ipl(isw)

c     Newmark-Beta updates

      if(noi.eq.1) then

        call dyna01(u(1,1,ip),urate,nneq,ndf,ndfp,ndfo,isw)

c     Backward Euler update

      elseif(noi.eq.2) then

        call dyna02(u(1,1,ip),urate,nneq,ndf,ndfp,ndfo,isw)

c     Conserving HHT update

      elseif(noi.eq.3) then

        call dyna03(u(1,1,ip),urate,nneq,ndf,ndfp,ndfo,isw)

c     Newmark explicit update

      elseif(noi.eq.4) then

        call dyna04(u(1,1,ip),urate,nneq,ndf,ndfp,ndfo,isw)

c     Three-parameter algorithm in Conservation form

      elseif(noi.eq.5) then

        call dyna05(u(1,1,ip),urate,nneq,ndf,ndfp,ndfo,isw)

c     STATICS:  Generalized Mid-point configuration

      elseif(noi.eq.6) then

        call dyna06(u(1,1,ip),urate,nneq,ndf,ndfp,ndfo,isw)

c     FIRST ORDER :  Generalized Mid-point configuration

      elseif(noi.eq.7) then

        call dyna07(u(1,1,ip),urate,nneq,ndf,ndfp,ndfo,isw)

c     CENTRAL DIFFERENCE: Explicit algorithm

      elseif(noi.eq.8) then

        call dyna08(u(1,1,ip),urate,nneq,ndf,ndfp,ndfo,isw)

c     BACKWARD DIFFERENCE FORMULA (BDF2): Implicit algorithm

      elseif(noi.eq.9) then

        call dyna09(u(1,1,ip),urate,nneq,ndf,ndfp,ndfo,isw)

c     Euler implicit for 2nd order problems

      elseif(noi.eq.10) then

        call dyna10(u(1,1,ip),urate,nneq,ndf,ndfp,ndfo,isw)

c     USER:  User generated routine

      elseif(noi.eq.-1) then

        call udynam(u(1,1,ip),u,urate,nneq,ndf,ndfp,ndfo,npart,isw)

      endif

      end
