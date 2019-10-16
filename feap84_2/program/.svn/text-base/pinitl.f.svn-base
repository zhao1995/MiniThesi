c$Id:$
      subroutine pinitl(lct,ctl,err)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Add 'init,mate' option                           16/01/2007
c       2. Add 'init,regi' option                           17/01/2007
c       3. Add 'init,acce' option                           24/04/2011
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Establish initial conditions for solution vectors
c               [init,disp] - set initial displacements
c               [init,rate] - set initial rates
c               [init,acce] - set initial acceleration
c               [init,spin,omg1,omg2,omg3] - set initial spins
c               [init,mate] - set initial velocity for material sets
c               [init,regi] - set initial velocity for regions

c      Inputs:
c        lct             - Parameter: disp,rate,spin
c        ctl(3)          - Values for omg1,omg2,omg3 for spins

c      Outputs:
c        hr(np(40))      - Initial displacements
c        hr(np(42))      - Initial rates
c        hr(np(42)+nneq) - Initial accelerations
c        err             - Set to .true. if error occurs
c-----[--.----+----.----+----.-----------------------------------------]
      implicit   none

      include   'cdata.h'
      include   'ddata.h'
      include   'iofile.h'
      include   'prflag.h'
      include   'print.h'
      include   'sdata.h'
      include   'pointer.h'
      include   'comblk.h'

      logical    err, pcomp, palloc,setvar
      character  lct*15
      real*8     ctl(3)

      save

      err = .false.

c     Rate initial conditions

      if(pcomp(lct,'rate',4) .or.
     &   pcomp(lct,'velo',4) .or.
     &   pcomp(lct,'acce',4) .or.
     &   pcomp(lct,'spin',4)) then

c       Check that memory allocated for rates

        if(np(42).eq.0) then
          setvar = palloc( 42,'VEL  ',nneq*nrt ,2)
        endif

        ivelfl = .true.

c       Rate nodal conditions

        if(pcomp(lct,'rate',4)) then
          call genvec(ndf,ndf,hr(np(42)),' rate ',prt,prth,err,.false.)
          call genclr(ndf,hr(np(42)), mr(np(190)), numnp)

c       Acceleration nodal conditions

        elseif(pcomp(lct,'acce',4)) then
          call genvec(ndf,ndf,hr(np(42)+nneq),' accel',
     &                prt,prth,err,.false.)
          call genclr(ndf,hr(np(42)+nneq), mr(np(190)), numnp)

c       Spin nodal conditions

        else
          setvar = palloc(111,'TEMP1',numnp,1)
          call genspi(hr(np(42)),mr(np(111)),mr(np(33)),hr(np(43)),
     &                nen,nen1,ndf,ndm,numel,numnp,ctl,prt,prth)
          setvar = palloc(111,'TEMP1',0,1)
        end if

c     Set initial conditions by material number

      elseif(pcomp(lct,'mate',4)) then

        setvar = palloc(111,'TEMP1',numnp,1)
        call pinitlm(hr(np(40)),hr(np(42)),mr(np(33)),mr(np(111)),
     &               nen,nen1,ndf,numel,numnp,prt,prth,1)
        setvar = palloc(111,'TEMP1',0,1)

c     Set initial conditions by region number

        setvar = palloc(111,'TEMP1',numnp,1)
      elseif(pcomp(lct,'regi',4)) then

        setvar = palloc(111,'TEMP1',numnp,1)
        call pinitlm(hr(np(40)),hr(np(42)),mr(np(33)),mr(np(111)),
     &               nen,nen1,ndf,numel,numnp,prt,prth,2)
        setvar = palloc(111,'TEMP1',0,1)

c     Displacement initial conditions

      else
        idisfl = .true.
        call genvec(ndf,ndf,hr(np(40)),' displ. ',prt,prth,err,.false.)
      endif

c     Check for errors

      if(err .and. ior.gt.0 ) then
        write(ilg,3000)
        write(iow,3000)
        call plstop()
      endif

c     Format

3000  format(/'  *ERROR* PINITL: Bad initial conditions')

      end
