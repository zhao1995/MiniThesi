c$Id:$
      subroutine plstop()

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Check of active problem before calling ptimpl.   13/05/2009
c       2. Remove check on 'pfeap_on' for parstop call      13/06/2009
c       3. Move pdelfl after tplot options                  15/09/2009
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Close any open plot windows and stop execution

c      Inputs:
c         none

c      Outputs:
c         none
c-----[--.----+----.----+----.-----------------------------------------]
      implicit   none

      include   'allotd.h'
      include   'allotn.h'
      include   'comsav.h'
      include   'counts.h'
      include   'debugs.h'
      include   'endata.h'
      include   'iofile.h'
      include   'pdatps.h'
      include   'plflag.h'
      include   'rdata.h'
      include   'rdat1.h'
      include   'tdata.h'
      include   'x11f.h'

      include   'setups.h'

      logical    lopen, setval,palloc
      character  dnam*5
      integer    m, mm, dnum, dpre
      real*4     etime,tary(2)

      save

c     Close PostScript file if open

      if (hdcpy) call fpplcl()

c     X11 device

      if(everon) call gdx11(6,xx,yy)

c     Write last log record for transient problems

      inquire(unit = ilg, opened=lopen)
      if(lopen) then
        if(niter.gt.1) then
          write(ilg,3001) nstep,niter,nform,ttim,dt,rel0,rnorm,rnmax,
     &                    aengy,etime(tary)
        else
          write(ilg,3002) nstep,niter,nform,ttim,dt,rel0,rnmax,
     &                    etime(tary)
        endif
        titer = titer + niter
        tform = tform + nform
        write(ilg,3003) titer,tform
        close(unit=ilg)
      endif

c     Do last tplot outputs if required

      if(.not.prob_on) then
        rfl = .false.
        call ptimpl()
        prob_on = .true.   ! Prevent time plots from causing recurrsion
      endif

c     Close open tplot files (keep)

      do m = 30,75
        inquire(unit = m, opened = lopen)
        if(lopen) then
          close(m)
        endif
      end do ! m

c     Delete existing files

      call pdelfl()

c     Delete memory use

      mm = ndict
      do m = mm,1,-1
        dnum   = dlist(m)
        dnam   = dict(m)
        dpre   = iprec(m)
        if(debug) then
          write(  *,*) ' RANK:',rank,dnum,dnam,dpre
          write(iow,*) ' RANK:',rank,dnum,dnam,dpre
        endif
        setval = palloc(dnum, dnam, 0, dpre)
      end do ! m

c     Close parallel arrays

      call parstop()

c     Close last output file

      close(unit=iow)

      stop

c     Format

3001  format(2i6,i5,1p,1e11.3,1p,5e10.2,    0p,1f10.2)
3002  format(2i6,i5,1p,1e11.3,1p,2e10.2,10x,1p,1e10.2,
     &       10x,0p,1f10.2)
3003  format(/' Total',i6,i5/'  ',34('-'),' END OF FEAP LOG ',35('-'))

      end
