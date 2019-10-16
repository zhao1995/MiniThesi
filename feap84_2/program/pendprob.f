c$Id:$
      subroutine pendprob

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    12/05/2009
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: End a problem

c      Inputs:
c        None

c      Outputs:
c        Close tplot files and clean up memory
c-----[--.----+----.----+----.-----------------------------------------]
      implicit   none

      include   'allotd.h'
      include   'allotn.h'
      include   'cdata.h'
      include   'comfil.h'
      include   'comsav.h'
      include   'contrl.h'
      include   'counts.h'
      include   'debugs.h'
      include   'endata.h'
      include   'iofile.h'
      include   'pdatps.h'
      include   'rdata.h'
      include   'rdat1.h'
      include   'sdata.h'
      include   'tdata.h'

      logical    errs
      real*4     etime, tary(2)

      save

c     Set last data in tplots for multiple problem case

      if(prob_on) then

        inquire(unit=ior,name=fnamp,exist=errs)

        if(errs) then

c         Clear plot files, delete scratch files and close output file

          if(hdcpy) call fpplcl()
          if(niter.gt.1) then
            write(ilg,3004) nstep,niter,nform,ttim,dt,rel0,rnorm,rnmax,
     &                      aengy,etime(tary)
          else
            write(ilg,3005) nstep,niter,nform,ttim,dt,rel0,rnorm,rnmax,
     &                      etime(tary)
          endif
          titer = titer + niter
          tform = tform + nform
          write(ilg,3006) titer,tform
          rfl = .false.
          call ptimpl()

c         Restore master output file number and name

          close(iow, status = 'keep')
          iow     = iow_sav
          fout    = fout_sav

        endif
      endif

c     Formats

3004  format(2i6,i5,1p,1e11.3,1p,5e10.2,    0p,1f10.2)
3005  format(2i6,i5,1p,1e11.3,1p,4e10.2,10x,0p,1f10.2)

3006  format(/'Total',i7,i5)

      end
