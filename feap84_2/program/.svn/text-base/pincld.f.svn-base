c$Id:$
      subroutine pincld(fname)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Print missing include file to screen             06/05/2012
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Control file for I/O from include files

c      Inputs:
c         fname    - Name of file to read include data from
c                    N.B. If fname = 'end', closes include file

c      Outputs:
c         none
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'comfil.h'
      include  'iodata.h'
      include  'ioincl.h'
      include  'iofile.h'
      include  'iosave.h'

      logical   errc,lopen,pcomp
      character fname*(*),dnam*15, fnamr*21, fext*8

      save

c     Perform inputs from an include file

      dnam = fname

      if(pcomp(dnam,'end',3)) then
        inquire(file=fnamr,opened=lopen,exist=errc)
        if(errc.and.lopen) then
          close(icf)
        else
          inquire(unit=icf,opened=lopen,exist=errc)
          if(errc.and.lopen) then
            close(icf)
          endif
        endif
        if(isf.gt.1) then
          fincld(isf) = ' '
        endif
        isf   = max(1,isf-1)
        ior   = isfile(isf)
        icf   = max(icl,abs(ior))
        lread = .false.
      else
        inquire(file = dnam, exist = errc)
        if(errc) then
          icf   = max(icl,abs(ior))
  1       inquire(unit=icf,opened=lopen)
          if(lopen) then
            icf = icf + 1
            go to 1
          endif
          fnamr =  dnam
          fext  =  dnam(1:8)
          call opnfil(fext,fnamr,-2,icf,lread)
          if(.not.lread) then
            if(ior.lt.0) return
            call plstop()
          endif
          isfile(isf) = ior
          ior         = icf
          isf         = isf + 1
          fincld(isf) = fnamr
          irecrd(isf) = 0
        else
          write(ilg,3000) dnam
          write(iow,3000) dnam
          write(  *,3000) dnam
          call plstop()
        endif
      endif

c     Format

3000  format(/' **ERROR** File ',a,' in INCLUDE statement does not',
     &        ' exist'/)
      end
