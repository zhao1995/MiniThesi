c$Id:$
      subroutine writer(ctc,b,nneq)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Dimension 'b(nneq,2)' for 'n+1','n' displacement 29/03/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Save nodal displacement and stress values for later
c               use

c      Inputs:
c         ctc       - Name of array, file rewind, or file name
c         b(*)      - Array to write
c         numnp     - Length of array to write

c      Outputs:
c         none
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'cdata.h'
      include  'chdata.h'
      include  'fdata.h'
      include  'iodata.h'
      include  'iofile.h'
      include  'pdata3.h'
      include  'prstrs.h'
      include  'strnum.h'
      include  'tdata.h'
      include  'comblk.h'

      include  'p_int.h'
      include  'p_point.h'

      logical   lflg,pcomp
      character ctc*(*),ct*4,fname*20,y*1
      integer   nneq
      real*8    b(nneq,2)

      save

      data      lflg/.false./

c     Set name

      ct = ctc

c     Save current eigenpairs

      if(pcomp(ct,'eige',4)) then

        if(lflg) then
          call peigsv(fname,1)
        else
          go to 920
        endif

c     Save current displacement state

      elseif(pcomp(ct,'disp',4)) then
        if(lflg) then
          write(iwd,err=910) ct,ttim
          write(iwd,err=910) b
        else
          go to 920
        endif

c     Save current nodal stress state

      elseif(pcomp(ct,'stre',4)) then
        if(lflg) then
          if(fl(11)) then
            write(iwd,err=910) ct,ttim
            fp(2) = nph -1 + npstr*numnp
            write(iwd,err=910) istv,(hr(point),point=nph,fp(2))
          else
            write(iow,2004)
            write(ilg,2004)
            if(ior.lt.0) write(*,2004)
          endif
        else
          go to 920
        endif

c     Rewind file

      elseif(pcomp(ct,'wind',4)) then
        if(lflg) then
          rewind iwd
        endif

c     Close file

      elseif(pcomp(ct,'clos',4)) then
        close(iwd)
        lflg = .false.

c     Set filename

      else
        fname = ct
        call prwext(fname)
        inquire(file=fname,exist=lflg)

c       This is a old file name

        if(lflg) then
          if(ior.lt.0) then
            write(xxx,2002) fname
            call pprint(xxx)
10          read (*,1000,err=11,end=12) y
            goto  13
11          call  errclr ('WRITER')
            goto  10
12          call  endclr ('WRITER',y)
13          if(y.ne.'y' .or. y.ne.'Y') return
          else
            write(iow,2003) fname
            write(ilg,2003) fname
            call plstop()
          endif
          open(iwd,file=fname,status='old',form='unformatted')
          rewind iwd

c       This is a new file name

        else
          if(ior.lt.0) write(*,2001) fname
          write(iow,2001) fname
          open(iwd,file=fname,status='new',form='unformatted')
        endif
        lflg = .true.
      endif
      return

c     Error messages

910   if(ior.gt.0) then
        write(iow,3001) ct
        write(ilg,3001) ct
        call plstop()
      else
        write(*,3001)
      endif
      return

920   if(ior.gt.0) then
        write(iow,3002)
        write(ilg,3002)
        call plstop()
      else
        write(*,3002)
      endif

c     Formats

1000  format(a1)

2001  format('   Output file for write operations is named ',a)

2002  format(' *WARNING* File ',a,' exists. Erase? (y or n) >')

2003  format(' *ERROR* WRITER: File ',a,' exists.')

2004  format(' *ERROR* WRITER: Nodal stresses do not exist for write')

3001  format(' *ERROR* WRITER: On a write command for ',a4)

3002  format(' *ERROR* WRITER: No write file is open.')

      end
