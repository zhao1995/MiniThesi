c$Id:$
      subroutine ptmplt(ftyp, ttim, tpl,ntplts, ntstep, iunit, oflag)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Output information to time history files

c      Inputs:
c         ftyp(*)   - Type of data
c         ttim      - Solution time for data
c         tpl(*)    - Time history data
c         ntplts    - Number of time history data items
c         ntstep    - Indicator for first time step
c         iunit     - Unit number (30-75)
c         oflag     - Close after write if true

c      Outputs:
c         none      - Data saved to disk
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'comfil.h'
      include  'iodata.h'

      logical   oflag, exst, lopen
      character ftyp*(*), fnamr*132, fext*8
      integer   ntplts,ntstep, iunit, n,nn,ntp,ntl,nunit, m, ipos
      real*8    ttim, tpl(ntplts), tdum

      save

c     Set file name for output of time history data

      fnamr = fplt

c     Extract name for file and strip all data after any period ('.')

      n = index(fnamr,'.')
      if(n.gt.0) then
        fnamr(n:128) = ' '
      endif

c     Locate character where added letter 'A' to 'J' to be added

      n     = ipos(fnamr,128)
      nn    = 0
      nunit = iunit
      do ntp = 1,ntplts,20
        nn             = nn + 1
        fnamr(n+1:n+1) = char(96+nn)
        ntl            = min(ntp+19,ntplts)

c       Add extender

        if(nn.eq.1) then
          fext  =  ftyp
          call addext(fnamr,fext,128,8)
        endif

c       Check if file exists, if it does delete it

        inquire(file=fnamr,exist=exst)
        if(exst.and.ntstep.eq.1) then
          open(unit=nunit, file = fnamr, form = 'formatted',
     &         access = 'sequential', status = 'unknown')
          close(nunit,status='delete')
        end if

c       Open file and find end

        inquire(file=fnamr,opened=lopen)
        if(.not.lopen) then
          open(unit=nunit, file = fnamr, form = 'formatted',
     &         access = 'sequential', status = 'unknown')
        endif
        if(oflag) then
  10      read(nunit,2000,end=11) tdum
          go to 10
        endif

c       Add line of data

  11    write(nunit,2000) ttim,(tpl(m),m=ntp,ntl)

c       Close file if necessary

        if(oflag) then
          close(nunit)
        endif

        nunit = nunit + 1
      end do ! ntp

c     Format

2000  format(1p,21e12.4)

      end
