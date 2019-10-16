c$Id:$
      subroutine pstart()

c      * * F E A P * * A Finite Element Analysis Program
c                        -      -       -        -
c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Replace MSFLIB/PORTLIB by DFLIB,DFPORT           08/12/2006
c       2. Add 'setups.h' and set ntasks = 1; rank = 0      24/02/2009
c-----[--+---------+---------+---------+---------+---------+---------+-]
c      Purpose: Start graphical outputs: Windows version with 'open'

c      Inputs: none

c      Outputs: none
c-----[--+---------+---------+---------+---------+---------+---------+-]

      use dflib
      use dfwin

      implicit  none

      include  'codat.h'
      include  'comfil.h'
      include  'iodata.h'
      include  'machnc.h'
      include  'memuse.h'
      include  'pfeapb.h'
      include  'prmptd.h'
      include  'setups.h'

      integer      status, ilen, i
      character    file_spec*512, fildata*512, pathname*512
      logical*4    ldrv
      character*(*), parameter:: filter_spec = "Text Files"C // "I*"C
     &             // "All Files"C // "*.*"C // ""C
      character*(*), parameter:: DLGTITLE = "Choose Input File"C

      type(T_OPENFILENAME) :: ofn

      interface
        logical(4)  function initialsettings
        end function
      end interface

c     Set flags for serial version

      pfeap_on   = .false.
      pfeap_gnod = .false.
      ntasks     = 1
      rank       = 0

c     Set maximum memory use: 0 = unlimited;
c                           > 0 maximum array size allocated

      maxuse = 0

c     Get

      ofn%lStructSize       = sizeof(ofn)
      ofn%hwndOwner         = GetForegroundWindow()
      ofn%hInstance         = NULL      ! Set hInstance if desired
      ofn%lpstrFilter       = loc(filter_spec)
      ofn%lpstrCustomFilter = NULL
      ofn%nMaxCustFilter    = 0
      ofn%nFilterIndex      = 1         !Specifies initial filter value
      ofn%lpstrFile         = loc(file_spec)
      ofn%nMaxFile          = sizeof(file_spec)
      ofn%nMaxFileTitle     = 0
      ofn%lpstrInitialDir   = NULL      !Open current directory
      ofn%lpstrTitle        = loc(DLGTITLE) !Give title to dialog box
      ofn%Flags             = OFN_PATHMUSTEXIST
      ofn%lpstrDefExt       = loc("txt"C)
      ofn%lpfnHook          = NULL
      ofn%lpTemplateName    = NULL

      if(ciflg) then
        status = GetOpenFileName(ofn)
        fileck = .false.  ! File checking at startup is off

c       Extract the file name
        if(status.eq.0) then
          call plstop()
        else
          ilen = index(file_spec,char(0)) - 1
          fildata = file_spec(1:ilen)
          do i = ilen,1,-1
            if(fildata(i:i).eq.'\') then
              pathname = fildata(1:i-1)
              finp     = fildata(i+1:ilen)
              exit
            endif
          end do ! i
          fout      = finp
          fres      = finp
          fsav      = finp
          fplt      = finp
          fout(1:1) = 'O'
          fres(1:1) = 'R'
          fsav(1:1) = 'S'
          fplt(1:1) = 'P'
          status    = CHANGEDIRQQ(pathname)
          open(ios,file='feapname',status='unknown')
          call getcon(epmac)
          write(ios,1000) finp,fout,fres,fsav,fplt,epmac
          close(ios,status='keep')
        endif

c     Files read from filnam inputs

      else
        fileck = .true.   ! File checking at startup is on
      endif

c     Set restart flag to false (reset by command line inputs only)

      reflg = .false.

c     Graphics Driver number for PC version

      call pdriver()

c     Initialize memory

      call pinitm()

c     Start Windows

      call pwopn()

c     Check user installation options

      call pinstall()

c     Set initial file names

      call filnam()

c     Format

1000  format(a/a/a/a/a/1p,1e25.15)

      end

      logical(4)   function initialsettings ( )

      use          DFLIB
      use          DFPORT

      implicit     none

      type(qwinfo) winfo
      logical      status, pcomp, dflag
      integer      istatus, j,l,len, ifile, lun, lmenu
      character    record*80,str*80

      save

c     Maximize Frame

      winfo.type = qwin$max
      status     = setwsizeqq(qwin$framewindow,winfo)

      initialsettings = .true.

      end function
