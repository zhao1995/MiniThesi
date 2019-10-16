c$Id:$
      subroutine usetio(ni,u,h,lenu,lenh,flgh,iounit, isw)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: I/O routines for saving displacement and history data.

c      Inputs:
c         ni        - Point number
c         u(lenu)   - Displacments (if isw = 1)
c         h(lenh)   - History data (if isw = 1)
c         lenu      - Number entries in 'u' array
c         lenh      - Number entries in 'h' array
c         flgh      - History data exists if true
c         iounit    - I/O logical unit number
c         isw       - Switch parameter: = 1 for read, = 2 for write

c      Outputs:
c         u(lenu)   - Displacments (if isw = 2)
c         h(lenh)   - History data (if isw = 2)
c-----[--.----+----.----+----.-----------------------------------------]
      implicit   none

      include   'debugs.h'
      include   'iofile.h'
      include   'rdata.h'
      include   'rdat1.h'
      include   'setups.h'

      logical    flgh
      integer    ni, lenu,lenh, iounit, isw, ri
      real*8     u(lenu),h(lenh)

      save

      if(debug) then
        call udebug(' usetio',isw)
      endif

c     Input displacement and history data

      if(isw.eq.1) then
        read (iounit,end=999) rnmax
        read (iounit,end=999) ri,u
        if(flgh) read (iounit,end=999) h
        if(ri.eq.0) then
          write(iow,*) ' WARNING: READ:',rank,' NI =',ni,' RI =',ri
          write(  *,*) ' WARNING: READ:',rank,' NI =',ni,' RI =',ri
        else
          ni = ri
        endif
c       if(debug) then
c         call mprint(u,1,lenu,1,'U_histin')
c         if(flgh) call mprint(h,1,lenh,1,'H_histin')
c       endif

c     Output displacement and history data

      elseif(isw.eq.2) then
        write(iounit) rnmax
        write(iounit) ni,u
        if(flgh) write(iounit) h
        call pflush(iounit)
c       if(debug) then
c         call mprint(u,1,lenu,1,'U_histout')
c         if(flgh) call mprint(h,1,lenh,1,'H_histout')
c       endif
      endif

      return

c     Error message

999   write(iow,*) ' ERROR: READ =',rank,' NI =',ni,' IOUNIT =',iounit,
     &             ' LENU =',lenu
      write(  *,*) ' ERROR: READ =',rank,' NI =',ni,' IOUNIT =',iounit,
     &             ' LENU =',lenu

      end
