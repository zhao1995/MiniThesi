c$Id:$
      logical function gettxtd (tx,nt,td,nn,ferrc)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Coded/modified by:                           Date:          rel.:
c               Giorgio Zavarise           April 10, 1996            1.0

c      Acronym: GIVe TX and TD data (text and numerical)

c      Purpose: Get numerical data from input file with error control

c      Inputs:
c         nt      - Amount of text data to be read
c         nn      - Amount of real data to be read
c         ferrc   - Flag for ERRor Control (skip,show,stop,back)
c                   skip -> search for the first valid line
c                   show -> search and print skipped lines
c                   stop -> stop if line not found immediately
c                   back -> give back the error control

c      Outputs:
c         tx      - Text data
c         td      - Real data
c         gettd   - Flag, returns true if error occurs
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'chdata.h'

      logical   errck,tinput, pcomp
      character tx(2)*15,ferrc*(*)
      integer   nt,nn
      real*8    td(*)

      save

c      call cdebug0 ('      - gettxtd',-1)

      errck = tinput(tx,nt,td,nn)

      if (errck) then
        if(pcomp(ferrc,'skip',4)) then
          do while(errck)
            errck = tinput(tx,nt,td,nn)
          end do

c         WARNING "show" to be completed

        elseif(pcomp(ferrc,'show',4)) then
          do while(errck)
c            call prtx(xxx,80
            errck = tinput(tx,nt,td,nn)
          end do

c         WARNING "stop" to be completed

        elseif(pcomp(ferrc,'stop',4)) then
          call plstop()

c         Back

        elseif(pcomp(ferrc,'back',4)) then
          continue
        endif
      endif

      gettxtd = errck

      end
