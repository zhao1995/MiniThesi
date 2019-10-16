c$Id:$
      subroutine skipcdat ()

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Assign length to cc of 4 characters              17/04/2007
c-----[--.----+----.----+----.-----------------------------------------]
c      Coded/modified by:                           Date:          rel.:
c               Giorgio Zavarise           April 10, 1996            1.0

c      Acronym: Equivalence with pnums

c      Purpose: Determine maximum contact material set, maximum
c               contact surface set, maximum contactpairs

c      Inputs:
c                 - From read of data file(s)

c      Outputs:
c         numcs   - # of contact surfaces
c         numcm   - # of contact material laws
c         numcp   - # of contact pairs
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'c_0.h'
      include  'c_dict.h'
      include  'chdata.h'
      include  'iofile.h'

      logical   errck,pcomp,gettxtd,whfl,whf2
      character cc*4,tx(2)*15
      integer   ksc
      real*8    td(15)

      save

      call cdebug0 ('    skipcdat',-1)

c     Skip command data

      errck = gettxtd(tx,2,td,0,'skip')
      cc   = tx(1)(1:4)
      whfl = .true.
      do while (whfl)
        do while (.not.pcomp(cc,'    ',4))
          errck = gettxtd(tx,2,td,0,'skip')
          cc    = tx(1)(1:4)
        end do
        errck = gettxtd(tx,2,td,0,'skip')
        cc    = tx(1)(1:4)

c       Check if blank line divide a subcommand or a new command

        ksc  = 0
        whf2 = .true.
        do while (whf2)
          ksc = ksc+1
          if (ksc.le.c_ncc) then
            if (pcomp(cc,cwd(ksc),4)) then
              whfl = .false.
              whf2 = .false.
              backspace ior
            endif
          else
            whf2 = .false.
          endif
        end do
      end do

      end
