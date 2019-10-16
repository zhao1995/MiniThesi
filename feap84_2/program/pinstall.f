c$Id:$
      subroutine pinstall()

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Remove set of 'incred'; for memory adjustment    09/07/2012
c       2. Remove set of manual files                       08/02/2013
c-----[--.----+----.----+----.-----------------------------------------]
c     Purpose: Set installation parameters

c     Inputs:  None

c     Outputs: Written in file: feap.ins
c-----[--.----+----.----+----.-----------------------------------------]
      implicit   none

      save

      call pinset(1)
      call pinset(2)

      end

      subroutine pinset(isw)

      implicit  none

      include  'chdata.h'
      include  'comfil.h'
      include  'codat.h'
      include  'hlpdat.h'
      include  'iofile.h'
      include  'iosave.h'
      include  'prmptd.h'
      include  'psize.h'

      character text(3)*15
      logical   pcomp,exst, flag, tinput
      integer   isw
      real*8    td(1)

      save

c-----[--.----+----.----+----.-----------------------------------------]

      if(isw.eq.1) then
        inquire(file='feap.ins',exist=exst)
        if(exst) then
          open( unit = ior, file = 'feap.ins', status = 'old',
     &          access = 'sequential')
        endif
      else
        inquire(file='./.feap.ins',exist=exst)
        if(exst) then
          open( unit = ior, file = './.feap.ins', status = 'old',
     &          access = 'sequential')
        endif
      endif
      if(exst) then
        eofile = .true.
        flag   = .false.
        do while(.not.flag)

          flag = tinput(text,3,td,0)

          if(pcomp(text(1),'noparse',7)) then

            coflg = .false. ! Numerical input mode, noparsing

          elseif(pcomp(text(1),'parse  ',7)) then

            coflg = .true.  ! Parse all input as expressions

c         Set graphics default options

          elseif(pcomp(text(1),'graphic',7)) then
            if    (pcomp(text(2),'prompt ',7)) then
              if(pcomp(text(3),'off',3)) then
                prompt = .false.
              else
                prompt = .true.
              endif
            elseif(pcomp(text(2),'default',7)) then
              if(pcomp(text(3),'off',3)) then
                defalt = .false.
              else
                defalt = .true.
              endif
            endif

c         Set PostScript default mode

          elseif(pcomp(text(1),'postscr',7)) then
            if    (pcomp(text(2),'color  ',7)) then
              pscolr = .true. ! Color PostScript
              if(pcomp(text(3),'reverse',7)) then
                psrevs = .true. ! Color order is reversed
              else
                psrevs = .false.! Color order is normal
              endif
            else
              pscolr = .false.! Grayscale PostScript
            endif

c         Set help display level:

          elseif(pcomp(text(1),'helplev',7)) then
            if    (pcomp(text(2),'basic  ',7)) then
              hlplev = 0      ! Basic
            elseif(pcomp(text(2),'interme',7)) then
              hlplev = 1      ! Intermediate
            elseif(pcomp(text(2),'advance',7)) then
              hlplev = 2      ! Advanced
            else
              hlplev = 3      ! Expert
            endif

c         Set file check on startup

          elseif(pcomp(text(1),'fileche',7)) then

            if(pcomp(text(2),'off',3)) then
              fileck = .false.
            else
              fileck = .true.
            endif

          endif
        end do ! while
        close(ior)
      endif

c-----[--.----+----.----+----.-----------------------------------------]

      end
