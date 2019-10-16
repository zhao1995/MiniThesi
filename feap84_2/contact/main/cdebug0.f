c$Id:$
      subroutine cdebug0 (strin,isw)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Coded/modified by:                           Date:          rel.:
c               Giorgio Zavarise           April 10, 1996            1.0

c      Acronym: Contact DEBUG basic

c      Purpose: Debug subroutine to print all contact subroutine called

c      Inputs :
c         strin   - Subroutine name
c         isw     - status flag

c      Outputs:
c                 - In file 'Cdebug0'
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'c_0.h'
      include  'c_contac.h'
      include  'counts.h'
      include  'tdata.h'

      character strin*(*)
      integer   isw, debf0

      save

c     Set file unit for list of contact calls

      if (strin.eq.'debf') then
        if (isw.ne.-1) then
          debf0 = isw
        else
          debf0 = 99
        endif
        open(unit = debf0, file ='Cdebug0', status = 'unknown')
      endif

c     Perform listing only if debug mode is active

      if (ifdb) then

c       Write flag only if different from zero

        if (isw.eq.-1) then
          write (debf0,*) strin
        else
          write (debf0,*) strin,isw,niter,ttim
        endif
      endif

      end
