c$Id:$
      subroutine uscreen(isw,iow)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c     Purpose: Display Specific Comments About Program

c     Input:
c       isw:   1 - prints before normal banner
c              2 - prints after  normal banner
c-----[--.----+----.----+----.-----------------------------------------]
      implicit none

      integer  isw,iow

      if(isw.eq.1) then
        if(iow.gt.0) then
c         write(iow,2001)  ! Message to file
        else
c         write(  *,2001)  ! Message to screen
        endif
      else
        if(iow.gt.0) then
c         write(iow,2002)  ! Message to output file
        else
c         write(  *,2002)  ! Message to screen
        endif
      endif

c     Format

c2001 format(/9x,'Message before FEAP banner'/)
c2002 format(/9x,'Message after  FEAP banner'/)

      end
