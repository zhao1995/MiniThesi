c$Id:$
      subroutine punits()

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c     Purpose: Set unit multipliers for length, force and time

c     Inputs:

c     Outputs:
c-----[--.----+----.----+----.-----------------------------------------]
      implicit   none

      include   'cdata.h'
      include   'iofile.h'
      include   'pglob1.h'
      include   'sdata.h'

      include   'pointer.h'
      include   'comblk.h'

      character  text*15
      logical    errck, pcomp, tinput
      real*8     td(1)

      save

      text = 'start'
      do while (.not.pcomp(text,'    ',4))

        errck = tinput(text,1,td,1)

c       Set length multipler

        if(pcomp(text,'leng',4)) then
          units(1) = td(1)

c       Set force multipler

        elseif(pcomp(text,'forc',4)) then
          units(2) = td(1)

c       Set time multipler

        elseif(pcomp(text,'time',4)) then
          units(3) = td(1)

        endif

      end do ! while

c     Output multiplier values

      write(iow,2000) units
      if(ior.lt.0) then
        write(*,2000) units
      endif

c     Apply multipliers

      if(units(1).ne.1.d0) then
        call punitm(hr(np(43)),ndm*numnp,units(1))
      endif

      if(units(2).ne.1.d0) then
        call punitm(hr(np(27)),ndf*numnp,units(2))
      endif

c     Format

2000  format(/5x,'U n i t    M u l t i p l e r s'//
     &       10x,'Length = ',1p,1e12.4/
     &       10x,'Force  = ',1p,1e12.4/
     &       10x,'Time   = ',1p,1e12.4/)

      end
