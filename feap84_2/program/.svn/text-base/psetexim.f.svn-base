c$Id:$
      subroutine psetexim(ix,nen,nen1,numel)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    18/05/2012
c      1. Change ix(nen+7 to ix(nen+6: implicit/explicit    01/07/2013
c      2. Correct call to tplot for td in 'elem' part       15/08/2013
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Set elements for implicit-explicit solutions
c               ix(nen+6,*) > 0 for explicit

c      Inputs:
c        ix(nen1,*) - Element data array
c        nen        - Number of nodes/element (max)
c        nen1       - Dimension of array
c        numel      - Number of elements

c      Outputs:
c        Set of explicit elements for solutions
c-----[--.----+----.----+----.-----------------------------------------]
      implicit   none

      include   'ddata.h'
      include   'iofile.h'
      include   'print.h'

      integer    nen,nen1,numel
      integer    ix(nen1,numel)

      logical    tinput, errck, pcomp
      character  type(2)*15
      integer    n,ma,it
      real*8     td(6)

c     Input of explict element sets

      type(1) = 'start'
      do while (.not.pcomp(type(1),'    ',4))
c       Input: mate ma it
        errck = tinput(type(1),1,td,2)
        if(pcomp(type(1),'mate',4)) then
          ma = nint(td(1))
          it = nint(td(2))
          do n = 1,numel
            if(ix(nen1,n).eq.ma) then
              ix(nen+6,n) = it + 1
              if(prt) write(iow,2001) n,it
            endif
          end do ! n
        elseif(pcomp(type(1),'elem',4)) then
          type(2) = 'start'
          do while(.not.pcomp(type(2),'    ',4))
            errck = tinput(type(2),0,td,2)
            n  = nint(td(1))
            it = nint(td(2))
            ix(nen+6,n) = it + 1
            if(prt) write(iow,2001) n,it
          end do ! while
        elseif(pcomp(type(1),'auto',4)) then
          imexfl = .true.
        elseif(pcomp(type(1),'off',4)) then
          imexfl = .false.
        endif
      end do ! while

c     Formats

2001  format(5x,'Implicit/Explicit Element =',i8,' Iterate =',i3)

      end
