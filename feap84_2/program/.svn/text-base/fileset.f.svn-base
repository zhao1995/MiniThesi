c$Id:$
      subroutine fileset(fout,flog,old, add)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Set first character of filename and destroy old one

c      Inputs:
c        fout    - Name of output file
c        old     - Old character
c        add     - New character

c      Outputs:
c        flog    - Name of log file
c-----[--.----+----.----+----.-----------------------------------------]
      implicit   none

      include   'iofile.h'

      character  fout*(*), flog*(*), old*1, add*1
      logical    initf
      integer    l,ll

      save

      flog = fout
      ll   = len(flog)
      do l = ll,2,-1
        if(flog(l:l).eq.old .and. flog(l-1:l-1).eq.'/') then
          flog(l:l) = add
          go to 100
        endif
      end do ! l
      flog(1:1) = add
100   inquire(file=flog,exist=initf)
      if(initf) then
        open (unit=ilg,file=flog,status='old')
        close(unit=ilg,          status='delete')
      endif

      end
