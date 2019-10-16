c$Id:$
        subroutine fparsop(fname,ftemp,nch)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Change 'include/finclude' to 'finclude'          23/01/2009
c       2. Remove common 'pfeapa' (values in 'setups')      05/02/2009
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Set PostScript filenames to receive plot data

c      Inputs:
c         none

c      Outputs:
c         fname     - Final name of file
c         ftemp     - Name to temporarily hold data
c         nch       - Number of characters in ftemp filename
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

#     include  "finclude/petsc.h"
#     include  "setups.h"

      integer   nch
      character fname*(*), ftemp*17

      save

c     Set initial file name

      fname = 'Feap_0001AAAA.eps'
      ftemp = 'temp_0001.eps'
      nch   =  13
      if(rank.lt.9) then
        write(fname(9:9),'(i1)') rank+1
        write(ftemp(9:9),'(i1)') rank+1
      elseif(rank.lt.99) then
        write(fname(8:9),'(i2)') rank+1
        write(ftemp(9:9),'(i2)') rank+1
      elseif(rank.lt.999) then
        write(fname(7:9),'(i3)') rank+1
        write(ftemp(9:9),'(i3)') rank+1
      elseif(rank.lt.9998) then
        write(fname(6:9),'(i4)') rank+1
        write(ftemp(9:9),'(i4)') rank+1
      endif

      end
