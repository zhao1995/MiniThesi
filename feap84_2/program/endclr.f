c$Id:$
      subroutine    endclr (subnam,chr)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: End-of-file clearing routine

c      Inputs:
c         subnam - Character array storing calling subroutine name

c      Outputs:
c         chr    - Blank character to clear error
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'iofile.h'

      character subnam*(*),chr*(*)

      save

      if (ior.gt.0)  then
        write(iow,3000) subnam,ior
        write(ilg,3000) subnam,ior
        call plstop()
      else
        write(*,3000) subnam,ior
        chr = ' '
      endif

c     Format

 3000 format (' *ERROR* ENDCLR: End of file encountered in ',a/
     &        '         Unit Number =',i4/)

      end
