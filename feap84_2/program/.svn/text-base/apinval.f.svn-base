c$Id:$
      subroutine apinval(xs,val,error)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose:  Moves character string into real value

c      Inputs:
c         xs(*)   - Character string

c      Outputs:
c         val     - Value extracted from character string
c         error   - Flag, true if error occurs
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      logical   error
      character xs*25
      real*8    val

      save

      read(xs,'(f25.0)',err=100) val
      return
100   error = .true.

      end
