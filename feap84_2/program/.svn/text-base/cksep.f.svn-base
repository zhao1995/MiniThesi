c$Id:$
      logical function cksep(x1)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Check for existence of separator characters in data.

c      Inputs:
c         x1  -  Character to check

c      Outputs:
c         cksep - True of character is a valid separator; else false.
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      character x1*1

      save

c     Input character separators are blank, comma, or equal sign

      cksep = (x1.eq.' ') .or. (x1.eq.',') .or. (x1.eq.'=')

      end
