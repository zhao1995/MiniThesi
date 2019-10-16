c$Id:$
      integer function ipos(file,nn)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Locate last character in character array

c      Inputs:
c         file(*) - Array to search
c         nn      - Length of array

c      Outputs:
c         ipos   - Position of last character
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      integer   n,nn
      character file(nn)*1

      save

      do n = nn,1,-1
        if(file(n).ne.' ') go to 100
      end do ! n
      n    = 0

100   ipos = n

      end
