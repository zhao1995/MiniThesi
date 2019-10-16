c$Id:$
      integer function cvtc2i(cdev)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Convert a character to an integer
c               ASCII character codes assumed

c      Inputs:
c         cdev    - Character array to convert

c      Outputs:
c         cvtc2i  - Integer value of character
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      character cdev(12)*1
      integer   i,n,nz

      save

c     Convert character to numerical

      nz= ichar('0')
      n = 0
      do i = 1,12
        if(cdev(i).eq.char(0) .or. cdev(i).eq.' ') go to 200
        n = 10*n + (ichar(cdev(i)) - nz)
      end do ! i

200   cvtc2i = n

      end
