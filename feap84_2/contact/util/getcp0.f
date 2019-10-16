c$Id:$
      subroutine getcp0(cp0,n1,n2,n3,value)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Coded/modified by:                           Date:          rel.:
c               R.L. Taylor              January 15, 2007            1.0

c      Acronym:  GET array CP0 values

c      Purpose: Input of contact type

c      Inputs :
c         n1           - first  index
c         n2           - second index
c         n3           - third  index
c         cp0(*,*,*)   - cp0 array

c      Outputs:
c         value        - cp0(n1,n2,n3)
c-----[--.----+----.----+----.-----------------------------------------]
      implicit   none

      include   'c_0.h'
      include   'c_comnd.h'

      integer    n1,n2,n3,value
      real*8     cp0(nr0,n0c3:nc03,*)

      value = nint(cp0(n1,n2,n3))

      end
