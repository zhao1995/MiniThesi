c$Id:$
      subroutine storec(ixc,ncen1,ixl,nnod,numcels)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Store local connection in global array

c      Inputs :

c      Outputs:
c-----[--.----+----.----+----.-----------------------------------------]
      implicit none

      integer  ncen1,nnod,numcels,i
      integer  ixc(ncen1,*),ixl(*)

      save

      do i = 1,nnod
        ixc(i,numcels) = ixl(i)
      end do

      end
