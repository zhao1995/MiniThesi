c$Id:$
        subroutine uftyp03(uptyp,nel,iel)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    31/08/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Set plot data for user element

c      Inputs:
c         uptyp   - User element topology
c         nel     - Number of element nodes
c         iel     - Element routine number

c      Output:
c         inord   - Number of plot
c-----[--.----+----.----+----.-----------------------------------------]
      implicit   none

      include   'pdata6.h'

      integer    uptyp,nel,iel

c     Plot for perspective contours: Do not change

      if(nel.eq.4) then

        inord(iel) = 5

        ipord( 1,iel) = 1
        ipord( 2,iel) = 2
        ipord( 3,iel) = 4
        ipord( 4,iel) = 3
        ipord( 5,iel) = 1

c     Plot for mesh: Set inord(iel) and ipord(i,iel)

      else


      endif

      end
