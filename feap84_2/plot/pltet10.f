c$Id:$
      subroutine pltet10(iel)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Set 3-D Plot Sequence for 10-node tetrahedra

c      Inputs:
c         iel       - Element number: > 0 for user    elements
c                                     < 0 for program elements

c      Outputs:
c         none      - Sequesnce returned in common /pdata6/
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'pdata5.h'
      include  'pdata6.h'

      integer   iel

      save

c     Set number of points

      if(iel.gt.0) then

        inord(iel)    = 15

c       Set plot sequence

        ipord( 1,iel) = 1
        ipord( 2,iel) = 5
        ipord( 3,iel) = 2
        ipord( 4,iel) = 6
        ipord( 5,iel) = 3
        ipord( 6,iel) = 7
        ipord( 7,iel) = 1
        ipord( 8,iel) = 8
        ipord( 9,iel) = 4
        ipord(10,iel) = 10
        ipord(11,iel) = 3
        ipord(12,iel) = 10
        ipord(13,iel) = 4
        ipord(14,iel) = 9
        ipord(15,iel) = 2

      elseif(iel.lt.0) then

        exord(-iel)    = 15

c       Set plot sequence

        epord( 1,-iel) = 1
        epord( 2,-iel) = 5
        epord( 3,-iel) = 2
        epord( 4,-iel) = 6
        epord( 5,-iel) = 3
        epord( 6,-iel) = 7
        epord( 7,-iel) = 1
        epord( 8,-iel) = 8
        epord( 9,-iel) = 4
        epord(10,-iel) = 10
        epord(11,-iel) = 3
        epord(12,-iel) = 10
        epord(13,-iel) = 4
        epord(14,-iel) = 9
        epord(15,-iel) = 2

      endif

      end
