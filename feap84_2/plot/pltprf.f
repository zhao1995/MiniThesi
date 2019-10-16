c$Id:$
      subroutine pltprf(jp,neq,lower)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c     Purpose: Display plot of profile layout

c     Inputs:
c        jp(*)   - Column pointers
c        neq     - Number of equations

c     Outputs:
c        none
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'eqsym.h'
      include  'pdata1.h'

      logical   lower
      integer   n,neq,jp(neq)
      real*8    x0,y0, x,y,c, pfact,rfact,ptone,ptnin,onept

      save

c     Set plot factors and start coordinate

      rfact = 0.5d0/(scaleg*fact)
      pfact = 0.8d0/dble(neq)*rfact

      x0    = 0.5d0*sx(1) - 0.5d0*rfact
      y0    = 0.5d0*sx(2) - 0.5d0*rfact

      ptone = 0.1d0*rfact
      onept = 1.0d0*rfact
      ptnin = 0.9d0*rfact

c     Upper part

      call pppcol(3,1)
      do n = 2,neq
        x = ptone + dble(n)*pfact
        y = onept - x + y0
        c = dble(jp(n) - jp(n-1))*pfact
        x = x + x0

        call plotl( x    , y    , 0.0d0, 3)
        call plotl( x    , y + c, 0.0d0, 2)

      end do ! n

c     Lower part

      if(lower) then
        call pppcol(4,1)
        do n = 2,neq
          x = ptone + dble(n)*pfact
          y = onept - x + y0
          c = dble(jp(n) - jp(n-1))*pfact
          x = x + x0

          call plotl( x - c, y    , 0.0d0, 3)
          call plotl( x    , y    , 0.0d0, 2)

        end do ! n
      end if

c     Diagonal

      call pppcol(2,1)
      call plotl( ptone + pfact + x0, ptnin - pfact + y0, 0.0d0, 3)
      call plotl( ptnin         + x0, ptone         + y0, 0.0d0, 2)

      end
