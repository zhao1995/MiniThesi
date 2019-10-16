c$Id:$
      subroutine defhvt2 (ww1,ww3)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Coded/modified by:                           Date:          rel.:
c               Robert L. Taylor             6 March 2003            1.0

c      Acronym: DEFine History VARiables for Tied 2d problems

c      Purpose: Define history variables for specific pair

c      Inputs :

c      Outputs:
c         ww1(*)   - Dictionary of variables for CH1 & CH2
c         ww3(*)   - Dictionary of variables for CH3
c-----[--.----+----.----+----.-----------------------------------------]
c     History variables for ch1 & ch2

c      Geometry
c         mastl   - MASTter Segment list
c         xi      - Natural coords of quadrature
c         g0      - Gap coordinates in reference state

c      Stiffness
c         lagm    - LAGrange Multiplier

c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      character ww1(*)*8,ww3(*)*8

      save

      call cdebug0 ('  defhvt2',0)

c     Variable definitions for the contact driver

c     ch1 & ch2 variables (dynamic, ch2 copied to ch1 each time step)

      ww1(1) = 'lagm'

c     ch3 (static, never change)

      ww3(1) = 'mastl'
      ww3(2) = 'g0'
      ww3(3) = 'xim'
      ww3(4) = 'xis'
      ww3(5) = 'area'
      ww3(7) = 'lint'

      end
