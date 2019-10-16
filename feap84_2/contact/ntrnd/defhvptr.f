c$Id:$
      subroutine defhvptr (w1,w3)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Coded/modified by:                           Date:          rel.:

c      Acronym: DEFine History Variables for contact element 1

c      Purpose: Define history variables for the specific pair
c               Variable definitions for CPTRND contact driver

c      Inputs :

c      Outputs:
c         w1(*)   - Dictionary of variables for CH1 & CH2
c         w3(*)   - Dictionary of variables for CH3
c-----[--.----+----.----+----.-----------------------------------------]
c     History variables for ch1 & ch2

c      Geometry
c         istgn   - Index of STatus for GN
c         gn      - Normal projection Gap  (+ if open)
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      character w1(*)*8,w3(*)*8

      save

      call cdebug0 ('  defhvptr',-1)


c     ADD NEW VARIABLES HERE - DO NOT CHANGE EXISTING DEFINITIONS

c     CH1 & CH2 VARIABLES (dynamic, CH2 copied in CH1 each time step)

c     Group of geometrical variables       (  1 -  50)

      w1(1)   = 'istgn'
      w1(2)   = 'gn'
      w1(21)  = 'lagmu'
      w1(22)  = 'lagmv'

c     Group of thermal variables           (101 - 120)

      w1(101)  = '     '

c     Group of thermal stiffness variables (121 - 150)

      w1(121)  = '     '

c     Group of augmentation variables      (151 - 170)

      w1(151)  = 'augfn'

c     Other groups of variables            (171 - 200)

      w1(171)  = '      '


c     CH3 VARIABLES (static, never automatically modified)

c     Group of geometrical variables       (  1 -  20)

      w3(1)   = '     '

c     Group of stiffness variables         ( 21 -  40)

      w3(21)  = 'lagmn'
      w3(22)  = 'lagmt'

c     Group of thermal variables           ( 41 -  60)

      w3(41)  = '      '

c     Group of thermal stiffness variables ( 61 -  80)

      w3(61)  = '      '

c     Other groups of variables            ( 91 - 100)

      end
