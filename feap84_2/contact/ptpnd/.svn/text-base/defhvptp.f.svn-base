c$Id:$
      subroutine defhvptp (w1,w3)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Coded/modified by:                           Date:          rel.:
c               Giorgio Zavarise & RLT      June 10, 1996            1.0

c      Acronym: DEFine History VARiables for PTP

c      Purpose: Define history variables for the specific pair

c      Inputs :

c      Outputs:
c         w1(*)   - Dictionary of variables for CH1 & CH2
c         w3(*)   - Dictionary of variables for CH3
c-----[--.----+----.----+----.-----------------------------------------]
c     History variables for ch1 & ch2

c      Geometry
c         masts   - MASTter Segment                              (nfacn)
c         lnc     - Local Node Closest                            (icls)
c         istgt   - Index of STatus for GT
c         istgn   - Index of STatus for GN
c         xi      - normalized projection slavenode on segment(xic,etac)
c         gn      - Normal projection Gap  (+ if open)            (gapn)
c         gt      - Tang.  projection
c         area    - Area associated to the slavenode
c         pslip   -      SLIP of tangent                         (pslip)
c         fdiss   - Friction DISSipation                         (fdiss)
c         wear    - WEAR on surface                               (wear)

c      Stiffness
c         fn      - Force Normal                                   (fn)
c         ft      - Tangential force (x,y,z - components)        (ftref)
c         istfr   - Index of STatus for Friction
c         lagmu   - LAGrange Multiplier (U=normal)               (lagmu)
c         lagmv   - LAGrange Multiplier (V=tangential)           (lagmv)

c     History variables for ch3
c         lagmn   - LAGrange Multiplier Normal                   (lagmn)
c         lagmt   - LAGrange Multiplier Tangential               (lagmt)

c      Geometry
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      character w1(*)*8,w3(*)*8

      save

      call cdebug0 ('  defhvptp',-1)

c     Variable definitions for the contact driver

c     ADD NEW VARIABLES - DO NOT CHANGE EXISTING DEFINITIONS
c     WARNING: Grouping is not necessary, but useful for clarity

c     CH1 & CH2 VARIABLES (dynamic, CH2 copied in CH1 at each time step)

c     group of geometrical variables       (  1 -  50)

      w1(1)   = 'masts'
      w1(2)   = 'lnc'
      w1(3)   = 'istgt'
      w1(4)   = 'istgn'
      w1(5)   = 'xi'
      w1(6)   = 'gn'
      w1(7)   = 'gt'
      w1(8)   = 'pslip'
      w1(9)   = 'fdiss'
      w1(10)  = 'wear'
      w1(21)  = 'lagmu'
      w1(22)  = 'lagmv'

c     group of stiffness variables         ( 51 - 100)

      w1(51)  = 'istfr'
      w1(52)  = 'fn'
      w1(53)  = 'ft'
      w1(54)  = 'kts'

c     group of thermal variables           (101 - 120)

      w1(101)  = '     '

c     group of thermal stiffness variables (121 - 150)

      w1(121)  = '     '

c     group of augmentation variables      (151 - 170)

      w1(151)  = 'augfn'
      w1(152)  = 'augft'

c     other groups of variables            (171 - 200)

      w1(171)  = '      '


c     CH3 VARIABLES (static, never automatically modified)

c     group of geometrical variables       (  1 -  20)

      w3(1)   = 'nfacn'
      w3(2)   = 'iclsn'
      w3(3)   = 'area'
      w3(4)   = 'gapn'
      w3(5)   = 'xin'

c     group of stiffness variables         ( 21 -  40)

      w3(21)  = 'lagmn'
      w3(22)  = 'lagmt'

c     group of thermal variables           ( 41 -  60)

      w3(41)  = '      '

c     group of thermal stiffness variables ( 61 -  80)

      w3(61)  = '      '

c     group of augmentation variables      ( 81 -  90)

      w3(81)  = '      '

c     other groups of variables            ( 91 - 100)

      w3(91)  = '      '

      end
