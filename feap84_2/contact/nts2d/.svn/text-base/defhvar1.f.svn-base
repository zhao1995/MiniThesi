c$Id:$
      subroutine defhvar1 (w1,w3)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Coded/modified by:                           Date:          rel.:
c               Giorgio Zavarise           April 10, 1996            1.0

c      Acronym: DEFine History VARiables # 1

c      Purpose: Define history variables for the specific pair
c               Variable definitions for NTS2D contact driver

c      Inputs :

c      Outputs:
c         w1(*)   - Dictionary of variables for CH1 & CH2
c         w3(*)   - Dictionary of variables for CH3
c-----[--.----+----.----+----.-----------------------------------------]
c     History variables for ch1 & ch2

c      Geometry
c         masts   - MASTter Segment
c         lnc     - Local Node Closest
c         istgt   - Index of STatus for GT
c         istgn   - Index of STatus for GN
c         s21     - Sinus    of unit vector tangent to master segm 12
c         c21     - Cosinus  of unit vector tangent to master segm 12
c         d21     - Distance of nodes 1 and 2 (length of segment)
c         csi     - normalized projection of slavenode on segment
c         gn      - Normal projection Gap  (+ if open)
c         gt      - Tang.  projection
c         area    - Area associated to the slavenode
c         dtd     - td increment from previous step
c         s21c    - Sinus of tangent vector for slavenode in Corner
c         c21c    - Cosinus of tangent vector for slavenode in Corner
c         csic    - csi for slavenode in Corner
c         td      - Tangential Displacement

c      Stiffness
c         fn      - Force Normal
c         dgnfn   - Derivative with respect to gn of fn
c         ft      - Tangential force
c         dtdft   - Derivative with respect to gt of ft
c         dgnft   - Derivative with respect to gn of ft
c         tde     - Stick part of dgt
c         tdp     - slip part of dgt
c         istfr   - Index of STatus for Friction

c      History variables for ch3

c      Geometry
c         d21oi   - d21 at the previous step to be used within iteration
c         csibb   - csi at the initial contact
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      character w1(*)*8,w3(*)*8

      save

      call cdebug0 ('  defhvar1',-1)


c     ADD NEW VARIABLES HERE - DO NOT CHANGE EXISTING DEFINITIONS

c     CH1 & CH2 VARIABLES (dynamic, CH2 copied in CH1 each time step)

c     Group of geometrical variables       (  1 -  50)

      w1(1)   = 'masts'
      w1(2)   = 'lnc'
      w1(3)   = 'istgt'
      w1(4)   = 'istgn'
      w1(5)   = 's21'
      w1(6)   = 'c21'
      w1(7)   = 'd21'
      w1(8)   = 'csi'
      w1(9)   = 'gn'
      w1(10)  = 'gt'
      w1(11)  = 'area'
      w1(12)  = 'dtd'
      w1(13)  = 's21c'
      w1(14)  = 'c21c'
      w1(15)  = 'csic'
      w1(16)  = 'td'
      w1(17)  = 'csibb'
      w1(21)  = 'lagmu'
      w1(22)  = 'lagmv'

c     Group of stiffness variables         ( 51 - 100)

      w1(51)  = 'fn'
      w1(52)  = 'dgnfn'
      w1(53)  = 'ft'
      w1(54)  = 'dtdft'
      w1(55)  = 'dgnft'
      w1(56)  = 'tde'
      w1(57)  = 'tdp'
      w1(58)  = 'istfr'

c     Group of thermal variables           (101 - 120)

      w1(101)  = '     '

c     Group of thermal stiffness variables (121 - 150)

      w1(121)  = '     '

c     Group of augmentation variables      (151 - 170)

      w1(151)  = 'augfn'
      w1(152)  = 'daugfn'
      w1(153)  = 'augft'

c     Other groups of variables            (171 - 200)

      w1(171)  = '      '


c     CH3 VARIABLES (static, never automatically modified)

c     Group of geometrical variables       (  1 -  20)

      w3(1)   = 'd21oi'
      w3(2)   = 'csioi'
      w3(3)   = 'fuoi'
      w3(4)   = 'fnmax'
      w3(5)   = 'fnresid'
      w3(6)   = 'fnaug'

c     Group of stiffness variables         ( 21 -  40)

      w3(21)  = 'lagmn'
      w3(22)  = 'lagmt'

c     Group of thermal variables           ( 41 -  60)

      w3(41)  = '      '

c     Group of thermal stiffness variables ( 61 -  80)

      w3(61)  = '      '

c     Group of augmentation variables      ( 81 -  90)

      w3(81)  = 'augft1'
      w3(82)  = 'augft2'
      w3(83)  = 'augstif1'
      w3(84)  = 'augstif2'
      w3(91)  = 'augfn1'
      w3(92)  = 'gn1'
      w3(93)  = 'augfn2'
      w3(94)  = 'gn2'

c     Other groups of variables            ( 91 - 100)

      end
