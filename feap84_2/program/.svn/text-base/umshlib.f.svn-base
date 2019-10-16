c$Id:$
      subroutine umshlib(i,tx,prt)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Add 'tx(*)' to argumant and umeshi calls         26/09/2011
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Interface for user mesh commands

c      Inputs:
c         i      - Command number
c         tx(*)  - Command line input data
c         prt    - Flag, output if true

c      Outputs:
c         None   - Users are responsible for providing outputs in
c                  umeshi routines
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      logical   prt
      character tx(*)*15
      integer   i

      save

      if(i.eq.1) then
        call umesh1(tx,prt)
      elseif(i.eq.2) then
        call umesh2(tx,prt)
      elseif(i.eq.3) then
        call umesh3(tx,prt)
      elseif(i.eq.4) then
        call umesh4(tx,prt)
      elseif(i.eq.5) then
        call umesh5(tx,prt)
      elseif(i.eq.6) then
        call umesh6(tx,prt)
      elseif(i.eq.7) then
        call umesh7(tx,prt)
      elseif(i.eq.8) then
        call umesh8(tx,prt)
      elseif(i.eq.9) then
        call umesh9(tx,prt)
      elseif(i.eq.10) then
        call umesh0(tx,prt)
      endif

      end
