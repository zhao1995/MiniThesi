c$Id:$
      subroutine umaclib(i,lct,ct,uprt)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Delete 'prt' from all calls to umacr's           09/07/2009
c       2. Add capability for umacr11 to umacr20, and       25/01/2012
c          change umacr0 to umacr10
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Interface for user command language instructions

c      Inputs:
c         i      - Command number
c         lct    - Character array describing option
c         ct(3)  - Command parameters
c         uprt   - Output flag set to value

c      Outputs:
c         N.B.  Users are responsible for generating command options
c               See programmer manual for example.
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'print.h'

      logical   uprt, oprt
      character lct*(*)
      integer   i
      real*8    ct(3)

      save

      oprt = prt
      prt  = uprt
      if(    i.eq.1) then
        call umacr1(lct,ct)
      elseif(i.eq.2) then
        call umacr2(lct,ct)
      elseif(i.eq.3) then
        call umacr3(lct,ct)
      elseif(i.eq.4) then
        call umacr4(lct,ct)
      elseif(i.eq.5) then
        call umacr5(lct,ct)
      elseif(i.eq.6) then
        call umacr6(lct,ct)
      elseif(i.eq.7) then
        call umacr7(lct,ct)
      elseif(i.eq.8) then
        call umacr8(lct,ct)
      elseif(i.eq.9) then
        call umacr9(lct,ct)
      elseif(i.eq.10) then
        call umacr10(lct,ct)
      elseif(i.eq.11) then
        call umacr11(lct,ct)
      elseif(i.eq.12) then
        call umacr12(lct,ct)
      elseif(i.eq.13) then
        call umacr13(lct,ct)
      elseif(i.eq.14) then
        call umacr14(lct,ct)
      elseif(i.eq.15) then
        call umacr15(lct,ct)
      elseif(i.eq.16) then
        call umacr16(lct,ct)
      elseif(i.eq.17) then
        call umacr17(lct,ct)
      elseif(i.eq.18) then
        call umacr18(lct,ct)
      elseif(i.eq.19) then
        call umacr19(lct,ct)
      elseif(i.eq.20) then
        call umacr20(lct,ct)
      endif
      prt = oprt

      end
