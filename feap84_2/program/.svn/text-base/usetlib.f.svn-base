c$Id:$
      subroutine usetlib(i)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Interface for user mesh manipulation set  commands

c      Inputs:
c         i      - Command number

c      Outputs:
c         None   - Users are responsible for providing outputs in
c                  umanii routines
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      integer   i

      save

      if(i.eq.1) then
        call umani1
      elseif(i.eq.2) then
        call umani2
      elseif(i.eq.3) then
        call umani3
      elseif(i.eq.4) then
        call umani4
      elseif(i.eq.5) then
        call umani5
      elseif(i.eq.6) then
        call umani6
      elseif(i.eq.7) then
        call umani7
      elseif(i.eq.8) then
        call umani8
      elseif(i.eq.9) then
        call umani9
      elseif(i.eq.10) then
        call umani0
      endif

      end
