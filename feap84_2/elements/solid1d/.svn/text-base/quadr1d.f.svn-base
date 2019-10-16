c$Id:$
      subroutine quadr1d(d)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/04/2009
c       1. Set 'npm' for mixed elements                     15/08/2009
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose:  1-D quadrature

c      Inputs:
c         d(*)        - Material set parameters

c      Outputs:
c         sg1(2,*)    - Quadrature points
c         shp1(2,*,*) - shape funcitons and first derivatives
c-----[--.----+----.----+----.-----------------------------------------]
      implicit none

      include  'eldata.h'
      include  'qudshp.h'

      real*8    d(*)

c     Set quadrature

      lint = min(5,nint(d(5)))
      if(lint.eq.0) then
        lint = nel
      endif
      if(nint(d(182)).eq.1) then
        call int1dn(lint,sg1)
      else
        call int1d(lint,sg1)
      endif
      quad = .true.
      npm  =  nel - 1

      end
