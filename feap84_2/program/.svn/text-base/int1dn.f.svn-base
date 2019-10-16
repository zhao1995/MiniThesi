c$Id:$
      subroutine int1dn(l,sw)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c     1. Use 'one3' from 'pconstant.h' instead of 'third'   14/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Nodal quadrature for 1-d element

c      Inputs:
c         l     - Number of points

c      Outputs:
c         sw(1,*) - Gauss points
c         sw(2,*) - Gauss weights
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'iofile.h'
      include  'pconstant.h'

      integer   l
      real*8    sw(2,*)

      save

c     Set end points

      sw(1,1) = -1.d0
      sw(1,2) =  1.d0

      if(l.eq.2) then

        sw(2,1) =  1.0d0
        sw(2,2) =  1.0d0

      elseif(l.eq.3) then

        sw(1,3) =  0.0d0
        sw(2,1) =  one3
        sw(2,2) =  one3
        sw(2,3) =  one3*4.d0

      elseif(l.eq.4) then

        sw(1,3) = -one3
        sw(1,4) =  one3
        sw(2,1) =  0.25d0
        sw(2,2) =  0.25d0
        sw(2,3) =  0.75d0
        sw(2,4) =  0.75d0

      else
        write(iow,2000) l
        write(ilg,2000) l
        if(ior.lt.0) then
          write(*,2000) l
        endif
        call plstop()

      endif

c     Format

2000  format(' *ERROR* INT1DN: Illegal quadrature order =',i4)

      end
