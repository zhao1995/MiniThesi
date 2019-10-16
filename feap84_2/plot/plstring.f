c$Id:$
      subroutine plstring(icol)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    05/02/2007
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Place string on plots

c      Inputs:
c         icol      - Color of text

c      Outputs:
c         none      - Plot outputs to screen/file
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'comfil.h'
      include  'pdatxt.h'

      integer   icol

      save

c     Draw line across bottom to box plot

      call plopen
      call pppcol(icol,1)
      call dplot( 0.00d0, 0.05d0, 3 )
      call dplot( 0.97d0, 0.05d0, 2 )

c     Place text on plot

      dtext = 0.0d0
      call pltsiz(2)
      call tplot(0.05d0, 0.02d0,record,60,0)
      call pltsiz(1)

      end
