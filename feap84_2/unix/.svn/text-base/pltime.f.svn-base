c$Id:$
      subroutine pltime()

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Place time on plot window

c      Inputs:  None

c      Outputs: To plot window
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'tdata.h'
      include  'pdatxt.h'

      character yy*15

      save

c     Display time for current view

      dtext = 0.00d0
      call pppcol(1,1)
      yy = ' '
      call tplot(1.12d0 , 0.135d0, yy, 15, 1)
      write(yy, '(a6,1p,1e9.2)' ) 'Time =',ttim
      call tplot(1.12d0 , 0.135d0, yy, 15, 1)

      end
