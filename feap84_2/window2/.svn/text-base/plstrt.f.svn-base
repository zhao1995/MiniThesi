c$Id:$
      subroutine plstrt()

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Open a graphics device to receive plot data

c      Inputs:
c         none

c      Outputs:
c         none      - Graphics window should appear after call
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'pdata2.h'
      include  'plflag.h'

      integer   i, status,vopnwk,vclrwk

      save

c     Open devices for Graphics output

      everon = .true.

      idx = 22000
      idy = 22000

      status = vopnwk()
      status = vclrwk()
      i = 1
      call pppcol(i,1)
      call plbord(i)

      end
