c$Id:$
      subroutine pmessage()

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Add include setups to check for main problem     08/04/2011
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Place message on start page

c      Inputs:

c      Outputs:
c-----[--.----+----.----+----.-----------------------------------------]
      implicit   none

      include   'setups.h'

c     Place any messages for display here

      if(rank.eq.0) write(*,2000)

2000  format(/10x,'--> Please report errors by e-mail to:'/
     &        10x,'    feap@ce.berkeley.edu '/)

      end
