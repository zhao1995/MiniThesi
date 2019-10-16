c$Id:$
      subroutine psweep(swang, nsinc,nxd,nxn)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--+---------+---------+---------+---------+---------+---------+-]
c      Purpose: Plot sweep of 2-d mesh through angle 'swang'

c      Inputs:
c        swang   - Sweep angle in degrees
c        nsinc   - Sweep increments
c        nxd     - Dimension of ix array
c        nxn     - Number of nodes on array

c      Outputs:
c        to screen
c-----[--+---------+---------+---------+---------+---------+---------+-]
      implicit   none

      include   'pointer.h'
      include   'comblk.h'

      integer    nsinc,nxd,nxn
      real*8     swang

      save

      call pppcol(1,1) ! White

      call pswsub(mr(plix), mr(np(78)), mr(np(62)), hr(np(53)),
     &            nxd,nxn, swang, nsinc )

      end
