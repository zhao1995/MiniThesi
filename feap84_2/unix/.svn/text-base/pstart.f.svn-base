c$Id:$
      subroutine pstart()

c      * * F E A P * * A Finite Element Analysis Program
c                        -      -       -        -
c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Set 'reflg' to false                             12/10/2007
c       2. Remove arguments from pinitm call                09/01/2012
c-----[--.----+----.----+----.-----------------------------------------]
c     Purpose:

c     Inputs:

c     Outputs:
c-----[--.----+----.----+----.-----------------------------------------]
      implicit      none

      include      'memuse.h'
      include      'pfeapb.h'
      include      'prmptd.h'
      include      'setups.h'
      include      'comblk.h'

      save

c     Set flags for serial execution

      pfeap_on   = .false.
      pfeap_gnod = .false.
      ntasks     =  1
      rank       =  0

c     Set maximum memory use

      maxuse = 0

c     Start for X11 graphics driver

      call pdriver()

c     Set for file checking at startup

      fileck = .true.

c     Set restart flag to false (reset on command line inputs only)

      reflg  = .false.

c     Initialize memory

      call pinitm()

c     Check user installation options

      call pinstall()

c     Set initial file names

      call filnam()

      end
