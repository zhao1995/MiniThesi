c$Id:$
      subroutine c3geoplt (ix1,ix2,pen1,pen2)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Coded/modified by:                           Date:          rel.:
c               Giorgio Zavarise           April 10, 1996            1.0

c      Acronym: Contact 2d GEOmetry PLot

c      Purpose: Plot geometry of 2D contact pair

c      Inputs:
c         x(*)    - Nodal coordinates in deformed position
c         ix1(*)  - Element nodal connection list for surface 1
c         ix2(*)  - Element nodal connection list for surface 2

c      Outputs:
c         none      - Plot outputs to screen/file
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'pointer.h'
      include  'comblk.h'

      logical   ifplt
      integer   ix1(*),ix2(*),pen1,pen2, ifsurf

      save

      call cdebug0 ('      c3geoplt',-1)

c     Get plot flag

      call setcplt (ifplt,ifsurf)

      if (ifplt) then
        call  c3dplot (hr(np(53)),ix1,ix2,pen1,pen2,ifsurf)
      endif

      end
