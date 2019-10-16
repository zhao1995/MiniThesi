c$Id:$
      subroutine suchkant (knoten,flag,ix2,xi1,xi2,nel2,nr,zweiter)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Coded/modified by:                           Date:          rel.:
c               Anna Haraldsson             February 1998            1.0

c      Acronym:

c      Purpose: liefert -> Knoten zu Kante, auf die xi zeigt,
c                          flag   fuer xi zeigt nicht direkt auf Kante,
c                          (Dann einfach linkes Element durchsuchen)

c      Inputs:
c         ix2(*)  - Surface facet nodes
c         xi1     - Natural coordinate of surface point
c         xi2     - Natural coordinate of surface point
c         nel2    - Facet number of ix2
c         nr      -
c         zweiter -

c      Outputs:
c         knoten  -  node of edge, on which xi points
c         zweiter -
c         flag    - .true. if xi doesn't point directly to this edge
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'c_0.h'
      include  'c_geom.h'
      include  'c_keyh.h'
      include  'c_pair.h'
      include  'c_tole.h'

      logical   flag
      integer   knoten,ix2(dnope2,*),im(4),nel2,nr,nod2,zweiter
      real*8    xi1,xi2

      save

      flag=.false.

      do nod2 = 1,4
        im(nod2) = ix2(nod2,nel2)
      end do ! nod2

      if (im(1).eq.nr) then

        if (xi1.le.-1.d0-tlouts) then
          knoten = im(4)
        elseif(xi2.le.-1.d0-tlouts) then
          knoten = im(2)
        else
c         write(*,*)'Hupps, xi zeigt nach draussen? xi=',xi1,xi2
        endif

        if ((xi1.le.-1.d0-tlouts).and.(xi2.le.-1.d0-tlouts)) then
          flag    = .true.
          knoten  =  im(4)
          zweiter =  im(2)
        endif

      elseif (im(2).eq.nr) then

        if (xi1.ge.1.d0+tlouts) then
          knoten = im(3)
        elseif(xi2.le.-1.d0-tlouts) then
          knoten = im(1)
        else
c         write(*,*)'Hupps, xi zeigt nach draussen? xi=',xi1,xi2
        endif

        if ((xi1.ge.1.d0+tlouts).and.(xi2.le.-1.d0-tlouts))  then
          flag    = .true.
          knoten  =  im(3)
          zweiter =  im(1)
        endif

      elseif (im(3).eq.nr) then

        if (xi1.ge.1.d0+tlouts) then
          knoten = im(2)
        elseif(xi2.ge.1.d0+tlouts) then
          knoten = im(4)
        else
c         write(*,*)'Hupps, xi zeigt nach draussen? xi=',xi1,xi2
        endif

        if ((xi1.ge.1.d0+tlouts).and.(xi2.ge.1.d0+tlouts)) then
          flag    = .true.
          knoten  =  im(2)
          zweiter =  im(4)
        endif

      elseif (im(4).eq.nr) then

        if (xi1.le.-1.d0-tlouts) then
          knoten = im(1)
        elseif(xi2.ge.1.d0+tlouts) then
          knoten = im(3)
        else
c         write(*,*)'Hupps, xi zeigt nach draussen? xi=',xi1,xi2
        endif

        if ((xi1.le.-1.d0-tlouts).and.(xi2.ge.1.d0+tlouts)) then
          flag    = .true.
          knoten  =  im(1)
          zweiter =  im(3)
        endif

      endif

      end
