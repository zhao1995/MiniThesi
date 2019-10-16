c$Id:$
      subroutine coutms(n,ntype,ics,dnope,nope,neps)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Correct list of output types. Output from ics.   26/01/2013
c          Get type from 'ntype'.
c-----[--.----+----.----+----.-----------------------------------------]
c     Purpose: Output surface data to file

c     Inputs:
c        n       - Surface number
c        ntype   - Type of surface
c        ics(*)  - List of node numbers for facets
c        dnope   - Dimension of ICS
c        nope    - Number of facet nodes
c        neps    - Number of facets`

c     Outputs:
c-----[--.----+----.----+----.-----------------------------------------]
      implicit   none

      include   'cdata.h'
      include   'iodata.h'

      character  type(9)*8
      integer    n, ntype, dnope,nope, neps
      integer    ics(dnope,neps), i

      data       type / 'LINE    ','TRIANGLE','QUAD    ','BEAM    ',
     &                  'POINT   ','RIGID   ','NURBS   ','TSPLINE ',
     &                  'PART    ' /

      write(ios,2000) n,type(ntype)
      do n = 1,neps
        write(ios,2001) n,' 0 ',(ics(i,n),i=1,nope)
      end do ! n

c     Formats

2000  format(/'SURFACE',i5/'  ',a/'    FACET')
2001  format(i8,a,8i8)

      end
