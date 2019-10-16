c$Id:$
      subroutine crel04(ics,ix,ip,regn,dnope,emax)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Coded/modified by:                           Date:          rel.:
c               Robert L. Taylor         October 10, 1996            1.0

c      Acronym: Contact SURfaces search for REGIons

c      Purpose: Set facets to left and right of current one

c      Inputs :
c         ics(dnope,*) - Facet node connections
c         ix(nen1,*)   - Element connection list
c         regn         - Region Number to use
c         dnope        - Dimension NOdes Per Element

c      Working array:
c         ip(*)        - Marker for region nodes

c      Outputs:
c         ics(dnope,*) - Point elements for region
c         emax         - No. Elements on Surface
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'cdata.h'
      include  'iofile.h'
      include  'sdata.h'

      integer   regn,dnope,emax,ix(nen1,*),ip(*),ics(dnope,*), i,n

      save

c     Set initial flags

      do n = 1,numnp
        ip(n) = 0
      end do ! n

c     Search elements for active region

      do n = 1,numel
        if(ix(nen1-1,n).eq.regn) then
          do i = 1,nen
            if(ix(i,n).gt.0) then
              ip(ix(i,n)) = 1
            endif
          end do ! i
        endif
      end do !

c     Set the surface nodes

      emax = 0
      do n = 1,numnp
        if(ip(n).gt.0) then
          emax        = emax + 1
          ics(1,emax) = n
        endif
      end do ! n

      end
