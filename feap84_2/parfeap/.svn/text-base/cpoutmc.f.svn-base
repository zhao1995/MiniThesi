c$Id:$
      subroutine cpoutmc(ics,dnope,nope,neps,
     &                   d, partn, nodx, revp, dflag, allfl)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    26/01/2013
c-----[--.----+----.----+----.-----------------------------------------]
c     Purpose: Output flags for active domain data

c     Inputs:
c        ics(*)   - List of node numbers for facets
c        dnope    - Dimension of ICS
c        nope     - Number of facet nodes
c        neps     - Number of facets`
c        d        - Current domain number
c        partn(*) - Partition numbers for nodes
c        nodx(*)  - Renumber list
c        revp(4,*)- Domain data

c     Outputs:
c        dflag    - Domain nodes on slave surface if true
c        allfl    - All surface in same domain
c-----[--.----+----.----+----.-----------------------------------------]
      implicit   none

      include   'cdata.h'
      include   'iodata.h'

      logical    dflag, eflag, allfl
      integer    dnope,nope, neps, d, n, e, ne
      integer    ics(dnope,neps), partn(*), nodx(0:*), revp(4,*)

c     Check if any facet is in current domain

      dflag = .false.
      ne    = 0
      do n = 1,neps
        eflag = .false.
        do e = 1,nope
          if(ics(e,n).gt.0) then
            if(partn(ics(e,n)) .eq. d         .and.
     &         nodx(ics(e,n))  .le. revp(1,d)) then
              dflag = .true.
              eflag = .true.
            endif
          endif
        end do ! e
        if(eflag) then
          ne = ne + 1
        endif
      end do ! n

c     Set all in domain

      allfl = neps.eq.ne

      end
