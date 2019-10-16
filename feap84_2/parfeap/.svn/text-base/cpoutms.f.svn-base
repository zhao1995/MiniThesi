c$Id:$
      subroutine cpoutms(nsurf,stype,ics,dnope,nope,neps,
     &                   d, partn,nodx,revp, dflag, allfl)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    26/01/2013
c-----[--.----+----.----+----.-----------------------------------------]
c     Purpose: Output surface data to domain file

c     Inputs:
c        nsurf    - Surface number
c        stype    - Type of surface
c        ics(*)   - List of node numbers for facets
c        dnope    - Dimension of ICS
c        nope     - Number of facet nodes
c        neps     - Number of facets`
c        d        - Current domain number
c        partn(*) - Partition numbers for nodes
c        nodx(*)  - Node renumber list
c        revp(4,*)- Domain data

c     Outputs:
c        dflag    - Domain nodes on slave surface if true
c-----[--.----+----.----+----.-----------------------------------------]
      implicit   none

      include   'cdata.h'
      include   'iodata.h'

      character  stype*8
      logical    dflag, eflag, allfl
      integer    nsurf, dnope,nope, neps, d, n, e,ee, ne
      integer    ics(dnope,neps), partn(*), nodx(0:*), revp(4,*)

c     Check if any facet is in current domain

      dflag = .false.
      do n = 1,neps
        do e = 1,nope
          if(ics(e,n).gt.0) then
            if(partn(ics(e,n)) .eq. d         .and.
     &         nodx(ics(e,n))  .le. revp(1,d)) then
              dflag = .true.
              go to 100
            endif
          endif
        end do ! e
      end do ! n

c     Output facet list

100   if(dflag .or. allfl) then
        write(ios,2000) nsurf,stype
        ne = 0
        ee = nope
        do n = 1,neps
          eflag = allfl
          do e = 1,nope
            if(ics(e,n).gt.0) then
              ee = e
              if(partn(ics(e,n)) .eq. d         .and.
     &           nodx(ics(e,n))  .le. revp(1,d)) then
                eflag = .true.
              endif
            endif ! ics(e,n) > 0
          end do ! e
          if(eflag) then
            ne = ne + 1
            write(ios,2001) ne,0,(nodx(ics(e,n)),e=1,ee)
     &                          ,(0,e=ee+1,nope)
          endif ! eflag
        end do ! n
      endif ! dflag

c     Formats

2000  format(/'SURFACE',i5/'  ',a/'    FACET')
2001  format(i8,i2,8i8:/(10i8:))

      end
