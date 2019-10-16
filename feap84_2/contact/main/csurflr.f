c$Id:$
      subroutine csurflr(ics,dnope,neps)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Coded/modified by:                           Date:          rel.:
c               Robert L. Taylor         October 10, 1996            1.0

c      Acronym: Contact SURfaces search Left & Right neighbors

c      Purpose: Set facets to left and right of current one

c      Inputs :
c         ics(dnope,*) - Facet node connections
c         dnope        - Dimension NOdes Per Element
c         neps         - No. Elements Per Surface

c      Outputs:
c         ics(dnope,*) - Left/Right Neighbors for Facet
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'iofile.h'

      logical   flag
      integer   dnope,neps,ics(dnope,neps), ii,n,m

      save

c     Loop over elements

      do n = 1,neps
        ics(dnope-1,n) = 0
        ics(dnope  ,n) = 0
      end do

c     Loop over elements

      do n = 1,neps
        ii = ics(1,n)
        flag = .false.
        do m = 1,neps
          if(ii.eq.ics(2,m)) then
            if(.not.flag) then
              ics(dnope-1,n) =  m
              ics(dnope  ,m) =  n
              flag          = .true.
            else
              write(iow,3000) m,n,ii
              write(ilg,3000) m,n,ii
            endif
          endif
        end do
      end do

c     Format

3000  format (/' *ERROR* CSURFLR: Facet,',i8,' and',i8,
     &         ' have same node,',i8/
     &         '         which also appeared previously')

      end
