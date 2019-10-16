c$Id:$
      subroutine crmat02 (cm)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Assign length to cc and c2 of 8 characters       17/04/2007
c-----[--.----+----.----+----.-----------------------------------------]
c      Coded/modified by:                           Date:          rel.:
c               Giorgio Zavarise           April 10, 1996            1.0

c      Acronym: Contact Read MATerial properties # 1

c      Purpose: Input of contact material properties

c      Inputs:
c         nsubc   - # of the SUB-Command
c         scdat   - Sub-Command DATa read
c                 - Read from file

c      Outputs:
c         cm(*)   - Contact materials data storage
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'iofile.h'
      include  'print.h'

      logical   errck,gettxtd,pcomp
      character tx(2)*15,cc*8,c2*8
      integer   i
      real*8    cm(*), td(16)

      save

      call cdebug0 ('      crmat01',-1)

c     Get command

      errck = gettxtd(tx,2,td,14,'skip')
      cc    = tx(1)(1:8)

      do while(.not.pcomp(cc,'    ',4))

c       Coulomb friction coefficient

        if(pcomp(cc,'fric',4))then
          c2 = tx(2)(1:8)
          if(pcomp(c2,'coul',4))then
            cm(1) = td(1)
          else
            write (  *,3000) cc,c2
            write (ilg,3000) cc,c2
            call coutcis(2)
            call plstop()
          endif

          if (prt) then
            write (iow,2001) c2,td(1)
          endif
        elseif(pcomp(cc,'nonl',4))then
          c2 = tx(2)(1:8)
          if(pcomp(c2,'coul',4))then
            cm(2) = td(1)
            cm(3) = td(2)
            cm(4) = td(3)
            cm(5) = td(4)
            cm(6) = td(5)
          else
            write (  *,3000) cc,c2
            write (ilg,3000) cc,c2
            call plstop()
          endif

          if (prt) then
            write (iow,2002) (td(i),i=1,5)
          endif
        else
          continue                ! other constants
        endif

        errck = gettxtd(tx,2,td,14,'skip')
        cc    = tx(1)(1:8)
      end do

c     Formats

2001  format (/5x,'FRIC: Data for Friction'/
     &         5x,'  Friction Model                  ',a/
     &         5x,'  Friction Coefficient            ',1p,e12.5)

2002  format (/5x,'NONL: Data for Friction'/
     &         5x,'  Non-linear Coefficient 1        ',1p,e12.5/
     &         5x,'  Non-linear Coefficient 2        ',1p,e12.5/
     &         5x,'  Non-linear Coefficient 3        ',1p,e12.5/
     &         5x,'  Non-linear Coefficient 4        ',1p,e12.5/
     &         5x,'  Non-linear Coefficient 5        ',1p,e12.5)

3000  format(' *ERROR* CRMAT02: Cont. mat. sub-option not defined ',
     &       a,' ',a)

      end
