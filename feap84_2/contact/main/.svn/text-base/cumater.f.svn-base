c$Id:$
      subroutine cumater (cm)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Assign length to cc of 8 characters              17/04/2007
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
      include  'chdata.h'

      logical   errck,gettxtd,pcomp
      character tx(2)*15,cc*8
      integer   i
      real*8    cm(*), td(14)

      save

      call cdebug0 ('      cumater',-1)

c     Get command

      errck = gettxtd(tx,2,td,14,'skip')
      cc    = tx(1)(1:8)

      do while(.not.pcomp(cc,'    ',4))

c       Coulomb friction coefficient

        if(pcomp(cc,'prop',4))then
          do i = 1,14
            cm(i) = td(i)
          end do ! i

          if (prt) then
            write (iow,2001) (cm(i),i=1,14)
          endif
        else
          write (ilg,3001) xxx
          write (iow,3001) xxx
          if(ior.lt.0) then
            write (*,3001) xxx
          endif
          call plstop()
        endif

        errck = gettxtd(tx,2,td,14,'skip')
        cc    = tx(1)(1:8)
      end do

c     Formats

2001  format (/5x,'Material Parameters for User Model'/
     &         5x,'  Parameter  1 = ',1p,e12.5/
     &         5x,'  Parameter  2 = ',1p,e12.5/
     &         5x,'  Parameter  3 = ',1p,e12.5/
     &         5x,'  Parameter  4 = ',1p,e12.5/
     &         5x,'  Parameter  5 = ',1p,e12.5/
     &         5x,'  Parameter  6 = ',1p,e12.5/
     &         5x,'  Parameter  7 = ',1p,e12.5/
     &         5x,'  Parameter  8 = ',1p,e12.5/
     &         5x,'  Parameter  9 = ',1p,e12.5/
     &         5x,'  Parameter 10 = ',1p,e12.5/
     &         5x,'  Parameter 11 = ',1p,e12.5/
     &         5x,'  Parameter 12 = ',1p,e12.5/
     &         5x,'  Parameter 13 = ',1p,e12.5/
     &         5x,'  Parameter 14 = ',1p,e12.5/)

3001  format (/5x,' *ERROR* CUMATER: Reading material type data.'/
     &         5x,'  Unrecognized string: '/
     &         1x,a)

      end
