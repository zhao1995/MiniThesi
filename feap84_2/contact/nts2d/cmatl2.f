c$Id:$
      subroutine cmatl2 (cp0,cm,ch1,ch2,ch3)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Coded/modified by:                           Date:          rel.:
c               Giorgio Zavarise           April 10, 1996            1.0

c      Acronym: Contact MATerial Law 2
c      Purpose: compute parameters related to contact material model

c      Inputs :
c         cp0(*)  - Contactpair control data
c         cm(*)   - Contact materials data storage
c         ch1(*)  - Contact history variables (old)
c         ch2(*)  - Contact history variables (current)

c      Outputs:
c         ch2(*)  - Contact history variables (current)
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'c_0.h'
      include  'c_comnd.h'
      include  'c_contac.h'
      include  'c_keyh.h'
      include  'c_mate.h'
      include  'c_pair.h'
      include  'counts.h'
      include  'iofile.h'

      real*8    cp0(nr0,n0c3:*),cm(*),ch1(*),ch2(*),ch3(*)

      save

      call cdebug0 ('      cmatl2',-1)

      write(*,*) ' Material Model 2 not Included'

      end
