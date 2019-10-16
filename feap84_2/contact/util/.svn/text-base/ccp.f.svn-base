c$Id:$
      function ccp(ncom)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Coded/modified by:                           Date:          rel.:
c               Giorgio Zavarise           April 10, 1996            1.0

c      Acronim: Contact Command Pointer

c      Purpose: Set appropriate pointer for command control table

c      Inputs :
c         ncom    - # of the command

c      Outputs:
c         ccp     - absolute pointer
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'c_0.h'
      include  'c_comnd.h'
      include  'c_ccp.h'
      include  'pointer.h'

      integer   ncom

      save

      ccp = np(131) + of0(ncom) - 1

      end
