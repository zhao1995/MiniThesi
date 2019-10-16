c$Id:$
      integer function scp (ncom,nsuc)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Coded/modified by:                           Date:          rel.:
c               Giorgio Zavarise           April 10, 1996            1.0

c      Acronim: Sub-Command Pointer

c      Purpose: Set appropriate pointer for SUB-COMMANDS in CIS array

c      Inputs :
c         ncom    - # of the command
c         nsuc    - # of the sub-command

c      Outputs:
c         scp     - absolute pointer
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'c_0.h'
      include  'c_comnd.h'
      include  'c_dicti.h'

      integer   ncom,nsuc, ofs,k

      save

      ofs = ofssc
      do k = 1,ncom-1
        ofs = ofs + nsc(k)
      end do
      scp = ofs + nsuc

      end
