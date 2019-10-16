c$Id:$
      integer function sop (ncom,nsuc,nsop)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Coded/modified by:                           Date:          rel.:
c               Giorgio Zavarise           April 10, 1996            1.0

c      Acronim: Sub-command Options Pointer

c      Purpose: Set appropriate pointer for sub-comm. OPTIONS in CIS
c               array

c      Inputs :
c         ncom    - # of the command
c         nsoc    - # of the sub-command option

c      Outputs:
c         sop     - absolute pointer
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'c_0.h'
      include  'c_comnd.h'
      include  'c_dicti.h'

      integer   ncom,nsuc,nsop, ofs,k

      save

      ofs = ofsso
      do k = 1,ncom-1
        ofs = ofs + nso(k)*nsc(k)
      end do
      sop = ofs + nso(ncom)*(nsuc-1) + nsop

      end
