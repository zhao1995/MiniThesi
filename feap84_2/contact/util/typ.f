c$Id:$
      integer function typ(ncom,ntyp)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Coded/modified by:                           Date:          rel.:
c               Giorgio Zavarise           April 10, 1996            1.0

c      Acronim: TYpe command Pointer

c      Purpose: Set the appropriate pointer for TYPES in CIS array

c      Inputs :
c         ncom    - # of the command
c         ntyp    - # of the type

c      Outputs:
c         typ     - absolute pointer
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'c_0.h'
      include  'c_comnd.h'
      include  'c_dicti.h'

      integer   ncom,ntyp, ofs,k

      save

      ofs = 0
      do k = 1,ncom-1
        ofs = ofs + nty(k)
      end do
      typ = ofs + ntyp

      end
