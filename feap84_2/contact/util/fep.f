c$Id:$
      integer function fep (ncom,nfea)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Coded/modified by:                           Date:          rel.:
c               Giorgio Zavarise           April 10, 1996            1.0

c      Acronim: FEature command Pointer

c      Purpose: Set the appropriate pointer for FEATURES in CIS array

c      Inputs :
c         ncom    - # of the command
c         nfea    - # of the feature

c      Outputs:
c         fep     - absolute pointer
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'c_0.h'
      include  'c_comnd.h'
      include  'c_dicti.h'

      integer   ncom,nfea, ofs,k

      save

      ofs = ofsfe
      do k = 1,ncom-1
        ofs = ofs + nfe(k)
      end do
      fep = ofs + nfea

      end
