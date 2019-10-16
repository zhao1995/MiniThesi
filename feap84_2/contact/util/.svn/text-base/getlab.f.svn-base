c$Id:$
      character*70 function getlab(ncs,ll)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Coded/modified by:                           Date:          rel.:
c               Giorgio Zavarise           April 10, 1996            1.0

c      Acronym: GET LABel of a command

c      Purpose: Get a comment string from command line

c      Inputs:
c         ncs     - # of Commas to be Skipped

c      Outputs:
c         ll      - Label Length
c         getlab  - command LABel
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'chdata.h'

      integer   ncs,ll, ip,kc,lp

      save

c      call cdebug0 ('      - getlab',-1)

      ip=3
      do kc=1,ncs
        do while(xxx(ip:ip).ne.',')
          ip=ip+1
          if(ip.eq.80)then
            ll=0
            getlab=' '
            return
          endif
        end do
        ip=ip+1
      end do

      lp=80
      do while(xxx(lp:lp).eq.' ')
        lp=lp-1
      end do

      getlab=xxx(ip:lp)
      ll    = lp - ip + 1

      end
