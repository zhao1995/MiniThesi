c$Id:$
      subroutine basmat(phib,mb,cmass,lmass,t, neq)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c     Purpose: Compute coupling term

c     Inputs:

c     Outputs:

c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

c#     include   "sdata.h"
      include   'fdata.h'
      include   'pfeapb.h'
      include   'pointer.h'
      include   'comblk.h'

      integer    neq, i
      real*8     t(*),cmass(*),lmass(*),phib(*),mb(*)

      save

c     Parallel solution

      if(pfeap_on) then

        call parbmat(phib,mb,cmass,lmass,t, neq)

c     Serial solution

      else

c       Consistent mass

        if(fl(1)) then
          call caprod(cmass(1),cmass(neq+1),phib,t,
     &                mr(np(90)),mr(np(91)),neq)

        do i = 1,neq
          t(i) = t(i) - mb(i)
        end do ! i

c       Lumped mass

        else
          do i = 1,neq
            t(i) = lmass(i)*phib(i) - mb(i)
          end do ! i
        endif
      endif

      end
