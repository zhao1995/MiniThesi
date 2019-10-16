c$Id:$
      subroutine mshcksn(is,iblend,numsd,numsn,numbd,errs)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Perform check on supernode values of blending function

c      Inputs:
c         is(16,*)       - Entries on blend sides
c         iblend(20,*)   - Blending function description
c         numsd          - Maximum number of sides
c         numsn          - Maximum number of supernodes
c         numbd          - Number of blending functions

c      Outputs:
c         errs           - Flag, true if errors detected
c-----[--.----+----.----+----.-----------------------------------------]

      implicit none

      include 'iofile.h'

      logical  errs
      integer  numsd,numsn,numbd,n,i,inc,is(16,numsd),iblend(20,numbd)

      save

c     Loop over SIDE nodes to check if any greater than number SNODes

      do n = 1,numsd
        if(is(1,n).eq.2) then
          inc = 2
        else
          inc = 1
        endif
        do i = 2,16,inc
          if(is(i,n).gt.numsn) then
            write(ilg,2000) n, numsn, i
            write(iow,2000) n, numsn, i
            errs = .true.
          endif
        end do ! i
      end do ! n

c     Loop over BLENd nodes to check if any greater than number SNODes

      do n = 1,numbd
        do i = 11,18
          if(iblend(i,n).gt.numsn) then
            write(ilg,2001) n, numsn, i
            write(iow,2001) n, numsn, i
            errs = .true.
          endif
        end do ! i
      end do ! n

c     Formats

2000  format(' *ERROR* SIDE',i5,' has value greater than maximum SNODE'
     &      ,' (',i5,') at entry',i5/)

2001  format(' *ERROR* BLENd',i5,' has value greater than maximum SNODE'
     &      ,' (',i5,') at entry',i5/)

      end
