c$Id:$
      subroutine mkside(n,iface,is,isd)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Form sides for edges of 3-d Block

c      Inputs:
c        n           - Face number of block
c        iface(4)    - Side numbers of face
c        isd         - Dimension of is array

c      Outputs:
c        is(isd,*)   - Side list
c-----[--.----+----.----+----.-----------------------------------------]
      implicit none

      include 'iofile.h'

      integer  i,j,k, i1,j1, ie,je, n,isd,is(isd,*),iface(4)

      save

c     Check first side for correct direction

      i1 = abs(iface(1))
      j1 = abs(iface(2))
      if(is(1,i1).eq.2) then
        do k = 3,isd
          if(is(k,i1).ne.0) ie = k
        end do ! k
      else
        ie = 3
      endif
      if(is(1,j1).eq.2) then
        do k = 3,isd
          if(is(k,j1).ne.0) je = k
        end do ! k
      else
        je = 3
      endif
      if(is(ie,i1).eq.is(2,j1)    .or. is(ie,i1).eq.is(je,j1)) then
        iface(1) =  i1
      elseif(is(2,i1).eq.is(2,j1) .or. is(2,i1).eq.is(je,j1)) then
        iface(1) = -i1
      else
        write(iow,2000) n,1,2
      endif

c     Check remaining directions

      do i = 1,3
        j  = i + 1
        i1 = abs(iface(i))
        j1 = abs(iface(j))
        if(is(1,i1).eq.2) then
          do k = 3,isd
            if(is(k,i1).ne.0) ie = k
          end do ! k
        else
          ie = 3
        endif
        if(is(1,j1).eq.2) then
          do k = 3,isd
            if(is(k,j1).ne.0) je = k
          end do ! k
        else
          je = 3
        endif
        if(iface(i).gt.0) then
          if(is(ie,i1).eq.is(2,j1)) then
            iface(j) =  j1
          elseif(is(ie,i1).eq.is(je,j1)) then
            iface(j) = -j1
          else
            write(iow,2000) n,i,j
          endif
        else
          if(is(2,i1).eq.is(2,j1)) then
            iface(j) =  j1
          elseif(is(2,i1).eq.is(je,j1)) then
            iface(j) = -j1
          else
            write(ilg,2000) n,i,j
            write(iow,2000) n,i,j
          endif
        endif
      end do ! i

c     Formats

2000  format(' *ERROR* MKSIDE: Face ',i2,' No match between sides',i2,
     &       ' and',i2)

      end
