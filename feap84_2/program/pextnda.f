c$Id:$
      subroutine pextnda(ix,ip, iplast)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Count number of elements connected to each node

c      Inputs:
c         ix(nen1,*)  - Element connection list

c      Outputs:
c         ip(*)       - List of repeated elements
c         iplast      - Size of repeated list array
c-----[--.----+----.----+----.-----------------------------------------]
      implicit   none

      include   'cdata.h'
      include   'sdata.h'

      integer    ix(nen1,numel), ip(numnp)
      integer    i,ii, n, iplast, nen0
      integer    nume

      save

c     Set the limits for searches

      if(ndm.eq.1) then
        nen0 = min(nen,2)
      elseif(ndm.eq.2) then
        nen0 = min(nen,4)
      elseif(ndm.eq.3) then
        nen0 = min(nen,8)
      endif

c     Count elements connected to array

      do n = 1,numnp
        ip(n) = 0
      end do ! n
      do n = 1,numel
        nume = 0
        do i = 1,nen0
          if(ix(i,n).gt.0) then
            nume = nume + 1
          endif
        end do ! i
        if(nume.eq.nen0) then
          do i = 1,nen0
            ii = ix(i,n)
            if(ii.gt.0) then
              ip(ii) = ip(ii) + 1
            endif
          end do ! i
        endif
      end do ! n

c     Convert 'ip' to pointer form

      do n = 2,numnp
        ip(n) = ip(n) + ip(n-1)
      end do ! n
      iplast = ip(numnp)

      end
