c$Id:$
      subroutine pextndb(ix,ip, ic, iplast)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Compute external nodes on mesh

c      Inputs:
c         ix(nen1,*)  - Element connection list
c         ip(*)       - Pointer array for elements

c      Outputs:
c         ic(*)       - List of repeated elements
c         iplast      - Size of ic list array
c-----[--.----+----.----+----.-----------------------------------------]
      implicit   none

      include   'cdata.h'
      include   'sdata.h'

      integer    ix(nen1,numel), ip(numnp), ic(*)
      integer    i,ii, j,jmin, n, iplast, nen0
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

c     Fill 'ic' array with

      do n = 1,iplast
        ic(n) = 0
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
              if(ii.eq.1) then
                jmin = 1
              else
                jmin = ip(ii-1) + 1
              endif
              do j = jmin,ip(ii)
                if(ic(j).eq.n) then
                  go to 100
                elseif(ic(j).eq.0) then
                  ic(j) = n
                  go to 100
                endif
              end do ! j
100           continue
            endif
          end do ! i
        endif
      end do ! n

      end
