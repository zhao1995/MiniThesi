c$Id:$
      subroutine pelbld(ix,id,ixg,nen,nen1,ndf,numel,elist)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Build list of elements connected to nodes

c      Inputs:

c      Outputs:
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'eqslv.h'
      include  'pointer.h'
      include  'comblk.h'

      logical   flag
      integer   i,ii,j,n, elist
      integer   nen,nen1,ndf,numel
      integer   ix(nen1,*),id(ndf,*),ixg(*)

      save

c     Determine number of restrained slave nodes

      neqg  = 0
      do n = 0,numg - 1
        do i = 1,ndf
          if(id(i,mr(np(205)+n)).le.0) then
            neqg = neqg + 1
          endif
        end do ! i
      end do ! n

c     Build list of elements needed for constructing the stiffness
c     coefficients and profile information

      elist = 0
      do n = 1,numel
        if(ix(nen1-1,n).ge.0) then
          flag = .true.
          do i = 1,nen
            ii = ix(i,n)
            if(ii.gt.0) then
              do j = 1,numg

c               Save element number if first occurrance

                if (flag) then
                  elist = elist + 1
                  ixg(elist) = n
                  flag = .false.
                endif

              end do ! j
            endif
          end do ! i
        endif
      end do ! n

      end
