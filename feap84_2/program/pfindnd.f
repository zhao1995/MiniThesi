c$Id:$
      subroutine pfindnd(ix,ib,ir)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c     Purpose:

c     Inputs:

c     Outputs:
c-----[--.----+----.----+----.-----------------------------------------]
      implicit   none

      include   'cdata.h'
      include   'sdata.h'

      integer    i,ii, m,mm,mp,m0, n,ne,nm,np,n0, nel
      integer    ib(*),ir(3,*), ix(nen1,*)

      save

      ne = 0
      do n = 1,numel
        n0  = 0
        nel = 0
        do i = 1,nen
          ii = ix(i,n)
          if(ii.gt.0) then
            nel = nel + 1
            if(ib(ii).gt.0) then
              n0 = n0 + 1
              np = i
            endif
          endif
        end do ! i
        if(n0.eq.1 .and. nel.gt.2) then
          ne       = ne + 1
          ir(1,ne) = n
          ir(2,ne) = np
          ir(3,ne) = nel
        endif
      end do ! n

c     'ne' now is the number of possible dangling surface nodes


      if(ne.gt.0) then
        do n = 1,ne-1
          n0 = ir(2,n)
          nm = n0 - 1
          if(nm.eq.0) nm = ir(3,n)
          np = mod(n0,ir(3,n)) + 1
          nm = ix(nm,ir(1,n))
          np = ix(np,ir(1,n))
          do m = n+1,ne
            m0 = ir(2,m)
            if(ix(m0,ir(1,m)).ne.ix(n0,ir(1,n))) then
              mm = m0 - 1
              if(mm.eq.0) mm = ir(3,m)
              mp = mod(m0,ir(3,m)) + 1
              mm = ix(mm,ir(1,m))
              mp = ix(mp,ir(1,m))
              if(mp.eq.nm) then
                ib(mp) = 1
              elseif(mm.eq.np) then
                ib(mm) = 1
              endif
            endif
          end do ! m
        end do ! n
      endif

      end
