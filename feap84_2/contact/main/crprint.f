c$Id:$
      subroutine crprint(dr,initfl,nfl)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose:  Computes range and profile of plot value requested

c      Inputs:
c         dr(*)  - Values for plot (N.B. checks dr(i) values)

c      Outputs:
c         nfl    - Returns zero if all values are zero
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'c_0.h'
      include  'c_keyh.h'
      include  'fdata.h'
      include  'iofile.h'
      include  'pdata2.h'
      include  'prange.h'
      include  'prmptd.h'
      include  'rpdata.h'

      logical   initfl
      integer   nfl, i,n
      real*8    dr(*), pr(9),drv,rs

      save

c     Compute profile of values

      if(rangfl) then
        rmn = rangmn
        rmx = rangmx
        if (initfl) then
          call pzero (pr,9)
          drv = 100.d0/nset
        endif
      else
        if (initfl) then
          call pzero (pr,9)
          drv = 100.d0/nset
          rmx = dr(1)
          rmn = dr(1)
        endif
        do n = 1,nset
          rmx = max(rmx,dr(n))
          rmn = min(rmn,dr(n))
        end do
      endif

      if(abs(rmx-rmn).gt.1.d-5*max(abs(rmx),abs(rmn))) then
        do n = 1,nset
          rs = (dr(n) - rmn)/(rmx - rmn)
          do i = 1,9
            if(rs.ge.0.1d0*i) pr(i) = pr(i) + drv
          end do
        end do
        if(prompt) then
          if(pfr) write(iow,2000) rmn,rmx,pr
          if(ior.lt.0 .and. idev.ne.2 .and. .not.defalt) then
            write(*,2000) rmn,rmx,pr
          endif
        endif
      else
        if(prompt) then
          if(pfr) write(iow,2000) rmn,rmx
          if(ior.lt.0 .and. .not.defalt) then
            write(*,2000) rmn,rmx
          endif
        endif
        if(abs(rmx).gt.1.d0) then
          rmn = rmn - 0.1*abs(rmx)
          rmx = rmx + 0.1*abs(rmx)
        else
          rmn = rmn - 1.0d0
          rmx = rmx + 1.0d0
        endif
        nfl = 0
      endif

2000  format('    Minimum is ',1p,e10.2,' Maximum is ',1p,e10.2:/
     &  22x,'10%   20%   30%   40%   50%   60%   70%   80%   90%'/
     &       '    Profile above is:',9f6.1)

      end
