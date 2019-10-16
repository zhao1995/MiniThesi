c$Id:$
      subroutine rprint(ix,nef,nen1,nume, dr,ndf,nfl)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Remove check on dr(1,nn).eq.0 in non0fl          14/05/2010
c       2. Add send/receive or ranges for parallel          09/03/2013
c       3. Reverse order on rscale, remove min/max on       23/05/2013
c          receive value.
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose:  Computes range and profile of plot value requested

c      Inputs:
c         ix(nen1,*)- Element connections
c         nef       - Maximum number of nodes/element
c         nen1      - Dimension of 'ix' array
c         nume      - Number of elements/faces
c         dr(ndf,*) - Values for plot (N.B. checks dr(1,i) values)
c         ndf       - Dimension of dr-array

c      Outputs:
c         nfl       - Returns -nfl if all values are < 1.0d-08
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'setups.h'

      include  'cdata.h'
      include  'fdata.h'
      include  'iofile.h'
      include  'pbody.h'
      include  'pdata2.h'
      include  'pointer.h'
      include  'prmptd.h'
      include  'rpdata.h'
      include  'comblk.h'

      include  'p_int.h'

      logical   non0fl
      integer   nef,nen1,nume, ndf,nfl, i,n,nn, ma
      integer   ix(nen1,*)
      real*8    drv,rs, dr(ndf,*),pr(9), rscale(2), tdatabuf(2,2)

      save

c     Compute minimum and maximum profile of values

      non0fl = .true.
      fp(1) = npty - 1
      do n = 1,nume
        ma = ix(nen1,n)
        if((ix(nen1-1,n).ge.nreg1 .and. ix(nen1-1,n).le.nreg2) .and.
     &     ((ma.gt.0 .and. maplt.eq.0) .or. ma.eq.maplt) ) then
          do i = 1,nef
            nn = ix(i,n)
            if(nn.gt.0) then
              if(mr(fp(1)+nn) .ge. 0) then
                if(non0fl) then
                  rmx    = dr(1,nn)
                  rmn    = dr(1,nn)
                  non0fl = .false.
                else
                  rmx = max(rmx,dr(1,nn))
                  rmn = min(rmn,dr(1,nn))
                endif
              endif
            endif
          end do ! i
        endif
      end do ! n

c     Do communication for parallel solutions

      rscale(1) = rmx
      rscale(2) = rmn
      call pfeapsr(rscale,tdatabuf, 2, .false.)
      rmx = rscale(1)
      rmn = rscale(2)

c     Check for zero condition

      if(non0fl) then
        rmx = 0.0d0
        rmn = 0.0d0
      endif

c     Check range for contour outputs

      do n = 1,9
        pr(n) = 0.0d0
      end do ! n
      drv   = 100.d0/numnp
      if(abs(rmx-rmn).gt.1.d-5*max(abs(rmx),abs(rmn))) then
        do n = 1,numnp
          if(mr(fp(1)+n) .ge. 0) then
            rs = (dr(1,n) - rmn)/(rmx - rmn)
            do i = 1,9
              if(rs.ge.0.1d0*i) pr(i) = pr(i) + drv
            end do ! i
          endif
        end do ! n
      else
        if(max(abs(rmx),abs(rmn)).gt.1.d-8) then
          rmn = rmn - 0.0999*max(abs(rmx),abs(rmn))
          rmx = rmx + 0.1001*max(abs(rmx),abs(rmn))
        else
          rmn = -1.0d-08 + 1.0d-12
          rmx =  1.0d-08 + 1.0d-12
          nfl = -abs(nfl)
        endif
      endif

c     Output range for min/max

      if(prompt) then
        if(pfr) write(iow,2000) rmn,rmx
        if(ior.lt.0 .and. .not.defalt) then
          write(*,2000) rmn,rmx
        endif
      endif

2000  format('    Minimum is ',1p,e10.2,' Maximum is ',1p,e10.2:/
     &  22x,'10%   20%   30%   40%   50%   60%   70%   80%   90%'/
     &       '    Profile above is:',9f6.1)

      end
