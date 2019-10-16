c$Id:$
      subroutine elhispl(htab, hdat, ma,nhv,l)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Revise to have 'aflg' option for inertial effect 28/12/2007
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Project history plot data

c      Inputs:
c        htab(:,:) - Table of plot components to plot
c        hdat(:)   - Values of history variables to plot
c        ma        - Material set for data
c        nhv       - Maxium number of available history variables
c        l         - Pt plotted

c      Outputs:
c        pl(:,:)   - Projection value
c-----[--.----+----.----+----.-----------------------------------------]
      implicit   none

      include   'cdata.h'
      include   'eldatp.h'

      integer    ma,nhv,l,htab(10,*)
      real*8     hdat(*)

      integer    i,ih

      do i = 1,10
        plhis(i,l) = 0.0d0
      end do ! i

      do i = 1,10
        ih = htab(i,ma)
        if(ih.gt.0 .and. ih.le.nhv) then
          plhis(i,l) = hdat(ih)
          hpltfl     = .true.
          plhmax     = max(plhmax,i)
        endif
      end do ! i

      end
