c$Id:$
      subroutine tcritnd(d,xl,ul,ndm,ndf,nel)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    03/02/2007
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Element critical time step estimate

c      Inputs:
c        d(184)    - Maximum wave speed
c        xl(ndm,*) - Element undeformed coordinates
c        ul(ndf,*) - Element displacements
c        ndm       - Mesh dimension
c        ndf       - Dof's/node
c      Outputs:
c        dtcr      - Through /tdata/
c-----[--.----+----.----+----.-----------------------------------------]
      implicit   none

      include   'tdata.h'

      integer    ndm,ndf,nel, a,b, i
      real*8     xlmin,xlen
      real*8     d(*), xl(ndm,*),ul(ndf,*)

      xlmin = 0.0d0
      do i = 1,ndm
        xlmin = xlmin + (xl(i,2) - xl(i,1) + ul(i,2) - ul(i,1))**2
      end do ! i
      do a = 1,nel-1
        do b = a+1,nel
          xlen = 0.0d0
          do i = 1,ndm
            xlen = xlen + (xl(i,b) - xl(i,a) + ul(i,b) - ul(i,a))**2
          end do ! i
          if(xlen.gt.0.0d0) then
            xlmin = min(xlmin,xlen)
          endif
        end do ! b
      end do ! a
      if(d(184).gt.0.0d0) then
        if(dtcr.eq.0.0d0) then
          dtcr = sqrt(xlmin)/d(184)
        else
          dtcr = min(dtcr,sqrt(xlmin)/d(184))
        endif
      endif

      end
