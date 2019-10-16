c$Id:$
      subroutine udnorm(urate)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Compute norms of transient derivatives.

c      Inputs:
c        urate(ndf,*,*) - Rate array

c      Outputs:
c-----[--.----+----.----+----.-----------------------------------------]
      implicit   none

      include   'cdata.h'
      include   'iofile.h'
      include   'part0.h'
      include   'print.h'
      include   'sdata.h'

      integer    i,n
      real*8     urate(ndf,numnp,*), vnorm(2)

      save

c     Compute and output velocity and acceleration norms

      vnorm(1) = 0.0d0
      vnorm(2) = 0.0d0
      do i = 1,ndf
        if(ndfp(i).eq.npart .and. ndfo(i).ge.1) then
          do n = 1,numnp
            vnorm(1) = vnorm(1) + urate(i,n,1)*urate(i,n,1)
          end do ! n
        endif
        if(ndfp(i).eq.npart .and. ndfo(i).ge.2) then
          do n = 1,numnp
            vnorm(2) = vnorm(2) + urate(i,n,2)*urate(i,n,2)
          end do ! n
        endif
      end do ! i
      vnorm(1) = sqrt(vnorm(1))
      if(vnorm(2).gt.0.0d0) then
        vnorm(2) = sqrt(vnorm(2))
        if(prnt.and.ior.lt.0) then
          write(*,2000) vnorm
        endif
        write(iow,2000) vnorm
      else
        if(prnt.and.ior.lt.0) then
          write(*,2000) vnorm(1)
        endif
        write(iow,2000) vnorm(1)
      endif

c     Formats

2000  format('   N o r m s   f o r   D y n a m i c s'/
     &   10x,'Velocity:',e13.5:,' Acceleration:',e13.5)

      end
