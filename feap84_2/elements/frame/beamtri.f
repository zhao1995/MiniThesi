c$Id:$
      subroutine beamtri(xl,ndm,d,xs,pt)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c     Purpose: Compute three dimensional frame transformation matrix

c     Inputs:

c     Outputs:
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'iofile.h'

      integer i,ndm, lref

      real*8  xl(ndm,*),d(*),pt(3,4),t2,dx1,dx2,theta,xs

c     Compute center (0,0) displacement

      do i = 1,3
        pt(i,1) = 0.5d0*(xl(i,1) + xl(i,2))
      end do !i

c     Compute normal 1 to element

      lref = nint(d(96))

      if(lref .eq. 1) then

        pt(1,3) = d(97) - xl(1,1)
        pt(2,3) = d(98) - xl(2,1)
        pt(3,3) = d(99) - xl(3,1)

      elseif(lref .eq. 2) then

        pt(1,3) = d(97)
        pt(2,3) = d(98)
        pt(3,3) = d(99)

      elseif(lref .eq. 3) then
        dx1 = 0.5d0*(xl(1,1) + xl(1,2))
        dx2 = 0.5d0*(xl(2,1) + xl(2,2))
        theta  = atan2(dx2,dx1)
        pt(1,3) = cos(theta)
        pt(2,3) = sin(theta)
        pt(3,3) = 0.0d0
      else
        write(ilg,*) ' NO REFERENCE VECTOR DEFINED'
        call plstop()
      endif

c     Compute tangent to element

      t2 = 0.d0
      do i = 1,3
        t2 = t2 + (xl(i,2) - xl(i,1))**2
      end do !i
      xs = 0.5d0*sqrt(t2)
      t2 = 0.5d0/xs
      do i = 1,3
        pt(i,4) = (xl(i,2) - xl(i,1))*t2
      end do !i

c     Compute normal 2 to element

      call vecp (pt(1,3),pt(1,4),pt(1,2))

      end
