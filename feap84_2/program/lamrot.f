c$Id:$
      subroutine lamrot(ua, dl)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Compute rotation matrix associated to rotation
c               vector ua(3).

c      Inputs:
c         ua(3)   - Rotation vector

c      Outputs:
c         dl(3,3) - Rotation matrix
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      real*8    uanrm, fac1, fac2, fac3, ua(3), dl(3,3)

      save

c     Compute norm of ua

      uanrm = sqrt(ua(1)*ua(1) + ua(2)*ua(2) + ua(3)*ua(3))
      fac1  = cos(uanrm)

c     Check For asymptotic aolution if ||ua|| near 0.0:

      if (uanrm.lt.1.d-04) then
        uanrm = uanrm**2
        fac2  = 1.0d0 - uanrm*(1.d0/6.d0
     &                - uanrm*(1.d0/120.d0
     &                - uanrm/5040.d0))

        fac3  = 0.5d0 - uanrm*(1.d0/24.d0
     &                - uanrm*(1.d0/720.d0
     &                - uanrm/40320.d0))
      else
        fac2  = sin(uanrm)/uanrm
        uanrm = uanrm * 0.5d0
        fac3  = 0.5d0 * (sin(uanrm)/uanrm)**2
      endif

c     Assemble rotation matrix

      dl(1,1)  =        fac1 + fac3*ua(1)*ua(1)
      dl(1,2)  = -ua(3)*fac2 + fac3*ua(1)*ua(2)
      dl(1,3)  =  ua(2)*fac2 + fac3*ua(1)*ua(3)

      dl(2,1)  =  ua(3)*fac2 + fac3*ua(2)*ua(1)
      dl(2,2)  =        fac1 + fac3*ua(2)*ua(2)
      dl(2,3)  = -ua(1)*fac2 + fac3*ua(2)*ua(3)

      dl(3,1)  = -ua(2)*fac2 + fac3*ua(3)*ua(1)
      dl(3,2)  =  ua(1)*fac2 + fac3*ua(3)*ua(2)
      dl(3,3)  =        fac1 + fac3*ua(3)*ua(3)

      end
