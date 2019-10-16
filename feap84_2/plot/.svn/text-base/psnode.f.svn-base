c$Id:$
      subroutine psnode(x,ndm,numnp, nz)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Get number for a node using mouse (2-d only)

c      Inputs:
c         x(ndm,*)  - Nodal coordinates in deformed state
c         numnp     - Number of nodal points

c      Outputs:
c         none      - Selected node numbers shown in text screen
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'pdata1.h'

      logical   noerr
      character button*1
      integer   n,ndm,numnp, nz
      real*8    x1,y1, xm,ym, x(ndm,numnp)

      save

c     Pick point from screen

      write(*,*) 'Use  LEFT  Button to Get NODE Number'
      write(*,*) 'Use MIDDLE or RIGHT Button to END'

      button = 'l'
1     call gin(x1,y1,noerr,button)

      if(button .ne.'l') return

      x1 = 0.5d0*(sx(1) + (x1 - s0(1))/scale)
      y1 = 0.5d0*(sx(2) + (y1 - s0(2))/scale)

c     Find closest node

      xm   = abs(x(1,1) - x1)**2 + (x(2,1) - y1)**2
      nz = 1
      do n = 2,numnp
        ym = (x(1,n)-x1)**2+(x(2,n)-y1)**2
        if(ym.lt.xm) then
          xm = ym
          nz = n
        endif
      end do ! n

      write(*,*) 'NODE = ',nz

      go to 1

      end
