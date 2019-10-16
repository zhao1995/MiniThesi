c$Id:$
      subroutine pltqud(iquad,ndm)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Set quadrant numbers to recieve plot data

c      Inputs:
c         ndm       - Spatial dimension of mesh

c      Outputs:
c         iquad(*)  - Quadrant indicators
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'iofile.h'
      include  'prmptd.h'

      logical   errck, pinput, qflg
      integer   ix,iy,iz,ndm, iquad(2,2,2)
      real*8    td(3)

      save

      call pzeroi(iquad,8)

      qflg = .true.

1     if(defalt) then
        call pzero(td,3)
      else
        if(ior.lt.0) then
          call pprint(' Input: x-quad, y-quad, z-quad <CR to exit>: ')
        endif
        errck = pinput(td,3)
        if(errck) go to 1
      endif

      if(abs(td(1)) + abs(td(2)) + abs(td(3)) .gt. 0.0d0) then
        if(qflg) then
          do ix = 1,8
            iquad(ix,1,1) = 1
          end do ! ix
          qflg = .false.
        end if

        if(td(1).ge.0.0d0) then
          ix = 1
        else
          ix = 2
        end if

        if(td(2).ge.0.0d0 .or. ndm.lt.2) then
          iy = 1
        else
          iy = 2
        end if

        if(td(3).ge.0.0d0 .or. ndm.lt.3) then
          iz = 1
        else
          iz = 2
        end if

        iquad(ix,iy,iz) = 0

      else
        return
      end if

      go to 1

      end
