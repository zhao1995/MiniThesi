c$Id:$
      function prop2p(lunit,l,tv,itime)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Proportional load table type 2 input/computation

c                 prop = vm(i)+(vm(i+1)-vm(i))*(t-tm(i))/(tm(i+1)-tm(i)
c                        tm(i) < t < tm(i+1)
c      Inputs:
c         l         - Number of data input pairs to input/record
c                     Compute proportional load if zero.

c      Outputs:
c         prop2     - Value of total proportional load type 2
c         tv(2,*)   - Table of times and values:
c                       tm(*) = tv(1,*); vm(*) = tv(2,*)
c         itime     - Activation indicator
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'comfil.h'
      include  'iodata.h'
      include  'iofile.h'
      include  'print.h'

      logical   errck, pinput, done
      integer   l, m, ilast, itime, lunit, iosav
      real*8    prop2p, tv(2,*), td(16)

      save

c     Input table of proportional loads

      ilast = 0
      itime = 1

c     Rewind file to permit setting the tv(2,*) array

      rewind lunit
      iosav = ior
      ior   = lunit

c     Start read

      done = .false.
      do while(.not.done)
102     errck = pinput(td,2*l)
        if(errck) go to 102
        do m = 1,l
          if(abs(td(2*m-1))+abs(td(2*m)).ne.0.0d0
     &                         .or. ilast.eq.0) then
            ilast       = ilast + 1
            tv(1,ilast) = td(2*m-1)
            tv(2,ilast) = td(2*m)
          else
            done        = .true.
          endif
        enddo ! m
      end do ! while

      close(unit = lunit, status = 'delete')
      ior = iosav

      prop2p = 0.0d0

      end
