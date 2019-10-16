c$Id:$
      real function etime(tarry)

c      * * F E A P * * A Finite Element Analysis Program
c                        -      -       -        -
c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--+---------+---------+---------+---------+---------+---------+-]
c      Purpose: System clock value -- Fortran 90 compatible routine

c      Inputs: none

c      Outputs: tarry(1) = cpu time
c               tarry(2) = system time (Not returned for this routine)
c-----[--+---------+---------+---------+---------+---------+---------+-]
      implicit  none

      include  'etime1.h'

      logical   vinput, errck
      character ctime*10
      real      tarry(2)
      real*8    hr,min,sec

c     Call standard system routine to get time

      call date_and_time(time= ctime)

c     Extract times from character array

      errck = vinput(ctime(1:2) ,2,hr  ,1)
      errck = vinput(ctime(3:4) ,2,min ,1)
      errck = vinput(ctime(5:10),5,sec ,1)

c     Convert to elapsed time in seconds

      tarry(1) = 60.d0*(min + 60.d0*hr) + sec - tim0
      tarry(2) = 0.0   ! No System timing is available

c     Adjust if clock passed mid-night

      if(tarry(1).lt.0.0) then
        tarry(1) = tarry(1) + 86400.0
        tim0     = tim0     - 86400.0
      endif

c     Return elapsed total time (CPU + System)

      etime    = tarry(1) + tarry(2)

      end
