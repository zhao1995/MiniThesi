c$Id:$
      subroutine umati8(type,vv, d, ud, n1,n3)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Dummy user material model routine

c      Inputs:
c         type   - Name of material model
c         vv(5)  - Command line real data
c         d(*)   - Program material parameter data

c      Outputs:
c         ud(*)  - Material parameter data for model
c         n1     - Number of history items/point (time   dependent)
c         n3     - Number of history items/point (time independent)
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      logical   pcomp
      character type*15
      integer   n1,n3
      real*8    vv(5),d(*),ud(*)

c     Set command name

      if(pcomp(type,'mat8',4)) then     ! Default  form DO NOT CHANGE
c       type = 'name'                   ! Specify 'name'

c     Input user data and save in ud(*) array

      else                              ! Perform input for user data

      endif

      end
