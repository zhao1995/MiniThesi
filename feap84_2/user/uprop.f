c$Id:$
      subroutine uprop(unumb,type,tmin,tmax,ap, t,upropld, isw)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: User proportional load routine

c      Inputs:
c        unumb   - User function number
c        type    - User type definition
c        tmin    - Minimum time applicable
c        tmax    - Maximum time applicable
c        ap(5)   - User function parameters
c        t       - Current time for evaluation point

c      Output:
c        upropld - Value of proportional load
c-----[--.----+----.----+----.-----------------------------------------]
      implicit   none

      include   'iofile.h'

      integer    unumb,type, isw
      real*8     tmin, tmax, ap(5), t, upropld

c     Print output of user data

      if(isw.eq.1) then

c     Compute function value for 'propld'

      elseif (isw.eq.2) then

        upropld = 0.0d0 ! Default return value

c       Compute load for Type 1

        if(unumb.eq.1) then

        endif

      endif

      end
