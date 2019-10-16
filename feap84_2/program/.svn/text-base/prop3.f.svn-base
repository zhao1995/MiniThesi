c$Id:$
      function prop3(i,t,ap,x,f,fpro)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c     1. Use 'pi' from 'pconstant.h'                        14/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Proportional load type 3 for sliding/rolling

c      Inputs:
c         i           - Proportional load for lookup
c         t           - Time for lookup
c         ap(5)       - 1-vel,2-vel,r-vel,x-0,y-0
c         fpro(ndf,*) - Nodal dof proportional load numbers

c      Outputs:
c         f(ndf,*,2)  - Displacement values at nodes
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'cdata.h'
      include  'iofile.h'
      include  'pconstant.h'
      include  'sdata.h'

      integer   i, n, fpro(ndf,numnp)
      real*8    prop3, t,u0,v0,r0, dx,dy,radius,theta
      real*8    ap(*),x(ndm,numnp),f(ndf,numnp,2)

      save

      u0 = ap(1)*t
      v0 = ap(2)*t
      r0 = ap(3)*t*pi/180.0d0

      do n = 1,numnp
        if(fpro(1,n).eq.i .and. fpro(2,n).eq.i) then
          dx    = x(1,n) - ap(4)
          dy    = x(2,n) - ap(5)
          theta = atan2(dy,dx)
          radius= sqrt(dx*dx + dy*dy)
          f(1,n,2)    = u0 + radius*(cos(theta+r0) - cos(theta))
          f(2,n,2)    = v0 + radius*(sin(theta+r0) - sin(theta))
        endif
      end do ! n

c     Return unit value for the proportional load

      prop3 = 1.d0

      end
