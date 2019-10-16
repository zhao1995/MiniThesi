c$Id:$
      subroutine setfor(f,f0,id,prop,nn, dr)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Add call to pldadd for load group forces         05/03/2009
c       2. Add 'id' array to argument of routine and call   29/05/2009
c          to 'pldadd'.
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Set current value of nodal forces for plots

c      Inputs:
c         f(*)    - Value of force controlled by prop
c         f0(*)   - Base value of force
c         id(*)   - Boundary codes
c         prop    - Proportional load factor
c         nn      - Number of dof

c      Outputs:
c         dr(*)   - Values of nodal forces
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'pointer.h'
      include  'prld1.h'
      include  'comblk.h'

      integer   i,n,nn,id(*)
      real*8    prop
      real*8    f(*),f0(*),dr(*)

      save

c     Check for proportional loading: f = current load,
c         f0(2*nn) = base load, f0 = user supplied loads

      do n = 1,nn
        i = mr(np(29)+n-1)
        if(i.le.0) then
          dr(n) = dr(n) + f(n)*prop     + f0(n) + f0(2*nn+n)
        else
          dr(n) = dr(n) + f(n)*prldv(i) + f0(n) + f0(2*nn+n)
        endif
      end do ! n

c     Add forces from 'load' groups

      call pldadd(mr(np(265)),mr(np(266)),hr(np(267)), dr,id,1,1, prop)

      end
