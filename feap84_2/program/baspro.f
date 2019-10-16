c$Id:$
      subroutine baspro(t,phi,pmass,n,mf,neq,vneq,mp)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c     Purpose: Compute coupling term

c     Inputs:

c     Outputs:

c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      integer   m,n,mf,neq,vneq,mp
      real*8    dot,t(*),phi(vneq,*),pmass(mf,mp)

      save

      do m= 1,mf
        pmass(m,n) = dot(phi(1,m), t,neq)
      end do ! m

      end
