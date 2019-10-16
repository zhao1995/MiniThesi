c$Id:$
      subroutine orthf_apl(a)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    12/12/2012
c      Purpose: Compute orthogonal vectors for plasticity functions

c      Inputs:
c        a(3,2)  - Input plastic vectors

c      Outputs:
c        a(3,2)  - Orthogonal plastic vectors
c-----[--.----+----.----+----.-----------------------------------------]
      implicit   none

      include   'iofile.h'

      real*8     a(3,2), xx1, tol
      integer    i

c     Tolerance set

      data       tol / 1.e-8 /

c     Compute first unit rotation vector

      xx1    = 1.0d0/sqrt(a(1,1)**2 + a(2,1)**2 + a(3,1)**2)
      if(xx1.lt.tol) then
        write(iow,3000) 1,(a(i,2),i=1,3)
        write(ilg,3000) 1,(a(i,2),i=1,3)
        call plstop()
      endif
      a(:,1) = a(:,1)*xx1

c     Compute second orthogonal rotation vector

      xx1    = sqrt(a(1,2)**2 + a(2,2)**2 + a(3,2)**2)
      if(xx1.lt.tol) then
        write(iow,3000) 2,(a(i,2),i=1,3)
        write(ilg,3000) 2,(a(i,2),i=1,3)
        call plstop()
      endif

c     Compute second orthogonal rotation vector

      xx1  = a(1,2)*a(1,1)+ a(2,2)*a(2,1) + a(3,2)*a(3,1)
      a(:,2) = a(:,2) - xx1*a(:,1)
      xx1    = sqrt(a(1,2)**2 + a(2,2)**2 + a(3,2)**2)
      if(xx1.lt.tol) then
        write(iow,3001) (a(i,2),i=1,3),(a(i,1),i=1,3)
        write(ilg,3001) (a(i,2),i=1,3),(a(i,1),i=1,3)
        call plstop()
      endif

c     Normalize to unit length

      xx1    = 1.0d0/xx1
      a(:,2) = a(:,2)*xx1

3000  format(' ** ERROR ** Vector',i2,'=',1p,3e12.4,' has zero length')
3001  format(' ** ERROR ** Vector 1 =',1p,3e12.4,' parallel to'/
     &       '             Vector 2 =',1p,3e12.4)

      end
