c$Id:$
      subroutine cupdat(id,u,du,nneq,ndf)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Static complex imaginary UPDATE

c      Inputs:
c         id (ndf,numnp)  - ID-array (vector storage)
c         u (3*nneq)      - Displacement vectors
c         du(nneq)        - Displacement increment from SOLVER
c         nneq            - numnp * ndf
c         ndf             - Number of DOF/node

c      Outputs:
c         u (3*nneq)      - Displacement vectors
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'part0.h'
      include  'tdatb.h'

      integer   i,j, n, nneq,ndf,nneq2, id(*)
      real*8    u(*),du(*)

      save

c     Update imaginary displacement and increments within step.

      nneq2 = nneq + nneq

      do i = 1,ndf
        if(ndfp(i).eq.npart) then
          do n = i,nneq,ndf
            j = id(n)

c           For active degrees-of-freedom compute values from solution

            if (j.gt.0) then
              u(n)       = u(n)      + cc1*du(j)
              u(n+nneq)  = u(n+nneq) + cc2*du(j)
              u(n+nneq2) =                 du(j)
            else
              u(n)       = 0.0d0
              u(n+nneq)  = 0.0d0
              u(n+nneq2) = 0.0d0
            endif
          end do ! n
        endif
      end do ! i

      end
