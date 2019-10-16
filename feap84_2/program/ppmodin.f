c$Id:$
      subroutine ppmodin(phi, y, u,v,tu,nv,tneq,mtyp)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Initial Conditions for modal solution of
c               Parallel interface: Dummy routine

c      Inputs:
c         phi(vneq,*) - Mass orthonormal eigenvectors
c         u(*)        - Initial displacements
c         v(*)        - Initial velocities
c         nv          - Number of eigenpairs (converged)
c         tneq        - Number of d.o.f. in model
c         mtyp        - Mass type (true for consistent,false for lumped)

c      Scratch:
c         tu(*,2)     - Temporary vector

c      Outputs:
c         y(nv,3)     - Eigensolution at time t_0
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'pfeapb.h'

      logical   mtyp
      integer   nv,tneq
      real*8    phi(vneq,*), y(nv,3), u(*),v(*), tu(tneq,*)

      save

      end
