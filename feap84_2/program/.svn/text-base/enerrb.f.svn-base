c$Id:$
      subroutine enerrb(mass,pi1,w1,rcg,momt,rangm,rener)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Compute energy for rigid body

c      Inputs:
c         mass   - Mass of rigid body
c         pi1    - Angular momenta vector
c         w1     - Angular velocity vector
c         rcg    - Array of translational velocities, etc.

c      Outputs:
c         momt   - Linear  momentum
c         rangm  - Angular momentum
c         rener  - Kinetic energy
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      integer   i
      real*8    mass,rener, rcg(3,11),pi1(3),w1(3),momt(3),rangm(3)

      save

c     Compute linear momentum

      do i=1,3
        momt(i) = mass*rcg(i,5)
      end do ! i

c     Total energy

      rener = 0.5d0*(w1(1)*pi1(1) + w1(2)*pi1(2) + w1(3)*pi1(3)
     &      + rcg(1,5)*momt(1) + rcg(2,5)*momt(2) + rcg(3,5)*momt(3))

c     Resultant angular momentum: pi + r x momt

      rangm(1) = pi1(1) + rcg(2,2)*momt(3)  - rcg(3,2)*momt(2)
      rangm(2) = pi1(2) + rcg(3,2)*momt(1)  - rcg(1,2)*momt(3)
      rangm(3) = pi1(3) + rcg(1,2)*momt(2)  - rcg(2,2)*momt(1)

      end
