c$Id:$
      subroutine usetci()

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c     Purpose: Set parameters to perform updates for given time step

c     Inputs :
c       dt       - Current time increment  ! TDATA.H
c       theta(3) - Command line parameters ! DDATA.H

c     Outputs:
c       gtan(3)  - Parmeters for: 1 - stiffness; 2 - damping; 3 - mass
c                  e.g., S = K*gtan(1) + C*gtan(2) + M*gtan(3)
c                                                          ! GLTRAN.H
c       cc1      - Parameter to update solution  at t_n+1  ! TDATB.H
c       cc2      - Parameter to update increment at t_n+1
c       cc3      - Parameter to update specified boundary displacement
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      save

      end
