c$Id:$
      subroutine uidset(id,ndf,numnp)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    18/12/2012
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Reset boundary codes to activate equations

c      Inputs:
c         id(ndf,*)  - Boundary conditions
c         ndf        - DOF/node
c         numnp      - Number nodes

c      Outputs:
c         id(ndf,*)  - Modified bBoundary conditions
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      integer   ndf,numnp
      integer   id(ndf,numnp)

c     Users to provide any modifications necessary

      end
