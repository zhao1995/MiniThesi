c$Id:$
      subroutine ptay_bc(numnp,ndf,id)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    15/05/2013
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Restrain all dof's in problem: Taylor boundary condition

c      Inputs:
c         numnp   - Number of nodes
c         ndf     - Number dofs/node

c      Outputs:
c         id(:,:) - Boundary condition array
c-----[--.----+----.----+----.-----------------------------------------]
      implicit   none

      integer    numnp,ndf
      integer    id(ndf,numnp)

c     Restrain all dofs

      id = 1

      end
