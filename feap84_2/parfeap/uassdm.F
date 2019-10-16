c$Id:$
      subroutine uassdm(s,eq,consist,lumped)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    03/22/2013
c-----[--.----+----.----+----.-----------------------------------------]
c     Purpose:  Parallel assembly for nodal stiffness

c     Inputs:
c       s(numpn*ndf)  - Diagonal mass
c       eq(numnp*ndf) - global equations numbers
c       consist       - Consistent assembly if true
c       lumped        - Lumped assembly if true

c     Outputs:
c       Mmat          - Global mass
c-----[--.----+----.----+----.-----------------------------------------]
      implicit   none

#     include   "finclude/petsc.h"
#     include   "pfeapb.h"
#     include   "pfeapc.h"          ! Contains Mmat
#     include   "sdata.h"

      PetscInt   ierr

      logical    consist, lumped
      integer    i, e
      integer    eq(*)
      real*8     s(*)

c     Assemble matrix

      if(consist) then
        do i = 1,numpn*ndf
          e = eq(i)
          if (e.gt.0 .and. s(i).ne.0.0d0) then
            call MatSetValue( Mmat, e-1, e-1, s(i), ADD_VALUES, ierr )
          endif
        end do ! i loop
      endif

c     Assemble lumped mass vector

      if(lumped) then
        do i = 1,numpn*ndf
          e = eq(i)
          if (e.gt.0 .and. s(i).ne.0.0d0) then
            call VecSetValue( Mdiag, e-1, s(i), ADD_VALUES, ierr )
          endif
        end do ! i loop
      endif

      end
