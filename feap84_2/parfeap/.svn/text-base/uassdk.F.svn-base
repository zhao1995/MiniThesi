c$Id:$
      subroutine uassdk(s,eq)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    03/22/2013
c-----[--.----+----.----+----.-----------------------------------------]
c     Purpose:  Parallel assembly for nodal stiffness

c     Inputs:
c       s(numpn*ndf)  - Diagonal stiffness
c       eq(numnp*ndf) - global equations numbers

c     Outputs:
c       Kmat          - Global stiffness
c-----[--.----+----.----+----.-----------------------------------------]
      implicit   none

#     include   "finclude/petsc.h"
#     include   "pfeapc.h"          ! Contains Kmat

#     include   "compas.h"
#     include   "rdat1.h"
#     include   "pfeapb.h"

#     include   "sdata.h"
#     include   "iofile.h"

      PetscInt   ierr

      integer    i, e
      integer    eq(*)
      real*8     s(*)

c     Assemble matrix

      do i = 1,numpn*ndf
        e = eq(i)
        if (e.gt.0 .and. s(i).ne.0.0d0) then
          call MatSetValue( Kmat, e-1, e-1, s(i), ADD_VALUES, ierr )
        endif
      end do ! i loop

c     Set flag indicating new tangent formed
      pfeap_tang = .true.

      end
