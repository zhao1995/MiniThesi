c$Id:$
      subroutine parform(isw)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    24/03/2013
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Set any last steps for parallel solutions

c      Inputs:
c         isw       - Switch: solution type for solves

c      Outputs:
c-----[--.----+----.----+----.-----------------------------------------]
      implicit    none

#     include    "finclude/petsc.h"
#     include    "fdata.h"
#     include    "pfeapb.h"
#     include    "pfeapc.h"

      PetscInt    ierr
      PetscViewer M_view

      integer     isw

c     Check for parallel solution

      if(pfeap_on) then

c       Assemble mass for parallel form

        if(isw.eq.5) then

          if(fl(1)) then             ! Consistent mass
            call MatAssemblyBegin (Mmat, MAT_FINAL_ASSEMBLY, ierr)
            call MatAssemblyEnd   (Mmat, MAT_FINAL_ASSEMBLY, ierr)
          else                       ! Lumped mass
            call VecAssemblyBegin (Mdiag, ierr)
            call VecAssemblyEnd   (Mdiag, ierr)
          endif

          if(pfeap_view) then
            call PetscViewerASCIIOpen(PETSC_COMM_WORLD,"mass.m",M_view,
     &                                ierr)
            call PetscViewerSetFormat(M_view,PETSC_VIEWER_ASCII_MATLAB,
     &                                ierr)
            if(fl(1)) then
              call MatView           (Mmat,M_view,ierr)
            else
              call VecView           (Mdiag,M_view,ierr)
            endif
            call PetscViewerDestroy  (M_view , ierr)
          endif
        endif ! isw

      endif ! pfeap_on

      end
