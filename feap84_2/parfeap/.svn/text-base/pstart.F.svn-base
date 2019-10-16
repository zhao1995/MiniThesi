c$Id:$
      subroutine pstart()

c      * * F E A P * * A Finite Element Analysis Program
c                        -      -       -        -
c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Change 'include/finclude' to 'finclude'          23/01/2009
c       2. Remove common 'pfeapa' (values in 'setups')      05/02/2009
c       3. Remove arguments from call to pinitm             05/11/2012
c       4. Completely removed pfeapa.h                      05/01/2013
c       5. Update for loss of VecValid, MatValid            05/01/2013
c       6. Remove processor                                 07/01/2013
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Set startup values

c      Inputs:
c         None

c      Outputs:
c-----[--.----+----.----+----.-----------------------------------------]

      implicit   none

#     include   "finclude/petsc.h"
#     include   "pfeapb.h"
#     include   "pfeapc.h"

#     include   "comfil.h"
#     include   "memuse.h"
#     include   "prmptd.h"
#     include   "setups.h"
#     include   "comblk.h"

      character*128  sfinp,sfout,sfres,sfsav,sfplt
      integer        ii, msg_stat(MPI_STATUS_SIZE)
      PetscErrorCode ierr

      save

c     Set maximum memory use

      maxuse = 0

c     Start for X11 graphics driver

      call pdriver()

c     Set for file checking at startup

      fileck = .true.

c     Initialize memory

      call pinitm()

c     Initialize the PETSc environment and get rank and size of job

      call PetscInitialize(PETSC_NULL_CHARACTER,   ierr)
      call MPI_Comm_rank(PETSC_COMM_WORLD, rank,   ierr)
      call MPI_Comm_size(PETSC_COMM_WORLD, ntasks, ierr)

      pfeap_on   = .false.
      pfeap_view = .false.

      rhs   = 0
      sol   = 0
      xvec  = 0
      yvec  = 0
      zvec  = 0
      Mdiag = 0
      Msqrt = 0
      Kmat  = 0
      Mmat  = 0
      Pmat  = 0

c     Check user installation options

      call pinstall()

c     Set initial file names

      if (rank.eq.0) then
        call filnam()

        do ii = 2,ntasks
          call adomnam(finp,sfinp,ii)
          call adomnam(fout,sfout,ii)
          call adomnam(fres,sfres,ii)
          call adomnam(fsav,sfsav,ii)
          call adomnam(fplt,sfplt,ii)

          call MPI_SEND(sfinp, 128, MPI_CHARACTER,
     &                 ii-1, 1, PETSC_COMM_WORLD, ierr)
          call MPI_SEND(sfout, 128, MPI_CHARACTER,
     &                 ii-1, 2, PETSC_COMM_WORLD, ierr)
          call MPI_SEND(sfres, 128, MPI_CHARACTER,
     &                 ii-1, 3, PETSC_COMM_WORLD, ierr)
          call MPI_SEND(sfsav, 128, MPI_CHARACTER,
     &                 ii-1, 4, PETSC_COMM_WORLD, ierr)
          call MPI_SEND(sfplt, 128, MPI_CHARACTER,
     &                 ii-1, 5, PETSC_COMM_WORLD, ierr)
        end do
      else
         call MPI_RECV(finp, 128, MPI_CHARACTER, 0, 1,
     &                PETSC_COMM_WORLD, msg_stat,ierr)
         call MPI_RECV(fout, 128, MPI_CHARACTER, 0, 2,
     &                PETSC_COMM_WORLD, msg_stat,ierr)
         call MPI_RECV(fres, 128, MPI_CHARACTER, 0, 3,
     &                PETSC_COMM_WORLD, msg_stat,ierr)
         call MPI_RECV(fsav, 128, MPI_CHARACTER, 0, 4,
     &                PETSC_COMM_WORLD, msg_stat,ierr)
         call MPI_RECV(fplt, 128, MPI_CHARACTER, 0, 5,
     &                PETSC_COMM_WORLD, msg_stat,ierr)
      endif

      end
