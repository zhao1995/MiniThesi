c$Id:$
      subroutine pfeapmi(intsr)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Change to integer send/receive and add barrier   17/11/2006
c       2. Fixed incorrect cast on tdatabuf                 30/04/2007
c       3. Change 'include/finclude' to 'finclude'          23/01/2009
c       4. Remove common 'pfeapa' (values in 'setups')      05/02/2009
c-----[--.----+----.----+----.-----------------------------------------]
c     Purpose:  Send and receive data using MPI: PETSc version

c     Inputs:
c        intsr - Integer value from current processor

c     Outputs:
c        intsr - Maximum value from all processors
c-----[--.----+----.----+----.-----------------------------------------]
      implicit   none

#     include   "finclude/petsc.h"
#     include   "setups.h"

      integer    intsr, intsro, i
      integer    tdatabuf

      integer    usr_msg, req
      parameter (usr_msg=11)

c     Petsc values

      PetscErrorCode ierr

      integer        msg_stat(MPI_STATUS_SIZE)

      save

c     Sending Data Asynchronously

      intsro = intsr

      do i = 0,ntasks-1
        if(rank.ne.i) then
          call MPI_Isend( intsro, 1, MPI_INTEGER, i,
     &                    usr_msg, MPI_COMM_WORLD, req, ierr )
        endif
      end do ! i

c     Receiving Data in blocking mode

      do i = 0,ntasks-1
        if(rank.ne.i) then
          call MPI_Recv( tdatabuf, 1, MPI_INTEGER, i,
     &                     usr_msg, MPI_COMM_WORLD, msg_stat, ierr )
          intsr = max(intsr,tdatabuf)
        endif
      end do ! i

      call MPI_BARRIER(MPI_COMM_WORLD, ierr)

      end
