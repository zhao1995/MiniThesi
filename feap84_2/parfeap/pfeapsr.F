c$Id:$
      subroutine pfeapsr(array,tdatabuf,pmax,srflag)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1.  Add error reporting for send/receive            16/03/2007
c       2. Change 'include/finclude' to 'finclude'          23/01/2009
c       3. Remove common 'pfeapa' (values in 'setups')      05/02/2009
c       4. Add 'srflag' to control return mode              09/03/2013
c       5. Fix min/max for srflag false                     23/05/2013
c       6. Fix srflag false to compute min                  06/01/2014
c-----[--.----+----.----+----.-----------------------------------------]
c     Purpose:  Send and receive real data using MPI: PETSc version

c     Inputs:
c        array(pmax) - Real array from this processor
c        pmax        - Number items in array
c        srflag      - Flag to control accumulation

c     Outputs:
c        array(pmax) -  Real array accumulated from all processors
c-----[--.----+----.----+----.-----------------------------------------]
      implicit   none

#     include   "finclude/petsc.h"
#     include   "iofile.h"
#     include   "setups.h"

      logical    srflag
      integer    i, j, pmax
      real*8     array(pmax),tdatabuf(pmax,2)

      integer    usr_msg, req
      parameter (usr_msg=12)

c     Petsc values

      PetscErrorCode ierr

      integer        msg_stat(MPI_STATUS_SIZE)

      save

c     Sending Data Asynchronously

      do j = 1, pmax
        tdatabuf(j,1) = array(j)
      end do ! j

      do i = 0,ntasks-1
        if(rank.ne.i) then
          call MPI_Isend( tdatabuf(1,1), pmax, MPI_DOUBLE_PRECISION, i,
     &                    usr_msg, MPI_COMM_WORLD, req, ierr )
c         Report send error

          if(ierr.ne.0) then
            write(iow,*) ' -->> Send Error[',rank+1,'->',i,']',ierr
            write(  *,*) ' -->> Send Error[',rank+1,'->',i,']',ierr
          endif

        endif
      end do ! i

c     Receiving Data in blocking mode

      do i = 0,ntasks-1
        if(rank.ne.i) then
          call MPI_Recv( tdatabuf(1,2), pmax, MPI_DOUBLE_PRECISION, i,
     &                   usr_msg, MPI_COMM_WORLD, msg_stat, ierr )
          if(ierr.ne.0) then
             write(iow,*) 'Recv Error[',i,'->',rank+1,']',ierr
             write(  *,*) 'Recv Error[',i,'->',rank+1,']',ierr
          endif
          if(srflag) then
            do j = 1, pmax
              array(j) = array(j) + tdatabuf(j,2)
            end do ! j
          else
            if(pmax.eq.1) then
              array(1) = max(array(1),tdatabuf(1,2))
            elseif(pmax.eq.2) then
              array(1) = max(array(1),tdatabuf(1,2))
              array(2) = min(array(2),tdatabuf(2,2))
            else
             write(iow,*) 'Error in pfeapsr pmax = ',pmax
             write(  *,*) 'Error in pfeapsr pmax = ',pmax
            endif
          endif
        endif
      end do ! i

      call MPI_BARRIER(MPI_COMM_WORLD, ierr)

      end
