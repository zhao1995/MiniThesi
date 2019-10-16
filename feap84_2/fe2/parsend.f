c$Id:$
      subroutine parsend()

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    12/21/2007
c       1. Revised form for multiple RVE mesh types         12/05/2010
c       2. Add sbuf(8) increase send size                   13/12/2010
c       3. Increase storage for rbuf to 24                  20/07/2012
c       4. Replace 'omacr1.h' by 'elpers.h'                 21/05/2013
c-----[--.----+----.----+----.-----------------------------------------]
c     Purpose: MPI Time update Exchanges

c     Input:
c       None

c     Output:
c       None
c-----[--.----+----.----+----.-----------------------------------------]
      implicit   none

      include   'debugs.h'
      include   'counts.h'
      include   'elpers.h'
      include   'hdatam.h'
      include   'iofile.h'
      include   'setups.h'
      include   'tdata.h'

      include   'mpif.h'

      integer    mm,nsbuf,nrbuf
      integer    usr_msg, msg_stat(MPI_STATUS_SIZE), ierr
      integer    nproce

      real*8     sbuf(24), rbuf(91)

c     Compute time update

      save

c     Set values

c     Exit for single task problems

      if(ntasks.le.1) then

        return

c     Send: Time update from Rank 0 model to each task

      elseif(rank.eq.0) then

c       Set size of buffers

        nsbuf = dsend + 7
        nrbuf = drecv

        sbuf(1) = 0           ! Update number
        sbuf(2) = nstep       ! From /counts/
        sbuf(3) = niter       ! From /counts/
        sbuf(4) = dt          ! From /tdata/
        sbuf(5) = 12          ! ISW element switch value
        sbuf(6) = 1.d0        ! (1 = Fortran true)
        sbuf(7) = 0.0d0
        sbuf(8) = 0.0d0

c       Send deformation gradient to each processor

        usr_msg = 12
        do mm = 1,ntasks-1

c         Assign processor number

          nproce  = mm

c         Send message

          if(debug) then
            call udebug('     MPI_SSend:Tgrad',usr_msg)
          endif
          call MPI_SSend( sbuf, nsbuf, MPI_DOUBLE_PRECISION, nproce,
     &                   usr_msg,  MPI_COMM_WORLD, ierr)
          if(ierr.ne.0) then
            write(*,*) ' IERR_send =',ierr
          endif
        end do ! mm

c       Receive message from processors

        usr_msg = 13
        do mm = 1,ntasks-1

c         Assign processor number

          nproce  = mm

c         Receive time update reply

          if(debug) then
            call udebug('     MPI_Recv:Flux',usr_msg)
          endif
          call MPI_Recv( rbuf, nrbuf, MPI_DOUBLE_PRECISION, nproce,
     &                   usr_msg,  MPI_COMM_WORLD, msg_stat, ierr)

        end do ! mm

      endif

      end
