c$Id:$
      subroutine parstop()

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Increase rbuf to 20 words                        20/01/2011
c       2. Set nrbuf to dsend + 7                           08/05/2012
c       3. Increase rbuf to 24 words                        20/)7/2012
c       4. Replace 'omacr1.h' by 'elpers.h'                 21/05/2013
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Close any open parallel array and delete memory use
c               Dummy routine in serial version.

c      Inputs:
c         none

c      Outputs:
c         none
c-----[--.----+----.----+----.-----------------------------------------]
      implicit   none

      include   'debugs.h'
      include   'elpers.h'
      include   'setups.h'
      include   'mpif.h'

      integer    i, ierr, nrbuf, usr_msg
      real*8     rbuf(24)

      save

      data       rbuf / -999.0d0, 23*0.0d0 /
      data       usr_msg / 12 /

      if(debug) then
        call udebug(' parstop',0)
      endif

c     Close parallel arrays

      if(rank.eq.0) then
        nrbuf = dsend + 7
        do i = 1,ntasks - 1
          if(debug) then
            call udebug('     MPI_Send:Stop',usr_msg)
          endif
          call MPI_SSend(rbuf, nrbuf, MPI_DOUBLE_PRECISION, i, usr_msg,
     &                  MPI_COMM_WORLD, ierr)
        end do ! i
      endif

      call MPI_Finalize( ierr )

      end
