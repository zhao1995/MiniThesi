c$Id:$
      subroutine pstart()

c      * * F E A P * * A Finite Element Analysis Program
c                        -      -       -        -
c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Set 'reflg' to false                             12/10/2007
c       2. Correct call to pinitm, remove comblk.h          01/05/2012
c       3. Replace 'omacr1.h' by 'elpers.h'                 21/05/2013
c-----[--.----+----.----+----.-----------------------------------------]
c     Purpose: Initialize and start all RVE processes

c     Inputs:

c     Outputs:
c-----[--.----+----.----+----.-----------------------------------------]
      implicit      none

      include      'comfil.h'
      include      'elpers.h'
      include      'memuse.h'
      include      'pfeapb.h'
      include      'prmptd.h'
      include      'setups.h'

      include      'mpif.h'           ! OpenMPI common block

      character*128 sfinp,sfout,sfres,sfsav,sfplt, fext*2
      integer       ierr, ii, msg_stat(MPI_STATUS_SIZE)

      save

c     Set flags for serial execution

      pfeap_on   = .false.
      pfeap_gnod = .false.
      ntasks     =  1
      rank       =  0

c     Set maximum memory use

      maxuse = 0

c     Initialize size of communication arrays

      dsend = 0
      drecv = 0

c     Start for X11 graphics driver

      call pdriver()

c     Set for file checking at startup

      fileck = .true.

c     Set restart flag to false (reset on command line inputs only)

      reflg  = .false.

c     Initialize memory

      call pinitm()

c     Initialize MPI

      call MPI_Init( ierr )

c     Get rank id for each process and total number of processes

      call MPI_Comm_rank ( MPI_COMM_WORLD, rank  , ierr)
      call MPI_Comm_size ( MPI_COMM_WORLD, ntasks, ierr)

c     write(*,*) ' Process:',rank,' of ', ntasks,' is alive'

c     Check user installation options

      call pinstall()

c     Set to broadcast interactive commands

      brcast = .true.

c     Set initial file names

      if(rank.eq.0) then

        call filnam()

        do ii = 1,ntasks-1
          fext = '00'
          if(ii.lt.10) then
            write(fext(2:2),'(i1)') ii
          else
            write(fext(1:2),'(i2)') ii
          endif
          sfinp = finp
          call addext(sfinp,fext,128,2)
          sfout = fout
          call addext(sfout,fext,128,2)
          sfres = fres
          call addext(sfres,fext,128,2)
          sfsav = fsav
          call addext(sfsav,fext,128,2)
          sfplt = fplt
          call addext(sfplt,fext,128,2)
          call  MPI_Ssend(sfinp, 128, MPI_CHARACTER,
     &                    ii, 1, MPI_COMM_WORLD, ierr)
          call  MPI_Ssend(sfout, 128, MPI_CHARACTER,
     &                    ii, 2, MPI_COMM_WORLD, ierr)
          call  MPI_Ssend(sfres, 128, MPI_CHARACTER,
     &                    ii, 3, MPI_COMM_WORLD, ierr)
          call  MPI_Ssend(sfsav, 128, MPI_CHARACTER,
     &                    ii, 4, MPI_COMM_WORLD, ierr)
          call  MPI_Ssend(sfplt, 128, MPI_CHARACTER,
     &                    ii, 5, MPI_COMM_WORLD, ierr)
        end do ! ii
      else
        call  MPI_Recv(finp, 128, MPI_CHARACTER,
     &                 0, 1, MPI_COMM_WORLD, msg_stat, ierr)
        call  MPI_Recv(fout, 128, MPI_CHARACTER,
     &                 0, 2, MPI_COMM_WORLD, msg_stat, ierr)
        call  MPI_Recv(fres, 128, MPI_CHARACTER,
     &                 0, 3, MPI_COMM_WORLD, msg_stat, ierr)
        call  MPI_Recv(fsav, 128, MPI_CHARACTER,
     &                 0, 4, MPI_COMM_WORLD, msg_stat, ierr)
        call  MPI_Recv(fplt, 128, MPI_CHARACTER,
     &                 0, 5, MPI_COMM_WORLD, msg_stat, ierr)

      endif

      end
