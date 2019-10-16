c$Id:$
      subroutine uprob(titl,prt)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    13/04/2009
c       1. Add 'titl' to argument list                      08/06/2009
c       2. Add input of elinks if necessary                 30/11/2009
c       3. Strip blanks from input line for 'elink'         02/12/2010
c       4. Change 'eln' to 'elin' for edge links            17/03/2011
c       5. Add call to setext for 'elin'                    27/02/2012
c       6. Extract prtype and finflg from end of filename   08/05/2012
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: MPI - Receive command to start mesh on micro problems
c               Allows for different RVE's on each processor.

c      Inputs:
c         prt    - Flag, output results if true

c      Outputs:
c         none   - Users are responsible for generating problems
c-----[--.----+----.----+----.-----------------------------------------]
      implicit   none

      include   'cblktr.h'
      include   'cdata.h'
      include   'cdat2.h'
      include   'comfil.h'
      include   'debugs.h'
      include   'elpers.h'
      include   'iofile.h'
      include   'mxsiz.h'
      include   'sdata.h'
      include   'setups.h'
      include   'pointer.h'
      include   'comblk.h'

      include   'mpif.h'

      logical    prt, setvar, palloc, pcomp, vinput
      character  titl*80,sbuf*132,filenam*128, text*256,xyz*256,fext*8
      character  value*15
      integer    i, ierr, ndum, usr_msg
      integer    msg_stat(MPI_STATUS_SIZE)
      real*8     vv(2)

      save

      if(debug) then
        call udebug(' uprob',0)
      endif

c     Set command

      if(rank.gt.0) then

c       Receive file name from macro program

        usr_msg = 14

        if(debug) then
          call udebug('     MPI_Recv:NEW',usr_msg)
        endif
        call MPI_Recv(sbuf, 132, MPI_CHARACTER, 0, usr_msg,
     &                  MPI_COMM_WORLD, msg_stat, ierr)
        filenam = sbuf(1:128)
        value   = ' '
        value(14:15) = sbuf(129:130)
        setvar = vinput(value,15,vv(1),1)
        value   = ' '
        value(14:15) = sbuf(131:132)
        setvar = vinput(value,15,vv(2),1)

c       Set findex to use for file offsets

        findex = 1
        do i = len_trim(filenam),1,-1
          if(pcomp(filenam(i:i),char(47),1) .or.       ! char(47) = '/'
     &       pcomp(filenam(i:i),char(92),1)) go to 100 ! char(92) = '\'
        end do ! i
        i = 0
100     findex = i + 1

c       Set inputs from file specified in second field

        call pincld(filenam)

c       Start problem

        nio = 0
        neo = 0
        mao = 0
        call pnewprob(0)

c       Check for edge links

1       read(ior,1000,end=2) xyz
        call pstrip(text,xyz,1)

        if(pcomp(text,'elin',4)) then
          ndum = 0
          call setext('elin', ndum, fext, .false.)
          call plinka(fext,'set','   ')
          leflg = .true.
          if(np(257).eq.0) then
            setvar = palloc(257,'ELINK',numnp*ndm,  1)
            do i = 1,ndm*numnp
              mr(np(257)+i-1) = 0
            end do ! i
            do i = 1,3
              dxlnk(i) = 0.0d0
            end do ! i
          endif
        else
          go to 1
        endif

2       call pincld('end')

c       Export problem type and deformation mode

        prtype = nint(vv(1))
        finflg = vv(2).le.0.0d0

      endif

c     Format

1000  format(a)

      end
