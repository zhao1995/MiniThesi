c$Id:$
      subroutine umacr11(lct,ctl)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    11/03/2009
c       1. Revise to permit different RVE models on each    13/04/2009
c          processor.
c       2. Set prtype and finflg to end of filename sent to 08/05/2012
c          each rve processor. (sbuf = 132 characters)
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: MPI - Manipulation command to select mesh for micro
c                     RVE problems.  Sets problem type.

c      Use    : Solution command name: 'rve point n_pnt'
c             point - 'point': Activate output to files
c             n_pnt -  RVE number for output

c      N.B. Currently restricted to 64 processors by dimensions in
c           include file 'oumatl.h'

c      Inputs:

c      Outputs:
c         none   - Users are responsible for generating outputs
c                  through common blocks, etc.  See programmer
c                  manual for example.
c-----[--.----+----.----+----.-----------------------------------------]
      implicit   none

      include   'iofile.h'
      include   'oumatl.h'
      include   'chdata.h'
      include   'setups.h'
      include   'umac1.h'

      include   'pointer.h'
      include   'comblk.h'

      include   'mpif.h'

      logical    pcomp
      character  lct*15, sbuf*132
      integer    i,ii,n,nsmax, ierr,usr_msg
      real*8     ctl(3)

      save

      if (pcomp(uct,'ma11',4)) then
          uct = 'rve '
      else

        if(rank.eq.0) then

c         Set point for outputs

          if(pcomp(lct,'poin',4)) then
            n_pnt = nint(ctl(1))
            write(iow,2000) n_pnt
            if(ior.lt.0) then
              write(*,2000) n_pnt
            endif
          else
            n_pnt = 0
          endif

          call uprocset(mr(np(269)), nsmax)
          usr_msg = 14
          ii = 0
          do n = 1,nsmax
            if(unproc(n).gt.0) then
              do i = 1,unproc(n)
                ii = ii + 1

c               Set problem file, type and deformation state

                sbuf(1:128) = matfile(n)
                write(sbuf(129:130),'(i2)') prtyp(1,n)
                write(sbuf(131:132),'(i2)') prtyp(2,n)

c               Send a file name

                if(sbuf(1:1).ne.' ') then
                  call MPI_SSend(sbuf, 132, MPI_CHARACTER, ii, usr_msg,
     &                        MPI_COMM_WORLD, ierr)
                endif
              end do ! i
            endif
          end do ! n

        endif

      end if

c     Format

2000  format(/5x,'RVE number for outputs =',i8)

      end

      subroutine uprocset(rvema, nsmax)

      implicit   none

      include   'debugs.h'
      include   'iofile.h'
      include   'setups.h'
      include   'oumatl.h'
      include   'pointer.h'
      include   'comblk.h'

      logical    err, setval, palloc
      integer    n,nsmax, maxn,maxp, minn,minp, totp, totsend
      integer    rvema(*)

      save

      if(debug) call iprint(rvema,1,nsend,1,'RVEMA')
      do n = 1,64
        umproc(n) = 0
      end do ! n

      err = .false.
      nsmax = 0
      do n = 1,nsend
        if(rvema(n).gt.64) then
          err = .true.
          exit
        elseif(rvema(n).gt.0) then
          umproc(rvema(n)) = umproc(rvema(n)) + 1
          nsmax            = max(nsmax,rvema(n))
        endif
      end do ! n

c     Error on number of processors

      if(err) then
        write(ilg,3000)
        call plstop()
      endif

c     Check for correct number of sends

      totsend = 0
      do n = 1,nsmax
        totsend = totsend + umproc(n)
      end do ! n

      if(totsend.ne.nsend) then
        write(ilg,3001) totsend,nsend
        call plstop()
      endif

c     Allocate number of processors for each material

      minp = 64
      minn = 0
      maxp = 0
      maxn = 0
      totp = 0
      do n = 1,nsmax
        if(umproc(n).gt.0) then
          unproc(n) = nint((ntasks-1)*dble(umproc(n))/dble(totsend))
          unproc(n) = max(1,unproc(n))
          totp      = totp + unproc(n)
          if(unproc(n).gt.maxp) then
            maxp = max(maxp,unproc(n))
            maxn = n
          elseif(unproc(n).lt.minp) then
            minp = max(minp,unproc(n))
            minn = n
          endif
        endif
      end do ! n

c     Check that number of processors is correct

      do while (totp.ne.ntasks - 1)
        if(totp .gt. ntasks - 1) then
          unproc(maxn) = unproc(maxn) - 1
          totp         = totp - 1
        elseif(totp.lt.ntasks - 1) then
          unproc(minn) = unproc(minn) + 1
          totp         = totp + 1
        endif
        minp = 64
        minn = 0
        maxp = 0
        maxn = 0
        do n = 1,nsmax
          if(unproc(n).gt.maxp) then
            maxp = max(maxp,unproc(n))
            maxn = n
          elseif(unproc(n).lt.minp) then
            minp = max(minp,unproc(n))
            minn = n
          endif
        end do ! n
      end do ! while

      if(debug) then
        call iprint(umproc,1,nsmax,1,'UMPROC')
        call iprint(unproc,1,nsmax,1,'UNPROC')
      endif

c     Set up order to send data to each processor

      nrow = 0
      do n = 1,nsmax
        if(unproc(n).gt.0) then
          nrow = max(nrow,(umproc(n) + unproc(n) - 1)/unproc(n))
        endif
      end do ! n
      ncol = ntasks - 1

c     Allocate array to store send list

      setval = palloc( 270,'RVESD',nrow*ncol, 1)

      call usetrvsd(mr(np(269)),mr(np(270)),nsmax)

c     Formats

3000  format(' *ERROR* More user materials than processors.'/
     &       '         Increase number of processors.')

3001  format(' *ERROR* Incorrect number of sends.'/
     &       '         Total sends =',i8,' Should be =',i8)

      end

      subroutine usetrvsd(rvema,ums,nsmax)

      implicit   none

      include   'debugs.h'
      include   'oumatl.h'
      include   'setups.h'

      integer    nsmax, nc,nr, i,j,n
      integer    rvema(*),ums(ncol,nrow)

      save

      j  = 0
      nc = 0
      do n = 1,nsmax
        if(unproc(n).gt.0) then
          nr = 0
          do i = 1,nsend
            if(rvema(i).eq.n) then
              j = mod(i-1,unproc(n)) + 1
              if(j.eq.1) then
                nr = nr + 1
              endif
              ums(nc+j,nr) = i
            endif
          end do ! i
          nc = nc + unproc(n)
        endif
      end do ! n

c     Tag last send to each processor

      do n = 1,ncol
        do nr = nrow,1,-1
          if(ums(n,nr).gt.0) then
            ums(n,nr) = -ums(n,nr)
            exit
          endif
        end do ! nr
      end do ! n

      if(debug) then
        call iprint(ums,ncol,nrow,nc,'UMS')
      endif

      end
