c$Id:$
      subroutine pdomain(prt)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Change 'include/finclude' to 'finclude'          23/01/2009
c       2. Remove common 'pfeapa' (values in 'setups')      05/02/2009
c       3. Separate blocking from bc-eqs being in           05/01/2013
c       4. Perform compatibiltiy check on matrix type       05/01/2013
c          with solver choice
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Input of domain data for parallel solutions

c      Inputs:
c         prt    - Flag, output results if true

c      Outputs:
c         none   - Users are responsible for generating outputs
c                  through common blocks, etc.  See programmer
c                  manual for example.
c-----[--.----+----.----+----.-----------------------------------------]
      implicit   none

#     include   "cdata.h"
#     include   "chdata.h"
#     include   "debugs.h"
#     include   "iofile.h"
#     include   "pfeapb.h"
#     include   "sdata.h"
#     include   "setups.h"

#     include   "pointer.h"
#     include   "comblk.h"

#     include   "finclude/petsc.h"

      PetscErrorCode ierr
      PetscBool      chk

      logical    prt,pcomp,errck, pinput, tinput, palloc, matno, noeq
      character  cmnd*15, lct*15, optstr*20
      integer    ii, jj, ndomains,nvals, leqn,geqn
      real*8     td(10)
      real*4     etime, tary(2)

c     Default is unblocked matrix storage with bc-eqs out

      pfeap_blk  = .false.
      pfeap_bcin = .false.
      pfeap_view = .false.
      if(ntasks.gt.1) then
        pfeap_on = .true.
        solver   = .false.
      else
        pfeap_on = .false.
        solver   = .true.
      endif
      pfeap_dstr = .false.
      noeq = .true.

c     Input of domain data

      cmnd = 'start'
      do while(.not.pcomp(cmnd,'end',3))
        errck = tinput(cmnd,1,td,0)

c       Local to Global numbers

        if(pcomp(cmnd,'loca',4)) then

          if(prt .and. debug) then
            write(*,*) 'rank[',rank,']',yyy(1:45),' ',etime(tary)
          endif
          if( np(244).eq.0 ) then
            if(pfeap_bcin) then
              errck = palloc(244,'GN',numnp*2,1)
            else
              errck = palloc(244,'GN',numnp  ,1)
            endif
          endif

          if(prt) then
            if(pfeap_bcin) then
              write(iow,2000) 'Equation Data'
            else
              write(iow,2000)
            endif
          endif
          do ii = 1,numnp
            errck = pinput(td,3)
            mr(np(244)+nint(td(1))-1) = nint(td(2))
            if(pfeap_bcin) then
              mr(np(244)+nint(td(1))-1+numnp) = nint(td(3))
              if(prt) then
                write(iow,2001) nint(td(1)),nint(td(2)),nint(td(3))
              endif
            else
              if(prt) then
                write(iow,2001) nint(td(1)),nint(td(2))
              endif
            endif
          end do ! ii

c       Equation numbers

        elseif(pcomp(cmnd,'equa',4)) then

          noeq = .false.
          if(prt .and. debug) then
            write(*,*) 'rank[',rank,']',yyy(1:45),' ',etime(tary)
          endif
          if( np(245).eq.0 ) then
            errck = palloc(245,'EQ',ndf*numnp,1)
          endif

          if(prt) write(iow,2002)
          do ii = 1,numnp
            errck = pinput(td,ndf+1)
            leqn  = ndf*(nint(td(1))-1)
            if(prt) write(iow,2003) nint(td(1))
            do jj = 0, ndf-1
              mr(np(245)+ leqn + jj) = nint(td(jj+2))
              if(prt) write(iow,2004) nint(td(jj+2))
            end do ! jj
            if(prt) write(iow,2005)
          end do ! ii

c       Matrix storage descriptors

        elseif(pcomp(cmnd,'matr',4)) then

          if(prt .and. debug) then
            write(*,*) 'rank[',rank,']',yyy(1:45),' ',etime(tary)
          endif
          errck = palloc(246,'DNNZ',ndf*numpn,1)
          errck = palloc(247,'ONNZ',ndf*numpn,1)

          if(prt) write(iow,2006)
          matno = .true.
          numpeq = 0
          do while (matno)
            errck = pinput(td,3)
            if (nint(td(1)).eq.0) then
              matno = .false.
            else
              if(prt) write(iow,2001) (nint(td(ii)),ii=1,3)
              mr( np(246)+ numpeq ) = nint(td(2))
              mr( np(247)+ numpeq ) = nint(td(3))
              numpeq = numpeq + 1
            endif
          end do ! while
          if(pfeap_blk) then
            numpeq = numpeq*nsbk
          endif
          write(iow,*) '          Number of partition equations:',numpeq

c       Get data for domains

        elseif(pcomp(cmnd,'getd',4)) then

          if(prt .and. debug) then
            write(*,*) 'rank[',rank,']',yyy(1:45),' ',etime(tary)
          endif
          lct = yyy(16:30)
          call setval (yyy(31:45),15, td(1))

          if(pcomp(lct,'poin',4)) then
            ndomains = nint(td(1))
            errck = palloc(248,'GETP',ndomains,1)

            if(prt) write(iow,2007) ndomains

            do ii = 1, ndomains
              errck = pinput(td(1),1)
              mr(np(248)+ii-1) = nint(td(1))
              if(prt) write(iow,2001) ii,nint(td(1))
            end do

          elseif(pcomp(lct,'valu',4)) then
            nvals = nint(td(1))
            errck = palloc(249,'GETV',nvals,1)

            if(prt) write(iow,2008)

            do ii = 1, nvals
              errck = pinput(td(1),1)
              mr(np(249)+ii-1) = nint(td(1))
            end do ! ii

            leqn = 0
            do ii = 1,ntasks
              do jj = 1,mr(np(248)+ii-1)
                if(prt) write(iow,2001) ii,mr(np(249)+leqn)
                leqn = leqn + 1
              end do ! jj
            end do ! ii
          endif

c       Send data for domains

        elseif(pcomp(cmnd,'send',4)) then

          if(prt .and. debug) then
            write(*,*) 'rank[',rank,']',yyy(1:45),' ',etime(tary)
          endif
          lct = yyy(16:30)
          call setval (yyy(31:45),15, td(1))

          if(pcomp(lct,'poin',4)) then
            ndomains = nint(td(1))
            errck = palloc(250,'SENP',ndomains,1)

            if(prt) write(iow,2009) ndomains

            do ii = 1, ndomains
              errck = pinput(td(1),1)
              mr(np(250)+ii-1) = nint(td(1))
              if(prt) write(iow,2001) ii,nint(td(1))
            end do

          elseif(pcomp(lct,'valu',4)) then
            nvals = nint(td(1))
            errck = palloc(251,'SENV',nvals,1)

            if(prt) write(iow,2010)

            do ii = 1, nvals
              errck = pinput(td(1),1)
              mr(np(251)+ii-1) = nint(td(1))
            end do ! ii

            leqn = 0
            do ii = 1,ntasks
             do jj = 1,mr(np(250)+ii-1)
              if(prt) write(iow,2001) ii,mr(np(251)+leqn)
              leqn = leqn + 1
             end do ! jj
            end do ! ii
          endif

c       matrix format for storage

        elseif(pcomp(cmnd,'form',4)) then

          lct = yyy(16:30)
          if(pcomp(lct,'baij',4)) then
            pfeap_blk = .true.
            if(prt) write(iow,2013)
            if(debug) then
              write(*,2013)
            endif
          else
            if(prt) write(iow,2014)
            if(debug) then
              write(*,2014)
            endif
          endif

        elseif(pcomp(cmnd,'bcin',4)) then

          pfeap_bcin = .true.
          if(prt) write(iow,2012)
          if(debug) then
            write(*,2012)
          endif

          call setval (yyy(31:45),15, td(1))
          nsbk = nint(td(1))
          if(prt) write(iow,2011) nsbk
          if(debug) then
            write(*,2011) nsbk
          endif

        endif

      end do ! while

c     Check to assign equations

      if(noeq) then
        if( np(245).eq.0 ) then
          errck = palloc(245,'EQ',ndf*numnp,1)
        endif

        if(prt) write(iow,2002)
        do ii = 1,numnp
          leqn  = ndf*(ii - 1) - 1
          geqn  = ndf*(mr(np(244)+ii-1+numnp) - 1)
          if(prt) write(iow,2003) ii
          do jj = 1, ndf
            mr(np(245)+ leqn + jj) = geqn + jj
            if(prt) write(iow,2004) geqn+jj
          end do ! jj
          if(prt) write(iow,2005)
        end do ! ii
      endif

      if(prt .and. debug) then
        write(*,*) 'rank[',rank,']',yyy(1:45),' ',etime(tary)
      endif

c     Check on compatibility of solver selection with matrix type

c     Formats
      call PetscOptionsGetString(PETSC_NULL_CHARACTER,'-pc_type',
     &     optstr,chk,ierr)
      if(chk .eqv. PETSC_TRUE) then

        if((optstr.eq.'hypre' .or. optstr.eq.'ml') .and. pfeap_blk) then
           write(iow,2015) optstr
           call plstop
        else if(optstr.eq.'gamg') then
           if(pfeap_blk) then
             write(iow,2015) optstr
             call plstop
           else if(.not.pfeap_bcin) then
             write(iow,2019) optstr
             call plstop
           endif
        else if(optstr.eq.'prometheus') then
           if(.not.pfeap_blk) then
             write(iow,2017) optstr
             call plstop
           else if(ndm.ne.3) then
             write(iow,2018) optstr
             call plstop
           endif
        endif

      endif

      call PetscOptionsGetString(PETSC_NULL_CHARACTER,
     &     '-pc_factor_mat_solver_package',optstr,chk,ierr)
      if(chk .eqv. PETSC_TRUE) then

        if((optstr.eq.'spooles' .or. optstr.eq.'mumps' .or.
     &      optstr.eq.'superlu_dist') .and. pfeap_blk)  then
           write(iow,2016) optstr
           call plstop
        endif
      endif

2000  format(//5x,'Local to global mapping',
     &        /5x,'Local node',10x,'Global node':,8x,a)

2001  format(i15,2i21)

2002  format(//5x,'Local node to global eq mapping',
     &        /5x,'Local node',10x,'Global eqs')

2003  format(5x,i10,5x,$)

2004  format(i10,$)

2005  format(/,$)

2006  format(//5x,'MATRIX non-zero storage data',
     &        /5x,'Global Eq',5x,
     &       'Diag non-zero',5x,'Off Diag non-zero')

2007  format(//5x,'Get data pointer list ',i5,' domains'
     &        /5x,'Domain',10x,'Num Data Items')

2008  format(//5x,'Get data values from',
     &        /5x,'Domain',10x,'Node')

2009  format(//5x,'Send data pointer list ',i5,' domains'
     &        /5x,'Domain',10x,'Num Data Items')

2010  format(//5x,'Send data values to',
     &        /5x,'Domain',10x,'Node')

2011  format(/5x,'Matrix block size  =',i4)

2012  format(/5x,'BC Eqns in Assembly')

2013  format(//5x,'Matrix format BAIJ',
     &        /5x,'Storage information given in blocked format')

2014  format(//5x,'Matrix format AIJ',
     &        /5x,'Storage information given in unblocked format')

2015  format(//5x,'**ERROR** PC TYPE ',a,' requires AIJ matrix format')

2016  format(//5x,'**ERROR* SOLVER PACKAGE ',a
     &           ,' requires AIJ matrix format')

2017  format(//5x,'**ERROR** PC TYPE ',a,' requires BAIJ matrix format')

2018  format(//5x,'**ERROR** PC TYPE ',a,' requires NDM = 3')

2019  format(//5x,'**ERROR** PC TYPE ',a,' requires BCIN')

      end
