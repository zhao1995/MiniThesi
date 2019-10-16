c$Id:$
      subroutine umacr1(lct,ctl)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2013: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Remove 'prt' from argument list                  09/07/2009
c       2. Added MKL Pardiso to USSL                        12/05/2013
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose:  Set user solver from available library

c      Inputs:
c         lct       - Command character parameters
c         ctl(3)    - Command numerical parameters

c      Outputs:
c         N.B.  Users are responsible for command actions.  See
c               programmers manual for example.
c-----[--.----+----.----+----.-----------------------------------------]
!User Sparse Solver Library
!
!USSL usage:
!
!USSL,MKLP,v1,v2
!  Turns on MKL Pardiso solver
!  v1=1 will print statistics to screen
!  v2=1 sets the symmetric solver to symmetric positive definite (default is indefinite)

      implicit  none

      include  'fdata.h'
      include  'iofile.h'
      include  'part0.h'
      include  'setups.h'
      include  'umac1.h'

      include  'user_solvers.h'

      logical   pcomp,prt, setvar,palloc
      character lct*15,tname*5
      real*8    ctl(3)

      save

c     Install command MKLP (for MKL Pardiso Solver)

      if(pcomp(uct,'mac1',4)) then      ! Usual    form

        write(*,2000)
        uct    = 'ussl'

      elseif(urest.eq.1) then           ! Read  restart data

      elseif(urest.eq.2) then           ! Write restart data

c     Set solver type

      else

        if(pcomp(lct,'off',3)) then ! clear solver

          solver = .true.
          setvar = palloc(225,'USOL1', 0,1)
          setvar = palloc(226,'USOL2', 0,1)
          write(tname,'(4hTANG,i1)') npart
          setvar = palloc(npart,tname, 0,2)


        elseif(pcomp(lct,'mklp',4)) then
          ! MKL Pardiso
          solver = .false.
          usol_pack = 1
          if(pfr) then
            write(iow,2001)
            if(ior.lt.0) then
              write(*,2001)
            endif
          endif
          if(nint(ctl(1)).eq.1) then
            msglvl = 1 ! print statistical information
            write(iow,2002)
            if(ior.lt.0) then
              write(*,2002)
            endif
          else
            msglvl = 0
            write(iow,2003)
            if(ior.lt.0) then
              write(*,2003)
            endif
          endif
          if(nint(ctl(2)).eq.1) then
            mtype = 2 ! symmetric, positive definite
            write(iow,2004)
            if(ior.lt.0) then
              write(*,2004)
            endif
          else
            mtype = -2 ! symmetric, indefinite
            write(iow,2005)
            if(ior.lt.0) then
              write(*,2005)
            endif
          endif


        else ! clear solver
          solver = .true.
          setvar = palloc(225,'USOL1', 0,1)
          setvar = palloc(226,'USOL2', 0,1)
          write(tname,'(4hTANG,i1)') npart
          setvar = palloc(npart,tname, 0,2)
        endif

      endif

2000  format(/' --> User solver library installed'/)
2001  format( '  MKL Pardiso Sparse Solver Activated'/)
2002  format( '  Print Pardiso statistics: On'/)
2003  format( '  Print Pardiso statistics: Off'/)
2004  format( '  TANG factorization: SPD matrix'/)
2005  format( '  TANG factorization: indefinite matrix'/)

      end
