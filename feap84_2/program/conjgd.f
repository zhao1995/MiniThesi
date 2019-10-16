c$Id:$
      subroutine conjgd(ad,ac,adr,aur,jpr,x,jc,ir,r,z,p,b,
     &                  neq,jmax,tol,rn,rn0,ittyp)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Remove interupt tests (intf and lintr)           17/11/2006
c       2. Remove unused format 1000                        13/11/2006
c       3. Define last 'tt' in terms of 'tary(1)'           08/05/2008
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Preconditioned conjugate gradient with preconditioning

c         See:

c         R.M. Ferencz, "EBE Preconditioning Techniques for Large
c                        Scale, Vectorized Finite Element Analysis
c                        in Nonlinear Solid and Structural Mechanics,"
c                        Box 3.2.2 Computational implementation of the
c                        modified preconditioned conjugate gradient
c                        algorithm, Ph.D. Dissertation, Stanford,
c                        March 1989.
c      Inputs:
c         ad(*)    - Diagonal entries of coefficient matrix
c         ac(*)    - Compressed part of symmetric coefficient matrix
c         jpr(*)   - Pointer array for columns/rows of block precond.
c         x(*)     - Residual
c         jc(*)    - Pointer array to locate column/rows in ir array
c         ir(*)    - Location of non-zero terms in symmetric array.
c         neq      - Number of equations
c         jmax     - Maximum number of iterations
c         tol      - Solution residual tolerance
c         rn0      - Initial residual
c         ittyp    - Preconditioner type

c      Outputs:
c         x(*)     - Solution
c         rn       - Final residual

c      Scratch:
c         adr(*)   - Storage for reciprocal diagonals
c         aur(*)   - Storage for symmetric block preconditioners
c         r(*)     - Scratch arrays
c         z(*)     - Scratch arrays
c         p(*)     - Scratch arrays
c         b(*)     - Scratch arrays
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'debugs.h'
      include  'iofile.h'
      include  'print.h'
      include  'sdata.h'

      logical   noconv, captim
      integer   neq,nc,ni,nj,jmax,ittyp, jpr(*),jc(*),ir(*)
      integer   megops, numcap,numnoc, numadd,nopcap,nopcgd,totops
      real*8    mflops
      real*4    tt0, tt, etime, tary(2)
      real*8    alp,bet,gam,rn,rn0,tol,tolr0,dot
      real*8    ad(neq),ac(*),adr(neq),aur(*),x(*)
      real*8    r(*),z(*),p(*),b(*)

      save

c     Initialize

      do ni = 1,neq
        r(ni) = x(ni)
        b(ni) = x(ni)
        x(ni) = 0.0d0
      end do ! ni

      if(rn0.eq.0.0d0) then
        rn0 = dot(r,r,neq)
      endif
      tolr0 = rn0*tol**2

c     Set-up initial solution vectors 'z' and 'p'

c     1. Diagonal preconditioning

      if(ittyp.eq.1) then

        do ni = 1,neq
          if(ad(ni).ne.0.0d0) then
            adr(ni) = 1.d0/ad(ni)
          else
            adr(ni) = 0.0d0
          end if
        end do ! ni

        do ni = 1,neq
          z(ni) = adr(ni)*r(ni)
          p(ni) = z(ni)
        end do ! ni

c     2. Block preconditioning

      else

c       Assemble upper block array

        alp = 1.0d0
        if(ittyp.gt.2) alp = 0.5d0
        call pzero(aur,jpr(neq))
        do ni = 2,neq
          do nj = jpr(ni-1)+1,jpr(ni)
            do nc = jc(ni-1)+1,jc(ni)
              if(ni-ir(nc).eq.jpr(ni)+1-nj) then
                aur(nj) = ac(nc)*alp
              endif
            end do ! nc
          end do ! nj
        end do ! ni

c       Solve preconditioning equations

        do ni = 1,neq
          adr(ni) = ad(ni)
          z(ni)   = r(ni)
        end do ! ni

        call datri(aur,aur,adr,jpr, neq, neq)
        call dasol(aur,aur,adr,z,jpr,neq,neq,bet,.false.)

        do ni = 1,neq
          p(ni) = z(ni)
        end do ! ni

      end if


      gam    = dot(r,z,neq)
      noconv = .true.
      captim = .true.
      nj     = 1
      tt0    = etime(tary)   ! Timing start
      tt0    = tary(1)
      do while (noconv)

c       Form A*p product for current tangent

        call caprod(ad,ac,p,z,jc,ir,neq)

        if(debug .and. captim) then
          captim = .false.
          tt     = etime(tary)
          tt     = tary(1) - tt0       ! CPU Time increment
          if(tt.gt.0.0d0) then
            nopcap = neq + 4*jc(neq) ! Number ops/A*p products
            if(nopcap.ge.1000000) then
              megops = (nopcap/1000000)
              nopcap =  mod(nopcap,1000000)
            endif
            mflops   = (dble(megops) + dble(totops)*1.d-06)/tt
            write(iow,2001) mflops,tt
            if(ior.lt.0) then
              write(*,2001) mflops,tt
            endif
          endif
        endif

c       Compute update vector

        alp = gam/dot(p,z,neq)

        do ni = 1,neq
          x(ni) = x(ni) + alp*p(ni)
        end do ! ni

        if(mod(nj,50).eq.0) then
          call caprod(ad,ac,x,z,jc,ir,neq)
          do ni = 1,neq
            r(ni) = b(ni) - z(ni)
          end do ! ni
        else
          do ni = 1,neq
            r(ni) = r(ni) - alp*z(ni)
          end do ! ni
        endif

c       Check convergence

        rn     = dot(r,r,neq)
        if(mod(nj,20).eq.0) then
          if(prnt) then
            tt = etime(tary)
            if(ior.lt.0) then
              write(*,2000) nj,sqrt(rn),tary,
     &                      sqrt(rn/rn0),tol
            endif
            write(iow,2000) nj,sqrt(rn),tary,
     &                      sqrt(rn/rn0),tol
          endif

        endif
        noconv = (rn.gt.tolr0) .and. ( nj.lt.jmax)
        nj     = nj + 1

c       Diagonal preconditioning

        if(ittyp.eq.1) then
          do ni = 1,neq
            z(ni) = adr(ni)*r(ni)
          end do ! ni

c       Block preconditioning

        else
          do ni = 1,neq
            z(ni) = r(ni)
          end do ! ni
          call dasol(aur,aur,adr,z,jpr,neq,neq,bet,.false.)
        endif

c       Compute final update vector for step

        bet = gam
        gam = dot(r,z,neq)
        bet = gam/bet
        do ni = 1,neq
          p(ni) = z(ni) + bet*p(ni)
        end do ! ni

      end do ! while

c     Output solution data

      if(mod(nj-1,20).ne.0 .and. prnt) then
        tt = etime(tary)
        if(ior.lt.0) then
          write(*,2000) nj-1,sqrt(rn),tary,
     &                  sqrt(rn/rn0),tol
        endif
        write(iow,2000) nj-1,sqrt(rn),tary,
     &                  sqrt(rn/rn0),tol
      endif

c     Produce flop rate

      if(prnt) then
        nj     = nj - 1             ! Number of actual iterations
        tt     = tary(1) - tt0      ! Increment of time
        if(tt.gt.1.0d-1) then
          numadd = (nj-1)/50          ! Additional A*p products
            numcap = nj + numadd        ! Total number A*p products
          numnoc = nj - numadd
          nopcap = neq + 4*jc(neq)  ! Number ops/A*p products
          if(nopcap.ge.1000000) then
            megops = (nopcap/1000000)*numcap
            nopcap =  mod(nopcap,1000000)
          endif
          nopcgd = neq*11            ! CG ops
          if(nopcgd.ge.1000000) then
            megops = megops + (nopcgd/1000000)*nj
            nopcgd = mod(nopcgd,1000000)
          endif
          totops = nopcap*numcap + neq*numadd + nopcgd*nj + neq*2*numnoc
          if(totops.ge.1000000) then
            megops = megops + totops/1000000
            totops = mod(totops,1000000)
          endif
          mflops   = (dble(megops) + dble(totops)*1.d-06)/tt
          write(iow,2002) mflops,tt
          if(ior.lt.0) then
            write(*,2002) mflops,tt
          endif
        endif
      endif

c     Format

 2000 format('  CG: Iter. =',i5,' Resid. =',1p,1e10.3,
     &        ' Cpu. =',1p,1e10.3,' Sys. =',1p,1e10.3/
     &                  14x,'Rel. Resid. =',1p,1e10.3,
     &           '     :    Solution tol =',1p,1e10.3)
 2001 format(/'  CG: CAPROD at ',f10.2,' Mflops. Time = ',f12.2)
 2002 format(/'  CG: Solve at ',f10.2,' Mflops. Time = ',f12.2)

      end
