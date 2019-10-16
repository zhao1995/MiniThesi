c$Id:$
      subroutine arfeaps(mode,d,v,msqrt,resid,workl,workd,nev,ncv,
     &                   maxitr,sigma,stol)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014 Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c     Author  : D.S. Bindel  6/21/2005
c     Modified: R.L. Taylor 12/14/2005
c       Original version                                    01/11/2006
c       1. Added output by peigout                          03/09/2007
c          Add computation of residuals.
c       2. Allow for computation of all eigenpairs          07/09/2007
c       3. Include 'sigma' in d(j,2) computation            09/09/2007
c       4. Add computation of residual for mode = 1.        13/09/2007
c          Normalize by v^T*M*v instead of d(1,j)
c       5. Remove prt argument                              21/04/2011
c       6. Remove vneq from argument                        16/02/2013
c-----[--.----+----.----+----.-----------------------------------------]
c     Purpose: Driver for symmetric eigenproblem (dynamics).
c              Modified from the symmetric driver in the ARPACK
c              distribution.

c              Uses mode 1 (Regular solve for lumped mass)
c              Uses mode 3 (shift-invert) with shift sigma

c      Parameter from common block
c         neq    - Number of equations

c      Inputs:
c         mode   - Solution mode: 1 for lump; 3 for other
c         msqrt  - Reciprocal square root of mass (mode = 1 only)
c         resid  - Working space
c         workl  - Working space
c         workd  - Working space
c         nev    - Number of eigenvalues
c         ncv    - Number of vectors
c         maxitr - Maximum number of iterations
c         sigma  - Shift
c         stol   - Tolerance on eigenvalues

c      Outputs:
c         d      - Eigenvalues
c         v      - Eigenvectors
c-----[--.----+----.----+----.-----------------------------------------]
      implicit   none

      include   'cdata.h'
      include   'compas.h'
      include   'debugs.h'
      include   'sdata.h'
      include   'iofile.h'
      include   'pointer.h'
      include   'comblk.h'

      integer    nev, ncv, maxitr, it, nnv
      real*8     sigma, stol
      real*8     d(ncv,2), v(neq,*), resid(*), workl(*), workd(*)
      real*8     msqrt(*)

c     %---------------%
c     | Local Arrays  |
c     %---------------%

      integer    maxncv
      parameter (maxncv = 3000)
      logical    select(maxncv)
      integer    iparam(11), ipntr(11)

c     %---------------%
c     | Local Scalars |
c     %---------------%

      character  bmat*1, which*2
      integer    ido, lworkl, info, i,j, ierr
      integer    nconv, ishfts, mode
      real*8     dnrm2
      logical    rvec

c     %------------%
c     | Parameters |
c     %------------%

      real*8     zero
      parameter (zero = 0.0D+0)

      save

c     %-----------------------%
c     | Executable statements |
c     %-----------------------%

      if(mode.eq.1) then
        bmat   = 'I'
        which  = 'SM'
      elseif(mode.eq.2) then
        bmat   = 'G'
        which  = 'SM'
      else
        bmat   = 'G'
        which  = 'LM'
      endif
      lworkl = ncv*(ncv+8)
      ido    = 0
      info   = 0

c     %---------------------------------------------------%
c     | This program uses exact shifts with respect to    |
c     | the current Hessenberg matrix (IPARAM(1) = 1).    |
c     | IPARAM(3) specifies the maximum number of Arnoldi |
c     | iterations allowed.  Mode 3 specified in the      |
c     | documentation of DSAUPD is used (IPARAM(7) = 3).  |
c     | All these options may be changed by the user.     |
c     | For details, see the documentation in DSAUPD.     |
c     %---------------------------------------------------%

      ishfts    =  1

      iparam(1) = ishfts
      iparam(3) = maxitr
      iparam(7) = mode

c     %-------------------------------------------%
c     | M A I N   L O O P (Reverse communication) |
c     %-------------------------------------------%

 10   continue

c       %---------------------------------------------%
c       | Repeatedly call the routine DSAUPD and take |
c       | actions indicated by parameter IDO until    |
c       | either convergence is indicated or maxitr   |
c       | has been exceeded.                          |
c       %---------------------------------------------%

        call dsaupd ( ido, bmat, neq, which, nev, stol, resid,
     &                ncv, v, neq, iparam, ipntr, workd, workl,
     &                lworkl, info )

        if (ido .eq. -1) then

c         %--------------------------------------------%
c         | Perform  y <--- OP*x = inv[A-SIGMA*M]*M*x  |
c         | to force the starting vector into the      |
c         | range of OP.  The user should supply       |
c         | his/her own matrix vector multiplication   |
c         | routine and a linear system solver here.   |
c         | The matrix vector multiplication routine   |
c         | takes workd(ipntr(1)) as the input vector. |
c         | The final result is returned to            |
c         | workd(ipntr(2)).                           |
c         %--------------------------------------------%

          if(mode.eq.1 .or. mode.eq.2) then
            call aropk(neq, workd(ipntr(1)), workd(ipntr(2)),
     &                      workd(ipntr(3)), mode)
          else
            call aropm(neq, workd(ipntr(1)), workd(ipntr(2)))
            call aropk(neq, workd(ipntr(2)), workd(ipntr(2)),
     &                      workd(ipntr(3)), mode)
          endif

c         %-----------------------------------------%
c         | L O O P   B A C K to call DSAUPD again. |
c         %-----------------------------------------%

          go to 10

        else if (ido .eq. 1) then

c         %-----------------------------------------%
c         | Perform y <-- OP*x = inv[A-sigma*M]*M*x |
c         | M*x has been saved in workd(ipntr(3)).  |
c         | the user only needs the linear system   |
c         | solver here that takes workd(ipntr(3)   |
c         | as input, and returns the result to     |
c         | workd(ipntr(2)).                        |
c         %-----------------------------------------%

          if(mode.eq.1 .or. mode.eq.2) then
            call aropk(neq, workd(ipntr(1)), workd(ipntr(2)),
     &                      workd(ipntr(3)), mode)
          else
            call aropk(neq, workd(ipntr(3)), workd(ipntr(2)),
     &                      workd(ipntr(3)), mode)
          endif

c         %-----------------------------------------%
c         | L O O P   B A C K to call DSAUPD again. |
c         %-----------------------------------------%

          go to 10

        else if (ido .eq. 2) then

c         %-----------------------------------------%
c         |          Perform  y <--- M*x            |
c         | Need the matrix vector multiplication   |
c         | routine here that takes workd(ipntr(1)) |
c         | as the input and returns the result to  |
c         | workd(ipntr(2)).                        |
c         %-----------------------------------------%

          if(mode.eq.2 .or. mode.eq.3) then
            call aropm(neq, workd(ipntr(1)), workd(ipntr(2)))
          endif

c         %-----------------------------------------%
c         | L O O P   B A C K to call DSAUPD again. |
c         %-----------------------------------------%

          go to 10

        end if

c     %-----------------------------------------%
c     | Either we have convergence, or there is |
c     | an error.                               |
c     %-----------------------------------------%

      if ( info .lt. 0 ) then

c       %--------------------------%
c       | Error message, check the |
c       | documentation in DSAUPD. |
c       %--------------------------%

        write(*,2000) info

      else

c       %-------------------------------------------%
c       | No fatal errors occurred.                 |
c       | Post-Process using DSEUPD.                |
c       |                                           |
c       | Computed eigenvalues may be extracted.    |
c       |                                           |
c       | Eigenvectors may also be computed now if  |
c       | desired.  (indicated by rvec = .true.)    |
c       %-------------------------------------------%

        rvec = .true.

        call dseupd ( rvec, 'All', select, d, v, neq, sigma,
     &       bmat, neq, which, nev, stol, resid, ncv, v, neq,
     &       iparam, ipntr, workd, workl, lworkl, ierr )

c       %----------------------------------------------%
c       | Eigenvalues are returned in the first column |
c       | of the two dimensional array D and the       |
c       | corresponding eigenvectors are returned in   |
c       | the first NEV columns of the two dimensional |
c       | array V if requested.  Otherwise, an         |
c       | orthogonal basis for the invariant subspace  |
c       | corresponding to the eigenvalues in D is     |
c       | returned in V.                               |
c       %----------------------------------------------%

        if ( ierr .ne. 0 ) then

c         %------------------------------------%
c         | Error condition:                   |
c         | Check the documentation of DSEUPD. |
c         %------------------------------------%

          write(*,2000) ierr

        else

          nnv   =  2*neq
          nconv =  iparam(5)
          if(ncv-nev.eq.1) then
            nconv = nconv + 1
          endif
          do j=1, nconv

c           %---------------------------%
c           | Compute the residual norm |
c           |                           |
c           |   ||  A*x - lambda*x ||   |
c           |                           |
c           | for the NCONV accurately  |
c           | computed eigenvalues and  |
c           | eigenvectors.  (iparam(5) |
c           | indicates how many are    |
c           | accurate to the requested |
c           | tolerance)                |
c           %---------------------------%

            if(mode.eq.1) then

              do i = 1,neq
                workd(i+nnv) = v(i,j)
              end do ! i
              call aropk(neq, workd(nnv+1), workd,  ! N.B. Compute A*v
     &                   workd(neq+1),  1)
              call daxpy (neq, -d(j,1)+sigma, workd(nnv+1), 1, workd, 1)
              d(j,2) = dnrm2(neq, workd, 1)/ dnrm2(neq,workd(nnv+1), 1)

c             Scale eigenvectors by reciprocal mass square root

              do i = 1,neq
                v(i,j) = v(i,j)*msqrt(i)
              end do ! i

            elseif(mode.eq.3 .and. ittyp.eq.-3) then ! Profile solver
              do i = 1,neq
                workd(i+nnv) = v(i,j)
              end do ! i
              call aropk(neq, workd(nnv+1), workd,  ! N.B. Compute A*v
     &                   workd, -1)                 !      from factors
              call aropm(neq, v(1,j), workd(neq+1))
              call daxpy (neq, -d(j,1)+sigma, workd(neq+1), 1, workd, 1)
              d(j,2) = dnrm2(neq, workd, 1)/ dnrm2(neq, workd(neq+1), 1)
            else
              d(j,2) = zero
            endif
          end do ! j

c         Output results

          it = iparam(9)
c         do j = 1,nev
          do j = 1,nconv
            workl(j) = sqrt(abs(d(j,1)))
          end do ! j
          call peigout('ARFEAP',d,workl,workl(1+nconv),d(1,2),
     &                  stol,nconv,it)

        end if

c       %------------------------------------------%
c       | Print additional convergence information |
c       %------------------------------------------%

        if ( info .eq. 1) then
          write(*,2002)
        else if ( info .eq. 3) then
          write(*,2003)
        end if

        if(ior.lt.0 .and. debug) then
          write(  *,2004) neq,nev,ncv,which,nconv,
     &                    iparam(3),iparam(9),stol
        elseif(debug) then
          write(iow,2004) neq,nev,ncv,which,nconv,
     &                    iparam(3),iparam(9),stol
        end if

      end if

c     Formats

2000  format(/10x,'Error with dsaupd, info = ',i8/
     &        10x,'Check documentation of dsaupd '/)


2002  format(/10x,'Maximum number of iterations reached.'/)

2003  format(/10x,'No shifts could be applied during implicit',/,
     &        10x,'Arnoldi update, try increasing NCV.'/)

2004  format(/15x,'ARFEAP',/,15x,' ====== '//
     &        10x,'Size of the matrix. . . . . . . . . . . . =',i8/
     &        10x,'Number of Ritz values requested . . . . . =',i8/
     &        10x,'Number of Arnoldi vectors generated (NCV) =',i8/
     &        10x,'Portion of the spectrum: ',a/
     &        10x,'Number of converged Ritz values . . . . . =',i8/
     &        10x,'Number Implicit Arnoldi update iterations =',i8/
     &        10x,'Number of OP*x. . . . . . . . . . . . . . =',i8/
     &        10x,'Convergence criterion . . . . . . . . . . =',
     &             1p,1e12.5/)

      end
