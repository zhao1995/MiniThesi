c$Id:$
      subroutine parfeaps(mode,d,v,resid,workl,workd,
     &                    numpeq,vdim,nev,ncv,maxitr,sigma,stol,prt)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014 Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c     Author:   D.S. Bindel  6/21/2005
c     Modified: R.L. Taylor 12/14/2005
c       Original version                                    01/11/2006
c       1. added 'idptr.h'                                  05/01/2013
c       2. updated psetb call to new format                 05/01/2013
c       3. Change pneq to numpeq                            07/01/2013
c       4. Check for dividve by zero of d(j,1)              11/04/2013
c          Fix output of nconv (iparam(5))
c-----[--.----+----.----+----.-----------------------------------------]
c     Purpose:  Driver for symmetric generalized eigenproblem (dynamics)
c               Modified from the symmetric driver in the ARPACK
c               distribution. Uses an interface to PETSc to perform
c               parallel solutions.
c
c               Uses mode 1 Regular solve for lumped mass
c               Uses mode 3 (shift-invert) with shift sigma

c      Inputs:
c         mode   - Solution mode: 1 for lump; 3 for other
c         resid  - Working space
c         workl  - Working space
c         workd  - Working space
c         numpeq - Number of equations in partition
c         vdim   - Leading dimension of vectors 'v'
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

#     include   "finclude/petsc.h"

#     include   "sdata.h"
#     include   "setups.h"
#     include   "iofile.h"
#     include   "idptr.h"
#     include   "pointer.h"
#     include   "comblk.h"

      logical    prt, setvar, palloc
      integer    numpeq, vdim, nev, ncv, maxitr, it, pmaxr,pmaxs

      real*8     sigma, stol
      real*8     d(ncv,2), v(vdim,*), resid(*), workl(*), workd(*)

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
      integer    ido, lworkl, info, j, ierr
      integer    nconv, ishfts, mode
      logical    rvec
      real*8     pdnrm2

c     %------------%
c     | Parameters |
c     %------------%

      real*8     zero
      parameter (zero = 0.0D+0)

c     save

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

      ishfts    = 1
c     mode      = 3

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

        call pdsaupd ( ido, bmat, numpeq, which, nev, stol, resid,
     &                 ncv, v, vdim, iparam, ipntr, workd, workl,
     &                 lworkl, info )


        if (ido .eq. -1) then

c         %--------------------------------------------%
c         | Mode 1:                                    |
c         |  Perform  y <--- OP*x = L*A*L*x            |
c         |                         L = 1/sqrt(M)      |
c         | Mode 2:                                    |
c         |  Perform  y <--- OP*x = inv[M]*A*x         |
c         | Mode 3:                                    |
c         |  Perform  y <--- OP*x = inv[A-SIGMA*M]*M*x |
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
            call paropk(numpeq, workd(ipntr(1)), workd(ipntr(2)))
          else
            call paropm(numpeq, workd(ipntr(1)), workd(ipntr(2)))
            call paropk(numpeq, workd(ipntr(2)), workd(ipntr(2)))
          endif

c         %-----------------------------------------%
c         | L O O P   B A C K to call DSAUPD again. |
c         %-----------------------------------------%

          go to 10

        else if (ido .eq. 1) then

c         %--------------------------------------------%
c         | Mode 1:                                    |
c         |  Perform  y <--- OP*x = L*A*L*x            |
c         |                         L = 1/sqrt(M)      |
c         | Mode 2:                                    |
c         |  Perform  y <--- OP*x = inv[M]*A*x         |
c         | Mode 3:                                    |
c         |  Perform  y <--- OP*x = inv[A-SIGMA*M]*M*x |
c         | M*x has been saved in workd(ipntr(3)).     |
c         | the user only needs the linear system      |
c         | solver here that takes workd(ipntr(3)      |
c         | as input, and returns the result to        |
c         | workd(ipntr(2)).                           |
c         %--------------------------------------------%

          if(mode.eq.1 .or. mode.eq.2) then
            call paropk(numpeq, workd(ipntr(1)), workd(ipntr(2)))
          else
            call paropk(numpeq, workd(ipntr(3)), workd(ipntr(2)))
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
            call paropm(numpeq, workd(ipntr(1)), workd(ipntr(2)))
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

        write(*,2001) info

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

        call pdseupd ( rvec, 'A', select, d, v, vdim, sigma,
     &       bmat, numpeq, which, nev, stol, resid, ncv, v, vdim,
     &       iparam, ipntr, workd, workl, lworkl, ierr )


c       Compute max message buffer size for ghost communication

        pmaxr = 0
        pmaxs = 0
        do j = 0, ntasks-1
          pmaxr = max( pmaxr,mr(np(248)+j) )
          pmaxs = pmaxs + mr(np(250)+j)
        end do ! j
        pmaxr = pmaxr*ndf
        pmaxs = pmaxs*ndf

c       Perform communication and set of ghost values in eigen-solution

        setvar = palloc( 120, 'TEMP0', max(1,pmaxr+pmaxs), 2)
        do j = 1,ncv
          call psetb(v(1,j),mr(np(248)),mr(np(249)),mr(np(250)),
     &         mr(np(251)),mr(id31), ndf, hr(np(120)),hr(np(120)+pmaxr))
          if(mode.eq.1) then
            call parmv(numpeq,v(1,j))
            call psetb(v(1,j),mr(np(248)),mr(np(249)),mr(np(250)),
     &         mr(np(251)),mr(id31), ndf, hr(np(120)),hr(np(120)+pmaxr))
          endif
        end do ! j
        setvar = palloc( 120, 'TEMP0',    0, 2)

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

          write(*,2001) ierr

        else

          nconv =  iparam(5)
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

            call parkv (numpeq,  v(1,j), workd)
            call paropm(numpeq,  v(1,j), workd(vdim+1))
            call daxpy (numpeq, -d(j,1), workd(vdim+1), 1, workd, 1)
            d(j,2) = pdnrm2(numpeq, workd, 1)
            if(abs(d(j,1)).gt.0.0d0) then
              d(j,2) = d(j,2) / abs(d(j,1))
            else
              d(j,2) = 1.0d0 ! zero
            endif
          end do ! j

        end if

c       Output computed eigenvalues

        it = iparam(9)
        if(ior.gt.0) then
           write(iow,2005) it,(d(j,1),j=1,nev)
           write(iow,2006) it,(d(j,2),j=1,nev)
        elseif(prt) then
           write(  *,2005) it,(d(j,1),j=1,nev)
           write(  *,2006) it,(d(j,2),j=1,nev)
        endif

c       %------------------------------------------%
c       | Print additional convergence information |
c       %------------------------------------------%

        if ( info .eq. 1) then
          write(*,2002)
        else if ( info .eq. 3) then
          write(*,2003)
        end if
      end if

      if(prt) then
        write(*,2004) numpeq,nev,ncv,which,iparam(5),iparam(3),
     &                iparam(9),stol
      endif
      write(iow,2004) numpeq,nev,ncv,which,iparam(5),iparam(3),
     &                iparam(9),stol

c     Formats

2001  format(/10x,'Error with dsaupd, info = ',i8/
     &        10x,'Check documentation of dsaupd '/)

2002  format(/10x,'Maximum number of iterations reached.'/)

2003  format(/10x,'No shifts could be applied during implicit',/,
     &        10x,'Arnoldi update, try increasing NCV.'/)

2004  format(/15x,'PARFEAP',/,15x,' ====== '//
     &        10x,'Size of the partion . . . . . . . . . . . =',i8/
     &        10x,'Number of Ritz values requested . . . . . =',i8/
     &        10x,'Number of Arnoldi vectors generated (NCV) =',i8/
     &        10x,'Portion of the spectrum: ',a/
     &        10x,'Number of converged Ritz values . . . . . =',i8/
     &        10x,'Number Implicit Arnoldi update iterations =',i8/
     &        10x,'Number of OP*x. . . . . . . . . . . . . . =',i8/
     &        10x,'Convergence criterion . . . . . . . . . . =',
     &             1p,1e12.5/)

2005  format(/'  PARFEAP: Current eigenvalues, iteration',i7/
     &        (5x,1p,4d17.8))

2006  format( '  PARFEAP: Current residuals,   iteration',i7/
     &        (5x,1p,4d17.8))

      end
