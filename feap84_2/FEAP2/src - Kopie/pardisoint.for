c----------------------------------------------------------------------+
c     Interface  to use PARDISO parallel solver in FEAP
c     INTEL 32/64Bit
c
c     W. Wagner BS KIT 01/15 Check for INTEL COMPOSER 15
c----------------------------------------------------------------------+
c
      subroutine mkptr4(isymcsr)
c-----------------------------------------------------------------------
c
c     Purpose: Define arrays for PARDISO direct solver
c              pointer with CSR storage in mkptr_csr
c
c     Inputs: -
c
c     Output: isymcsr   1 = symmetric 2 = unsymmetric matrix
c
c     W. Wagner BS UKA 09/10
c-----------------------------------------------------------------------
c
      USE cdata
      USE iofile
      USE pardi
      USE pdata2
      USE soltyp
      USE doalloc
      implicit double precision (a-h,o-z)

#ifndef _BATCH_
      if(idev.ne.3) stop 'solver only possible for idev=3(INTEL)'
#endif

c.... set parameter
      isymcsr = ctis(1)
      nproc   = ctis(2)

c.... default values
      if(isymcsr.ne.2) isymcsr = 1  ! 1= symmetric matrix 2= unsymmetric

c.... no of processors
      call get_number_proc(nproc)

      if(ior.lt.0) write(*,*) 'PARDISO: Number of processors used',nproc
      write(iow,*) 'PARDISO: Number of processors used',nproc

c.... Separate need of x and b for Ax=b and R=x^T*b
      call ralloc(drpar,neq,'Pardiso-x',lrhp)

c.... set Pardiso parameter
c     ipt(64)     Internal solver memory pointer
      call pzeroi(ipt,64)

c     iparm(64)
      call pzeroi(iparm,64)

c     dparm(64)
      call pzero(dparm,64)

c     maxfct Max no of factors with identical nonzero sparsity structure
      maxfct = 1

c     mnum        matrix to factorize, def=1
      mnum = 1

c     nrhs        number of right-hand sides that need to be solved for
      nrhs = 1

c     mtype       defines the matrix type
c                 1 real and structurally symmetric
c                 2 real and symmetric positive definite
c                -2 real and symmetric indefinite
c                11 real and nonsymmetric
      if(isymcsr.eq.1) mtype = -2
      if(isymcsr.eq.2) mtype = 11

C     PARDISO solver
      solver = 0 ! use sparse    direct solver

c
c     iphase     11 Analysis
c                12 Analysis, numerical factorization
c                13 Analysis, numerical factorization,
c                             solve, iterative refinement
c
c                22 Numerical factorization
c                23 Numerical factorization, solve, iterative refinement
c
c                33 Solve, iterative refinement
c
c                -1 Release all internal memory for all matrices
      iphase = 0

c     msglvl     Message level information.
c                0 no output
c                1 solver prints statistical information to the screen
      msglvl    = 1
      ierror    = 0 ! initialize error flag 

c.... Set up PARDISO control parameter

      iparm(1) = 1   ! no solver default
cww   iparm(2) = 2   ! fill-in reordering from METIS
      iparm(2) = 3   ! fill-in reordering from METIS-parallel (OpenMP)version
      iparm(3) = nproc ! no. of processors, but not used by INTEL-Pardiso
      iparm(4) = 0   ! no iterative-direct algorithm
      iparm(5) = 0   ! no user fill-in reducing permutation
      iparm(6) = 0   ! =0 AX=B   B not overwritten 
      iparm(7) = 0   ! Output: No. of iterative refinement steps 
      iparm(8) = 9   ! numbers of iterative refinement steps
      iparm(9) = 0   ! not in use
      iparm(10) = 13 ! perturb the pivot elements with 1E-13
      iparm(11) = 1  ! use nonsymmetric permutation and scaling MPS
      iparm(12) = 0  ! Solve AX=B
      iparm(13) = 1  ! maximum weighted matching algorithm is switched-on (default for non-symmetric) Try iparm(13) = 1 in case of inappropriate accuracy
      iparm(14) = 0  ! Output: number of perturbed pivots
      iparm(15) = 0  ! Output: Peak memory symbolic factorization in KBytes-Phase 1
      iparm(16) = 0  ! Output: Permanent memory symbolic factorization in KBytes-Phase 1
      iparm(17) = 0  ! Output: Memory numerical factorization and solution in KBytes-Phase 2
      iparm(18) = -1 ! Output: number of nonzeros in the factor LU
      iparm(19) = -1 ! Output: Mflops for LU factorization
      iparm(20) = 0  ! Output: number of CG Iterations
      iparm(21) = 1  ! Bunch-Kaufman pivoting
      iparm(22) = 0  ! Output: number of positive eigenvalues
      iparm(23) = 0  ! Output: number of negative eigenvalues
      iparm(24) = 1  ! Faster Parallel Numerical Factorization
      iparm(25) = 0  ! Parallel Forward/Backward Solve
      iparm(30) = 1  ! Output: number of zero or negative pivots
cww   iparm(33) = 1  ! compute determinant of real symmetric indefinite matrix, only BASEL

      
c     iparm(60) = 0  ! Default IC InCore
cww   iparm(60) = 1  ! IC/OOC InCore/Out Of Core
cww   iparm(60) = 2  ! OOC


      return
      end
c
c----------------------------------------------------------------------+
c
      subroutine mkptr8(isymcsr)
c-----------------------------------------------------------------------
c
c     Purpose: Define arrays for PARDISO solver
c              pointer with CSR storage in mkptr_csr
c
c     Inputs: -
c
c     Output: isymcsr   1 = symmetric 2 = unsymmetric matrix
c
c     W. Wagner BS UKA 06/09
c-----------------------------------------------------------------------
c
      USE cdata
      USE iofile
      USE pardi
      USE pdata2
      USE soltyp
      USE doalloc
      implicit double precision (a-h,o-z)
      logical ldummy

      if(idev.ne.3) stop 'solver only possible for idev=3(INTEL)'

c.... set parameter

c.... default values
      isymcsr = 1 ! only symmetric matrix allowed
        nproc = 1 ! only sequentiell allowed

c.... no of processors
      call get_number_proc(nproc)

      if(ior.lt.0) write(*,*) 'PARDISO: Number of processors used',nproc
      write(iow,*) 'PARDISO: Number of processors used',nproc

c.... Separate need of x and b for Ax=b and R=x^T*b
      call ralloc(drpar,neq,'Pardiso-x',lrhp)

c.... set Pardiso parameter
c     ipt(64)     Internal solver memory pointer
      call pzeroi(ipt,64)

c     iparm(64)
      call pzeroi(iparm,64)

c     dparm(64)
      call pzero(dparm,64)

c     maxfct Max no of factors with identical nonzero sparsity structure
      maxfct = 1

c     mnum        matrix to factorize, def=1
      mnum = 1

c     nrhs        number of right-hand sides that need to be solved for
      nrhs = 1

c     mtype       defines the matrix type
c                 1 real and structurally symmetric
c                 2 real and symmetric positive definite
c                -2 real and symmetric indefinite
c                11 real and nonsymmetric
      if(isymcsr.eq.1) mtype = -2
      if(isymcsr.eq.2) mtype = 11

C     PARDISO solver
      solver = 1 ! use iterative multilevel solver


      if(ierror.ne.0) call erropardi(ierror)

c
c     iphase     11 Analysis, reordering and symbolic factorization
c                12 Analysis, numerical factorization
c                13 Analysis, numerical factorization,
c                             solve, iterative refinement
c
c                22 Numerical factorization
c                23 Numerical factorization, solve, iterative refinement
c
c                33 Solve, iterative refinement
c
c                 0 Release all internal memory for L and U
c                -1 Release all internal memory for all matrices
      iphase = 0

c     msglvl     Message level information.
c                0 no output
c                1 solver prints statistical information to the screen
c                2 solver prints statistical information to pardiso-ml.out
      msglvl = 2

      iparm(3) = nproc

cww   iparm(23) = nneg ! Inertia: Number of negative eigenvalues is returned in  iparm(23)

      iparm(32) = 1    ! use multirecursive iterative solver
      iparm(33) = 1    ! compute determinant

      niter     = ctis(1)
      dtol      = ctis(2)
      niterimp  = ctis(3)


      dparm(1) = niter
      if(niter.eq.0)
     +dparm(1) = 300   !  Max Iteration in Krylov-Subspace Iteration
      dparm(2) = dtol
      if(dtol.eq.0.d0)
     +dparm(2) = 1.d-6 !  Relative Residual Tolerance for Convergencd
      dparm(3) = 5000  !  Size of Coarse Grid Matrix
      dparm(4) = 10    !  Maximum Number of Levels in Grid Hierachy
      dparm(5) = 1.e-2 !  Dropping Threshold for Incomplete Factor
      dparm(6) = 5.e-3 !  Dropping Threshold for Schurcomplement
      dparm(7) = 10000 !  Max number of Fill-in for each column in the factor (def=10)
      dparm(8) = 5000  !  Bound for the Norm of the Inverse of L (def=500)
      dparm(9) = niterimp
      if(niterimp.eq.0)
     +dparm(9) = 25    !  No of Iterations without improvement of solution

      return
      end
c

c----------------------------------------------------------------------
c
      subroutine datri4(a,ia,neq,nneg)
c-----------------------------------------------------------------------
c
c.... Purpose: factorization  using PARDISO solver
c
c-----------------------------------------------------------------------
c      Inputs:
c         a(*)       - Non factored terms of A
c         ia(*)      - Pointer to row    CSR sparse technique
c         ja(*)      - Pointer to column CSR=csrja
c         neq        - Number of equations
c
c      Outputs:
c         nneg       - Number of negative diagonal elements
c
c      Open:
c
c     W. Wagner BS UKA 04/09
c----------------------------------------------------------------------
      USE fdata
      USE iofile
      USE iscsr
      USE pardi
      implicit double precision (a-h,o-z)
      real*4  tary,tary1
      dimension a(*),ia(*)

      if(iphase.eq.0) then ! only for the first time
c....   Reordering and Symbolic Factorization, This step also allocates
c       all memory that is necessary for the factorization

        iphase    = 11      ! only reordering and symbolic factorization
c       idum, ddum: dummys
        call etimef(tary)
        if(pfr)                write(iow,2014) tary
        if(pfr .and. ior.lt.0) write(*  ,2014) tary

        call pardiso (ipt,maxfct,mnum,mtype,iphase,neq,a,ia,csrja,
     1                idum,nrhs,iparm,msglvl,ddum,ddum,ierror,dparm)

        call etimef(tary1)
        if(pfr)                write(iow,2015) tary1-tary
        if(pfr .and. ior.lt.0) write(*  ,2015) tary1-tary

        if(ierror.ne.0) call erropardi(ierror)

        if(solver.eq.0) then
          if(ior.lt.0) then
      write(*,*)'  Peak Memory[KB] for  symb. factorization: ',iparm(15)
      write(*,*)'  Permanent Memory[KB] symb. factorization: ',iparm(16)
      write(*,*)'  Number of nonzeros in factors           : ',iparm(18)
      write(*,*)'  Number of factorization MFLOPS          : ',iparm(19)
cww   write(*,*)'  OOC Memory[KB] for   symb. factorization: ',iparm(63)
        end if
          end if
      write(iow,*)'  Peak Memory[KB] for  symb. factorization: ',
     +iparm(15)
      write(iow,*)'  Permanent Memory[KB] symb. factorization: ',
     +iparm(16)
      write(iow,*)'  Number of nonzeros in factors           : ',
     +iparm(18)
      write(iow,*)'  Number of factorization MFLOPS          : ',
     +iparm(19)
      end if

c.... Factorization

      iphase     = 22  ! only factorization

      call pardiso (ipt,maxfct,mnum,mtype,iphase,neq,a,ia,csrja,
     1              idum,nrhs,iparm,msglvl,ddum,ddum,ierror,dparm)

      if(ierror.ne.0) call erropardi(ierror)

c     Number of neg. eigenvalues for symmetric indef. matrices(mtype=-2)
      if(mtype.eq.-2) nneg=iparm(23)

2014  format('   Begin    symbolic factorization',31x,'t=',0pf9.2)
2015  format('   Time for symbolic factorization',31x,'t=',0pf9.4)
      return
      end
c
c----------------------------------------------------------------------
c
      subroutine dasol4(a,b,ia,neq,energy)
c----------------------------------------------------------------------
c
c      Purpose: Solution of the problem Ax=b  for PARDISO-solver
c
c      Inputs:
c         a(*)        - Non factored terms of A
c         b(*)=drpar  - right hand side vector b
c         ia(*)       - Pointer to row    CSR sparse technique
c         ja(*)       - Pointer to column CSR =csrja
c         neq         - Number of equations
c
c      Outputs:
c         b(*)        - Solution vector x
c         energy      - Energy residual
c
c
c     W. Wagner BS UKA 04/09
c----------------------------------------------------------------------
      USE iofile
      USE iscsr
      USE pardi
      implicit double precision (a-h,o-z)

      dimension a(*),b(*),ia(*)

c.... test: copy a
c      write(iow,*) 'a'
c      do i = 1,549
c       write(iow,*) i,a(i)
c      end do

c.... test: copy b
c      write(iow,*) 'b'
c      do i = 1,neq
c       write(iow,*) i, b(i)
c      end do


c.... test ia,ja,a
c      call ioprof(ia,csrja,neq,1)


c...  copy rhs to drpar
      call matcop(b,neq,1,drpar)

c..   Back substitution and iterative refinement
c##   iparm(8)  = 1   ! max numbers of iterative refinement steps
      iparm(8)  = 2   ! max numbers of iterative refinement steps
      iphase    = 33  ! only factorization

      call pardiso (ipt,maxfct,mnum,mtype,iphase,neq,a,ia,csrja,
     1              idum,nrhs,iparm,msglvl,drpar,b,ierror,dparm)

      if(ierror.ne.0) call erropardi(ierror)

c.... control for iterative multilevel solver
      if(solver.eq.1) then
        write(  *,2000) dparm(35),dparm(34)
        write(iow,2000) dparm(35),dparm(34)
      end if

c.... energy
      energy = ddot(neq,b,1,drpar,1)

c     iphase=-1 (clear) here not used! arrays will be overwritten by Pardiso
c     but used separately in case of RINP!

c     formats
2000  format(1x,'No. of. Iterations: ',f5.0,' Residuum',e12.5)

      return
      end
c
c----------------------------------------------------------------------
c
      subroutine erropardi(ierror)
c-----------------------------------------------------------------------
c
c.... Purpose: print errors PARDISO solver
c
c-----------------------------------------------------------------------
      USE iofile
      implicit double precision (a-h,o-z)
      character*63 errorind(11)
      data errorind /
     +' -1  Input inconsistent.                                       ',
     +' -2  Not enough memory.                                        ',
     +' -3  Reordering problem.                                       ',
     +' -4  Zero pivot, num. fact. or iterative refinement problem.   ',
     +' -5  Unclassified (internal) error.                            ',
     +' -6  Preordering failed (matrix types 11, 13 only).            ',
     +' -7  Diagonal matrix problem.                                  ',
     +' -8  32-bit integer overflow problem                           ',
     +' -9  not enough memory for OOC (from INTEL)                    ',
     +' -10 problems with opening OOC temporary files                 ',
     +' -11 read/write problems with the OOC data file                '/
c
      ierrori =iabs(ierror)
                   write(iow,1000) errorind(ierrori)
      if(ior.lt.0) write(  *,1000) errorind(ierrori)

1000  format(5x,'Error in solver Pardiso',/,a62)

      return
      end
c
c----------------------------------------------------------------------
c
      subroutine clear12(neq)
c----------------------------------------------------------------------
c
c      Purpose: release internal memory used by PARDISO-solver
c
c      Inputs: -
c
c      Comments:
c        iphase=-1 (clear) has to be used only in case of RINP
c        otherwise: arrays will be overwritten by Pardiso
c
c
c     W. Wagner BS UKA 04/09
c----------------------------------------------------------------------
      USE pardi
      implicit double precision (a-h,o-z)

c..   release internal memory
      iphase    = -1

      call pardiso (ipt,maxfct,mnum,mtype,iphase,neq,ddum,idum,idum,
     1              idum,nrhs,iparm,msglvl,ddum,ddum,ierror,dparm)

      return
      end
c
c----------------------------------------------------------------------
c
      subroutine detkt4(nneg,det1)
c----------------------------------------------------------------------
c     Purpose: calculate determinant of stiffness matrix
c              for Pardiso solver: value det1
c      Inputs:
c
c      Outputs:
c        nneg - no. of negative diagonal elements
c        det1 - value of exp(determinant)
c
c----------------------------------------------------------------------
c
      USE iofile
      USE pardi
      implicit double precision (a-h,o-z)

c     Number of neg. eigenvalues for symmetric indef. matrices(mtype=-2)
      nneg=0
      if(mtype.eq.-2) nneg=iparm(23)

c     Determinant
      det1 = 0.d0
      if(iparm(33).eq.1)  then
        det1 = abs(dparm(33))  ! only absolute value needed
      end if

c     if(iparm(33).eq.1)  then  ! original pardiso (including sign)
c       if(dparm(33).le.0) then
c         det1 = -abs(dparm(33))
c       else
c         det1 =      dparm(33)
c       end if
c     end if

      return
      end
c














