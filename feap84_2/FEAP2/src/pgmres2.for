      subroutine mkptr7(ia,isymcsr)
c-----------------------------------------------------------------------
c
c     Purpose: Define arrays for PGMRES2 solver
c              pointer with CSR storage in mkptr_csr
c
c     Inputs:
c
c     Output:
c
c-----------------------------------------------------------------------
      USE cdata
      USE ispgmr
      USE soltyp
      USE doalloc
      implicit double precision (a-h,o-z)
      logical ldummy
      dimension ia(neq+1)

      ldummy = .true.

c.... set parameter for pgmres2
      iterpgmr = ctis(1)
       tolpgmr = ctis(2)
          imp  = ctis(3)

c.... default values
      if(iterpgmr.eq.0 )  iterpgmr = 150    ! No. of Iterations
      if(tolpgmr.eq.0.d0) tolpgmr  = 1.e-08 ! Tolerance
      if(imp.eq.0)        imp      = 50     ! Size of krylov subspace, typical 50, arnoldi size, LÖhnert: 100 here 50=fixed
      isymcsr  = 2                          ! only unsymmetric matrix


      if (ldummy) then
c....   pgmres  
        call ralloc( rmpgmrx, neq,       'GMRES2- mpgmrx',ldummy)
        call ralloc(rmpgmrvv,(imp+1)*neq,'GMRES2-mpgmrvv',ldummy)
        call ralloc(rmpgmrrs,(imp+1),    'GMRES2-mpgmrrs',ldummy)
        call ralloc( rmpgmrc, imp,       'GMRES2- mpgmrc',ldummy)
        call ralloc( rmpgmrs, imp,       'GMRES2- mpgmrs',ldummy)
      end if

      return
      end

c-----------------------------------------------------------------------


      subroutine dasol7 (a,b,ia,neq,energy)
c----------------------------------------------------------------------
c
c      Purpose: Solution of the problem Ax=b  for the PGMRES2 solver
c
c      Inputs:
c         a(*)       - Non factored terms of A
c         b(*)       - right hand side vector b
c         ia(*)      - Pointer to row    CSR sparse technique
c         ja(*)      - Pointer to column CSR
c         neq        - Number of equations
c
c      Outputs:
c         b(*)       - Solution vector x
c         energy     - Energy residual
c
c       Comment:
c         Coefficient matrix is not(!) decomposed into its triangular
c         factors using datri before using dasol.
c
c         Iteration parameters can be set by macro pgmr,iter
c
c       Open:
c
c----------------------------------------------------------------------
      USE iofile
      USE fdata
      USE iscsr
      USE ispgmr
      USE isprec
      implicit double precision (a-h,o-z)
      dimension a(*),b(*),ia(*)

      real*4 t1,t2

      call pzero(rmpgmrx,neq)

c.... here entry for guess of x (mgmrx) is possible
c
c     generate random number as initial guess.
c
c     do  j = 1, n
c        x(j) = rand(j)
c     end do

      call etimef(t1)

c.... preconditioner
      call precond(a,ia,csrja,csrka)

c.....solver
      iout = iow
      iout = 1
      call pgmres (neq,imp,b,rmpgmrx,rmpgmrvv,tolpgmr,iterpgmr,iout,
     +             a,csrja,ia,rmpcalu,impcjlu,impcju,ierr,
     +             rmpgmrc,rmpgmrs,rmpgmrrs,its,ro)

cww   write(*,*) 'error',ierr

cww   if(ierr.eq. 0) write(*,*) 'PGMRES2 successful return.'
cww   if(ierr.eq. 1) write(*,*) 'PGMRES2 no convergence achieved.'
cww   if(ierr.eq.-1) write(*,*) 'PGMRES2 starting residual is zero.'

      if(pfr)                write(iow,2000) its,ro
      if(pfr .and. ior.lt.0) write(*  ,2000) its,ro

cww   call etimef(t2)
cww   if(pfr)                write(iow,2004) t2-t1
cww   if(pfr .and. ior.lt.0) write(*  ,2004) t2-t1

c.... energy
      energy = ddot(neq,b,1,rmpgmrx,1)

c.... copy x->b for further FEAP processing
      call matcop(rmpgmrx,neq,1,b)

      return
2000  format(3x,'PGMRES2: Iteration',i4,' Norm  ',g15.5)
2004  format(3x,'Time for PGMRES2-Solution',37x,'t=',0pf9.4)
      end
c


      subroutine pgmres ( n, im, rhs, sol, vv, eps, maxits, iout,
     +  aa, ja, ia, alu, jlu, ju, ierr, c,s,rs, its,ro )

      !*******************************************************************************
      !
      !! PGMRES is an ILUT - Preconditioned GMRES solver.
      !
      !  Discussion:
      !
      !    This is a simple version of the ILUT preconditioned GMRES algorithm.
      !    The ILUT preconditioner uses a dual strategy for dropping elements
      !    instead  of the usual level of-fill-in approach. See details in ILUT
      !    subroutine documentation. PGMRES uses the L and U matrices generated
      !    from the subroutine ILUT to precondition the GMRES algorithm.
      !    The preconditioning is applied to the right. The stopping criterion
      !    utilized is based simply on reducing the residual norm by epsilon.
      !    This preconditioning is more reliable than ilu0 but requires more
      !    storage. It seems to be much less prone to difficulties related to
      !    strong nonsymmetries in the matrix. We recommend using a nonzero tol
      !    (tol=.005 or .001 usually give good results) in ILUT. Use a large
      !    lfil whenever possible (e.g. lfil = 5 to 10). The higher lfil the
      !    more reliable the code is. Efficiency may also be much improved.
      !    Note that lfil=n and tol=0.0 in ILUT  will yield the same factors as
      !    Gaussian elimination without pivoting.
      !
      !    ILU(0) and MILU(0) are also provided for comparison purposes
      !    USAGE: first call ILUT or ILU0 or MILU0 to set up preconditioner and
      !    then call pgmres.
      !
      !  Modified:
      !
      !    07 January 2004
      !
      !  Author:
      !
      !    Youcef Saad
      !
      !  Parameters:
      !
      !    Input, integer N, the order of the matrix.
      !
      !    Input, integer IM, the size of the Krylov subspace.  IM should not
      !    exceed 50 in this version.  This restriction can be reset by changing
      !    the parameter command for KMAX below.
      !    changed by ww as variable input
      !
      ! rhs   == real vector of length n containing the right hand side.
      !          Destroyed on return.
      !
      ! sol   == real vector of length n containing an initial guess to the
      !          solution on input. approximate solution on output
      !
      ! eps   == tolerance for stopping criterion. process is stopped
      !          as soon as ( ||.|| is the euclidean norm):
      !          || current residual||/||initial residual|| <= eps
      !
      ! maxits== maximum number of iterations allowed
      !
      ! iout  == output unit number number for printing intermediate results
      !          if (iout <= 0) nothing is printed out.
      !
      !    Input, real AA(*), integer JA(*), IA(N+1), the matrix in CSR
      !    Compressed Sparse Row format.
      !
      !
      ! alu,jlu== A matrix stored in Modified Sparse Row format containing
      !           the L and U factors, as computed by routine ilut.
      !
      ! ju     == integer array of length n containing the pointers to
      !           the beginning of each row of U in alu, jlu as computed
      !           by routine ILUT.
      !
      ! on return:
      !
      ! sol   == contains an approximate solution (upon successful return).
      ! ierr  == integer. Error message with the following meaning.
      !          ierr = 0 --> successful return.
      !          ierr = 1 --> convergence not achieved in itmax iterations.
      !          ierr =-1 --> the initial guess seems to be the exact
      !                       solution (initial residual computed was zero)
      !
      ! work arrays:
      !
      ! vv    == work array of length  n x (im+1) (used to store the Arnoli
      !          basis)

cww   !  c,s  == work arrays of length im
cww   !  rs   == work arrays of length im+1


        implicit none

cww     integer, parameter :: kmax = 50
        integer n

        real*8 aa(*)
        real*8 alu(*)
cww     real*8 c(kmax)
        real*8 c(im)
        real*8 ddot
        real*8 eps
        real*8 eps1
        real*8, parameter :: epsmac = 1.0D-16
        real*8 gam
cww        real*8 hh(kmax+1,kmax)
        real*8 hh(im+1,im)
        integer i
        integer i1
        integer ia(n+1)
        integer ierr
        integer ii
        integer im
        integer iout
        integer its
        integer j
        integer ja(*)
        integer jj
        integer jlu(*)
        integer ju(*)
        integer k
        integer k1
        integer maxits
cww        integer n1
        real*8 rhs(n)
        real*8 ro
cww    real*8 rs(kmax+1)
        real*8 rs(im+1)
cww     real*8 s(kmax)
        real*8 s(im)
        real*8 sol(n)
        real*8 t
        real*8 vv(n,*)
      !
      !  Arnoldi size should not exceed KMAX=50 in this version. modified ww variable input
      !  To reset modify parameter KMAX accordingly.
      !


cww     n1 = n + 1
        its = 0
      !
      !  Outer loop starts here.
      !  Compute initial residual vector.
      !
        call ope ( n, sol, vv, aa, ja, ia )

        vv(1:n,1) = rhs(1:n) - vv(1:n,1)

20      continue

        ro = sqrt(ddot(n,vv,1,vv,1))

        if ( 0 < iout .and. its == 0 ) then
          write(iout, 199) its, ro
        end if

        if ( ro == 0.0D+00 ) then
          ierr = -1
          return
        end if

        t = 1.0D+00 / ro
        vv(1:n,1) = vv(1:n,1) * t

        if ( its == 0 ) then
          eps1 = eps * ro
          if(eps1.lt.epsmac) eps1=epsmac  ! added ww otherwise no convergence in last netwon step
        end if
      !
      !  Initialize first term of RHS of Hessenberg system.
      !
         rs(1) = ro
         i = 0

4       continue

         i = i + 1
         its = its + 1
         i1 = i + 1
         call lusol ( n, vv(1,i), rhs, alu, jlu, ju )
         call ope ( n, rhs, vv(1,i1), aa, ja, ia )
      !
      !  Modified Gram - Schmidt.
      !
         do j = 1, i
            t=ddot(n,vv(1,j),1,vv(1,i1),1)
            hh(j,i) = t
            call daxpty ( n, -t, vv(1,j), vv(1,i1))
         end do

         t = sqrt(ddot(n,vv(1,i1),1,vv(1,i1),1))
         hh(i1,i) = t

         if ( t /= 0.0D+00 ) then

           t = 1.0D+00 / t
           vv(1:n,i1) = vv(1:n,i1) * t

        end if
      !
      !  Update factorization of HH.
      !
        if ( i == 1 ) then
          go to 121
        end if
      !
      !  Perform previous transformations on I-th column of H.
      !
         do k = 2, i
           k1 = k-1
           t = hh(k1,i)
           hh(k1,i) = c(k1) * t + s(k1) * hh(k,i)
           hh(k,i) = -s(k1) * t + c(k1) * hh(k,i)
         end do

121     continue

        gam = sqrt ( hh(i,i)**2 + hh(i1,i)**2 )
      !
      !  If GAMMA is zero then any small value will do.
      !  It will affect only residual estimate.
      !
         if ( gam == 0.0D+00 ) then
           gam = epsmac
         end if
      !
      !  Get the next plane rotation.
      !
         c(i) = hh(i,i) / gam
         s(i) = hh(i1,i) / gam
         rs(i1) = -s(i) * rs(i)
         rs(i) = c(i) * rs(i)
      !
      !  Determine residual norm and test for convergence.
      !
         hh(i,i) = c(i) * hh(i,i) + s(i) * hh(i1,i)
         ro = abs ( rs(i1) )
cww131      format(1h ,2e14.4)

        if ( 0 < iout ) then
          write(iout, 199) its, ro
        end if

        if ( i < im .and. eps1 < ro ) then
          go to 4
        end if
      !
      !  Now compute solution.  First solve upper triangular system.
      !
         rs(i) = rs(i) / hh(i,i)

         do ii = 2, i
            k = i - ii + 1
            k1 = k + 1
            t = rs(k)
            do j = k1, i
               t = t - hh(k,j) * rs(j)
            end do
            rs(k) = t / hh(k,k)
        end do
      !
      !  Form linear combination of V(*,i)'s to get solution.
      !
         t = rs(1)
         rhs(1:n) = vv(1:n,1) * t

         do j = 2, i
            t = rs(j)
            rhs(1:n) = rhs(1:n) + t * vv(1:n,j)
        end do
      !
      !  Call preconditioner.
      !
         call lusol ( n, rhs, rhs, alu, jlu, ju )

         sol(1:n) = sol(1:n) + rhs(1:n)
      !
      !  Restart outer loop when necessary.
      !
         if ( ro <= eps1 ) then
           ierr = 0
           return
         end if

         if ( maxits < its ) then
           ierr = 1
           return
         end if
      !
      !  Else compute residual vector and continue.
      !
         do j = 1, i
           jj = i1 - j + 1
           rs(jj-1) = -s(jj-1) * rs(jj)
           rs(jj) = c(jj-1) * rs(jj)
         end do

         do j = 1, i1
            t = rs(j)
            if ( j == 1 ) then
              t = t - 1.0D+00
            end if
            call daxpty ( n, t, vv(1,j),  vv)
        end do

  199   format(' its =', i4, ' res. norm =', G14.6)
      !
      !  Restart outer loop.
      !
         go to 20
      end


      subroutine ope ( n, x, y, a, ja, ia )
      !*******************************************************************************
      !
      !! OPE sparse matrix * vector multiplication
      !
      !  Modified:
      !
      !    07 January 2004
      !
      !  Author:
      !
      !    Youcef Saad
      !
      !  Parameters:
      !
      !    Input, integer N, the order of the matrix.
      !
      !    Input, real X(N), the vector to be multiplied.
      !
      !    Output, real Y(N), the product A * X.
      !
      !    Input, real A(*), integer JA(*), IA(N+1), the matrix in CSR
      !    Compressed Sparse Row format.
      !
        implicit none

        integer n

        real*8 a(*)
        integer i
        integer ia(n+1)
        integer ja(*)
        integer k
        integer k1
        integer k2
        real*8 x(n)
        real*8 y(n)

        do i = 1, n
          k1 = ia(i)
          k2 = ia(i+1) -1
          y(i) = 0.0D+00
          do k = k1, k2
            y(i) = y(i) + a(k) * x(ja(k))
          end do
        end do

        return
      end


