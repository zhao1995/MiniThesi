      subroutine mkptr5(isymcsr)
c-----------------------------------------------------------------------
c
c     Purpose: Define arrays for PBCG solver
c              pointer with CSR storage in mkptr_csr
c     Inputs:
c
c     Output:
c
c-----------------------------------------------------------------------
      USE cdata
      USE ispcg
      USE soltyp
      USE doalloc
      implicit double precision (a-h,o-z)
      logical ldummy

      ldummy = .true.

c.....pcg
      if (ldummy) then 
        call ralloc(amcgz ,neq,'PBCG-mcgz ',ldummy)
        call ralloc(amcgzz,neq,'PBCG-mcgzz',ldummy)
        call ralloc(amcgr ,neq,'PBCG-mcgr ',ldummy)
        call ralloc(amcgrr,neq,'PBCG-mcgrr',ldummy)
        call ralloc(amcgp ,neq,'PBCG-mcgp ',ldummy)
        call ralloc(amcgpp,neq,'PBCG-mcgpp',ldummy)
        call ralloc(amcgx ,neq,'PBCG-mcgx ',ldummy)
      end if

c.... set iteration parameter, modify with [pbcg,iter]
      itercg  = ctis(1)
       tolcg  = ctis(2)
      itolcg  = ctis(3)
      isymcsr = ctis(4)

c.... default values
      if(itercg.eq.0   ) itercg = 150    ! No. of Iterations
      if( tolcg.eq.0.d0)  tolcg = 1.e-08 ! Tolerance
      if(itolcg.eq.0)    itolcg = 1      ! Type of Tolerance Norm
      if(isymcsr.ne.2)   isymcsr = 1     ! 1= symmetric matrix 2= unsymmetric

      return
      end


      subroutine dasol5 (a,b,ia,neq,energy)
c----------------------------------------------------------------------
c
c      Purpose: Solution of the problem Ax=b  for the PBCG solver
c
c      Inputs:
c         a(*)        - Non factored terms of A
c         b(*)        - right hand side vector b
c         ia(*)       - Pointer to row    CSR sparse technique
c         ja=csrja  - Pointer to column CSR
c         neq         - Number of equations
c
c      Outputs:
c         b(*)        - Solution vector x
c         energy      - Energy residual
c
c       Comment:
c         Coefficient matrix is not(!) decomposed into its triangular
c         factors using datri before using dasol.
c
c         Iteration parameters can be set by macro pbcg,iter
c
c----------------------------------------------------------------------
      USE fdata
      USE iofile
      USE iscsr
      USE ispcg
      USE isprec
      implicit double precision (a-h,o-z)
      dimension a(*),b(*),ia(*)
      real tary,tary1

      call etimef(tary)

c.... preconditioner
      call precond(a,ia,csrja,csrka)

      call pzero(amcgx,neq)

c.... here entry for guess of x (mcgx) is possible

c.... PBCG-algorithm
      call pbcg(a,ia,csrja,b,amcgx,amcgz,amcgzz,amcgr,
     +          amcgrr,amcgp,amcgpp,rmpcalu,impcjlu,impcju,
     +          neq,itercg,itolcg,
     +          tolcg,its,rof,pfr)

      if(pfr .and. ior.lt.0)  write(*  ,2000) its,rof
      if(pfr)                 write(iow,2000) its,rof

c.... energy
      energy = ddot(neq,b,1,amcgx,1)

c.... copy x->b for further FEAP processing
      call pmove(amcgx,b,neq)

cww   call etimef(tary1)
cww   if(pfr)                write(iow,2001) tary1-tary,
cww   if(pfr .and. ior.lt.0) write(*  ,2001) tary1-tary

      return

2000  format(3x,'PBCG: Number of iterations',i4,' Norm(Ku-b)  ',g15.5)
c2001  format(3x,'Time for PBCG-Solution',40x,'t=',0pf9.4)
      end
c


      subroutine pbcg(a,ia,ja,b,x,z,zq,r,rq,p,pq,alu,jlu,ju,
     +               neq,niter,itol,tol,its,rnorm,pfr)
c-----------------------------------------------------------------------
c
c      Purpose: solve a linear set of equations  Ku=f with the
c               PBCG-method (preconditioned-Bi-Conjugate-Gradient method
c               Ref: see e.g. Numerical Recipes: SR linbcg
c
c      Inputs:
c      b     - RHS (neq) f
c      x     - solution(neq) u  initial guess
c      a     - matrix K
c      ia    - Pointer array to row      of K
c      ja    - Pointer array to column   of K
c      z     - iteration vector (neq):  M  *z =r , K  *p  = z
c      zq    - iteration vector (neq),  M^T*zq=rq, K^T*pq = zq
c      r     - iteration vector residual: r = K  *x-f (neq)
c      rq    - iteration vector residual: rq= K^T*x-f (neq)
c      p     - direction vector for update of x (neq)
c      pq    - direction vector
c      alu   - Matrix ALU in MSR storage
c      jlu   - pointer to non-diagonals of ALU
c      ju    - pointer to column of ALU
c      niter - max. iterations
c      itol  - typ of norm
c      tol   - tolerance for convergence
c
c      Outputs:
c      x     - Solution vector x=u
c
c      Open:
c
c      test of tolerance criteria
c
c     (c) WW 5/2006
c-----------------------------------------------------------------------
      USE iofile
      USE iscsr
      USE pdata8
      USE soltyp
      implicit double precision(a-h,o-z)
      logical fa,tr,pfr
      dimension b(*),ia(*),ja(*),a(*),x(*),alu(*),jlu(*),ju(*),
     1          z(*),zq(*),r(*),rq(*),p(*),pq(*),ait(1)
      data fa,tr/.false.,.true./
c
      eps=1.d-14
      itmax = niter

c.... initial residual, if x .ne.0    K*u = r
      call pzero(r,neq)
      call promul_csr(a,x,r,ia,ja,neq,.true.,isymcsr)

c.... b --> r
      do j = 1,neq
         r(j) = b(j)-r(j)
        rq(j) = r(j)
      end do

c.... set norms
      if(itol.eq.1) then
        bnrm = sqrt(ddot(neq,b,1,b,1))
        call lusol (neq,r,z,alu,jlu,ju ) ! M*z = r_0

      else if(itol.eq.2) then
        call lusol (neq,b,z,alu,jlu,ju ) ! M*z = b
        bnrm = sqrt(ddot(neq,z,1,z,1))
        call lusol (neq,r,z,alu,jlu,ju ) ! M*z = r_0

      else if (itol.eq.3.or.itol.eq.4) then
        call lusol (neq,b,z,alu,jlu,ju ) ! M*z = b
        if(itol.eq.3) bnrm = sqrt(ddot(neq,z,1,z,1))
        if(itol.eq.4) bnrm = valmax(z,neq)
        call lusol (neq,r,z,alu,jlu,ju ) ! M*z = r_0
        if(itol.eq.3) znrm = sqrt(ddot(neq,z,1,z,1))
        if(itol.eq.4) znrm = valmax(z,neq)
      end if

c.... PBCG iteration - main loop
      i = 0

100   i = i+1

c.....Diagonal pre-conditioning with M^T (!)
      call lutsol (neq,rq,zq,alu,jlu,ju ) ! M^T*zq = rq

c.... Coefficient beta_k and direction vectors p,pq
      bknum = ddot(neq,z,1,rq,1)

      if (i.eq.1) then
        call pmove( z, p,neq)
        call pmove(zq,pq,neq)
      else
        bk = bknum/bkden
        do j = 1,neq
           p(j) =  z(j) + bk* p(j)
          pq(j) = zq(j) + bk*pq(j)
        end do
      end if
      bkden = bknum

c.... coefficient alpha_k, new iterate u and new residuals r, rq
      call pzero(z,neq)
      call promul_csr(a,p,z,ia,ja,neq,.true.,isymcsr)

      akden = ddot(neq,z,1,pq,1)
      ak    = bknum/akden

      if(isymcsr.eq.2) then
        call atmux (neq,pq,zq,a,ja,ia)                    ! K^T*pq = zq
      else
        call pzero(zq,neq)
        call promul_csr(a,pq,zq,ia,ja,neq,.true.,isymcsr) !  K^T*pq = zq  sym K^T = K!
      end if


      do j = 1,neq
         x(j) =  x(j) + ak* p(j)
         r(j) =  r(j) - ak* z(j)
        rq(j) = rq(j) - ak*zq(j)
      end do

c.... solve M*z=r and check stopping criterion
      call lusol (neq,r,z,alu,jlu,ju ) ! M*z = r

      if(itol.eq.1) then
         rnorm = sqrt(ddot(neq,r,1,r,1))/bnrm
      else if(itol.eq.2) then
         rnorm = sqrt(ddot(neq,z,1,z,1))/bnrm
      else if(itol.eq.3.or.itol.eq.4) then
        zm1nrm = znrm
        if(itol.eq.3) znrm = sqrt(ddot(neq,z,1,z,1))
        if(itol.eq.4) znrm = valmax(z,neq)
        if(abs(zm1nrm-znrm).gt.eps*znrm) then
          if(itol.eq.3) wnrm = sqrt(ddot(neq,p,1,p,1))
          if(itol.eq.4) wnrm = valmax(p,neq)
          dxnrm=abs(ak)*wnrm
          rnorm=znrm/abs(zm1nrm-znrm)*dxnrm
        else
          rnorm=znrm/bnrm
          goto 101
        end if
        if(itol.eq.3) xnrm = sqrt(ddot(neq,x,1,x,1))
        if(itol.eq.4) xnrm = valmax(x,neq)
        if(rnorm.le.0.5d0*xnrm) then
          rnorm=rnorm/xnrm
        else
          rnorm=znrm/bnrm
          goto 101
        end if
      end if

101   write(*  ,2000) i,rnorm  ! for testing

      if(rnorm.gt.tol.and.i.le.itmax) go to 100
      if(rnorm.gt.tol.and.i.gt.itmax) then
c....   more equations only for interactive mode
        if(pfr)                write(iow,2001) i,rnorm
        if(pfr .and. ior.lt.0) write(*  ,2001) i,rnorm
        if(iplot.ne.0) then
          write(*  ,2002)
          call dinput(ait,1)
          miter = ait(1)
          if (miter.eq.0) then
            go to 102
          else
            itmax = itmax + miter
            goto 100
          end if
         end if
      end if
c.... solution
102   its  = i

      return

2000  format(3x,'PBCG: Iteration',i4,'  Norm ',g15.5)
2001  format(3x,'PBCG: Max. Numb. of It. exeeded',i4,'  Norm ',g15.5)
2002  format(3x,'Give number of additional iterations (0 to stop) : ',$)

      end

c---------------------------------------------------------------------
      double precision function valmax(a,n)
c--------------------------------------------------------------
c
c     Purpose: find abs maximum of array
c
c--------------------------------------------------------------
      implicit double precision (a-h,o-z)
      dimension a(n)
       imax=1
      do i=1,n
        if(abs(a(i)).gt.abs(a(imax))) imax=i
      end do
      valmax = abs(a(imax))
      return
      end
c

