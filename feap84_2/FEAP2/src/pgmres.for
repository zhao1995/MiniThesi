      subroutine mkptr6(isymcsr)
c-----------------------------------------------------------------------
c
c     Purpose: Define arrays for PGMRES solver 
c              pointer with CSR storage in mkptr_csr 
c
c     Inputs:
c
c     Output:
c              
c-----------------------------------------------------------------------
      USE cdata
      USE isgmr
      USE soltyp
      USE doalloc
      implicit double precision (a-h,o-z)
      logical ldummy

      ldummy = .true. 
      
c.... set iteration parameter, modify with [gmr,iter] 
      itergmr = ctis(1)
       tolgmr = ctis(2)
          im  = ctis(3)
      isymcsr = ctis(4)

c.... default values
      if(itergmr.eq.0 )   itergmr = 150    ! No. of Iterations   
      if(tolgmr.eq.0.d0)  tolgmr  = 1.e-08 ! Tolerance 
      if(im.eq.0)         im      = 50     ! Size of krylov subspace, typical 50, arnoldi size, LÖhnert: 100
      if(isymcsr.ne.2)    isymcsr = 1      ! 1= symmetric matrix 2= unsymmetric 

c.... gmres
      if (ldummy) then
        call ralloc( rmgmrx,neq,       'GMRES- mgmrx',ldummy)
        call ralloc(rmgmrss,(im+1)*neq,'GMRES-mgmrss',ldummy)
        call ralloc(rmgmrhh,(im+1)*im, 'GMRES-mgmrhh',ldummy)
        call ralloc(rmgmrrs,(im+1),    'GMRES-mgmrrs',ldummy)
        call ralloc( rmgmrc,im,        'GMRES- mgmrc',ldummy)
        call ralloc( rmgmrs,im,        'GMRES- mgmrs',ldummy)
      end if

      return
      end  
c
      subroutine dasol6 (a,b,ia,neq,energy)
c----------------------------------------------------------------------
c
c      Purpose: Solution of the problem Ax=b  for the PGMRES solver
c
c      Inputs:
c         a(*)       - Non factored terms of A
c         b(*)       - right hand side vector b
c         ia(*)      - Pointer to row    CSR sparse technique
c         ja=csrja - Pointer to column CSR
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
c----------------------------------------------------------------------
      USE fdata
      USE iofile
      USE iscsr
      USE isgmr
      USE isprec
      implicit double precision (a-h,o-z)
      dimension a(*),b(*),ia(*)
      real tary,tary1

      call etimef(tary)

c.... preconditioner 
      call precond(a,ia,csrja,csrka)

      call pzero(rmgmrx,neq)    
c.... here entry for guess of x (mgmrx) is possible      
  
c.... GMRES-algorithm
      call gmres(a,ia,csrja,b,rmgmrx,rmgmrss,rmgmrhh,
     +     rmgmrrs,rmgmrc,rmgmrs,rmpcalu,impcjlu,impcju,
     +     neq,im,tolgmr,itergmr,its,rof,pfr)


      if(pfr .and. ior.lt.0)  write(*  ,2000) its,rof
      if(pfr)                 write(iow,2000) its,rof

c.... energy
      energy = ddot(neq,b,1,rmgmrx,1)

c.... copy x->b for further FEAP processing     
      call pmove(rmgmrx,b,neq)

cww   call etimef(tary1)
cww   if(pfr)                write(iow,2001) tary1-tary,
cww   if(pfr .and. ior.lt.0) write(*  ,2001) tary1-tary

      return

2000  format(3x,'PGMRES: Number of iterations',i4,' Norm(Ku-b)',g15.5)
c2001  format(3x,'Time for PGMRES-Solution',38x,'t=',0pf9.4)

      end
c
       subroutine gmres (a,ia,ja,b,x,ss,hh,rs,c,s,alu,jlu,ju,
     +                   neq,im,tol,niter,its,rof,pfr)
c-----------------------------------------------------------------------
c
c      Purpose: solve a linear set of equations  Ax=b with the
c               generalized minimum residual method (GMRES) including
c               preconditioning
c               
c
c      Inputs: 
c      a(neq )      - terms of K
c      b(neq)       - RHS
c      ia(*)        - Pointer array to row      of K              
c      ja(*)        - Pointer array to column   of K
c      alu(*)       - Matrix ALU in MSR storage
c      jlu(*)       - pointer to non-diagonals of ALU
c      ju(*)        - pointer to column of ALU 
c      x(neq)       - 0 on input, approximate solution on output
c      ss(neq,im+1) - work space of size n x (im+1)
c      hh (im+1,im) - Hessenberg matrix
c      neq          - size of problem
c      im           - size of krylov subspace, typical 50, arnoldi size
c      tol          - tolerance for stopping criterion
c                     ||current residual||/||initial residual|| <= tol
c      niter        - max. iterations
c      pfr          - print option
c
c      Outputs:
c      x            - Solution vector x
c      its          - Number of iterations used
c      rof          - residuum k*x-b
c
c
c      Comments:    - comments in notation of Löhnert (Diss F04/2 Hannover) added
c                   - preconditioner added K_ii see Löhnert
c                     (c) based on version 05,23,1985, Y. Saad

c      Open:        - convergence check Löhnert is different (on global residual!)
c                   - im = 50 Änderungs option
c                   - Neustart mit x_0 .ne.0 
c                   - Householder instead of Gram-Schmidt
c                   - Preconditioner
c
c     (c) ww 5/2006 based on a code of SAAD
c-----------------------------------------------------------------------
      USE iofile
      USE iscsr
      USE pdata8
      implicit double precision(a-h,o-z)

      logical pfr

      dimension ss(neq,*),x(*),hh(im+1,im),c(im),s(im),rs(im+1)
      dimension a(*),ia(*),ja(*),b(*),ait(1),alu(*),jlu(*),ju(*) 
      dimension wr(im+1),wi(im+1),hh1(im+1,im)


      data tolmac/1.d-16/
      its = 0
      itmax = niter 

      call pzero(ss,neq*(im+1))
      call pzero(hh,im*(im+1))
      call pzero(c,im)
      call pzero(s,im)
      call pzero(rs,im+1)

c.... outer loop
10    continue
      call pzero(hh1,im*im)

      
c.... compute initial residual vector
c.... r_0 = p - K u_0
      call pzero(ss(1,1),neq)                     
      call promul_csr(a,x,ss(1,1),ia,ja,neq,.true.,isymcsr)

      do j=1,neq
        ss(j,1) = b(j) - ss(j,1)
      end do


c.....preconditioner M*c = r_0 
      call lusol (neq,ss(1,1),ss(1,1),alu,jlu,ju )


c.... beta
      ro = dsqrt(ddot(neq,ss(1,1),1,ss(1,1),1))
      if (ro.eq.0.0d0) return
        
c.... v_1 = c/beta         
      do j=1, neq
        ss(j,1) = ss(j,1)/ro
      end do

      if (its .eq. 0) tol1=tol*ro

      if(tol1.lt.tolmac) tol1=tolmac  ! added ww otherwise no convergence in last step  

      if(pfr)                write(iow,2000) its,ro
      if(pfr .and. ior.lt.0) write(*  ,2000) its,ro

c.... initialize 1-st term  of rhs of hessenberg system..
      rs(1) = ro
      i     = 0
20    i     = i+1
      its   = its + 1
      i1    = i + 1

c.... w_j = K*v_j
      call pzero(ss(1,i1),neq)                         
      call promul_csr(a,ss(1,i),ss(1,i1),ia,ja,neq,.true.,isymcsr)

c.....preconditioner M*w_j = K*v_j 
      call lusol (neq,ss(1,i1),ss(1,i1),alu,jlu,ju )

c.... modified Gram-Schmidt orthogonalization
      do j=1, i

c....   h_ij = w_j * v_i
        t = ddot(neq,ss(1,j),1,ss(1,i1),1)
        hh(j,i)  = t
        hh1(j,i) = t

c....   w_j = w_j-h_ij*v_i
        call daxpty(neq,-t,ss(1,j),ss(1,i1))
      end do   

c.... h_j+1,j=||w_j||  
      t = dsqrt(ddot(neq,ss(1,i1),1,ss(1,i1),1))
      hh(i1,i)  = t
      hh1(i1,i) = t
      if (t.eq.0.0d0) go to 30

c.... v_j+1 = w_j/h_j+1,j
      t = 1.0d0 / t
      do k=1,neq
        ss(k,i1) = ss(k,i1)*t
      end do

c.... done with modified Gram-Schmidt and Arnoldi step

c.... now  update factorization of hh
30    if (i .eq. 1) goto 40
c.... perform previous transformations  on i-th column of h
      do k=2,i
        k1 = k-1
        t = hh(k1,i)
        hh(k1,i) = c(k1)*t + s(k1)*hh(k,i)
        hh(k,i) = -s(k1)*t + c(k1)*hh(k,i)
      end do     

40    gam = dsqrt(hh(i,i)**2 + hh(i1,i)**2)
      if (gam.eq.0.0d0) gam = tolmac

c.... determine next plane rotation
      c(i)   = hh(i,i)/gam
      s(i)   = hh(i1,i)/gam
      rs(i1) = -s(i)*rs(i)
      rs(i)  =  c(i)*rs(i)

c...  determine residual norm and test for convergence
      hh(i,i)  = c(i)*hh(i,i) + s(i)*hh(i1,i)

      ro = dabs(rs(i1))
      if(pfr)                write(iow,2000) its,ro
      if(pfr .and. ior.lt.0) write(*  ,2000) its,ro

      if (i .lt. im .and. (ro .gt. tol1))  go to 20
     
c.... now compute solution. first solve upper triangular system.

c.....y_m      

      rs(i) = rs(i)/hh(i,i)
      do ii=2,i
        k=i-ii+1
        k1 = k+1
        t=rs(k)
        do j=k1,i
          t = t-hh(k,j)*rs(j)
        end do
        rs(k) = t/hh(k,k)
      end do    

c.... done with back substitution

c.... now form linear combination to get solution
c.... u_m = u_0 + V_m *y_m
      do j=1, i
        t = rs(j)
        call daxpty(neq,t,ss(1,j),x)
      end do    

c.... differs from Löhnert!!
c.... restart outer loop  when necessary
      if (ro .gt. tol1 .and. its .le. itmax) goto 10
      if (ro .le. tol1 .and. its .le. itmax) goto 50
             
c.... more equations only for interactive mode
cww      if(iplot.ne.0) then
cww        write(*  ,2000) its,ro
cww        write(*  ,2001) 
cww        call dinput(ait,1)
cww        miter = ait(1)
cww        if (miter.eq.0) then
cww          goto 50 
cww        else 
cww          itmax = itmax + miter
cww          goto 10
cww        end if
cww     end if 

      write(iow,2004) itmax
      write(*  ,2004) itmax

c.... determine final residual norm : test with matrix of problem
50    call pzero(ss(1,1),neq)                      
      call promul_csr(a,x,ss(1,1),ia,ja,neq,.true.,isymcsr)
 
      do i=1,neq
        ss(i,1) = ss(i,1)-b(i)
      end do    

      rof =  dsqrt( ddot(neq,ss(1,1),1,ss(1,1),1))

c...  calculate sign of diagonals, not really correct in all situations, dependes on precondition 
cfg      nn = ii-1
cfg      call hqr(im+1,nn,1,nn,hh1,wr,wi,ierr)
cfg      do i = 1,nn
cfg        if(wr(i).lt.0.and.wi(i).ge.0) write(  *,2002)
cfg        if(wr(i).lt.0.and.wi(i).lt.0) write(  *,2003) 
cfg      end do 
      return
      
2000  format(3x,'PGMRES: Iteration',i4,' Norm      ',g15.5)
2001  format(3x,'Give no. of additional iterations outer loop ',
     +          '(0 to stop) : ',$)
2002  format('**PGMRES warning: sign of diagonal changed   ')
2003  format('**PGMRES warning: sign of diagonal changed ? ')
2004  format('**PGMRES warning: no convergence after ',i5,' iterations')

      end
