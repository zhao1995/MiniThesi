      subroutine feast(a,b,ia,ja,ib,jb,e,x,res,neq,imas,emin,emax,m0,
     +                 info,iorf,iow,prt)
c ----------------------------------------------------------------------
c.... Purpose: interface to eigenvalue solver FEAST                                     
c
c....  Inputs:
c      a          - matrix A
c      b          - matrix B    
c      ia(*)      - Pointer to row    CSR sparse technique
c      ja(*)      - Pointer to column CSR=csrja
c      ib(*)      - Pointer to row    CSR only for lumped B=B_l or B=1   
c      jb(*)      - Pointer to column CSR only for lumped B=B_l or B=1
c      neq        - size of problem
c      imas       - type of problem  1=general Ax=Bx 2=special Ax=lx
c      nits       - no of iterations
c      m0         - initial guess for subspace dimension to be used
c      iorf       - input  channel
c      iow        - output channel 
c      prt        - Flag, output iteration arrays if true
c           
c      ifpm(128)        - pass various parameters to FEAST 
c
c      Outputs:
c      epsout     - relative error on the trace
c      loop       - number of refinement loops executed
c      m1         - eigenvalues found
c      e(m0)      - array of eigenvalues   with m1 entries
c      x(neq,m0)  - array of eigenvectors  with m1 entries
c      res(m0)    - residuals              with m1 entries
c      info       - result of execution

c     Ref.: Eric Polizzi: 
c           Density-matrix-based algorithm for solving eigenvalue problems
c           PHYSICAL REVIEW B 79, 115112_1-6 2009
c           http://www.ecs.umass.edu/~polizzi/feast/
c
c     (c) WW IBS KIT 08/2013
c
c     Comments:
c     -integrated MKL-version INTEL
c     -special problem not used, only valid for Ax=1x   
c     -Ax=Bx with B=B_l solved by general problem!  
c     -Ax=1x with 1     solved by general problem!  
c     -B should be positive definite 
c      
c     OPEN
c     m=>=1.5*m1, but m1 is not known
c
c ----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)

      include 'mkl.fi' ! INTEL-MKL

      logical prt
      dimension a(*),b(*),e(m0),x(neq,m0),res(m0) ! with length m0>m
      dimension ifpm(128),ia(*),ja(*),ib(*),jb(*)
      character*1 uplo

      if(m0.lt.2) then 
        write(*,*) '2 values for gerschgorin circles necessary'
        info = 1
        return
      end if

c.... Abschaetzung der Eigenwerte von A -Gerschgorin-Kreise
      call asumr (x,neq,a,ia,ja,m0,emin1,emax1) 

c.... set values 
      if(emin.eq.0.d0.and.emax.eq.0.d0) then
        emin=emin1
        emax=emax1
      end if
      
c.... number of fixpoints within limit values
      nfix=0 
      do i=1,neq
        if(x(i,1).ge.emin.and.x(i,1).le.emax) nfix=nfix+1
      end do

      write(*,*) 'nfix,emin,emax',nfix,emin,emax

      
c...  solv EV-Problem
      call pzeroi(ifpm,128)

      call pzero(e,m0) 
      call pzero(res,m0) 
      call pzero(x,m0*neq)
       
      nits=25
      m01=m0 ! m0 wird in dfeast auf Zahl der gefundenen EV gesetzt.
      
c...  set ifpm         
      ifpm( 1)=1    ! print runtime status to the screen. 
      ifpm( 2)=8    ! number of contour points Ne = 8
      ifpm( 4)=nits ! Maximum number of refinement loops allowed.def=25
      ifpm( 5)=1    ! user supplied initial subspace is used. 
      ifpm(64)=1    ! use Pardiso with definitions defined in pardisoinit ! only for Pardiso
      
      call feastinit(ifpm)

c.... type of matrix
c     uplo       - 'U' upper triangular parts of A and B
c                  'L' lower triangular parts of A and B
c                  'F'             full matrices A and B 
      uplo='U'

      if(imas.eq.1) then 
c....   FEAST general Ax=lBx,    A symm B symm   
        call dfeast_scsrgv(uplo,neq,a,ia,ja,b,ia,ja,ifpm, 
     +       epsout,loop,emin,emax,m01,e,x,m1,res,info)
      else if(imas.eq.2) then !  Ax=l*M_Lx
c....   FEAST special Ax=lBx,    A symm B diagonal   
        call dfeast_scsrgv(uplo,neq,a,ia,ja,b,ib,jb,ifpm, 
     +       epsout,loop,emin,emax,m01,e,x,m1,res,info)
c....   FEAST special Ax=lx,    A symm B=1    
c        call dfeast_scsrev(uplo,neq,a,ia,ja,ifpm, 
c     +       epsout,loop,emin,emax,m01,e,x,m1,res,info)
      end if

      if(prt) then
c....   info
        if(iorf.lt.0) write(*,*)   'info,loop,m1', info,loop,m1
                      write(iow,*) 'info,loop,m1', info,loop,m1

c....   eigenvalues+residuals

        if(iorf.lt.0) write(  *,*) 'i eigenvalue(i) residual(i)'
                      write(iow,*) 'i eigenvalue(i) residual(i)'

        do i=1,m0
          if(iorf.lt.0)write(  *,*) i,e(i),res(i)
                       write(iow,*) i,e(i),res(i)
        end do  
      end if 
  
c.... scale the vectors to have maximum element of 1.0
      do n = 1,m01
        if(n.le.m01) call scalev(x(1,n),neq)
      end do

c.... errormessages
      call errofeast(info) 

      return
      end
c
      subroutine errofeast(ierror)
c-----------------------------------------------------------------------
c
c.... Purpose: print errors FEAST EVsolver
c
c-----------------------------------------------------------------------
      USE iofile
      implicit double precision (a-h,o-z)
      character*63 errorind(12)
      character*8  erromess(12)
      data errorind /
     1'Successful exit                                                ', ! 0  
     3'Warning No eigenvalue found in the search interval.            ', ! 1  
     3'Warning No Convergence (no of iterat. loops >fpm(4)).          ', ! 2  
     4'Warning Size of the subspace m0 is too small (m0<m).           ', ! 3  
     5'Problem with emin,emax (emin>=emax).                           ', ! 200  
     6'Problem with size of initial subspace m0 (m0<=0).              ', ! 201  
     7'Problem with size of the system n (n<=0).                      ', ! 202  
     8'Internal error for allocation memory.                          ', ! -1  
     9'Internal error: not enough memory inner system solver.         ', ! -2  
     +'Internal error: matrix B may not be positive definite.         ', ! -3  
     1'Problem with -(100+i)-th argument of EVsolver interface.       ', ! -(100+i)
     2'Problem with(100+i)i-th value of input parameter (fpm(i)).     '/ !  (100+i) 

      data erromess /
     1'Message ', 
     2'Warning ', 
     3'Warning ', 
     4'Warning ', 
     5'Error   ', 
     6'Error   ', 
     7'Error   ', 
     8'Error   ', 
     9'Error   ', 
     +'Error   ', 
     1'Error   ', 
     2'Error   '/ 
     
      ind=0 
      if(ierror.eq.  0) ind=1 
      if(ierror.eq.  1) ind=2 
      if(ierror.eq.  2) ind=3 
      if(ierror.eq.  3) ind=4 
      if(ierror.eq.200) ind=5 
      if(ierror.eq.201) ind=6 
      if(ierror.eq.202) ind=7 
      if(ierror.eq. -1) ind=8 
      if(ierror.eq. -2) ind=9 
      if(ierror.eq. -3) ind=10 
      if(ind.eq.0) then
       if(ierror.gt.0)   ind=11
       if(ierror.lt.0)   ind=12 
      end if

      if(ind.eq.0) stop 'ErroFEAST'

c.... message
                   write(iow,1000) erromess(ind),ierror,errorind(ind)
      if(ior.lt.0) write(  *,1000) erromess(ind),ierror,errorind(ind)

1000  format(a,i8,' in EVsolver FEAST',/,a)

      return
      end
c
c--------------------------------------------------------------
c
      subroutine asumr (x,n,a,ia,ja,m0,emin,emax) 
c-----------------------------------------------------------------------
c     Purpose: Calculates 
c              1)D:   diagonal of A in x(n,1)
c              2)R:   R_i = sum_i=1^n |a_ij| i.ne.j of A in x(n,2)  
c              3)Emin=D-R   
c              4)Emax=D+R   
c
c              Matrix A is symmetric and stored in CSR format. 
c
c     Inputs:
c      n       - row dimension of A(ia,ja)
c      a(ia,ja)- matrix A in CSR format.
c      m0      - no of eigenvalues 
c
c     Outputs:
c      x(n,*)  - real array 
c      emin    - min. value of EV-range  
c      emax    - max. value of EV-range 
c
c     Theory: Gerschgorin circles 
c
c
c-----------------------------------------------------------------------

      implicit double precision (a-h,o-z) 
      real*8  x(n,*), a(*)
      integer ia(*),ja(*)

c.... loop over all rows
      do i=1,n,1

c....   diagonal part = Fixpoint        
        x(i,1) = a(ia(i))

c....   radius. r_i = sum_i=1^n |a_ij| i.ne.j 
c....   loop over all entries in row, except first=diagonal 
        do j=ia(i)+1,ia(i+1)-1 
c....     column index
          k=ja(j)
c....     add values
          x(i,2) = x(i,2) + abs(a(j))
          x(k,2) = x(k,2) + abs(a(j))
        end do  
      end do   

c.... limits Emin,Emax

c      write(*,*) 'first'
c      do i=1,n
c        write(*,100) 'i,a_ii,r,mim,max',i,
c     +  x(i,1),x(i,2),x(i,1)-x(i,2),x(i,1)+x(i,2)
c      end do

c.... sort in ascending order
      call sortfeast(x,n)

c      write(*,*) 'second'
c      do i=1,n
c        write(*,100) 'i,a_ii,r,mim,max',i,
c     +  x(i,1),x(i,2),x(i,1)-x(i,2),x(i,1)+x(i,2)
c      end do

c.... limit values  with emini>=0
      emin = 0.d0 
      emax = 0.d0 
      do i = 1,m0
        emini=x(i,1)-x(i,2)
        emaxi=x(i,1)+x(i,2)
        if(emini.lt.0.d0) emini=0.d0
        emin=min(emin,emini)
        emax=max(emax,emaxi)
      end do

100   format(a,i4,4e14.7)

      return
      end
c

      subroutine sortfeast(x,n)
c-----------------------------------------------------------------------
c
c     Purpose:  Sorts the entries in array x(n,2) 
c               into increasing order with respect to x(n,1).
c
c     Inputs:
c        x(n,2) - Array of integers to sort
c        n      - Number of entries to sort
c
c     Output:
c        x(n,2) - Sorted array
c
c-----------------------------------------------------------------------
c
      implicit double precision (a-h,o-z) 

      logical   flag
      real*8    x(n,*),xh(2)

c     Sort integer array using bubble sort algorithm

      flag = .true.
      do while( flag )
        flag = .false.
        do i = 1, n-1
          if (x(i,1).gt.x(i+1,1)) then
            flag   = .true.
c....       save i
            xh(1)  = x(i,1)
            xh(2)  = x(i,2)
c....       set i              
            x(i,1)   = x(i+1,1)
            x(i,2)   = x(i+1,2)
c....       set i+1              
            x(i+1,1) = xh(1)
            x(i+1,2) = xh(2)
          end if
        end do ! i
      end do ! while flag = true
      end
c
