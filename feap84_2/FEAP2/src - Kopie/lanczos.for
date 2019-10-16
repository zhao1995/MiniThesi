      subroutine lanczos(a,b,v,d,g,dp,dtol,p,t1,t2,t3,dh,gh,vh,jp,nf,nv,
     +                   neq,imas,shift,tol,prt,its)
c----------------------------------------------------------------------
c
c      Purpose: Lanczos iteration to extract lowest nf eigenpairs
c               Solves:  (A - shift*B)*V = B*V*d for V and d
c               Theory see e.g. Book of Bathe chap.11
c
c      Inputs:
c         a(*)      - Coefficient matrix (tangent stiffness)
c         b(*)      - Coefficient matrix (mass or geometric stiffness)
c         jp(*)     - Pointer array for row/columns of tangent
c         nf        - Number of pairs to converge
c         nv        - Size of work space
c         neq       - Size of A,B
c         imas      - Switch: =1 for consistent B; =2 for diagonal B
c         shift     - Value of shift
c         tol       - Tolerance to converge pairs
c         prt       - Flag, output iteration arrays if true
c         its       - Maximum number of lnczos iterations = 1
c
c      Scratch:
c         g(*)      - Subdiagonal elements in tridiagonal matrix T
c         dp(*)     - Previous iterated eigenvalues
c         dtol(*)   - Tolerance of eigenvalue iterations
c         p(nv,*)   - Eigenvectors of T*P = 1*P*d
c         t1(neq)   - Working vector
c         t2(neq)   - Working vector
c         t3(neq)   - Working vector
c         dh(nv)    - Working vector
c         gh(nv)    - Working vector
c         vh(nv)    - Working vector
c
c      Outputs:
c         v(neq,*)  - Eigenvectors
c         d(*)      - Eigenvalues
c
c----------------------------------------------------------------------
      USE evdata
      USE iofile
      USE iscsr
      USE pdata8
      USE soltyp
      implicit double precision (a-h,o-z)
      logical conv,prt
      dimension a(*),b(*),v(neq,*),d(*),g(*),dp(*),dtol(*),p(nv,nv),
     1          t1(*),t2(*),t3(*),dh(*),gh(*),vh(*),jp(*),ait(1)

c.... compute the initial iteration vectors
      call pzero(v,nv*neq)
      num  = 0
cww   nmas = 0
c.... on diagonal
c.... count number of values less than the current shift
c.... count the number of nonzero masses                    ! not used !
      if (istyp.eq.0) then  ! Standard solver
        do  n = 1,neq
          if(a(n).lt.0.0d0) num  = num  + 1
cww       if(b(n).ne.0.0d0) nmas = nmas + 1
        end do

      else if(istyp.eq.1.or.istyp.eq.2) then ! SM
        do n = 1,neq
          if( a(jp(n)).lt.0.0d0) num  = num  + 1
cww       if(imas.eq.1) then
cww         if( b(jp(n)).ne.0.0d0) nmas = nmas + 1
cww       else
cww         if(b(n).ne.0.0d0) nmas = nmas + 1
cww       end if
        end do

      else if(istyp.ge.3.and.istyp.le.8) then ! all other CSR solver
        call dnzero_csr(neq,a,csrka, num,2,1)

      end if
                   write(iow,2002) num
      if(ior.lt.0) write(  *,2002) num
c
      do n = 1,neq
        if (istyp.eq.0) then ! Standard solver
          dm = b(n)

        else if (istyp.eq.1.or.istyp.eq.2) then ! SM
          if(imas.eq.1) then
            dm = b(jp(n))
          else
            dm = b(n)
          end if

        else if (istyp.ge.3.and.istyp.le.8) then ! all other CSR solver
          call pick_csr(dm,n,b,csrka,imas)

        end if

        if(dm.ne.0.0d0) then
          call random_number(rnum)
          v(n,1) = rnum
        end if
      end do
      do i = 1,nv
        dtol(i) = 1.d0
      end do

      conv = .false.
      itlim = its
      itsta = 1

301   continue

      do 300 it = itsta, itlim
        call pzero(dh,nv)
        call pzero(gh,nv)
c
c....   normalize starting vector
        call pzero(t2,neq)
        call promul(b(neq+1),b(neq+1),b(1),v(1,1),t2,jp,neq,imas,1)

        dm   = 1.d0
        do n = 1,neq
c
          if (istyp.eq.0) then ! Standard Solver
            dm = b(n)

          else if (istyp.eq.1.or.istyp.eq.2) then ! SM
            if(imas.eq.1) then
              dm = b(jp(n))
            else
              dm = b(n)
            end if

          else if (istyp.ge.3.and.istyp.le.8) then ! all other CSR solver
            call pick_csr(dm,n,b,csrka,imas)

          end if
c
         if(dm.eq.0.0d0) t2(n) = 0.d0
       end do

        fact = ddot(neq,v(1,1),1,t2,1)
        if(fact.lt.0.d0) then
          write (*,*) 'Negative scaling factor in LANCZOS! Use SUBSPACE'
          return
        end if
        fact = 1.d0/dsqrt(fact)
        do n = 1,neq
          v(n,1) = fact*v(n,1)
        end do
c
c....   alpha_1 = q_1 * B * A^-1 * B * q_1
        call pzero(t2,neq)
        call promul(b(neq+1),b(neq+1),b(1),v(1,1),t2,jp,neq,imas,1)
        call pmove(t2,t1,neq)
        call dasol(a(neq+1),a(neq+1),a(1),t1,jp,neq,dummy)
        dh(1) = ddot(neq,t1,1,t2,1)  ! alpha
        gh(1) = 0.d0

c....   Lanczos iteration
        k  = 0
100     k  = k+1
        k1 = k + 1

        call pmove(dh,d,k)
        call pmove(gh,g,k)

c....   r_k = A^-1 * B * q_k - alpha_k*q_k - beta_k-1 * q_k-1
        call daxpty(neq,-d(k),v(1,k),t1)
        if(k.gt.1)call daxpty(neq,-g(k),v(1,k-1),t1)
        if(it.ne.1.and.dsqrt(ddot(neq,t1,1,t1,1)).lt.1.d-9*dabs(d(k)))
     +     go to 101

c....   New Lanczos vector
        call pzero(t3,neq)
        call promul(b(neq+1),b(neq+1),b(1),t1,t3,jp,neq,imas,1)

        dm   = 1.d0
        do n = 1,neq
c
          if (istyp.eq.0) then  ! Standard
            dm = b(n)

          else if (istyp.eq.1.or.istyp.eq.2) then ! SM
            if(imas.eq.1) then
              dm = b(jp(n))
            else
              dm = b(n)
            end if

          else if (istyp.ge.3.and.istyp.le.8) then ! all other CSR solver
            call pick_csr(dm,n,b,csrka,imas)

          end if
c
         if(dm.eq.0.0d0) t1(n) = 0.d0
        end do

        fact = ddot(neq,t1,1,t3,1)
        if(fact.lt.0.d0) then
          write (*,*) 'Negative scaling factor in LANCZOS! Use SUBSPACE'
          return
        end if
        fact = 1.d0/dsqrt(fact)

        do n = 1,neq
         v(n,k1) = fact*t1(n)
        end do

101     continue

c....   Complete reorthogonalization against the previous Lanczos vectors
        call pzero(t1,neq)
        call promul(b(neq+1),b(neq+1),b(1),v(1,k1),t1,jp,neq,imas,1)
        do i = 1,k
          fact = ddot(neq,v(1,i),1,t1,1)
          call daxpty(neq,-fact,v(1,i),v(1,k1))
        end do

c....   alpha_k+1 = q_k+1 * B * A^-1 * B * q_k+1,  beta_k = q_k * B * A^-1 * B * q_k+1
        call pmove(t2,t3,neq)
        call pzero(t2,neq)
        call promul(b(neq+1),b(neq+1),b(1),v(1,k1),t2,jp,neq,imas,1)
        call pmove(t2,t1,neq)
        call dasol(a(neq+1),a(neq+1),a(1),t1,jp,neq,dummy)

        d(k1) = ddot(neq,t2,1,t1,1)
        g(k1) = ddot(neq,t3,1,t1,1)

c....   solve reduced eigenvalue problem  T * p = d * p
        do i = 1,k1
         do j = 1,k1
          p(i,j) = 0.d0
         end do
         p(i,i) = 1.d0
        end do
        call pmove(d,dh,k1)
        call pmove(g,gh,k1)
        call tql21(nv,k1,d,g,p,ierr)
        if(ierr.ne.0) write(*,*)'No convergence in tql2'

c....   check for convergence
        do n = 1,k1
          if(d(n).ne.0.0d0) dtol(n) = dabs((d(n)-dp(n))/d(n))
          dp(n) = d(n)
        end do

        if(prt) then
          if(ior.gt.0) write(iow,2000) it,(d(n),n=1,nv)
          if(ior.lt.0) write(  *,2000) it,(d(n),n=1,nv)
          if(itlim.gt.1) then
            if(ior.gt.0) write(iow,2001) it,(dtol(n),n=1,nv)
            if(ior.lt.0) write(  *,2001) it,(dtol(n),n=1,nv)
          end if
        end if
        do n = 1,nf
          if(dtol(n).gt.tol) go to 220
        end do
        conv = .true.
c
220     continue

        if(conv) go to 305
        if(ior.lt.0) then
          write(*,2004)  (it-1)*(nv-1) + k
        end if


        if (k.lt.nv-1)  go to 100
c....   compute the new iteration vector
c
        do 230 i = 1,neq
          do l = 1,k1
            vh(l) = v(i,l)
          end do
          do 230 j = 1,k1
            v(i,j) = 0.0d0
            do l = 1,k1
              v(i,j) = v(i,j) + vh(l)*p(l,j)
            end do
230     continue

300   continue

      if(ior.lt.0) go to 306
c
c.... output converged solution
c
305   continue
      do 255 i = 1,neq
        do l = 1,k1
          vh(l) = v(i,l)
        end do
        do 255 j = 1,k1
          v(i,j) = 0.0d0
          do l = 1,k1
            v(i,j) = v(i,j) + vh(l)*p(l,j)
          end do
255   continue

      do n = 1,k1
        if(d(n).ne.0.0d0) d(n) = 1.d0/d(n) + shift
        dp(n)= dsqrt(dabs(d(n)))
        if(n.le.nv) call scalev(v(1,n),neq)
      end do

      write(iow,2010) (it-1)*(nv-1) + k -1
      do n = 1,min(2*nf,nv,neq)
        write(iow,2011) n,d(n),dtol(n),dp(n)
      end do
      if(ior.lt.0) then
        write(*,2010) (it-1)*(nv-1) + k -1
        do n = 1,min(2*nf,nv,neq)
          write(*,2011) n,d(n),dtol(n),dp(n)
        end do
      end if
      return
c.... further iterations if not converged and if needed
306   continue
      do n = 1,k1
        if(d(n).ne.0.0d0) d(n) = 1.d0/d(n) + shift
      end do
c      write(*,2000) itt,(d(n),n=1,nv)
c      write(*,2001) itt,(dtol(n),n=1,nv)
      write(*,2010) (it-2)*(nv-1) + k
      do n=1,min(2*nf,nv,neq)
        write(*,2011) n,d(n),dtol(n)
      end do
      do n = 1,nv
        if(d(n).ne.0.0d0) d(n) = 1.d0/(d(n) - shift)
      end do
      itsta = itlim+1
c.... more equations only for interactive mode
      if(iplot.ne.0) then
        write(*,2005)
        call dinput(ait,1)
        itadd = ait(1)
      else
        itadd = 0
      end if
      if(itadd.eq.0) goto 305
      itlim = itlim + itadd
      goto 301
2000  format(/5x,'current eigenvalues, iteration',i4/(4e20.8))
2001  format(5x,'current residuals,   iteration',i4/(4e20.8))
2002  format(5x,'There are',i4,' eigenvalues less than the shift'/)
c2003  format(/5x,'square root of eigenvalues'/(4e20.8))
2004  format(5x,'Completed Lanczos iteration',i4)
2005  format(5x,'Iteration not converged',/,
     1       5x,'Give number of additional iterations (0 to stop) : ',$)
2010  format(/5x,'Lanczos solution after iteration ',i3,/
     1        5x,'  No.','     Eigenvalue        ','  Residual      '
     2                 ,'      Square root of EV ')
2011  format( 4x,i5,3e20.8)
      end
c
      subroutine tql21(nm,n,d,e,z,ierr)
c---------------------------------------------------------------------+
c
c     Purpose:
c     this subroutine is a translation of the algol procedure tql2,
c     num. math. 11, 293-306(1968) by bowdler, martin, reinsch, and
c     wilkinson.
c     handbook for auto. comp., vol.ii-linear algebra, 227-240(1971).
c
c     this subroutine finds the eigenvalues and eigenvectors
c     of a symmetric tridiagonal matrix by the ql method.
c     the eigenvectors of a full symmetric matrix can also
c     be found if  tred2  has been used to reduce this
c     full matrix to tridiagonal form.
c
c     on input
c
c        nm must be set to the row dimension of two-dimensional
c          array parameters as declared in the calling program
c          dimension statement.
c        n is the order of the matrix.
c        d contains the diagonal elements of the input matrix.
c        e contains the subdiagonal elements of the input matrix
c          in its last n-1 positions.  e(1) is arbitrary.
c        z contains the transformation matrix produced in the
c          reduction by  tred2, if performed.  if the eigenvectors
c          of the tridiagonal matrix are desired, z must contain
c          the identity matrix.
c
c      on output
c        d contains the eigenvalues in ascending order.  if an
c          error exit is made, the eigenvalues are correct but
c          unordered for indices 1,2,...,ierr-1.
c        e has been destroyed.
c        z contains orthonormal eigenvectors of the symmetric
c          tridiagonal (or full) matrix.  if an error exit is made,
c          z contains the eigenvectors associated with the stored
c          eigenvalues.
c
c        ierr is set to
c          zero       for normal return,
c          j          if the j-th eigenvalue has not been
c                     determined after 30 iterations.
c
c     calls pythag for  dsqrt(a*a + b*b) .
c
c     questions and comments should be directed to burton s. garbow,
c     mathematics and computer science div, argonne national laboratory
c
c     this version dated august 1983.
c
c---------------------------------------------------------------------+
c
      integer i,j,k,l,m,n,ii,l1,l2,nm,mml,ierr
      double precision d(n),e(n),z(nm,n)
      double precision c,c2,c3,dl1,el1,f,g,h,p,r,s,s2,tst1,tst2,pythag
c
c
      ierr = 0
      if (n .eq. 1) go to 1001
c
      do 100 i = 2, n
  100 e(i-1) = e(i)
c
      f = 0.0d0
      tst1 = 0.0d0
      e(n) = 0.0d0
c
      do 240 l = 1, n
         j = 0
         h = dabs(d(l)) + dabs(e(l))
         if (tst1 .lt. h) tst1 = h
c     .......... look for small sub-diagonal element ..........
         do 110 m = l, n
            tst2 = tst1 + dabs(e(m))
            if (tst2 .eq. tst1) go to 120
c     .......... e(n) is always zero, so there is no exit
c                through the bottom of the loop ..........
  110    continue
c
  120    if (m .eq. l) go to 220
  130    if (j .eq. 30) go to 1000
         j = j + 1
c     .......... form shift ..........
         l1 = l + 1
         l2 = l1 + 1
         g = d(l)
         p = (d(l1) - g) / (2.0d0 * e(l))
         r = pythag(p,1.0d0)
         d(l) = e(l) / (p + dsign(r,p))
         d(l1) = e(l) * (p + dsign(r,p))
         dl1 = d(l1)
         h = g - d(l)
         if (l2 .gt. n) go to 145
c
         do 140 i = l2, n
  140    d(i) = d(i) - h
c
  145    f = f + h
c     .......... ql transformation ..........
         p = d(m)
         c = 1.0d0
         c2 = c
         el1 = e(l1)
         s = 0.0d0
         mml = m - l
c     .......... for i=m-1 step -1 until l do -- ..........
         do 200 ii = 1, mml
            c3 = c2
            c2 = c
            s2 = s
            i = m - ii
            g = c * e(i)
            h = c * p
            r = pythag(p,e(i))
            e(i+1) = s * r
            s = e(i) / r
            c = p / r
            p = c * d(i) - s * g
            d(i+1) = h + s * (c * g + s * d(i))
c     .......... form vector ..........
            do 180 k = 1, n
               h = z(k,i+1)
               z(k,i+1) = s * z(k,i) + c * h
               z(k,i) = c * z(k,i) - s * h
  180       continue
c
  200    continue
c
         p = -s * s2 * c3 * el1 * e(l) / dl1
         e(l) = s * p
         d(l) = c * p
         tst2 = tst1 + dabs(e(l))
         if (tst2 .gt. tst1) go to 130
  220    d(l) = d(l) + f
  240 continue
c     .......... order eigenvalues and eigenvectors ..........
      do 300 ii = 2, n
         i = ii - 1
         k = i
         p = d(i)
c
         do 260 j = ii, n
            if (dabs(d(j)) .le. dabs(p)) go to 260
            k = j
            p = d(j)
  260    continue
c
         if (k .eq. i) go to 300
         d(k) = d(i)
         d(i) = p
c
         do 280 j = 1, n
            p = z(j,i)
            z(j,i) = z(j,k)
            z(j,k) = p
  280    continue
c
  300 continue
c
      go to 1001
c     .......... set error -- no convergence to an
c                eigenvalue after 30 iterations ..........
 1000 ierr = l
 1001 return
      end
