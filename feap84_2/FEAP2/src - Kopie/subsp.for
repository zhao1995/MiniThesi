      subroutine subsp(a,b,v,t,g,h,d,dp,dtol,p,z,jp,nf,nv,neq,imas,
     1                 shift,tol,prt,its)
c----------------------------------------------------------------------
c
c      Purpose: Subspace iteration to extract lowest nf eigenpairs
c               Solves:  (A - shift*B)*V = B*V*d for V and d
c               Theory see e.g. Book of Hughes chap.10
c               terms added for sm-solver MAS+GEOM  WW 3/04
c
c      Inputs:
c         a(*)      - Coefficient matrix (tangent stiffness)
c         b(*)      - Coefficient matrix (mass or geometric stiffness)
c         jp(*)     - Pointer array for row/columns of tangent
c         nf        - Number of pairs to converge
c         nv        - Size of subspace problem > or = nf
c         neq       - Size of A,B
c         imas      - Switch: =1 for consistent B; =2 for diagonal B
c         shift     - Value of shift
c         tol       - Tolerance to converge pairs
c         prt       - Flag, output iteration arrays if true
c         its       - Maximum number of subspace iterations
c
c         z
c
c      Scratch:
c         t(neq)    - Working vector
c         g(*)      - Projection of A - shift*B matrix
c         h(*)      - Projection of B           matrix
c         dp(*)     - Previous iteration values
c         dtol(*)   - Tolerance of eigenvalue iterations
c         p(nv,*)   - Eigenvectors of G*P = H*P*d
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
      dimension a(*),b(*),v(neq,*),t(*),g(*),h(*),d(*),dp(*),dtol(*),
     1          p(nv,*),z(neq,*),jp(*),ait(1)

c.... compute the initial iteration vectors
      call pzero(v,nv*neq)
      num  = 0
      nmas = 0
c.... on diagonal
c.... count number of values less than the current shift
c.... count the number of nonzero masses
      if (istyp.eq.0) then ! Standard solver
        do  n = 1,neq
          if(a(n).lt.0.0d0) num  = num  + 1
          if(b(n).ne.0.0d0) nmas = nmas + 1
        end do

      else if(istyp.eq.1.or.istyp.eq.2) then ! SM
        do n = 1,neq
          if( a(jp(n)).lt.0.0d0) num  = num  + 1
          if(imas.eq.1) then
            if( b(jp(n)).ne.0.0d0) nmas = nmas + 1
          else
            if(b(n).ne.0.0d0) nmas = nmas + 1
          end if
        end do

      else if(istyp.ge.3.and.istyp.le.8) then ! all other CSR solver
        call dnzero_csr(neq,a,csrka, num,2,   1)
        call dnzero_csr(neq,b,csrka,nmas,1,imas)

      end if

                   write(iow,2002) num
      if(ior.lt.0) write(  *,2002) num
c
      nmas = nmas/nv
      i = 0
      j = 1

      do n = 1,neq
        if (istyp.eq.0) then ! Standard
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
          v(n,j) = dm
          i      = i + 1
          if(mod(i,nmas).eq.0) j = j + 1
          j      = min0(j,nv)

        end if
      end do

      do i = 1,nv
        dp(i) = 0.d0
        dtol(i) = 1.d0
        call scalev(v(1,i),neq)
      end do
c
c.... compute the new vectors and project 'a' onto 'g'
c
      conv = .false.
      itlim = its
      if(nv.eq.nf) itlim = 1
      itsta = 1
301   continue
      do 300 it = itsta,itlim
        itt = it
c
c.... project the 'b' matrix to form 'h' and compute 'z' vectors
c
        call sprojb(b,v,t,h,jp,neq,nv,imas)
c
c.... project the 'a' matrix to form 'g'
c
        call sproja(a,v,z,g,dtol,jp,neq,nv)
c
c.... solve the reduced eigenproblem
c
        if(imtyp.eq.2) then ! ->geom
c
c....     eigenproblem: 'h*p = g*p*1/d'
c
          call geig(h,g,d,p,t,nv,prt)
          do 130 n = 1,nv
            if(d(n).ne.0.0d0) d(n) = 1.d0/d(n)
130       continue
c....     resort the eigenvalues and vectors to ascending order
          num = nv + 1
          do 134 n = 1,nv/2
            dm       = d(n)
            d(n)     = d(num-n)
            d(num-n) = dm
            do 132 i = 1,nv
              dm         = p(i,n)
              p(i,n)     = p(i,num-n)
              p(i,num-n) = dm
132         continue
134       continue
        else ! -> lmas,iden,cmas
c
c....     eigenproblem: 'g*p = h*p*d'
c
          call geig(g,h,d,p,t,nv,prt)
        end if
c
c....   check for convergence
c
        do 200 n = 1,nv
          if(d(n).ne.0.0d0) dtol(n) = dabs((d(n)-dp(n))/d(n))
          dp(n) = d(n)
200     continue
        if(prt) then
          if(ior.gt.0) write(iow,2000) it,(d(n),n=1,nv)
          if(ior.lt.0) write(  *,2000) it,(d(n),n=1,nv)
          if(itlim.gt.1) then
            if(ior.gt.0) write(iow,2001) it,(dtol(n),n=1,nv)
            if(ior.lt.0) write(  *,2001) it,(dtol(n),n=1,nv)
          end if
        end if
        do 210 n = 1,nf
          if(dtol(n).gt.tol) go to 220
210     continue
        conv = .true.
c
c...  divide eigenvectors by eigenvalue to prevent overflows
c
220     do 235 i = 1,nv
          div = d(i)
          if(p(i,i).lt.-0.00001d0) div = -div
          do 230 j = 1,nv
            p(j,i) = p(j,i)/div
230       continue
235     continue
c
c.... compute the new iteration vector 'u' from 'z'
c
        do 255 i = 1,neq
        do 255 j = 1,nv
          v(i,j) = 0.0d0
          do 250 k = 1,nv
            v(i,j) = v(i,j) + z(i,k)*p(k,j)
250       continue
255     continue
        if(conv) go to 305
        if(ior.lt.0) then
          write(*,2004) it
        end if
300   continue
      if(ior.lt.0) go to 306
c
c.... scale the vectors to have maximum element of 1.0
c
305   do 310 n = 1,nv
        d(n) = 1.d0/d(n) + shift
        dp(n)= dsqrt(dabs(d(n)))
        if(n.le.nv) call scalev(v(1,n),neq)
310   continue
      write(iow,2010) itt
      do n = 1,nv
        write(iow,2011) n,d(n),dtol(n),dp(n)
      end do
      if(ior.lt.0) then
        write(*,2010) itt
        do n = 1,nv
          write(*,2011) n,d(n),dtol(n),dp(n)
        end do
      end if
      return
c.... further iterations if not converged and if needed
306   do 311 n = 1,nv
        d(n) = 1.d0/d(n) + shift
311   continue
c      write(*,2000) itt,(d(n),n=1,nv)
c      write(*,2001) itt,(dtol(n),n=1,nv)
      write(*,2010) itt
      do n=1,nv
        write(*,2011) n,d(n),dtol(n)
      end do
      do 312 n = 1,nv
        d(n) = 1.d0/(d(n) - shift)
312   continue
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
c2003  format(/5x,'square root of eigenvalues'/(4e 20.8))
2004  format(5x,'Completed subspace iteration',i4)
2005  format(5x,'Iteration not converged',/,
     1       5x,'Give number of additional iterations (0 to stop) : ',$)
2010  format(/5x,'Subspace solution after iteration ',i3,/
     1        5x,'  No.','     Eigenvalue        ','  Residual      '
     2                 ,'      Square root of EV ')
2011  format( 4x,i5,3e20.8)
      end
c
      subroutine sproja(a,v,z,g,dtol,jp,neq,nv)
c----------------------------------------------------------------------
c
c      Purpose: Compute subspace projection of 'a' to form 'g'
c
c      Inputs:     -
c         a(*)
c         v(neq,*) - Set of iteration vectors
c         z(neq,*) -
c         dtol(*)  -
c         jp(*)    - Pointer array for row/columns of tangent
c         neq      - Number of equations in A
c         nv       - Size of projected matrix

c      Scratch:

c      Outputs:
c         g(*)     - Projected matrix V_trans * A * V
c
c----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      dimension a(*),v(neq,*),z(neq,*),g(*),dtol(*),jp(*)

c.... forward reduce the eigenvector estimates
      do 400 i = 1,nv
c.... copy a vector 'v' into 'z'
      call pmove (v(1,i),z(1,i),neq)
c.... solve the equations
      call dasol(a(neq+1),a(neq+1),a(1),z(1,i),jp,neq,en)

400   continue
c.... compute the projection of the stiffness
      k = 0
      do 500 j = 1,nv
      do 500 i = 1,j
      k = k + 1
500   g(k) = ddot(neq,v(1,i),1,z(1,j),1)

      return
      end
c
      subroutine sprojb(b,v,t,h,jp,neq,nv,imas)
c----------------------------------------------------------------------
c
c      Purpose: Compute subspace projection of 'b' to form 'h'
c
c      Inputs:
c         b(*)     - Symmetric coefficient matrix for eigenproblem
c         v(neq,*) - Set of iteration vectors
c         jp(*)    - Pointer array for row/columns of tangent
c         neq      - Number of equations in B
c         nv       - Size of projected matrix
c         imas     - Mass type: 1 = consistent; 2 = diagonal.
c
c      Scratch:
c         t(neq)   - Working vector
c
c      Outputs:
c         h(*)     - Projected matrix V_trans * B * V
c
c----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      dimension b(*),v(neq,*),t(*),h(*),jp(*)
c.... compute 'z' and the 'b' projection to form 'h'
      do 110 j = 1,nv
c.... compute 'z'
cww      go to (1,2), imas
cwwc.... consistent mass
cww1     call pzero(t,neq)
cww      call promul(b(neq+1),b(neq+1),b(1),v(1,j),t,jp,neq,1,1)
cww      go to 3
cwwc.... lumped mass
cww2     do 200 i = 1,neq
cww200   t(i) = v(i,j)*b(i)

      call pzero(t,neq)
      call promul(b(neq+1),b(neq+1),b(1),v(1,j),t,jp,neq,imas,1)

c.... project the'z' and 'v' vectors to form 'h'
3     k = j*(j+1)/2
      do 100 i = j,nv
      h(k) = ddot(neq,t,1,v(1,i),1)
100   k = k + i
      do 110 i = 1,neq
110   v(i,j) = t(i)

      return
      end
c


      subroutine geig(g,h,d,p,t,nv,prt)
      implicit double precision (a-h,o-z)
c----------------------------------------------------------------------
c      Purpose: Solve general eigenproblem 'g*p = h*p*d'
c
c      Inputs:
c         g(*)   - Left hand projected array
c         h(*)   - Right hand projected array
c         nv     - Size of problem
c         prt    - Output computations if true
c
c      Outputs:
c         d(*)   - Eigenvalues of problem
c         p(*,*) - Eigenvectors of problem
c
c----------------------------------------------------------------------
c
      logical prt
      dimension g(*),h(*),d(*),p(nv,*),t(*)
c
c.... compute the choleski factors of 'h'
c
      if(prt) call wprojm(g,nv,1)
      if(prt) call wprojm(h,nv,2)

      call chlfac(h,nv)
c
c.... compute the standard eigenvalue problem matrix 'c'
c
      call chlfwd(h,g,p,nv)
c
c.... perform the eignfunction decomposition of 'c'
c
      call eisql(g,d,t,p,nv,ir)
c
c.... compute the vectors of the original problem
c
      call chlbac(h,p,nv)
      if(prt) call mprint(p,nv,nv,nv,'vectors p')
      return
      end
c
      subroutine wprojm(a,nn,ia)
c----------------------------------------------------------------------
c
c      Purpose: Outputs projected subspace arrays: G and H
c
c      Inputs:
c         a(*)        - Array to output
c         nn          - Number row/columns in array
c         ia          - Header to write (G or H)
c
c      Outputs:
c         none
c----------------------------------------------------------------------
      USE iofile
      implicit double precision (a-h,o-z)
      character ah(2)
      dimension a(*)
      data ah(1),ah(2) /'g','h'/
      if(ior.gt.0) write(iow,2000) ah(ia)
      if(ior.lt.0) write(  *,2000) ah(ia)
      i = 1
      do 100 n = 1,nn
      j = i + n - 1
      if(ior.gt.0) write(iow,2001) (a(k),k=i,j)
      if(ior.lt.0) write(  *,2001) (a(k),k=i,j)
100   i = i + n
      return
2000  format(' matrix ',a1)
2001  format(1p8d10.2)
      end
c
      subroutine chlbac(u,s,nn)
c----------------------------------------------------------------------
c
c      Purpose: Back substitution for Cholesky factors in eigen
c               solutions.
c
c      Inputs:
c        u(*)   - Unreduced array
c        s(*,*) - Factored array of matrix
c        nn     - Size of arrays

c      Outputs:
c        u(*)   - Solution after back substitution
c----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      dimension u(*),s(nn,nn)
      j = nn
      jd = nn*(nn+1)/2
      do 100 i = 1,nn
        s(nn,i) = s(nn,i)/u(jd)
100   continue
200   jd = jd - j
      j = j - 1
      if(j.le.0) return
      do 300 i = 1,nn
        call colbac(u(jd+1),s(1,i),u(jd),j)
300   continue
      go to 200
      end
c
      subroutine chlfac(a,nn)
c----------------------------------------------------------------------
c
c      Purpose: Factor a positive definite matrix using Cholesky
c               Method.  Array stored by columns for upper part.
c
c      Inputs:
c         a(*) - Unfactored array
c         nn     - Size of array
c
c      Outputs:
c         a(*,*) - Factored array
c
c----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      dimension a(*)

      if(a(1).lt.0.0d0) then
        j=1
        goto 10
      end if
      a(1) = dsqrt(a(1))
      if(nn.eq.1) return
      jd = 1
      do 200 j = 2,nn
        jm = j - 1
        id = 0
        do 100 i = 1,jm
          if(i-1.gt.0) a(jd+i) = a(jd+i) - ddot(i-1,a(id+1),1,a(jd+1),1)
          id = id + i
          a(jd+i) = a(jd+i)/a(id)
100     continue
        fact    = a(jd+j) - ddot(jm,a(jd+1),1,a(jd+1),1)
        if(fact.lt.0.0d0) goto 10
        a(jd+j) = dsqrt(fact)
        jd = jd + j
200   continue
      return
10    write(*,20) j
20    format(5x,'negative element in choleski factorization',/,
     1       5x,'equation number ',i5,' check mass matrix!')


      return
      end
c
      subroutine chlfwd(u,g,s,nn)
c----------------------------------------------------------------------
c
c      Purpose: Use Cholesky factors to project onto a standard
c                                                      eigenproblem
c
c      Inputs:
c         g(*)  - Symmetric projected matrix
c         u(*)  - Upper factor for projection
c         nn    - Size of arrays
c
c      Outputs:
c         s(*,*) - Projected array
c----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      dimension u(*),g(*),s(nn,nn)
      s(1,1) = g(1)/u(1)
      if(nn.eq.1) go to 300
      id = 1
      do 200 i = 2,nn
        s(1,i) = g(id+1)/u(1)
        im = i - 1
        jd = 0
        do 100 j = 1,im
         s(i,j) = (g(id+j)-ddot(im,u(id+1),1,s(1,j),1))/u(id+i)
         if(j.gt.1)s(j,i)=(g(id+j)-ddot(j-1,u(jd+1),1,s(1,i),1))/u(jd+j)
         jd = jd + j
100     continue
        id = id + i
        s(i,i) = (g(id) - ddot(im,u(id-im),1,s(1,i),1))/u(id)
200   continue
c
c.... complete projection
c
300   g(1) = s(1,1)/u(1)
      if(nn.eq.1) return
      jd = 2
      do 500 j = 2,nn
        g(jd) = s(j,1)/u(1)
        id = 2
        do 400 i = 2,j
          im = i - 1
          g(jd+im) = (s(j,i) - ddot(im,u(id),1,g(jd),1))/u(id+im)
          id = id + i
400     continue
        jd = jd + j
500   continue

      return
      end
c
      subroutine colbac(u,s,d,jj)
c----------------------------------------------------------------------
c
c      Purpose: Backsubstitution macro for eigen solution
c
c      Inputs:
c         s(*)  - Unreduced column
c         u(*)  - Column of upper array already reduced
c         d     - Solution value for 'u' column
c         jj    - Length to reduce
c
c      Outputs:
c         s(*)  - Reduced column
c----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      dimension u(*),s(*)
      dd = s(jj+1)
      do 100 j = 1,jj
        s(j) = s(j) - dd*u(j)
100   continue
      s(jj) = s(jj)/d
      return
      end
c
      subroutine eisql(a,d,e,z,n,ierr)
c----------------------------------------------------------------------
c
c      Purpose: Compute eigenvalues and eigenvectors for standard
c               eigenproblem

c      Inputs:
c         a(*)   - Matrix for wanted eigenvalues
c         n      - size of eigenproblem
c
c      Outputs:
c         d(n)   - Eigenvalues
c         z(n,n) - Eigenvectors
c         ierr   - Error indicator
c
c      Scratch:
c         e(*)   - Working vector
c
c      Comment: eispac ql algorithm adapted from 'tred2' and 'tql2'

c----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      dimension a(*),d(*),e(*),z(n,n)
      double precision machep
      data machep/0.13877788d-16/
      n2 = 0
      do 100 i = 1,n
      do 100 j = 1,i
      n2 = n2 + 1
100   z(i,j) = a(n2)
      if(n.eq.1) go to 320
      n2 = n + 2
      do 300 ii = 2,n
      i = n2 - ii
      l = i - 1
      h = 0.0d0
      scale = 0.0d0
      if(l.lt.2) go to 130
      do 120 k = 1,l
120   scale = scale + abs(z(i,k))
      if(scale.ne.0.0d0) go to 140
130   e(i) = z(i,l)
      go to 290
140   do 150 k = 1,l
      z(i,k) = z(i,k)/scale
150   h = h + z(i,k)*z(i,k)
      f = z(i,l)
      g = -dsign(sqrt(h),f)
      e(i) = scale*g
      h = h - f*g
      z(i,l) = f - g
      f = 0.0d0
      do 240 j = 1,l
      z(j,i) = z(i,j)/h
      g = 0.0d0
      do 180 k = 1,j
180   g = g + z(j,k)*z(i,k)
      jp1 = j + 1
      if(l.lt.jp1) go to 220
      do 200 k = jp1,l
200   g = g + z(k,j)*z(i,k)
220   e(j) = g/h
      f = f + e(j)*z(i,j)
240   continue
      hh = f/(h+h)
      do 260 j = 1,l
      f = z(i,j)
      g = e(j) - hh*f
      e(j) = g
      do 260 k = 1,j
260   z(j,k) = z(j,k) - f*e(k) - g*z(i,k)
290   d(i) = h
300   continue
c.... set transformation array for ql
320   d(1) = z(1,1)
      z(1,1) = 1.0d0
      e(1) = 0.0d0
      ierr = 0
      if(n.eq.1) go to 1001
      do 500 i = 2,n
      l = i - 1
      if(d(i).eq.0.0d0) go to 380
      do 360 j = 1,l
      g = 0.0d0
      do 340 k = 1,l
340   g = g + z(i,k)*z(k,j)
      do 360 k = 1,l
360   z(k,j) = z(k,j) - g*z(k,i)
380   d(i) = z(i,i)
      z(i,i) = 1.0d0
      do 400 j = 1,l
      z(i,j) = 0.0d0
400   z(j,i) = 0.0d0
500   continue
c.... begin 'ql' algorithm on tridagonal matrix now stored in 'd' and 'e
      do 600 i = 2,n
600   e(i-1) = e(i)
      f = 0.0d0
      b = 0.0d0
      e(n) = 0.0d0
      do 840 l = 1,n
      j = 0
      h = machep*(abs(d(l)) + abs(e(l)))
      if(b.lt.h) b = h
      do 710 m = l,n
      if(abs(e(m)).le.b) go to 720
710   continue
720   if(m.eq.l) go to 820
730   if(j.eq.30) go to 1000
      j = j + 1
      l1 = l + 1
      g = d(l)
      p = (d(l1)-g)/(e(l)+e(l))
      r = sqrt(p*p+1.0d0)
      d(l) = e(l)/(p+dsign(r,p))
      h = g - d(l)
      do 740 i = l1,n
740   d(i) = d(i) - h
      f = f + h
      p = d(m)
      c = 1.0d0
      s = 0.0d0
      mml = m - l
      do 800 ii = 1,mml
      i = m - ii
      g = c*e(i)
      h = c*p
      if(abs(p).lt.abs(e(i))) go to 750
      c = e(i)/p
      r = sqrt(c*c+1.0d0)
      e(i+1) = s*p*r
      s = c/r
      c = 1.0d0/r
      go to 760
750   c = p/e(i)
      r = sqrt(c*c+1.0d0)
      e(i+1) = s*e(i)*r
      s = 1.0d0/r
      c = c*s
760   p = c*d(i) - s*g
      d(i+1) = h + s*(c*g + s*d(i))
      do 780 k = 1,n
      h = z(k,i+1)
      z(k,i+1) = s*z(k,i) + c*h
780   z(k,i  ) = c*z(k,i) - s*h
800   continue
      e(l) = s*p
      d(l) = c*p
      if(abs(e(l)).gt.b) go to 730
820   d(l) = d(l) + f
840   continue
      do 900 ii = 2,n
      i = ii - 1
      k = i
      p = d(i)
      do 860 j = ii,n
      if(abs(d(j)).le.abs(p)) go to 860
      k = j
      p = d(j)
860   continue
      if(k.eq.i) go to 900
      d(k) = d(i)
      d(i) = p
      do 880 j = 1,n
      p = z(j,i)
      z(j,i) = z(j,k)
880   z(j,k) = p
900   continue
      go to 1001
1000  ierr = l
1001  return
      end
c
