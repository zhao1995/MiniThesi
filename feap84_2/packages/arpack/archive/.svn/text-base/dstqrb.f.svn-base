c-----------------------------------------------------------------------
c\BeginDoc

c  This source file has been adapted from the arpack distribution and is
c  distributed with feap under the following license conditions

c  Rice BSD Software License

c  Permits source and binary redistribution of the software
c  ARPACK and P_ARPACK  for both non-commercial and commercial use.

c  Copyright (c) 2001, Rice University
c  Developed by D.C. Sorensen, R.B. Lehoucq, C. Yang, and K. Maschhoff.
c  All rights reserved.

c  Redistribution and use in source and binary forms, with or without
c  modification, are permitted provided that the following conditions are met:

c  _ Redistributions of source code must retain the above copyright notice,
c    this list of conditions and the following disclaimer.
c  _ Redistributions in binary form must reproduce the above copyright notice,
c    this list of conditions and the following disclaimer in the documentation
c    and/or other materials provided with the distribution.
c  _ If you modify the source for these routines we ask that you change the
c    name of the routine and comment the changes made to the original.
c  _ Written notification is provided to the developers of  intent to use this
c    software.

c  Also, we ask that use of ARPACK is properly cited in any resulting
c  publications or software documentation.

c  _ Neither the name of Rice University (RICE) nor the names of its
c    contributors may be used to endorse or promote products derived from
c    this software without specific prior written permission.

c  THIS SOFTWARE IS PROVIDED BY RICE AND CONTRIBUTORS "AS IS" AND
c  ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
c  TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
c  A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL RICE OR
c  CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
c  SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT
c  NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
c  LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
c  CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,
c  STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
c  ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
c  ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

c\Name: dstqrb

c\Description:
c  Computes all eigenvalues and the last component of the eigenvectors
c  of a symmetric tridiagonal matrix using the implicit QL or QR method.

c  This is mostly a modification of the LAPACK routine dsteqr.
c  See Remarks.

c\Usage:
c  call dstqrb
c     ( N, D, E, Z, WORK, INFO )

c\Arguments
c  N       Integer.  (INPUT)
c          The number of rows and columns in the matrix.  N >= 0.

c  D       Double precision array, dimension (N).  (INPUT/OUTPUT)
c          On entry, D contains the diagonal elements of the
c          tridiagonal matrix.
c          On exit, D contains the eigenvalues, in ascending order.
c          If an error exit is made, the eigenvalues are correct
c          for indices 1,2,...,INFO-1, but they are unordered and
c          may not be the smallest eigenvalues of the matrix.

c  E       Double precision array, dimension (N-1).  (INPUT/OUTPUT)
c          On entry, E contains the subdiagonal elements of the
c          tridiagonal matrix in positions 1 through N-1.
c          On exit, E has been destroyed.

c  Z       Double precision array, dimension (N).  (OUTPUT)
c          On exit, Z contains the last row of the orthonormal
c          eigenvector matrix of the symmetric tridiagonal matrix.
c          If an error exit is made, Z contains the last row of the
c          eigenvector matrix associated with the stored eigenvalues.

c  WORK    Double precision array, dimension (max(1,2*N-2)).  (WORKSPACE)
c          Workspace used in accumulating the transformation for
c          computing the last components of the eigenvectors.

c  INFO    Integer.  (OUTPUT)
c          = 0:  normal return.
c          < 0:  if INFO = -i, the i-th argument had an illegal value.
c          > 0:  if INFO = +i, the i-th eigenvalue has not converged
c                              after a total of  30*N  iterations.

c\Remarks
c  1. None.

c-----------------------------------------------------------------------

c\BeginLib

c\Local variables:
c     xxxxxx  real

c\Routines called:
c     daxpy   Level 1 BLAS that computes a vector triad.
c     dcopy   Level 1 BLAS that copies one vector to another.
c     dswap   Level 1 BLAS that swaps the contents of two vectors.
c     lsame   LAPACK character comparison routine.
c     dlae2   LAPACK routine that computes the eigenvalues of a 2-by-2
c             symmetric matrix.
c     dlaev2  LAPACK routine that eigendecomposition of a 2-by-2 symmetric
c             matrix.
c     dlamch  LAPACK routine that determines machine constants.
c     dlanst  LAPACK routine that computes the norm of a matrix.
c     dlapy2  LAPACK routine to compute sqrt(x**2+y**2) carefully.
c     dlartg  LAPACK Givens rotation construction routine.
c     dlascl  LAPACK routine for careful scaling of a matrix.
c     dlaset  LAPACK matrix initialization routine.
c     dlasr   LAPACK routine that applies an orthogonal transformation to
c             a matrix.
c     dlasrt  LAPACK sorting routine.
c     dsteqr  LAPACK routine that computes eigenvalues and eigenvectors
c             of a symmetric tridiagonal matrix.
c     xerbla  LAPACK error handler routine.

c\Authors
c     Danny Sorensen               Phuong Vu
c     Richard Lehoucq              CRPC / Rice University
c     Dept. of Computational &     Houston, Texas
c     Applied Mathematics
c     Rice University
c     Houston, Texas

c\SCCS Information: @(#)
c FILE: stqrb.F   SID: 2.5   DATE OF SID: 8/27/96   RELEASE: 2

c\Remarks
c     1. Starting with version 2.5, this routine is a modified version
c        of LAPACK version 2.0 subroutine SSTEQR. No lines are deleted,
c        only commented out and new lines inserted.
c        All lines commented out have "c$$$" at the beginning.
c        Note that the LAPACK version 1.0 subroutine SSTEQR contained
c        bugs.

c\EndLib

c-----------------------------------------------------------------------

      subroutine dstqrb ( n, d, e, z, work, info )
      implicit   none

c     %------------------%
c     | Scalar Arguments |
c     %------------------%

      integer    info, n

c     %-----------------%
c     | Array Arguments |
c     %-----------------%

      Double precision
     &           d( n ), e( n-1 ), z( n ), work( 2*n-2 )

c     .. parameters ..
      Double precision
     &                   zero         , one
      parameter        ( zero = 0.0D+0, one   = 1.0D+0  )
      Double precision
     &                   two          , three
      parameter        ( two  = 2.0D+0, three = 3.0D+0 )
      integer            maxit
      parameter        ( maxit = 30 )
c     ..
c     .. local scalars ..
      integer            i, icompz, ii, iscale, j, jtot, k, l, l1, lend,
     &                   lendm1, lendp1, lendsv, lm1, lsv, m, mm, mm1,
     &                   nm1, nmaxit
      Double precision
     &                   anorm, b, c, eps, eps2, f, g, p, r, rt1, rt2,
     &                   s, safmax, safmin, ssfmax, ssfmin, tst
c     ..
c     .. external functions ..
      logical            lsame
      Double precision
     &                   dlamch, dlanst, dlapy2
      external           lsame, dlamch, dlanst, dlapy2
c     ..
c     .. external subroutines ..
      external           dlae2, dlaev2, dlartg, dlascl, dlaset, dlasr,
     &                   dlasrt, dswap, xerbla
c     ..
c     .. intrinsic functions ..
      intrinsic          abs, max, sign, sqrt
c     ..
c     .. executable statements ..

c     test the input parameters.

      info = 0

c$$$      IF( LSAME( COMPZ, 'N' ) ) THEN
c$$$         ICOMPZ = 0
c$$$      ELSE IF( LSAME( COMPZ, 'V' ) ) THEN
c$$$         ICOMPZ = 1
c$$$      ELSE IF( LSAME( COMPZ, 'I' ) ) THEN
c$$$         ICOMPZ = 2
c$$$      ELSE
c$$$         ICOMPZ = -1
c$$$      END IF
c$$$      IF( ICOMPZ.LT.0 ) THEN
c$$$         INFO = -1
c$$$      ELSE IF( N.LT.0 ) THEN
c$$$         INFO = -2
c$$$      ELSE IF( ( LDZ.LT.1 ) .OR. ( ICOMPZ.GT.0 .AND. LDZ.LT.MAX( 1,
c$$$     $         N ) ) ) THEN
c$$$         INFO = -6
c$$$      END IF
c$$$      IF( INFO.NE.0 ) THEN
c$$$         CALL XERBLA( 'SSTEQR', -INFO )
c$$$         RETURN
c$$$      END IF

c    *** New starting with version 2.5 ***

      icompz = 2
c    *************************************

c     quick return if possible

      if( n.eq.0 )
     $   return

      if( n.eq.1 ) then
         if( icompz.eq.2 )  z( 1 ) = one
         return
      end if

c     determine the unit roundoff and over/underflow thresholds.

      eps    = dlamch( 'e' )
      eps2   = eps**2
      safmin = dlamch( 's' )
      safmax = one / safmin
      ssfmax = sqrt( safmax ) / three
      ssfmin = sqrt( safmin ) / eps2

c     compute the eigenvalues and eigenvectors of the tridiagonal
c     matrix.

c$$      if( icompz.eq.2 )
c$$$     $   call dlaset( 'full', n, n, zero, one, z, ldz )

c     *** New starting with version 2.5 ***

      if ( icompz .eq. 2 ) then
         do j = 1, n-1
            z(j) = zero
         end do ! j
         z( n ) = one
      end if
c     *************************************

      nmaxit = n*maxit
      jtot = 0

c     determine where the matrix splits and choose ql or qr iteration
c     for each block, according to whether top or bottom diagonal
c     element is smaller.

      l1 = 1
      nm1 = n - 1

   10 continue
      if( l1.gt.n )
     $   go to 160
      if( l1.gt.1 )
     $   e( l1-1 ) = zero
      if( l1.le.nm1 ) then
         do m = l1, nm1
            tst = abs( e( m ) )
            if( tst.eq.zero )
     $         go to 30
            if( tst.le.( sqrt( abs( d( m ) ) )*sqrt( abs( d( m+
     $          1 ) ) ) )*eps ) then
               e( m ) = zero
               go to 30
            end if
         end do ! m
      end if
      m = n

   30 continue
      l      = l1
      lsv    = l
      lend   = m
      lendsv = lend
      l1     = m + 1
      if( lend.eq.l )
     $   go to 10

c     scale submatrix in rows and columns l to lend

      anorm  = dlanst( 'i', lend-l+1, d( l ), e( l ) )
      iscale = 0
      if( anorm.eq.zero )
     $   go to 10
      if( anorm.gt.ssfmax ) then
         iscale = 1
         call dlascl( 'g', 0, 0, anorm, ssfmax, lend-l+1, 1, d( l ), n,
     $                info )
         call dlascl( 'g', 0, 0, anorm, ssfmax, lend-l, 1, e( l ), n,
     $                info )
      else if( anorm.lt.ssfmin ) then
         iscale = 2
         call dlascl( 'g', 0, 0, anorm, ssfmin, lend-l+1, 1, d( l ), n,
     $                info )
         call dlascl( 'g', 0, 0, anorm, ssfmin, lend-l, 1, e( l ), n,
     $                info )
      end if

c     choose between ql and qr iteration

      if( abs( d( lend ) ).lt.abs( d( l ) ) ) then
         lend = lsv
         l = lendsv
      end if

      if( lend.gt.l ) then

c        ql iteration

c        look for small subdiagonal element.

   40    continue
         if( l.ne.lend ) then
            lendm1 = lend - 1
            do m = l, lendm1
               tst = abs( e( m ) )**2
               if( tst.le.( eps2*abs( d( m ) ) )*abs( d( m+1 ) )+
     $             safmin )go to 60
            end do ! m
         end if

         m = lend

   60    continue
         if( m.lt.lend )
     $      e( m ) = zero
         p = d( l )
         if( m.eq.l )
     $      go to 80

c        if remaining matrix is 2-by-2, use dlae2 or dlaev2
c        to compute its eigensystem.

         if( m.eq.l+1 ) then
            if( icompz.gt.0 ) then
               call dlaev2( d( l ), e( l ), d( l+1 ), rt1, rt2, c, s )
               work( l ) = c
               work( n-1+l ) = s
c$$$               call dlasr( 'r', 'v', 'b', n, 2, work( l ),
c$$$     $                     work( n-1+l ), z( 1, l ), ldz )

c              *** New starting with version 2.5 ***

               tst      = z(l+1)
               z(l+1) = c*tst - s*z(l)
               z(l)   = s*tst + c*z(l)
c              *************************************
            else
               call dlae2( d( l ), e( l ), d( l+1 ), rt1, rt2 )
            end if
            d( l   ) = rt1
            d( l+1 ) = rt2
            e( l   ) = zero
            l        = l + 2
            if( l.le.lend )
     $         go to 40
            go to 140
         end if

         if( jtot.eq.nmaxit )
     $      go to 140
         jtot = jtot + 1

c        form shift.

         g = ( d( l+1 )-p ) / ( two*e( l ) )
         r = dlapy2( g, one )
         g = d( m ) - p + ( e( l ) / ( g+sign( r, g ) ) )

         s = one
         c = one
         p = zero

c        inner loop

         mm1 = m - 1
         do i = mm1, l, -1
            f = s*e( i )
            b = c*e( i )
            call dlartg( g, f, c, s, r )
            if( i.ne.m-1 )
     $         e( i+1 ) = r
            g = d( i+1 ) - p
            r = ( d( i )-g )*s + two*c*b
            p = s*r
            d( i+1 ) = g + p
            g = c*r - b

c           if eigenvectors are desired, then save rotations.

            if( icompz.gt.0 ) then
               work( i ) = c
               work( n-1+i ) = -s
            end if

         end do ! i

c        if eigenvectors are desired, then apply saved rotations.

         if( icompz.gt.0 ) then
            mm = m - l + 1
c$$$            call dlasr( 'r', 'v', 'b', n, mm, work( l ), work( n-1+l ),
c$$$     $                  z( 1, l ), ldz )

c             *** New starting with version 2.5 ***

              call dlasr( 'r', 'v', 'b', 1, mm, work( l ),
     &                    work( n-1+l ), z( l ), 1 )
c             *************************************
         end if

         d( l ) = d( l ) - p
         e( l ) = g
         go to 40

c        eigenvalue found.

   80    continue
         d( l ) = p

         l = l + 1
         if( l.le.lend )
     $      go to 40
         go to 140

      else

c        qr iteration

c        look for small superdiagonal element.

   90    continue
         if( l.ne.lend ) then
            lendp1 = lend + 1
            do m = l, lendp1, -1
               tst = abs( e( m-1 ) )**2
               if( tst.le.( eps2*abs( d( m ) ) )*abs( d( m-1 ) )+
     $             safmin )go to 110
            end do ! m
         end if

         m = lend

  110    continue
         if( m.gt.lend )
     $      e( m-1 ) = zero
         p = d( l )
         if( m.eq.l )
     $      go to 130

c        if remaining matrix is 2-by-2, use dlae2 or dlaev2
c        to compute its eigensystem.

         if( m.eq.l-1 ) then
            if( icompz.gt.0 ) then
               call dlaev2( d( l-1 ), e( l-1 ), d( l ), rt1, rt2, c, s )
c$$$               work( m ) = c
c$$$               work( n-1+m ) = s
c$$$               call dlasr( 'r', 'v', 'f', n, 2, work( m ),
c$$$     $                     work( n-1+m ), z( 1, l-1 ), ldz )

c               *** New starting with version 2.5 ***

                tst      = z(l)
                z(l)   = c*tst - s*z(l-1)
                z(l-1) = s*tst + c*z(l-1)
c               *************************************
            else
               call dlae2( d( l-1 ), e( l-1 ), d( l ), rt1, rt2 )
            end if
            d( l-1 ) = rt1
            d( l ) = rt2
            e( l-1 ) = zero
            l = l - 2
            if( l.ge.lend )
     $         go to 90
            go to 140
         end if

         if( jtot.eq.nmaxit )
     $      go to 140
         jtot = jtot + 1

c        form shift.

         g = ( d( l-1 )-p ) / ( two*e( l-1 ) )
         r = dlapy2( g, one )
         g = d( m ) - p + ( e( l-1 ) / ( g+sign( r, g ) ) )

         s = one
         c = one
         p = zero

c        inner loop

         lm1 = l - 1
         do i = m, lm1
            f = s*e( i )
            b = c*e( i )
            call dlartg( g, f, c, s, r )
            if( i.ne.m )
     $         e( i-1 ) = r
            g = d( i ) - p
            r = ( d( i+1 )-g )*s + two*c*b
            p = s*r
            d( i ) = g + p
            g = c*r - b

c           if eigenvectors are desired, then save rotations.

            if( icompz.gt.0 ) then
               work( i ) = c
               work( n-1+i ) = s
            end if

         end do ! i

c        if eigenvectors are desired, then apply saved rotations.

         if( icompz.gt.0 ) then
            mm = l - m + 1
c$$$            call dlasr( 'r', 'v', 'f', n, mm, work( m ), work( n-1+m ),
c$$$     $                  z( 1, m ), ldz )

c           *** New starting with version 2.5 ***

            call dlasr( 'r', 'v', 'f', 1, mm, work( m ), work( n-1+m ),
     &                  z( m ), 1 )
c           *************************************
         end if

         d( l ) = d( l ) - p
         e( lm1 ) = g
         go to 90

c        eigenvalue found.

  130    continue
         d( l ) = p

         l = l - 1
         if( l.ge.lend )
     $      go to 90
         go to 140

      end if

c     undo scaling if necessary

  140 continue
      if( iscale.eq.1 ) then
         call dlascl( 'g', 0, 0, ssfmax, anorm, lendsv-lsv+1, 1,
     $                d( lsv ), n, info )
         call dlascl( 'g', 0, 0, ssfmax, anorm, lendsv-lsv, 1, e( lsv ),
     $                n, info )
      else if( iscale.eq.2 ) then
         call dlascl( 'g', 0, 0, ssfmin, anorm, lendsv-lsv+1, 1,
     $                d( lsv ), n, info )
         call dlascl( 'g', 0, 0, ssfmin, anorm, lendsv-lsv, 1, e( lsv ),
     $                n, info )
      end if

c     check for no convergence to an eigenvalue after a total
c     of n*maxit iterations.

      if( jtot.lt.nmaxit )
     $   go to 10
      do i = 1, n - 1
         if( e( i ).ne.zero )
     $      info = info + 1
      end do ! i
      go to 190

c     order eigenvalues and eigenvectors.

  160 continue
      if( icompz.eq.0 ) then

c        use quick sort

         call dlasrt( 'i', n, d, info )

      else

c        use selection sort to minimize swaps of eigenvectors

         do ii = 2, n
            i = ii - 1
            k = i
            p = d( i )
            do j = ii, n
               if( d( j ).lt.p ) then
                  k = j
                  p = d( j )
               end if
            end do ! j
            if( k.ne.i ) then
               d( k ) = d( i )
               d( i ) = p
c$$$               call dswap( n, z( 1, i ), 1, z( 1, k ), 1 )
c           *** New starting with version 2.5 ***

               p    = z(k)
               z(k) = z(i)
               z(i) = p
c           *************************************
            end if
         end do ! ii
      end if

  190 continue

c     %---------------%
c     | End of dstqrb |
c     %---------------%

      end
