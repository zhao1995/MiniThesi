      SUBROUTINE DLASET( UPLO, M, N, ALPHA, BETA, A, LDA )

c  -- LAPACK auxiliary routine (version 2.0) --
c     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
c     Courant Institute, Argonne National Lab, and Rice University
c     October 31, 1992
      implicit   none

c     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            LDA, M, N
      DOUBLE PRECISION   ALPHA, BETA
c     ..
c     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * )
c     ..

c  Purpose
c  =======

c  DLASET initializes an m-by-n matrix A to BETA on the diagonal and
c  ALPHA on the offdiagonals.

c  Arguments
c  =========

c  UPLO    (input) CHARACTER*1
c          Specifies the part of the matrix A to be set.
c          = 'U':      Upper triangular part is set; the strictly lower
c                      triangular part of A is not changed.
c          = 'L':      Lower triangular part is set; the strictly upper
c                      triangular part of A is not changed.
c          Otherwise:  All of the matrix A is set.

c  M       (input) INTEGER
c          The number of rows of the matrix A.  M >= 0.

c  N       (input) INTEGER
c          The number of columns of the matrix A.  N >= 0.

c  ALPHA   (input) DOUBLE PRECISION
c          The constant to which the offdiagonal elements are to be set.

c  BETA    (input) DOUBLE PRECISION
c          The constant to which the diagonal elements are to be set.

c  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
c          On exit, the leading m-by-n submatrix of A is set as follows:

c          if UPLO = 'U', A(i,j) = ALPHA, 1<=i<=j-1, 1<=j<=n,
c          if UPLO = 'L', A(i,j) = ALPHA, j+1<=i<=m, 1<=j<=n,
c          otherwise,     A(i,j) = ALPHA, 1<=i<=m, 1<=j<=n, i.ne.j,

c          and, for all UPLO, A(i,i) = BETA, 1<=i<=min(m,n).

c  LDA     (input) INTEGER
c          The leading dimension of the array A.  LDA >= max(1,M).

c =====================================================================

c     .. Local Scalars ..
      INTEGER            I, J
c     ..
c     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
c     ..
c     .. Intrinsic Functions ..
      INTRINSIC          MIN
c     ..
c     .. Executable Statements ..

      IF( LSAME( UPLO, 'U' ) ) THEN

c        Set the strictly upper triangular or trapezoidal part of the
c        array to ALPHA.

         DO J = 2, N
            DO I = 1, MIN( J-1, M )
               A( I, J ) = ALPHA
            END DO ! I
         END DO ! J

      ELSE IF( LSAME( UPLO, 'L' ) ) THEN

c        Set the strictly lower triangular or trapezoidal part of the
c        array to ALPHA.

         DO J = 1, MIN( M, N )
            DO I = J + 1, M
               A( I, J ) = ALPHA
            END DO ! I
         END DO ! J

      ELSE

c        Set the leading m-by-n submatrix to ALPHA.

         DO J = 1, N
            DO I = 1, M
               A( I, J ) = ALPHA
            END DO ! I
         END DO ! J
      END IF

c     Set the first min(M,N) diagonal elements to BETA.

      DO I = 1, MIN( M, N )
         A( I, I ) = BETA
      END DO ! I

c     End of DLASET

      END
