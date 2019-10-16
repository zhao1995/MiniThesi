      SUBROUTINE DLACPY( UPLO, M, N, A, LDA, B, LDB )

c  -- LAPACK auxiliary routine (version 2.0) --
c     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
c     Courant Institute, Argonne National Lab, and Rice University
c     February 29, 1992
      implicit   none

c     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            LDA, LDB, M, N
c     ..
c     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), B( LDB, * )
c     ..

c  Purpose
c  =======

c  DLACPY copies all or part of a two-dimensional matrix A to another
c  matrix B.

c  Arguments
c  =========

c  UPLO    (input) CHARACTER*1
c          Specifies the part of the matrix A to be copied to B.
c          = 'U':      Upper triangular part
c          = 'L':      Lower triangular part
c          Otherwise:  All of the matrix A

c  M       (input) INTEGER
c          The number of rows of the matrix A.  M >= 0.

c  N       (input) INTEGER
c          The number of columns of the matrix A.  N >= 0.

c  A       (input) DOUBLE PRECISION array, dimension (LDA,N)
c          The m by n matrix A.  If UPLO = 'U', only the upper triangle
c          or trapezoid is accessed; if UPLO = 'L', only the lower
c          triangle or trapezoid is accessed.

c  LDA     (input) INTEGER
c          The leading dimension of the array A.  LDA >= max(1,M).

c  B       (output) DOUBLE PRECISION array, dimension (LDB,N)
c          On exit, B = A in the locations specified by UPLO.

c  LDB     (input) INTEGER
c          The leading dimension of the array B.  LDB >= max(1,M).

c  =====================================================================

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
         DO J = 1, N
           DO I = 1, MIN( J, M )
             B( I, J ) = A( I, J )
           END DO ! I
         END DO ! J
      ELSE IF( LSAME( UPLO, 'L' ) ) THEN
         DO J = 1, N
           DO I = J, M
             B( I, J ) = A( I, J )
           END DO ! I
         END DO ! J
      ELSE
         DO J = 1, N
           DO I = 1, M
             B( I, J ) = A( I, J )
           END DO ! I
         END DO ! J
      END IF

c     End of DLACPY

      END
