      SUBROUTINE DGEQR2( M, N, A, LDA, TAU, WORK, INFO )

c  -- LAPACK routine (version 2.0) --
c     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
c     Courant Institute, Argonne National Lab, and Rice University
c     February 29, 1992
      implicit   none

c     .. Scalar Arguments ..
      INTEGER            INFO, LDA, M, N
c     ..
c     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), TAU( * ), WORK( * )
c     ..

c  Purpose
c  =======

c  DGEQR2 computes a QR factorization of a real m by n matrix A:
c  A = Q * R.

c  Arguments
c  =========

c  M       (input) INTEGER
c          The number of rows of the matrix A.  M >= 0.

c  N       (input) INTEGER
c          The number of columns of the matrix A.  N >= 0.

c  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
c          On entry, the m by n matrix A.
c          On exit, the elements on and above the diagonal of the array
c          contain the min(m,n) by n upper trapezoidal matrix R (R is
c          upper triangular if m >= n); the elements below the diagonal,
c          with the array TAU, represent the orthogonal matrix Q as a
c          product of elementary reflectors (see Further Details).

c  LDA     (input) INTEGER
c          The leading dimension of the array A.  LDA >= max(1,M).

c  TAU     (output) DOUBLE PRECISION array, dimension (min(M,N))
c          The scalar factors of the elementary reflectors (see Further
c          Details).

c  WORK    (workspace) DOUBLE PRECISION array, dimension (N)

c  INFO    (output) INTEGER
c          = 0: successful exit
c          < 0: if INFO = -i, the i-th argument had an illegal value

c  Further Details
c  ===============

c  The matrix Q is represented as a product of elementary reflectors

c     Q = H(1) H(2) . . . H(k), where k = min(m,n).

c  Each H(i) has the form

c     H(i) = I - tau * v * v'

c  where tau is a real scalar, and v is a real vector with
c  v(1:i-1) = 0 and v(i) = 1; v(i+1:m) is stored on exit in A(i+1:m,i),
c  and tau in TAU(i).

c  =====================================================================

c     .. Parameters ..
      DOUBLE PRECISION   ONE
      PARAMETER        ( ONE = 1.0D+0 )
c     ..
c     .. Local Scalars ..
      INTEGER            I, K
      DOUBLE PRECISION   AII
c     ..
c     .. External Subroutines ..
      EXTERNAL           DLARF, DLARFG, XERBLA
c     ..
c     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
c     ..
c     .. Executable Statements ..

c     Test the input arguments

      INFO = 0
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -4
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DGEQR2', -INFO )
         RETURN
      END IF

      K = MIN( M, N )

      DO I = 1, K

c        Generate elementary reflector H(i) to annihilate A(i+1:m,i)

         CALL DLARFG( M-I+1, A( I, I ), A( MIN( I+1, M ), I ), 1,
     $                TAU( I ) )
         IF( I.LT.N ) THEN

c           Apply H(i) to A(i:m,i+1:n) from the left

            AII       = A( I, I )
            A( I, I ) = ONE
            CALL DLARF( 'Left', M-I+1, N-I, A( I, I ), 1, TAU( I ),
     $                  A( I, I+1 ), LDA, WORK )
            A( I, I ) = AII
         END IF
      END DO ! I

c     End of DGEQR2

      END
