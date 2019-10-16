      SUBROUTINE DGER  ( M, N, ALPHA, X, INCX, Y, INCY, A, LDA )
      implicit   none
c     .. Scalar Arguments ..
      DOUBLE PRECISION   ALPHA
      INTEGER            INCX, INCY, LDA, M, N
c     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), X( * ), Y( * )
c     ..

c  Purpose
c  =======

c  DGER   performs the rank 1 operation

c     A := alpha*x*y' + A,

c  where alpha is a scalar, x is an m element vector, y is an n element
c  vector and A is an m by n matrix.

c  Parameters
c  ==========

c  M      - INTEGER.
c           On entry, M specifies the number of rows of the matrix A.
c           M must be at least zero.
c           Unchanged on exit.

c  N      - INTEGER.
c           On entry, N specifies the number of columns of the matrix A.
c           N must be at least zero.
c           Unchanged on exit.

c  ALPHA  - DOUBLE PRECISION.
c           On entry, ALPHA specifies the scalar alpha.
c           Unchanged on exit.

c  X      - DOUBLE PRECISION array of dimension at least
c           ( 1 + ( m - 1 )*abs( INCX ) ).
c           Before entry, the incremented array X must contain the m
c           element vector x.
c           Unchanged on exit.

c  INCX   - INTEGER.
c           On entry, INCX specifies the increment for the elements of
c           X. INCX must not be zero.
c           Unchanged on exit.

c  Y      - DOUBLE PRECISION array of dimension at least
c           ( 1 + ( n - 1 )*abs( INCY ) ).
c           Before entry, the incremented array Y must contain the n
c           element vector y.
c           Unchanged on exit.

c  INCY   - INTEGER.
c           On entry, INCY specifies the increment for the elements of
c           Y. INCY must not be zero.
c           Unchanged on exit.

c  A      - DOUBLE PRECISION array of DIMENSION ( LDA, n ).
c           Before entry, the leading m by n part of the array A must
c           contain the matrix of coefficients. On exit, A is
c           overwritten by the updated matrix.

c  LDA    - INTEGER.
c           On entry, LDA specifies the first dimension of A as declared
c           in the calling (sub) program. LDA must be at least
c           max( 1, m ).
c           Unchanged on exit.


c  Level 2 Blas routine.

c  -- Written on 22-October-1986.
c     Jack Dongarra, Argonne National Lab.
c     Jeremy Du Croz, Nag Central Office.
c     Sven Hammarling, Nag Central Office.
c     Richard Hanson, Sandia National Labs.


c     .. Parameters ..
      DOUBLE PRECISION   ZERO
      PARAMETER        ( ZERO = 0.0D+0 )
c     .. Local Scalars ..
      DOUBLE PRECISION   TEMP
      INTEGER            I, INFO, IX, J, JY, KX
c     .. External Subroutines ..
      EXTERNAL           XERBLA
c     .. Intrinsic Functions ..
      INTRINSIC          MAX
c     ..
c     .. Executable Statements ..

c     Test the input parameters.

      INFO = 0
      IF     ( M.LT.0 )THEN
         INFO = 1
      ELSE IF( N.LT.0 )THEN
         INFO = 2
      ELSE IF( INCX.EQ.0 )THEN
         INFO = 5
      ELSE IF( INCY.EQ.0 )THEN
         INFO = 7
      ELSE IF( LDA.LT.MAX( 1, M ) )THEN
         INFO = 9
      END IF
      IF( INFO.NE.0 )THEN
         CALL XERBLA( 'DGER  ', INFO )
         RETURN
      END IF

c     Quick return if possible.

      IF( ( M.EQ.0 ).OR.( N.EQ.0 ).OR.( ALPHA.EQ.ZERO ) )
     $   RETURN

c     Start the operations. In this version the elements of A are
c     accessed sequentially with one pass through A.

      IF( INCY.GT.0 )THEN
         JY = 1
      ELSE
         JY = 1 - ( N - 1 )*INCY
      END IF
      IF( INCX.EQ.1 )THEN
         DO J = 1, N
            IF( Y( JY ).NE.ZERO )THEN
               TEMP = ALPHA*Y( JY )
               DO I = 1, M
                  A( I, J ) = A( I, J ) + X( I )*TEMP
               END DO ! I
            END IF
            JY = JY + INCY
         END DO ! J
      ELSE
         IF( INCX.GT.0 )THEN
            KX = 1
         ELSE
            KX = 1 - ( M - 1 )*INCX
         END IF
         DO J = 1, N
            IF( Y( JY ).NE.ZERO )THEN
               TEMP = ALPHA*Y( JY )
               IX   = KX
               DO I = 1, M
                  A( I, J ) = A( I, J ) + X( IX )*TEMP
                  IX        = IX        + INCX
               END DO ! I
            END IF
            JY = JY + INCY
         END DO ! J
      END IF

c     End of DGER  .

      END
