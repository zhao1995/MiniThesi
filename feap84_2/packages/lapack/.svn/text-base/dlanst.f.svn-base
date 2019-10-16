      DOUBLE PRECISION FUNCTION DLANST( NORM, N, D, E )

c  -- LAPACK auxiliary routine (version 2.0) --
c     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
c     Courant Institute, Argonne National Lab, and Rice University
c     February 29, 1992
      implicit   none

c     .. Scalar Arguments ..
      CHARACTER          NORM
      INTEGER            N
c     ..
c     .. Array Arguments ..
      DOUBLE PRECISION   D( * ), E( * )
c     ..

c  Purpose
c  =======

c  DLANST  returns the value of the one norm,  or the Frobenius norm, or
c  the  infinity norm,  or the  element of  largest absolute value  of a
c  real symmetric tridiagonal matrix A.

c  Description
c  ===========

c  DLANST returns the value

c     DLANST = ( max(abs(A(i,j))), NORM = 'M' or 'm'
c              (
c              ( norm1(A),         NORM = '1', 'O' or 'o'
c              (
c              ( normI(A),         NORM = 'I' or 'i'
c              (
c              ( normF(A),         NORM = 'F', 'f', 'E' or 'e'

c  where  norm1  denotes the  one norm of a matrix (maximum column sum),
c  normI  denotes the  infinity norm  of a matrix  (maximum row sum) and
c  normF  denotes the  Frobenius norm of a matrix (square root of sum of
c  squares).  Note that  max(abs(A(i,j)))  is not a  matrix norm.

c  Arguments
c  =========

c  NORM    (input) CHARACTER*1
c          Specifies the value to be returned in DLANST as described
c          above.

c  N       (input) INTEGER
c          The order of the matrix A.  N >= 0.  When N = 0, DLANST is
c          set to zero.

c  D       (input) DOUBLE PRECISION array, dimension (N)
c          The diagonal elements of A.

c  E       (input) DOUBLE PRECISION array, dimension (N-1)
c          The (n-1) sub-diagonal or super-diagonal elements of A.

c  =====================================================================

c     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
c     ..
c     .. Local Scalars ..
      INTEGER            I
      DOUBLE PRECISION   ANORM, SCALE, SUM
c     ..
c     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
c     ..
c     .. External Subroutines ..
      EXTERNAL           DLASSQ
c     ..
c     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, SQRT
c     ..
c     .. Executable Statements ..

      IF( N.LE.0 ) THEN
         ANORM = ZERO
      ELSE IF( LSAME( NORM, 'M' ) ) THEN

c        Find max(abs(A(i,j))).

         ANORM = ABS( D( N ) )
         DO I = 1, N - 1
            ANORM = MAX( ANORM, ABS( D( I ) ) )
            ANORM = MAX( ANORM, ABS( E( I ) ) )
         END DO ! I
      ELSE IF( LSAME( NORM, 'O' ) .OR. NORM.EQ.'1' .OR.
     $         LSAME( NORM, 'I' ) ) THEN

c        Find norm1(A).

         IF( N.EQ.1 ) THEN
            ANORM = ABS( D( 1 ) )
         ELSE
            ANORM = MAX( ABS( D( 1 ) )+ABS( E( 1 ) ),
     $              ABS( E( N-1 ) )+ABS( D( N ) ) )
            DO I = 2, N - 1
               ANORM = MAX( ANORM, ABS( D( I ) )+ABS( E( I ) )+
     $                 ABS( E( I-1 ) ) )
            END DO ! I
         END IF
      ELSE IF( ( LSAME( NORM, 'F' ) ) .OR. ( LSAME( NORM, 'E' ) ) ) THEN

c        Find normF(A).

         SCALE = ZERO
         SUM = ONE
         IF( N.GT.1 ) THEN
            CALL DLASSQ( N-1, E, 1, SCALE, SUM )
            SUM = 2*SUM
         END IF
         CALL DLASSQ( N, D, 1, SCALE, SUM )
         ANORM = SCALE*SQRT( SUM )

      ELSE
         ANORM = 0.0d0
         CALL XERBLA( 'NORM  ', 1 )
      END IF

      DLANST = ANORM

c     End of DLANST

      END
