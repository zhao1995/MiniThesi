      SUBROUTINE DLASR( SIDE, PIVOT, DIRECT, M, N, C, S, A, LDA )

c  -- LAPACK auxiliary routine (version 2.0) --
c     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
c     Courant Institute, Argonne National Lab, and Rice University
c     October 31, 1992
      implicit   none

c     .. Scalar Arguments ..
      CHARACTER          DIRECT, PIVOT, SIDE
      INTEGER            LDA, M, N
c     ..
c     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), C( * ), S( * )
c     ..

c  Purpose
c  =======

c  DLASR   performs the transformation

c     A := P*A,   when SIDE = 'L' or 'l'  (  Left-hand side )

c     A := A*P',  when SIDE = 'R' or 'r'  ( Right-hand side )

c  where A is an m by n real matrix and P is an orthogonal matrix,
c  consisting of a sequence of plane rotations determined by the
c  parameters PIVOT and DIRECT as follows ( z = m when SIDE = 'L' or 'l'
c  and z = n when SIDE = 'R' or 'r' ):

c  When  DIRECT = 'F' or 'f'  ( Forward sequence ) then

c     P = P( z - 1 )*...*P( 2 )*P( 1 ),

c  and when DIRECT = 'B' or 'b'  ( Backward sequence ) then

c     P = P( 1 )*P( 2 )*...*P( z - 1 ),

c  where  P( k ) is a plane rotation matrix for the following planes:

c     when  PIVOT = 'V' or 'v'  ( Variable pivot ),
c        the plane ( k, k + 1 )

c     when  PIVOT = 'T' or 't'  ( Top pivot ),
c        the plane ( 1, k + 1 )

c     when  PIVOT = 'B' or 'b'  ( Bottom pivot ),
c        the plane ( k, z )

c  c( k ) and s( k )  must contain the  cosine and sine that define the
c  matrix  P( k ).  The two by two plane rotation part of the matrix
c  P( k ), R( k ), is assumed to be of the form

c     R( k ) = (  c( k )  s( k ) ).
c              ( -s( k )  c( k ) )

c  This version vectorises across rows of the array A when SIDE = 'L'.

c  Arguments
c  =========

c  SIDE    (input) CHARACTER*1
c          Specifies whether the plane rotation matrix P is applied to
c          A on the left or the right.
c          = 'L':  Left, compute A := P*A
c          = 'R':  Right, compute A:= A*P'

c  DIRECT  (input) CHARACTER*1
c          Specifies whether P is a forward or backward sequence of
c          plane rotations.
c          = 'F':  Forward, P = P( z - 1 )*...*P( 2 )*P( 1 )
c          = 'B':  Backward, P = P( 1 )*P( 2 )*...*P( z - 1 )

c  PIVOT   (input) CHARACTER*1
c          Specifies the plane for which P(k) is a plane rotation
c          matrix.
c          = 'V':  Variable pivot, the plane (k,k+1)
c          = 'T':  Top pivot, the plane (1,k+1)
c          = 'B':  Bottom pivot, the plane (k,z)

c  M       (input) INTEGER
c          The number of rows of the matrix A.  If m <= 1, an immediate
c          return is effected.

c  N       (input) INTEGER
c          The number of columns of the matrix A.  If n <= 1, an
c          immediate return is effected.

c  C, S    (input) DOUBLE PRECISION arrays, dimension
c                  (M-1) if SIDE = 'L'
c                  (N-1) if SIDE = 'R'
c          c(k) and s(k) contain the cosine and sine that define the
c          matrix P(k).  The two by two plane rotation part of the
c          matrix P(k), R(k), is assumed to be of the form
c          R( k ) = (  c( k )  s( k ) ).
c                   ( -s( k )  c( k ) )

c  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
c          The m by n matrix A.  On exit, A is overwritten by P*A if
c          SIDE = 'R' or by A*P' if SIDE = 'L'.

c  LDA     (input) INTEGER
c          The leading dimension of the array A.  LDA >= max(1,M).

c  =====================================================================

c     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
c     ..
c     .. Local Scalars ..
      INTEGER            I, INFO, J
      DOUBLE PRECISION   CTEMP, STEMP, TEMP
c     ..
c     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
c     ..
c     .. External Subroutines ..
      EXTERNAL           XERBLA
c     ..
c     .. Intrinsic Functions ..
      INTRINSIC          MAX
c     ..
c     .. Executable Statements ..

c     Test the input parameters

      INFO = 0
      IF( .NOT.( LSAME( SIDE, 'L' ) .OR. LSAME( SIDE, 'R' ) ) ) THEN
         INFO = 1
      ELSE IF( .NOT.( LSAME( PIVOT, 'V' ) .OR. LSAME( PIVOT,
     $         'T' ) .OR. LSAME( PIVOT, 'B' ) ) ) THEN
         INFO = 2
      ELSE IF( .NOT.( LSAME( DIRECT, 'F' ) .OR. LSAME( DIRECT, 'B' ) ) )
     $          THEN
         INFO = 3
      ELSE IF( M.LT.0 ) THEN
         INFO = 4
      ELSE IF( N.LT.0 ) THEN
         INFO = 5
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = 9
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DLASR ', INFO )
         RETURN
      END IF

c     Quick return if possible

      IF( ( M.EQ.0 ) .OR. ( N.EQ.0 ) )
     $   RETURN
      IF( LSAME( SIDE, 'L' ) ) THEN

c        Form  P * A

         IF( LSAME( PIVOT, 'V' ) ) THEN
            IF( LSAME( DIRECT, 'F' ) ) THEN
               DO J = 1, M - 1
                  CTEMP = C( J )
                  STEMP = S( J )
                  IF( ( CTEMP.NE.ONE ) .OR. ( STEMP.NE.ZERO ) ) THEN
                     DO I = 1, N
                        TEMP = A( J+1, I )
                        A( J+1, I ) = CTEMP*TEMP - STEMP*A( J, I )
                        A( J, I ) = STEMP*TEMP + CTEMP*A( J, I )
                     END DO ! I
                  END IF
               END DO ! J
            ELSE IF( LSAME( DIRECT, 'B' ) ) THEN
               DO J = M - 1, 1, -1
                  CTEMP = C( J )
                  STEMP = S( J )
                  IF( ( CTEMP.NE.ONE ) .OR. ( STEMP.NE.ZERO ) ) THEN
                     DO I = 1, N
                        TEMP = A( J+1, I )
                        A( J+1, I ) = CTEMP*TEMP - STEMP*A( J, I )
                        A( J, I ) = STEMP*TEMP + CTEMP*A( J, I )
                     END DO ! I
                  END IF
               END DO ! J
            END IF
         ELSE IF( LSAME( PIVOT, 'T' ) ) THEN
            IF( LSAME( DIRECT, 'F' ) ) THEN
               DO J = 2, M
                  CTEMP = C( J-1 )
                  STEMP = S( J-1 )
                  IF( ( CTEMP.NE.ONE ) .OR. ( STEMP.NE.ZERO ) ) THEN
                     DO I = 1, N
                        TEMP = A( J, I )
                        A( J, I ) = CTEMP*TEMP - STEMP*A( 1, I )
                        A( 1, I ) = STEMP*TEMP + CTEMP*A( 1, I )
                     END DO ! I
                  END IF
               END DO ! J
            ELSE IF( LSAME( DIRECT, 'B' ) ) THEN
               DO J = M, 2, -1
                  CTEMP = C( J-1 )
                  STEMP = S( J-1 )
                  IF( ( CTEMP.NE.ONE ) .OR. ( STEMP.NE.ZERO ) ) THEN
                     DO I = 1, N
                        TEMP = A( J, I )
                        A( J, I ) = CTEMP*TEMP - STEMP*A( 1, I )
                        A( 1, I ) = STEMP*TEMP + CTEMP*A( 1, I )
                     END DO ! I
                  END IF
               END DO ! J
            END IF
         ELSE IF( LSAME( PIVOT, 'B' ) ) THEN
            IF( LSAME( DIRECT, 'F' ) ) THEN
               DO J = 1, M - 1
                  CTEMP = C( J )
                  STEMP = S( J )
                  IF( ( CTEMP.NE.ONE ) .OR. ( STEMP.NE.ZERO ) ) THEN
                     DO I = 1, N
                        TEMP = A( J, I )
                        A( J, I ) = STEMP*A( M, I ) + CTEMP*TEMP
                        A( M, I ) = CTEMP*A( M, I ) - STEMP*TEMP
                     END DO ! I
                  END IF
               END DO ! J
            ELSE IF( LSAME( DIRECT, 'B' ) ) THEN
               DO J = M - 1, 1, -1
                  CTEMP = C( J )
                  STEMP = S( J )
                  IF( ( CTEMP.NE.ONE ) .OR. ( STEMP.NE.ZERO ) ) THEN
                     DO I = 1, N
                        TEMP = A( J, I )
                        A( J, I ) = STEMP*A( M, I ) + CTEMP*TEMP
                        A( M, I ) = CTEMP*A( M, I ) - STEMP*TEMP
                     END DO ! I
                  END IF
               END DO ! J
            END IF
         END IF
      ELSE IF( LSAME( SIDE, 'R' ) ) THEN

c        Form A * P'

         IF( LSAME( PIVOT, 'V' ) ) THEN
            IF( LSAME( DIRECT, 'F' ) ) THEN
               DO J = 1, N - 1
                  CTEMP = C( J )
                  STEMP = S( J )
                  IF( ( CTEMP.NE.ONE ) .OR. ( STEMP.NE.ZERO ) ) THEN
                     DO I = 1, M
                        TEMP = A( I, J+1 )
                        A( I, J+1 ) = CTEMP*TEMP - STEMP*A( I, J )
                        A( I, J ) = STEMP*TEMP + CTEMP*A( I, J )
                     END DO ! I
                  END IF
               END DO ! J
            ELSE IF( LSAME( DIRECT, 'B' ) ) THEN
               DO J = N - 1, 1, -1
                  CTEMP = C( J )
                  STEMP = S( J )
                  IF( ( CTEMP.NE.ONE ) .OR. ( STEMP.NE.ZERO ) ) THEN
                     DO I = 1, M
                        TEMP = A( I, J+1 )
                        A( I, J+1 ) = CTEMP*TEMP - STEMP*A( I, J )
                        A( I, J ) = STEMP*TEMP + CTEMP*A( I, J )
                     END DO ! I
                  END IF
               END DO ! J
            END IF
         ELSE IF( LSAME( PIVOT, 'T' ) ) THEN
            IF( LSAME( DIRECT, 'F' ) ) THEN
               DO J = 2, N
                  CTEMP = C( J-1 )
                  STEMP = S( J-1 )
                  IF( ( CTEMP.NE.ONE ) .OR. ( STEMP.NE.ZERO ) ) THEN
                     DO I = 1, M
                        TEMP = A( I, J )
                        A( I, J ) = CTEMP*TEMP - STEMP*A( I, 1 )
                        A( I, 1 ) = STEMP*TEMP + CTEMP*A( I, 1 )
                     END DO ! I
                  END IF
               END DO ! J
            ELSE IF( LSAME( DIRECT, 'B' ) ) THEN
               DO J = N, 2, -1
                  CTEMP = C( J-1 )
                  STEMP = S( J-1 )
                  IF( ( CTEMP.NE.ONE ) .OR. ( STEMP.NE.ZERO ) ) THEN
                     DO I = 1, M
                        TEMP = A( I, J )
                        A( I, J ) = CTEMP*TEMP - STEMP*A( I, 1 )
                        A( I, 1 ) = STEMP*TEMP + CTEMP*A( I, 1 )
                     END DO ! I
                  END IF
               END DO ! J
            END IF
         ELSE IF( LSAME( PIVOT, 'B' ) ) THEN
            IF( LSAME( DIRECT, 'F' ) ) THEN
               DO J = 1, N - 1
                  CTEMP = C( J )
                  STEMP = S( J )
                  IF( ( CTEMP.NE.ONE ) .OR. ( STEMP.NE.ZERO ) ) THEN
                     DO I = 1, M
                        TEMP = A( I, J )
                        A( I, J ) = STEMP*A( I, N ) + CTEMP*TEMP
                        A( I, N ) = CTEMP*A( I, N ) - STEMP*TEMP
                     END DO ! I
                  END IF
               END DO ! J
            ELSE IF( LSAME( DIRECT, 'B' ) ) THEN
               DO J = N - 1, 1, -1
                  CTEMP = C( J )
                  STEMP = S( J )
                  IF( ( CTEMP.NE.ONE ) .OR. ( STEMP.NE.ZERO ) ) THEN
                     DO I = 1, M
                        TEMP = A( I, J )
                        A( I, J ) = STEMP*A( I, N ) + CTEMP*TEMP
                        A( I, N ) = CTEMP*A( I, N ) - STEMP*TEMP
                     END DO ! I
                  END IF
               END DO ! J
            END IF
         END IF
      END IF

c     End of DLASR

      END
