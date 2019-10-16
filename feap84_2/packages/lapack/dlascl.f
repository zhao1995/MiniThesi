      SUBROUTINE DLASCL( TYPE, KL, KU, CFROM, CTO, M, N, A, LDA, INFO )

c  -- LAPACK auxiliary routine (version 2.0) --
c     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
c     Courant Institute, Argonne National Lab, and Rice University
c     February 29, 1992
      implicit   none

c     .. Scalar Arguments ..
      CHARACTER          TYPE
      INTEGER            INFO, KL, KU, LDA, M, N
      DOUBLE PRECISION   CFROM, CTO
c     ..
c     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * )
c     ..

c  Purpose
c  =======

c  DLASCL multiplies the M by N real matrix A by the real scalar
c  CTO/CFROM.  This is done without over/underflow as long as the final
c  result CTO*A(I,J)/CFROM does not over/underflow. TYPE specifies that
c  A may be full, upper triangular, lower triangular, upper Hessenberg,
c  or banded.

c  Arguments
c  =========

c  TYPE    (input) CHARACTER*1
c          TYPE indices the storage type of the input matrix.
c          = 'G':  A is a full matrix.
c          = 'L':  A is a lower triangular matrix.
c          = 'U':  A is an upper triangular matrix.
c          = 'H':  A is an upper Hessenberg matrix.
c          = 'B':  A is a symmetric band matrix with lower bandwidth KL
c                  and upper bandwidth KU and with the only the lower
c                  half stored.
c          = 'Q':  A is a symmetric band matrix with lower bandwidth KL
c                  and upper bandwidth KU and with the only the upper
c                  half stored.
c          = 'Z':  A is a band matrix with lower bandwidth KL and upper
c                  bandwidth KU.

c  KL      (input) INTEGER
c          The lower bandwidth of A.  Referenced only if TYPE = 'B',
c          'Q' or 'Z'.

c  KU      (input) INTEGER
c          The upper bandwidth of A.  Referenced only if TYPE = 'B',
c          'Q' or 'Z'.

c  CFROM   (input) DOUBLE PRECISION
c  CTO     (input) DOUBLE PRECISION
c          The matrix A is multiplied by CTO/CFROM. A(I,J) is computed
c          without over/underflow if the final result CTO*A(I,J)/CFROM
c          can be represented without over/underflow.  CFROM must be
c          nonzero.

c  M       (input) INTEGER
c          The number of rows of the matrix A.  M >= 0.

c  N       (input) INTEGER
c          The number of columns of the matrix A.  N >= 0.

c  A       (input/output) DOUBLE PRECISION array, dimension (LDA,M)
c          The matrix to be multiplied by CTO/CFROM.  See TYPE for the
c          storage type.

c  LDA     (input) INTEGER
c          The leading dimension of the array A.  LDA >= max(1,M).

c  INFO    (output) INTEGER
c          0  - successful exit
c          <0 - if INFO = -i, the i-th argument had an illegal value.

c  =====================================================================

c     .. Parameters ..
      DOUBLE PRECISION   ZERO        , ONE
      PARAMETER        ( ZERO = 0.0D0, ONE = 1.0D0 )
c     ..
c     .. Local Scalars ..
      LOGICAL            DONE
      INTEGER            I, ITYPE, J, K1, K2, K3, K4
      DOUBLE PRECISION   BIGNUM, CFROM1, CFROMC, CTO1, CTOC, MUL, SMLNUM
c     ..
c     .. External Functions ..
      LOGICAL            LSAME
      DOUBLE PRECISION   DLAMCH
      EXTERNAL           LSAME, DLAMCH
c     ..
c     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, MIN
c     ..
c     .. External Subroutines ..
      EXTERNAL           XERBLA
c     ..
c     .. Executable Statements ..

c     Test the input arguments

      INFO = 0

      IF(      LSAME( TYPE, 'G' ) ) THEN
         ITYPE = 0
      ELSE IF( LSAME( TYPE, 'L' ) ) THEN
         ITYPE = 1
      ELSE IF( LSAME( TYPE, 'U' ) ) THEN
         ITYPE = 2
      ELSE IF( LSAME( TYPE, 'H' ) ) THEN
         ITYPE = 3
      ELSE IF( LSAME( TYPE, 'B' ) ) THEN
         ITYPE = 4
      ELSE IF( LSAME( TYPE, 'Q' ) ) THEN
         ITYPE = 5
      ELSE IF( LSAME( TYPE, 'Z' ) ) THEN
         ITYPE = 6
      ELSE
         ITYPE = -1
      END IF

      IF( ITYPE.EQ.-1 ) THEN
         INFO = -1
      ELSE IF( CFROM.EQ.ZERO ) THEN
         INFO = -4
      ELSE IF( M.LT.0 ) THEN
         INFO = -6
      ELSE IF( N.LT.0 .OR. ( ITYPE.EQ.4 .AND. N.NE.M ) .OR.
     $         ( ITYPE.EQ.5 .AND. N.NE.M ) ) THEN
         INFO = -7
      ELSE IF( ITYPE.LE.3 .AND. LDA.LT.MAX( 1, M ) ) THEN
         INFO = -9
      ELSE IF( ITYPE.GE.4 ) THEN
         IF( KL.LT.0 .OR. KL.GT.MAX( M-1, 0 ) ) THEN
            INFO = -2
         ELSE IF( KU.LT.0 .OR. KU.GT.MAX( N-1, 0 ) .OR.
     $            ( ( ITYPE.EQ.4 .OR. ITYPE.EQ.5 ) .AND. KL.NE.KU ) )
     $             THEN
            INFO = -3
         ELSE IF( ( ITYPE.EQ.4 .AND. LDA.LT.KL+1 ) .OR.
     $            ( ITYPE.EQ.5 .AND. LDA.LT.KU+1 ) .OR.
     $            ( ITYPE.EQ.6 .AND. LDA.LT.2*KL+KU+1 ) ) THEN
            INFO = -9
         END IF
      END IF

      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DLASCL', -INFO )
         RETURN
      END IF

c     Quick return if possible

      IF( N.EQ.0 .OR. M.EQ.0 )
     $   RETURN

c     Get machine parameters

      SMLNUM = DLAMCH( 'S' )
      BIGNUM = ONE / SMLNUM

      CFROMC = CFROM
      CTOC = CTO

   10 CONTINUE
      CFROM1 = CFROMC*SMLNUM
      CTO1 = CTOC / BIGNUM
      IF( ABS( CFROM1 ).GT.ABS( CTOC ) .AND. CTOC.NE.ZERO ) THEN
         MUL = SMLNUM
         DONE = .FALSE.
         CFROMC = CFROM1
      ELSE IF( ABS( CTO1 ).GT.ABS( CFROMC ) ) THEN
         MUL = BIGNUM
         DONE = .FALSE.
         CTOC = CTO1
      ELSE
         MUL = CTOC / CFROMC
         DONE = .TRUE.
      END IF

      IF( ITYPE.EQ.0 ) THEN

c        Full matrix

         DO J = 1, N
            DO I = 1, M
               A( I, J ) = A( I, J )*MUL
            END DO ! I
         END DO ! J

      ELSE IF( ITYPE.EQ.1 ) THEN

c        Lower triangular matrix

         DO J = 1, N
            DO I = J, M
               A( I, J ) = A( I, J )*MUL
            END DO ! I
         END DO ! J

      ELSE IF( ITYPE.EQ.2 ) THEN

c        Upper triangular matrix

         DO J = 1, N
            DO I = 1, MIN( J, M )
               A( I, J ) = A( I, J )*MUL
            END DO ! I
         END DO ! J

      ELSE IF( ITYPE.EQ.3 ) THEN

c        Upper Hessenberg matrix

         DO J = 1, N
            DO I = 1, MIN( J+1, M )
               A( I, J ) = A( I, J )*MUL
            END DO ! I
         END DO ! J

      ELSE IF( ITYPE.EQ.4 ) THEN

c        Lower half of a symmetric band matrix

         K3 = KL + 1
         K4 = N + 1
         DO J = 1, N
            DO I = 1, MIN( K3, K4-J )
               A( I, J ) = A( I, J )*MUL
            END DO ! I
         END DO ! J

      ELSE IF( ITYPE.EQ.5 ) THEN

c        Upper half of a symmetric band matrix

         K1 = KU + 2
         K3 = KU + 1
         DO J = 1, N
            DO I = MAX( K1-J, 1 ), K3
               A( I, J ) = A( I, J )*MUL
            END DO ! I
         END DO ! J

      ELSE IF( ITYPE.EQ.6 ) THEN

c        Band matrix

         K1 = KL + KU + 2
         K2 = KL + 1
         K3 = 2*KL + KU + 1
         K4 = KL + KU + 1 + M
         DO J = 1, N
            DO I = MAX( K1-J, K2 ), MIN( K3, K4-J )
               A( I, J ) = A( I, J )*MUL
            END DO ! I
         END DO ! J

      END IF

      IF( .NOT.DONE )
     $   GO TO 10

c     End of DLASCL

      END
