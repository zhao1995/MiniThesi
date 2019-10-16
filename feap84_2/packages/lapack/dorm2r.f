      SUBROUTINE DORM2R( SIDE, TRANS, M, N, K, A, LDA, TAU, C, LDC,
     $                   WORK, INFO )

c  -- LAPACK routine (version 2.0) --
c     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
c     Courant Institute, Argonne National Lab, and Rice University
c     February 29, 1992
      implicit   none

c     .. Scalar Arguments ..
      CHARACTER          SIDE, TRANS
      INTEGER            INFO, K, LDA, LDC, M, N
c     ..
c     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), C( LDC, * ), TAU( * ), WORK( * )
c     ..

c  Purpose
c  =======

c  DORM2R overwrites the general real m by n matrix C with

c        Q * C  if SIDE = 'L' and TRANS = 'N', or

c        Q'* C  if SIDE = 'L' and TRANS = 'T', or

c        C * Q  if SIDE = 'R' and TRANS = 'N', or

c        C * Q' if SIDE = 'R' and TRANS = 'T',

c  where Q is a real orthogonal matrix defined as the product of k
c  elementary reflectors

c        Q = H(1) H(2) . . . H(k)

c  as returned by DGEQRF. Q is of order m if SIDE = 'L' and of order n
c  if SIDE = 'R'.

c  Arguments
c  =========

c  SIDE    (input) CHARACTER*1
c          = 'L': apply Q or Q' from the Left
c          = 'R': apply Q or Q' from the Right

c  TRANS   (input) CHARACTER*1
c          = 'N': apply Q  (No transpose)
c          = 'T': apply Q' (Transpose)

c  M       (input) INTEGER
c          The number of rows of the matrix C. M >= 0.

c  N       (input) INTEGER
c          The number of columns of the matrix C. N >= 0.

c  K       (input) INTEGER
c          The number of elementary reflectors whose product defines
c          the matrix Q.
c          If SIDE = 'L', M >= K >= 0;
c          if SIDE = 'R', N >= K >= 0.

c  A       (input) DOUBLE PRECISION array, dimension (LDA,K)
c          The i-th column must contain the vector which defines the
c          elementary reflector H(i), for i = 1,2,...,k, as returned by
c          DGEQRF in the first k columns of its array argument A.
c          A is modified by the routine but restored on exit.

c  LDA     (input) INTEGER
c          The leading dimension of the array A.
c          If SIDE = 'L', LDA >= max(1,M);
c          if SIDE = 'R', LDA >= max(1,N).

c  TAU     (input) DOUBLE PRECISION array, dimension (K)
c          TAU(i) must contain the scalar factor of the elementary
c          reflector H(i), as returned by DGEQRF.

c  C       (input/output) DOUBLE PRECISION array, dimension (LDC,N)
c          On entry, the m by n matrix C.
c          On exit, C is overwritten by Q*C or Q'*C or C*Q' or C*Q.

c  LDC     (input) INTEGER
c          The leading dimension of the array C. LDC >= max(1,M).

c  WORK    (workspace) DOUBLE PRECISION array, dimension
c                                   (N) if SIDE = 'L',
c                                   (M) if SIDE = 'R'

c  INFO    (output) INTEGER
c          = 0: successful exit
c          < 0: if INFO = -i, the i-th argument had an illegal value

c  =====================================================================

c     .. Parameters ..
      DOUBLE PRECISION   ONE
      PARAMETER        ( ONE = 1.0D+0 )
c     ..
c     .. Local Scalars ..
      LOGICAL            LEFT, NOTRAN
      INTEGER            I, I1, I2, I3, IC, JC, MI, NI, NQ
      DOUBLE PRECISION   AII
c     ..
c     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
c     ..
c     .. External Subroutines ..
      EXTERNAL           DLARF, XERBLA
c     ..
c     .. Intrinsic Functions ..
      INTRINSIC          MAX
c     ..
c     .. Executable Statements ..

c     Test the input arguments

      INFO = 0
      LEFT = LSAME( SIDE, 'L' )
      NOTRAN = LSAME( TRANS, 'N' )

c     NQ is the order of Q

      IF( LEFT ) THEN
         NQ = M
      ELSE
         NQ = N
      END IF
      IF( .NOT.LEFT .AND. .NOT.LSAME( SIDE, 'R' ) ) THEN
         INFO = -1
      ELSE IF( .NOT.NOTRAN .AND. .NOT.LSAME( TRANS, 'T' ) ) THEN
         INFO = -2
      ELSE IF( M.LT.0 ) THEN
         INFO = -3
      ELSE IF( N.LT.0 ) THEN
         INFO = -4
      ELSE IF( K.LT.0 .OR. K.GT.NQ ) THEN
         INFO = -5
      ELSE IF( LDA.LT.MAX( 1, NQ ) ) THEN
         INFO = -7
      ELSE IF( LDC.LT.MAX( 1, M ) ) THEN
         INFO = -10
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DORM2R', -INFO )
         RETURN
      END IF

c     Quick return if possible

      IF( M.EQ.0 .OR. N.EQ.0 .OR. K.EQ.0 )
     $   RETURN

      IF( ( LEFT .AND. .NOT.NOTRAN ) .OR. ( .NOT.LEFT .AND. NOTRAN ) )
     $     THEN
         I1 = 1
         I2 = K
         I3 = 1
      ELSE
         I1 = K
         I2 = 1
         I3 = -1
      END IF

      IF( LEFT ) THEN
         NI = N
         JC = 1
      ELSE
         MI = M
         IC = 1
      END IF

      DO I = I1, I2, I3
         IF( LEFT ) THEN

c           H(i) is applied to C(i:m,1:n)

            MI = M - I + 1
            IC = I
         ELSE

c           H(i) is applied to C(1:m,i:n)

            NI = N - I + 1
            JC = I
         END IF

c        Apply H(i)

         AII = A( I, I )
         A( I, I ) = ONE
         CALL DLARF( SIDE, MI, NI, A( I, I ), 1, TAU( I ), C( IC, JC ),
     $               LDC, WORK )
         A( I, I ) = AII
      END DO ! I

c     End of DORM2R

      END
