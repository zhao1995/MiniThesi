      SUBROUTINE DLARTG( F, G, CS, SN, R )

c  -- LAPACK auxiliary routine (version 2.0) --
c     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
c     Courant Institute, Argonne National Lab, and Rice University
c     September 30, 1994
      implicit   none

c     .. Scalar Arguments ..
      DOUBLE PRECISION   CS, F, G, R, SN
c     ..

c  Purpose
c  =======

c  DLARTG generate a plane rotation so that

c     [  CS  SN  ]  .  [ F ]  =  [ R ]   where CS**2 + SN**2 = 1.
c     [ -SN  CS  ]     [ G ]     [ 0 ]

c  This is a slower, more accurate version of the BLAS1 routine DROTG,
c  with the following other differences:
c     F and G are unchanged on return.
c     If G=0, then CS=1 and SN=0.
c     If F=0 and (G .ne. 0), then CS=0 and SN=1 without doing any
c        floating point operations (saves work in DBDSQR when
c        there are zeros on the diagonal).

c  If F exceeds G in magnitude, CS will be positive.

c  Arguments
c  =========

c  F       (input) DOUBLE PRECISION
c          The first component of vector to be rotated.

c  G       (input) DOUBLE PRECISION
c          The second component of vector to be rotated.

c  CS      (output) DOUBLE PRECISION
c          The cosine of the rotation.

c  SN      (output) DOUBLE PRECISION
c          The sine of the rotation.

c  R       (output) DOUBLE PRECISION
c          The nonzero component of the rotated vector.

c  =====================================================================

c     .. Parameters ..
      DOUBLE PRECISION   ZERO
      PARAMETER        ( ZERO = 0.0D0 )
      DOUBLE PRECISION   ONE
      PARAMETER        ( ONE  = 1.0D0 )
      DOUBLE PRECISION   TWO
      PARAMETER        ( TWO  = 2.0D0 )
c     ..
c     .. Local Scalars ..
      LOGICAL            FIRST
      INTEGER            COUNT, I
      DOUBLE PRECISION   EPS, F1, G1, SAFMIN, SAFMN2, SAFMX2, SCALE
c     ..
c     .. External Functions ..
      DOUBLE PRECISION   DLAMCH
      EXTERNAL           DLAMCH
c     ..
c     .. Intrinsic Functions ..
      INTRINSIC          ABS, INT, LOG, MAX, SQRT
c     ..
c     .. Save statement ..
      SAVE               FIRST, SAFMX2, SAFMIN, SAFMN2
c     ..
c     .. Data statements ..
      DATA               FIRST / .TRUE. /
c     ..
c     .. Executable Statements ..

      IF( FIRST ) THEN
         FIRST  = .FALSE.
         SAFMIN = DLAMCH( 'S' )
         EPS    = DLAMCH( 'E' )
         SAFMN2 = DLAMCH( 'B' )**INT( LOG( SAFMIN / EPS ) /
     $            LOG( DLAMCH( 'B' ) ) / TWO )
         SAFMX2 = ONE / SAFMN2
      END IF
      IF( G.EQ.ZERO ) THEN
         CS = ONE
         SN = ZERO
         R  = F
      ELSE IF( F.EQ.ZERO ) THEN
         CS = ZERO
         SN = ONE
         R  = G
      ELSE
         F1 = F
         G1 = G
         SCALE = MAX( ABS( F1 ), ABS( G1 ) )
         IF( SCALE.GE.SAFMX2 ) THEN
            COUNT = 0
   10       CONTINUE
            COUNT = COUNT + 1
            F1    = F1*SAFMN2
            G1    = G1*SAFMN2
            SCALE = MAX( ABS( F1 ), ABS( G1 ) )
            IF( SCALE.GE.SAFMX2 )
     $         GO TO 10
            R  = SQRT( F1**2+G1**2 )
            CS = F1 / R
            SN = G1 / R
            DO I = 1, COUNT
               R = R*SAFMX2
            END DO ! I
         ELSE IF( SCALE.LE.SAFMN2 ) THEN
            COUNT = 0
   30       CONTINUE
            COUNT = COUNT + 1
            F1    = F1*SAFMX2
            G1    = G1*SAFMX2
            SCALE = MAX( ABS( F1 ), ABS( G1 ) )
            IF( SCALE.LE.SAFMN2 )
     $         GO TO 30
            R  = SQRT( F1**2+G1**2 )
            CS = F1 / R
            SN = G1 / R
            DO I = 1, COUNT
               R = R*SAFMN2
            END DO ! I
         ELSE
            R  = SQRT( F1**2+G1**2 )
            CS = F1 / R
            SN = G1 / R
         END IF
         IF( ABS( F ).GT.ABS( G ) .AND. CS.LT.ZERO ) THEN
            CS = -CS
            SN = -SN
            R  = -R
         END IF
      END IF

c     End of DLARTG

      END
