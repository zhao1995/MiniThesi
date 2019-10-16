      DOUBLE PRECISION FUNCTION DLAPY2( X, Y )

c  -- LAPACK auxiliary routine (version 2.0) --
c     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
c     Courant Institute, Argonne National Lab, and Rice University
c     October 31, 1992
      implicit   none

c     .. Scalar Arguments ..
      DOUBLE PRECISION   X, Y
c     ..

c  Purpose
c  =======

c  DLAPY2 returns sqrt(x**2+y**2), taking care not to cause unnecessary
c  overflow.

c  Arguments
c  =========

c  X       (input) DOUBLE PRECISION
c  Y       (input) DOUBLE PRECISION
c          X and Y specify the values x and y.

c  =====================================================================

c     .. Parameters ..
      DOUBLE PRECISION   ZERO
      PARAMETER        ( ZERO = 0.0D0 )
      DOUBLE PRECISION   ONE
      PARAMETER        ( ONE  = 1.0D0 )
c     ..
c     .. Local Scalars ..
      DOUBLE PRECISION   W, XABS, YABS, Z
c     ..
c     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, MIN, SQRT
c     ..
c     .. Executable Statements ..

      XABS = ABS( X )
      YABS = ABS( Y )
      W    = MAX( XABS, YABS )
      Z    = MIN( XABS, YABS )
      IF( Z.EQ.ZERO ) THEN
         DLAPY2 = W
      ELSE
         DLAPY2 = W*SQRT( ONE+( Z / W )**2 )
      END IF

c     End of DLAPY2

      END
