      SUBROUTINE DLARFG( N, ALPHA, X, INCX, TAU )

c  -- LAPACK auxiliary routine (version 2.0) --
c     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
c     Courant Institute, Argonne National Lab, and Rice University
c     September 30, 1994
      implicit   none

c     .. Scalar Arguments ..
      INTEGER            INCX, N
      DOUBLE PRECISION   ALPHA, TAU
c     ..
c     .. Array Arguments ..
      DOUBLE PRECISION   X( * )
c     ..

c  Purpose
c  =======

c  DLARFG generates a real elementary reflector H of order n, such
c  that

c        H * ( alpha ) = ( beta ),   H' * H = I.
c            (   x   )   (   0  )

c  where alpha and beta are scalars, and x is an (n-1)-element real
c  vector. H is represented in the form

c        H = I - tau * ( 1 ) * ( 1 v' ) ,
c                      ( v )

c  where tau is a real scalar and v is a real (n-1)-element
c  vector.

c  If the elements of x are all zero, then tau = 0 and H is taken to be
c  the unit matrix.

c  Otherwise  1 <= tau <= 2.

c  Arguments
c  =========

c  N       (input) INTEGER
c          The order of the elementary reflector.

c  ALPHA   (input/output) DOUBLE PRECISION
c          On entry, the value alpha.
c          On exit, it is overwritten with the value beta.

c  X       (input/output) DOUBLE PRECISION array, dimension
c                         (1+(N-2)*abs(INCX))
c          On entry, the vector x.
c          On exit, it is overwritten with the vector v.

c  INCX    (input) INTEGER
c          The increment between elements of X. INCX > 0.

c  TAU     (output) DOUBLE PRECISION
c          The value tau.

c  =====================================================================

c     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
c     ..
c     .. Local Scalars ..
      INTEGER            J, KNT
      DOUBLE PRECISION   BETA, RSAFMN, SAFMIN, XNORM
c     ..
c     .. External Functions ..
      DOUBLE PRECISION   DLAMCH, DLAPY2, DNRM2
      EXTERNAL           DLAMCH, DLAPY2, DNRM2
c     ..
c     .. Intrinsic Functions ..
      INTRINSIC          ABS, SIGN
c     ..
c     .. External Subroutines ..
      EXTERNAL           DSCAL
c     ..
c     .. Executable Statements ..

      IF( N.LE.1 ) THEN
         TAU = ZERO
         RETURN
      END IF

      XNORM = DNRM2( N-1, X, INCX )

      IF( XNORM.EQ.ZERO ) THEN

c        H  =  I

         TAU = ZERO
      ELSE

c        general case

         BETA = -SIGN( DLAPY2( ALPHA, XNORM ), ALPHA )
         SAFMIN = DLAMCH( 'S' ) / DLAMCH( 'E' )
         IF( ABS( BETA ).LT.SAFMIN ) THEN

c           XNORM, BETA may be inaccurate; scale X and recompute them

            RSAFMN = ONE / SAFMIN
            KNT = 0
   10       CONTINUE
            KNT = KNT + 1
            CALL DSCAL( N-1, RSAFMN, X, INCX )
            BETA = BETA*RSAFMN
            ALPHA = ALPHA*RSAFMN
            IF( ABS( BETA ).LT.SAFMIN )
     $         GO TO 10

c           New BETA is at most 1, at least SAFMIN

            XNORM = DNRM2( N-1, X, INCX )
            BETA = -SIGN( DLAPY2( ALPHA, XNORM ), ALPHA )
            TAU = ( BETA-ALPHA ) / BETA
            CALL DSCAL( N-1, ONE / ( ALPHA-BETA ), X, INCX )

c           If ALPHA is subnormal, it may lose relative accuracy

            ALPHA = BETA
            DO J = 1, KNT
               ALPHA = ALPHA*SAFMIN
            END DO ! J
         ELSE
            TAU = ( BETA-ALPHA ) / BETA
            CALL DSCAL( N-1, ONE / ( ALPHA-BETA ), X, INCX )
            ALPHA = BETA
         END IF
      END IF

c     End of DLARFG

      END
