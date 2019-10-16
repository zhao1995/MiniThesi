      SUBROUTINE DLASSQ( N, X, INCX, SCALE, SUMSQ )

c  -- LAPACK auxiliary routine (version 2.0) --
c     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
c     Courant Institute, Argonne National Lab, and Rice University
c     October 31, 1992
      implicit   none

c     .. Scalar Arguments ..
      INTEGER            INCX, N
      DOUBLE PRECISION   SCALE, SUMSQ
c     ..
c     .. Array Arguments ..
      DOUBLE PRECISION   X( * )
c     ..

c  Purpose
c  =======

c  DLASSQ  returns the values  scl  and  smsq  such that

c     ( scl**2 )*smsq = x( 1 )**2 +...+ x( n )**2 + ( scale**2 )*sumsq,

c  where  x( i ) = X( 1 + ( i - 1 )*INCX ). The value of  sumsq  is
c  assumed to be non-negative and  scl  returns the value

c     scl = max( scale, abs( x( i ) ) ).

c  scale and sumsq must be supplied in SCALE and SUMSQ and
c  scl and smsq are overwritten on SCALE and SUMSQ respectively.

c  The routine makes only one pass through the vector x.

c  Arguments
c  =========

c  N       (input) INTEGER
c          The number of elements to be used from the vector X.

c  X       (input) DOUBLE PRECISION
c          The vector for which a scaled sum of squares is computed.
c             x( i )  = X( 1 + ( i - 1 )*INCX ), 1 <= i <= n.

c  INCX    (input) INTEGER
c          The increment between successive values of the vector X.
c          INCX > 0.

c  SCALE   (input/output) DOUBLE PRECISION
c          On entry, the value  scale  in the equation above.
c          On exit, SCALE is overwritten with  scl , the scaling factor
c          for the sum of squares.

c  SUMSQ   (input/output) DOUBLE PRECISION
c          On entry, the value  sumsq  in the equation above.
c          On exit, SUMSQ is overwritten with  smsq , the basic sum of
c          squares from which  scl  has been factored out.

c =====================================================================

c     .. Parameters ..
      DOUBLE PRECISION   ZERO
      PARAMETER        ( ZERO = 0.0D+0 )
c     ..
c     .. Local Scalars ..
      INTEGER            IX
      DOUBLE PRECISION   ABSXI
c     ..
c     .. Intrinsic Functions ..
      INTRINSIC          ABS
c     ..
c     .. Executable Statements ..

      IF( N.GT.0 ) THEN
         DO IX = 1, 1 + ( N-1 )*INCX, INCX
            IF( X( IX ).NE.ZERO ) THEN
               ABSXI = ABS( X( IX ) )
               IF( SCALE.LT.ABSXI ) THEN
                  SUMSQ = 1 + SUMSQ*( SCALE / ABSXI )**2
                  SCALE = ABSXI
               ELSE
                  SUMSQ = SUMSQ + ( ABSXI / SCALE )**2
               END IF
            END IF
         END DO ! IX
      END IF

c     End of DLASSQ

      END
