      SUBROUTINE DLAE2( A, B, C, RT1, RT2 )

c  -- LAPACK auxiliary routine (version 2.0) --
c     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
c     Courant Institute, Argonne National Lab, and Rice University
c     October 31, 1992
      implicit   none

c     .. Scalar Arguments ..
      DOUBLE PRECISION   A, B, C, RT1, RT2
c     ..

c  Purpose
c  =======

c  DLAE2  computes the eigenvalues of a 2-by-2 symmetric matrix
c     [  A   B  ]
c     [  B   C  ].
c  On return, RT1 is the eigenvalue of larger absolute value, and RT2
c  is the eigenvalue of smaller absolute value.

c  Arguments
c  =========

c  A       (input) DOUBLE PRECISION
c          The (1,1) element of the 2-by-2 matrix.

c  B       (input) DOUBLE PRECISION
c          The (1,2) and (2,1) elements of the 2-by-2 matrix.

c  C       (input) DOUBLE PRECISION
c          The (2,2) element of the 2-by-2 matrix.

c  RT1     (output) DOUBLE PRECISION
c          The eigenvalue of larger absolute value.

c  RT2     (output) DOUBLE PRECISION
c          The eigenvalue of smaller absolute value.

c  Further Details
c  ===============

c  RT1 is accurate to a few ulps barring over/underflow.

c  RT2 may be inaccurate if there is massive cancellation in the
c  determinant A*C-B*B; higher precision or correctly rounded or
c  correctly truncated arithmetic would be needed to compute RT2
c  accurately in all cases.

c  Overflow is possible only if RT1 is within a factor of 5 of overflow.
c  Underflow is harmless if the input data is 0 or exceeds
c     underflow_threshold / macheps.

c =====================================================================

c     .. Parameters ..
      DOUBLE PRECISION   ONE
      PARAMETER        ( ONE = 1.0D0 )
      DOUBLE PRECISION   TWO
      PARAMETER        ( TWO = 2.0D0 )
      DOUBLE PRECISION   ZERO
      PARAMETER        ( ZERO = 0.0D0 )
      DOUBLE PRECISION   HALF
      PARAMETER        ( HALF = 0.5D0 )
c     ..
c     .. Local Scalars ..
      DOUBLE PRECISION   AB, ACMN, ACMX, ADF, DF, RT, SM, TB
c     ..
c     .. Intrinsic Functions ..
      INTRINSIC          ABS, SQRT
c     ..
c     .. Executable Statements ..

c     Compute the eigenvalues

      SM = A + C
      DF = A - C
      ADF = ABS( DF )
      TB = B + B
      AB = ABS( TB )
      IF( ABS( A ).GT.ABS( C ) ) THEN
         ACMX = A
         ACMN = C
      ELSE
         ACMX = C
         ACMN = A
      END IF
      IF( ADF.GT.AB ) THEN
         RT = ADF*SQRT( ONE+( AB / ADF )**2 )
      ELSE IF( ADF.LT.AB ) THEN
         RT = AB*SQRT( ONE+( ADF / AB )**2 )
      ELSE

c        Includes case AB=ADF=0

         RT = AB*SQRT( TWO )
      END IF
      IF( SM.LT.ZERO ) THEN
         RT1 = HALF*( SM-RT )

c        Order of execution important.
c        To get fully accurate smaller eigenvalue,
c        next line needs to be executed in higher precision.

         RT2 = ( ACMX / RT1 )*ACMN - ( B / RT1 )*B
      ELSE IF( SM.GT.ZERO ) THEN
         RT1 = HALF*( SM+RT )

c        Order of execution important.
c        To get fully accurate smaller eigenvalue,
c        next line needs to be executed in higher precision.

         RT2 = ( ACMX / RT1 )*ACMN - ( B / RT1 )*B
      ELSE

c        Includes case RT1 = RT2 = 0

         RT1 = HALF*RT
         RT2 = -HALF*RT
      END IF

c     End of DLAE2

      END
