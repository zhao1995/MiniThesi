      SUBROUTINE DLAEV2( A, B, C, RT1, RT2, CS1, SN1 )

c  -- LAPACK auxiliary routine (version 2.0) --
c     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
c     Courant Institute, Argonne National Lab, and Rice University
c     October 31, 1992
      implicit   none

c     .. Scalar Arguments ..
      DOUBLE PRECISION   A, B, C, CS1, RT1, RT2, SN1
c     ..

c  Purpose
c  =======

c  DLAEV2 computes the eigendecomposition of a 2-by-2 symmetric matrix
c     [  A   B  ]
c     [  B   C  ].
c  On return, RT1 is the eigenvalue of larger absolute value, RT2 is the
c  eigenvalue of smaller absolute value, and (CS1,SN1) is the unit right
c  eigenvector for RT1, giving the decomposition

c     [ CS1  SN1 ] [  A   B  ] [ CS1 -SN1 ]  =  [ RT1  0  ]
c     [-SN1  CS1 ] [  B   C  ] [ SN1  CS1 ]     [  0  RT2 ].

c  Arguments
c  =========

c  A       (input) DOUBLE PRECISION
c          The (1,1) element of the 2-by-2 matrix.

c  B       (input) DOUBLE PRECISION
c          The (1,2) element and the conjugate of the (2,1) element of
c          the 2-by-2 matrix.

c  C       (input) DOUBLE PRECISION
c          The (2,2) element of the 2-by-2 matrix.

c  RT1     (output) DOUBLE PRECISION
c          The eigenvalue of larger absolute value.

c  RT2     (output) DOUBLE PRECISION
c          The eigenvalue of smaller absolute value.

c  CS1     (output) DOUBLE PRECISION
c  SN1     (output) DOUBLE PRECISION
c          The vector (CS1, SN1) is a unit right eigenvector for RT1.

c  Further Details
c  ===============

c  RT1 is accurate to a few ulps barring over/underflow.

c  RT2 may be inaccurate if there is massive cancellation in the
c  determinant A*C-B*B; higher precision or correctly rounded or
c  correctly truncated arithmetic would be needed to compute RT2
c  accurately in all cases.

c  CS1 and SN1 are accurate to a few ulps barring over/underflow.

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
      INTEGER            SGN1, SGN2
      DOUBLE PRECISION   AB, ACMN, ACMX, ACS, ADF, CS, CT, DF, RT, SM,
     $                   TB, TN
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
         SGN1 = -1

c        Order of execution important.
c        To get fully accurate smaller eigenvalue,
c        next line needs to be executed in higher precision.

         RT2 = ( ACMX / RT1 )*ACMN - ( B / RT1 )*B
      ELSE IF( SM.GT.ZERO ) THEN
         RT1 = HALF*( SM+RT )
         SGN1 = 1

c        Order of execution important.
c        To get fully accurate smaller eigenvalue,
c        next line needs to be executed in higher precision.

         RT2 = ( ACMX / RT1 )*ACMN - ( B / RT1 )*B
      ELSE

c        Includes case RT1 = RT2 = 0

         RT1 = HALF*RT
         RT2 = -HALF*RT
         SGN1 = 1
      END IF

c     Compute the eigenvector

      IF( DF.GE.ZERO ) THEN
         CS = DF + RT
         SGN2 = 1
      ELSE
         CS = DF - RT
         SGN2 = -1
      END IF
      ACS = ABS( CS )
      IF( ACS.GT.AB ) THEN
         CT = -TB / CS
         SN1 = ONE / SQRT( ONE+CT*CT )
         CS1 = CT*SN1
      ELSE
         IF( AB.EQ.ZERO ) THEN
            CS1 = ONE
            SN1 = ZERO
         ELSE
            TN = -CS / TB
            CS1 = ONE / SQRT( ONE+TN*TN )
            SN1 = TN*CS1
         END IF
      END IF
      IF( SGN1.EQ.SGN2 ) THEN
         TN = CS1
         CS1 = -SN1
         SN1 = TN
      END IF
      RETURN

c     End of DLAEV2

      END
