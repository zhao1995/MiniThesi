      SUBROUTINE DLARNV( IDIST, ISEED, N, X )

c  -- LAPACK auxiliary routine (version 2.0) --
c     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
c     Courant Institute, Argonne National Lab, and Rice University
c     September 30, 1994
      implicit   none

c     .. Scalar Arguments ..
      INTEGER            IDIST, N
c     ..
c     .. Array Arguments ..
      INTEGER            ISEED( 4 )
      DOUBLE PRECISION   X( * )
c     ..

c  Purpose
c  =======

c  DLARNV returns a vector of n random real numbers from a uniform or
c  normal distribution.

c  Arguments
c  =========

c  IDIST   (input) INTEGER
c          Specifies the distribution of the random numbers:
c          = 1:  uniform (0,1)
c          = 2:  uniform (-1,1)
c          = 3:  normal (0,1)

c  ISEED   (input/output) INTEGER array, dimension (4)
c          On entry, the seed of the random number generator; the array
c          elements must be between 0 and 4095, and ISEED(4) must be
c          odd.
c          On exit, the seed is updated.

c  N       (input) INTEGER
c          The number of random numbers to be generated.

c  X       (output) DOUBLE PRECISION array, dimension (N)
c          The generated random numbers.

c  Further Details
c  ===============

c  This routine calls the auxiliary routine DLARUV to generate random
c  real numbers from a uniform (0,1) distribution, in batches of up to
c  128 using vectorisable code. The Box-Muller method is used to
c  transform numbers from a uniform to a normal distribution.

c  =====================================================================

c     .. Parameters ..
      DOUBLE PRECISION   ONE, TWO
      PARAMETER        ( ONE = 1.0D+0, TWO = 2.0D+0 )
      INTEGER            LV
      PARAMETER        ( LV = 128 )
      DOUBLE PRECISION   TWOPI
      PARAMETER        ( TWOPI = 6.2831853071795864769252867663D+0 )
c     ..
c     .. Local Scalars ..
      INTEGER            I, IL, IL2, IV
c     ..
c     .. Local Arrays ..
      DOUBLE PRECISION   U( LV )
c     ..
c     .. Intrinsic Functions ..
      INTRINSIC          COS, LOG, MIN, SQRT
c     ..
c     .. External Subroutines ..
      EXTERNAL           DLARUV
c     ..
c     .. Executable Statements ..

      DO IV = 1, N, LV / 2
         IL = MIN( LV / 2, N-IV+1 )
         IF( IDIST.EQ.3 ) THEN
            IL2 = 2*IL
         ELSE
            IL2 = IL
         END IF

c        Call DLARUV to generate IL2 numbers from a uniform (0,1)
c        distribution (IL2 <= LV)

         CALL DLARUV( ISEED, IL2, U )

         IF( IDIST.EQ.1 ) THEN

c           Copy generated numbers

            DO I = 1, IL
               X( IV+I-1 ) = U( I )
            END DO ! I
         ELSE IF( IDIST.EQ.2 ) THEN

c           Convert generated numbers to uniform (-1,1) distribution

            DO I = 1, IL
               X( IV+I-1 ) = TWO*U( I ) - ONE
            END DO ! I
         ELSE IF( IDIST.EQ.3 ) THEN

c           Convert generated numbers to normal (0,1) distribution

            DO I = 1, IL
               X( IV+I-1 ) = SQRT( -TWO*LOG( U( 2*I-1 ) ) )*
     $                       COS( TWOPI*U( 2*I ) )
            END DO ! I
         END IF
      END DO ! IV

c     End of DLARNV

      END
