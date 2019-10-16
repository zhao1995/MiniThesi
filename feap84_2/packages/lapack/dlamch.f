      DOUBLE PRECISION FUNCTION DLAMCH( CMACH )

c  -- LAPACK auxiliary routine (version 2.0) --
c     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
c     Courant Institute, Argonne National Lab, and Rice University
c     October 31, 1992
      implicit   none

c     .. Scalar Arguments ..
      CHARACTER          CMACH
c     ..

c  Purpose
c  =======

c  DLAMCH determines double precision machine parameters.

c  Arguments
c  =========

c  CMACH   (input) CHARACTER*1
c          Specifies the value to be returned by DLAMCH:
c          = 'E' or 'e',   DLAMCH := eps
c          = 'S' or 's ,   DLAMCH := sfmin
c          = 'B' or 'b',   DLAMCH := base
c          = 'P' or 'p',   DLAMCH := eps*base
c          = 'N' or 'n',   DLAMCH := t
c          = 'R' or 'r',   DLAMCH := rnd
c          = 'M' or 'm',   DLAMCH := emin
c          = 'U' or 'u',   DLAMCH := rmin
c          = 'L' or 'l',   DLAMCH := emax
c          = 'O' or 'o',   DLAMCH := rmax

c          where

c          eps   = relative machine precision
c          sfmin = safe minimum, such that 1/sfmin does not overflow
c          base  = base of the machine
c          prec  = eps*base
c          t     = number of (base) digits in the mantissa
c          rnd   = 1.0 when rounding occurs in addition, 0.0 otherwise
c          emin  = minimum exponent before (gradual) underflow
c          rmin  = underflow threshold - base**(emin-1)
c          emax  = largest exponent before overflow
c          rmax  = overflow threshold  - (base**emax)*(1-eps)

c =====================================================================

c     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
c     ..
c     .. Local Scalars ..
      LOGICAL            FIRST, LRND
      INTEGER            BETA, IMAX, IMIN, IT
      DOUBLE PRECISION   BASE, EMAX, EMIN, EPS, PREC, RMACH, RMAX, RMIN,
     $                   RND, SFMIN, SMALL, T
c     ..
c     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
c     ..
c     .. External Subroutines ..
      EXTERNAL           DLAMC2
c     ..
c     .. Save statement ..
      SAVE               FIRST, EPS, SFMIN, BASE, T, RND, EMIN, RMIN,
     $                   EMAX, RMAX, PREC
c     ..
c     .. Data statements ..
      DATA               FIRST / .TRUE. /
c     ..
c     .. Executable Statements ..

      IF( FIRST ) THEN
         FIRST = .FALSE.
         CALL DLAMC2( BETA, IT, LRND, EPS, IMIN, RMIN, IMAX, RMAX )
         BASE = BETA
         T = IT
         IF( LRND ) THEN
            RND = ONE
            EPS = ( BASE**( 1-IT ) ) / 2
         ELSE
            RND = ZERO
            EPS = BASE**( 1-IT )
         END IF
         PREC = EPS*BASE
         EMIN = IMIN
         EMAX = IMAX
         SFMIN = RMIN
         SMALL = ONE / RMAX
         IF( SMALL.GE.SFMIN ) THEN

c           Use SMALL plus a bit, to avoid the possibility of rounding
c           causing overflow when computing  1/sfmin.

            SFMIN = SMALL*( ONE+EPS )
         END IF
      END IF

      IF( LSAME( CMACH, 'E' ) ) THEN
         RMACH = EPS
      ELSE IF( LSAME( CMACH, 'S' ) ) THEN
         RMACH = SFMIN
      ELSE IF( LSAME( CMACH, 'B' ) ) THEN
         RMACH = BASE
      ELSE IF( LSAME( CMACH, 'P' ) ) THEN
         RMACH = PREC
      ELSE IF( LSAME( CMACH, 'N' ) ) THEN
         RMACH = T
      ELSE IF( LSAME( CMACH, 'R' ) ) THEN
         RMACH = RND
      ELSE IF( LSAME( CMACH, 'M' ) ) THEN
         RMACH = EMIN
      ELSE IF( LSAME( CMACH, 'U' ) ) THEN
         RMACH = RMIN
      ELSE IF( LSAME( CMACH, 'L' ) ) THEN
         RMACH = EMAX
      ELSE IF( LSAME( CMACH, 'O' ) ) THEN
         RMACH = RMAX
      ELSE
         RMACH = 0.0d0
         CALL XERBLA( 'CMACH', 1 )
      END IF

      DLAMCH = RMACH

c     End of DLAMCH

      END

c***********************************************************************

      SUBROUTINE DLAMC1( BETA, T, RND, IEEE1 )

c  -- LAPACK auxiliary routine (version 2.0) --
c     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
c     Courant Institute, Argonne National Lab, and Rice University
c     October 31, 1992
      implicit   none

c     .. Scalar Arguments ..
      LOGICAL            IEEE1, RND
      INTEGER            BETA, T
c     ..

c  Purpose
c  =======

c  DLAMC1 determines the machine parameters given by BETA, T, RND, and
c  IEEE1.

c  Arguments
c  =========

c  BETA    (output) INTEGER
c          The base of the machine.

c  T       (output) INTEGER
c          The number of ( BETA ) digits in the mantissa.

c  RND     (output) LOGICAL
c          Specifies whether proper rounding  ( RND = .TRUE. )  or
c          chopping  ( RND = .FALSE. )  occurs in addition. This may not
c          be a reliable guide to the way in which the machine performs
c          its arithmetic.

c  IEEE1   (output) LOGICAL
c          Specifies whether rounding appears to be done in the IEEE
c          'round to nearest' style.

c  Further Details
c  ===============

c  The routine is based on the routine  ENVRON  by Malcolm and
c  incorporates suggestions by Gentleman and Marovich. See

c     Malcolm M. A. (1972) Algorithms to reveal properties of
c        floating-point arithmetic. Comms. of the ACM, 15, 949-951.

c     Gentleman W. M. and Marovich S. B. (1974) More on algorithms
c        that reveal properties of floating point arithmetic units.
c        Comms. of the ACM, 17, 276-277.

c =====================================================================

c     .. Local Scalars ..
      LOGICAL            FIRST, LIEEE1, LRND
      INTEGER            LBETA, LT
      DOUBLE PRECISION   A, B, C, F, ONE, QTR, SAVEC, T1, T2
c     ..
c     .. External Functions ..
      DOUBLE PRECISION   DLAMC3
      EXTERNAL           DLAMC3
c     ..
c     .. Save statement ..
      SAVE               FIRST, LIEEE1, LBETA, LRND, LT
c     ..
c     .. Data statements ..
      DATA               FIRST / .TRUE. /
c     ..
c     .. Executable Statements ..

      IF( FIRST ) THEN
         FIRST = .FALSE.
         ONE = 1

c        LBETA,  LIEEE1,  LT and  LRND  are the  local values  of  BETA,
c        IEEE1, T and RND.

c        Throughout this routine  we use the function  DLAMC3  to ensure
c        that relevant values are  stored and not held in registers,  or
c        are not affected by optimizers.

c        Compute  a = 2.0**m  with the  smallest positive integer m such
c        that

c           fl( a + 1.0 ) = a.

         A = 1
         C = 1

c+       WHILE( C.EQ.ONE )LOOP
   10    CONTINUE
         IF( C.EQ.ONE ) THEN
            A = 2*A
            C = DLAMC3( A, ONE )
            C = DLAMC3( C, -A )
            GO TO 10
         END IF
c+       END WHILE

c        Now compute  b = 2.0**m  with the smallest positive integer m
c        such that

c           fl( a + b ) .gt. a.

         B = 1
         C = DLAMC3( A, B )

c+       WHILE( C.EQ.A )LOOP
   20    CONTINUE
         IF( C.EQ.A ) THEN
            B = 2*B
            C = DLAMC3( A, B )
            GO TO 20
         END IF
c+       END WHILE

c        Now compute the base.  a and c  are neighbouring floating point
c        numbers  in the  interval  ( beta**t, beta**( t + 1 ) )  and so
c        their difference is beta. Adding 0.25 to c is to ensure that it
c        is truncated to beta and not ( beta - 1 ).

         QTR = ONE / 4
         SAVEC = C
         C = DLAMC3( C, -A )
         LBETA = nint(C + QTR)

c        Now determine whether rounding or chopping occurs,  by adding a
c        bit  less  than  beta/2  and a  bit  more  than  beta/2  to  a.

         B = LBETA
         F = DLAMC3( B / 2, -B / 100 )
         C = DLAMC3( F, A )
         IF( C.EQ.A ) THEN
            LRND = .TRUE.
         ELSE
            LRND = .FALSE.
         END IF
         F = DLAMC3( B / 2, B / 100 )
         C = DLAMC3( F, A )
         IF( ( LRND ) .AND. ( C.EQ.A ) )
     $      LRND = .FALSE.

c        Try and decide whether rounding is done in the  IEEE  'round to
c        nearest' style. B/2 is half a unit in the last place of the two
c        numbers A and SAVEC. Furthermore, A is even, i.e. has last  bit
c        zero, and SAVEC is odd. Thus adding B/2 to A should not  change
c        A, but adding B/2 to SAVEC should change SAVEC.

         T1 = DLAMC3( B / 2, A )
         T2 = DLAMC3( B / 2, SAVEC )
         LIEEE1 = ( T1.EQ.A ) .AND. ( T2.GT.SAVEC ) .AND. LRND

c        Now find  the  mantissa, t.  It should  be the  integer part of
c        log to the base beta of a,  however it is safer to determine  t
c        by powering.  So we find t as the smallest positive integer for
c        which

c           fl( beta**t + 1.0 ) = 1.0.

         LT = 0
         A = 1
         C = 1

c+       WHILE( C.EQ.ONE )LOOP
   30    CONTINUE
         IF( C.EQ.ONE ) THEN
            LT = LT + 1
            A = A*LBETA
            C = DLAMC3( A, ONE )
            C = DLAMC3( C, -A )
            GO TO 30
         END IF
c+       END WHILE

      END IF

      BETA = LBETA
      T = LT
      RND = LRND
      IEEE1 = LIEEE1
      RETURN

c     End of DLAMC1

      END

c***********************************************************************

      SUBROUTINE DLAMC2( BETA, T, RND, EPS, EMIN, RMIN, EMAX, RMAX )

c  -- LAPACK auxiliary routine (version 2.0) --
c     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
c     Courant Institute, Argonne National Lab, and Rice University
c     October 31, 1992
      implicit   none

c     .. Scalar Arguments ..
      LOGICAL            RND
      INTEGER            BETA, EMAX, EMIN, T
      DOUBLE PRECISION   EPS, RMAX, RMIN
c     ..

c  Purpose
c  =======

c  DLAMC2 determines the machine parameters specified in its argument
c  list.

c  Arguments
c  =========

c  BETA    (output) INTEGER
c          The base of the machine.

c  T       (output) INTEGER
c          The number of ( BETA ) digits in the mantissa.

c  RND     (output) LOGICAL
c          Specifies whether proper rounding  ( RND = .TRUE. )  or
c          chopping  ( RND = .FALSE. )  occurs in addition. This may not
c          be a reliable guide to the way in which the machine performs
c          its arithmetic.

c  EPS     (output) DOUBLE PRECISION
c          The smallest positive number such that

c             fl( 1.0 - EPS ) .LT. 1.0,

c          where fl denotes the computed value.

c  EMIN    (output) INTEGER
c          The minimum exponent before (gradual) underflow occurs.

c  RMIN    (output) DOUBLE PRECISION
c          The smallest normalized number for the machine, given by
c          BASE**( EMIN - 1 ), where  BASE  is the floating point value
c          of BETA.

c  EMAX    (output) INTEGER
c          The maximum exponent before overflow occurs.

c  RMAX    (output) DOUBLE PRECISION
c          The largest positive number for the machine, given by
c          BASE**EMAX * ( 1 - EPS ), where  BASE  is the floating point
c          value of BETA.

c  Further Details
c  ===============

c  The computation of  EPS  is based on a routine PARANOIA by
c  W. Kahan of the University of California at Berkeley.

c =====================================================================

c     .. Local Scalars ..
      LOGICAL            FIRST, IEEE, IWARN, LIEEE1, LRND
      INTEGER            GNMIN, GPMIN, I, LBETA, LEMAX, LEMIN, LT,
     $                   NGNMIN, NGPMIN
      DOUBLE PRECISION   A, B, C, HALF, LEPS, LRMAX, LRMIN, ONE, RBASE,
     $                   SIXTH, SMALL, THIRD, TWO, ZERO
c     ..
c     .. External Functions ..
      DOUBLE PRECISION   DLAMC3
      EXTERNAL           DLAMC3
c     ..
c     .. External Subroutines ..
      EXTERNAL           DLAMC1, DLAMC4, DLAMC5
c     ..
c     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, MIN
c     ..
c     .. Save statement ..
      SAVE               FIRST, IWARN, LBETA, LEMAX, LEMIN, LEPS, LRMAX,
     $                   LRMIN, LT
c     ..
c     .. Data statements ..
      DATA               FIRST / .TRUE. / , IWARN / .FALSE. /
c     ..
c     .. Executable Statements ..

      IF( FIRST ) THEN
         FIRST = .FALSE.
         ZERO = 0
         ONE = 1
         TWO = 2

c        LBETA, LT, LRND, LEPS, LEMIN and LRMIN  are the local values of
c        BETA, T, RND, EPS, EMIN and RMIN.

c        Throughout this routine  we use the function  DLAMC3  to ensure
c        that relevant values are stored  and not held in registers,  or
c        are not affected by optimizers.

c        DLAMC1 returns the parameters  LBETA, LT, LRND and LIEEE1.

         CALL DLAMC1( LBETA, LT, LRND, LIEEE1 )

c        Start to find EPS.

         B = LBETA
         A = B**( -LT )
         LEPS = A

c        Try some tricks to see whether or not this is the correct  EPS.

         B = TWO / 3
         HALF = ONE / 2
         SIXTH = DLAMC3( B, -HALF )
         THIRD = DLAMC3( SIXTH, SIXTH )
         B = DLAMC3( THIRD, -HALF )
         B = DLAMC3( B, SIXTH )
         B = ABS( B )
         IF( B.LT.LEPS )
     $      B = LEPS

         LEPS = 1

c+       WHILE( ( LEPS.GT.B ).AND.( B.GT.ZERO ) )LOOP
   10    CONTINUE
         IF( ( LEPS.GT.B ) .AND. ( B.GT.ZERO ) ) THEN
            LEPS = B
            C = DLAMC3( HALF*LEPS, ( TWO**5 )*( LEPS**2 ) )
            C = DLAMC3( HALF, -C )
            B = DLAMC3( HALF, C )
            C = DLAMC3( HALF, -B )
            B = DLAMC3( HALF, C )
            GO TO 10
         END IF
c+       END WHILE

         IF( A.LT.LEPS )
     $      LEPS = A

c        Computation of EPS complete.

c        Now find  EMIN.  Let A = + or - 1, and + or - (1 + BASE**(-3)).
c        Keep dividing  A by BETA until (gradual) underflow occurs. This
c        is detected when we cannot recover the previous A.

         RBASE = ONE / LBETA
         SMALL = ONE
         DO I = 1, 3
            SMALL = DLAMC3( SMALL*RBASE, ZERO )
         END DO ! I
         A = DLAMC3( ONE, SMALL )
         CALL DLAMC4( NGPMIN, ONE, LBETA )
         CALL DLAMC4( NGNMIN, -ONE, LBETA )
         CALL DLAMC4( GPMIN, A, LBETA )
         CALL DLAMC4( GNMIN, -A, LBETA )
         IEEE = .FALSE.

         IF( ( NGPMIN.EQ.NGNMIN ) .AND. ( GPMIN.EQ.GNMIN ) ) THEN
            IF( NGPMIN.EQ.GPMIN ) THEN
               LEMIN = NGPMIN
c            ( Non twos-complement machines, no gradual underflow;
c              e.g.,  VAX )
            ELSE IF( ( GPMIN-NGPMIN ).EQ.3 ) THEN
               LEMIN = NGPMIN - 1 + LT
               IEEE = .TRUE.
c            ( Non twos-complement machines, with gradual underflow;
c              e.g., IEEE standard followers )
            ELSE
               LEMIN = MIN( NGPMIN, GPMIN )
c            ( A guess; no known machine )
               IWARN = .TRUE.
            END IF

         ELSE IF( ( NGPMIN.EQ.GPMIN ) .AND. ( NGNMIN.EQ.GNMIN ) ) THEN
            IF( ABS( NGPMIN-NGNMIN ).EQ.1 ) THEN
               LEMIN = MAX( NGPMIN, NGNMIN )
c            ( Twos-complement machines, no gradual underflow;
c              e.g., CYBER 205 )
            ELSE
               LEMIN = MIN( NGPMIN, NGNMIN )
c            ( A guess; no known machine )
               IWARN = .TRUE.
            END IF

         ELSE IF( ( ABS( NGPMIN-NGNMIN ).EQ.1 ) .AND.
     $            ( GPMIN.EQ.GNMIN ) ) THEN
            IF( ( GPMIN-MIN( NGPMIN, NGNMIN ) ).EQ.3 ) THEN
               LEMIN = MAX( NGPMIN, NGNMIN ) - 1 + LT
c            ( Twos-complement machines with gradual underflow;
c              no known machine )
            ELSE
               LEMIN = MIN( NGPMIN, NGNMIN )
c            ( A guess; no known machine )
               IWARN = .TRUE.
            END IF

         ELSE
            LEMIN = MIN( NGPMIN, NGNMIN, GPMIN, GNMIN )
c         ( A guess; no known machine )
            IWARN = .TRUE.
         END IF
c**
c Comment out this if block if EMIN is ok
         IF( IWARN ) THEN
            FIRST = .TRUE.
            WRITE( 6, FMT = 9999 )LEMIN
         END IF
c**

c        Assume IEEE arithmetic if we found denormalised  numbers above,
c        or if arithmetic seems to round in the  IEEE style,  determined
c        in routine DLAMC1. A true IEEE machine should have both  things
c        true; however, faulty machines may have one or the other.

         IEEE = IEEE .OR. LIEEE1

c        Compute  RMIN by successive division by  BETA. We could compute
c        RMIN as BASE**( EMIN - 1 ),  but some machines underflow during
c        this computation.

         LRMIN = 1
         DO I = 1, 1 - LEMIN
            LRMIN = DLAMC3( LRMIN*RBASE, ZERO )
         END DO ! I

c        Finally, call DLAMC5 to compute EMAX and RMAX.

         CALL DLAMC5( LBETA, LT, LEMIN, IEEE, LEMAX, LRMAX )
      END IF

      BETA = LBETA
      T = LT
      RND = LRND
      EPS = LEPS
      EMIN = LEMIN
      RMIN = LRMIN
      EMAX = LEMAX
      RMAX = LRMAX

      RETURN

 9999 FORMAT( / / ' WARNING. The value EMIN may be incorrect:-',
     $      '  EMIN = ', I8, /
     $      ' If, after inspection, the value EMIN looks',
     $      ' acceptable please comment out ',
     $      / ' the IF block as marked within the code of routine',
     $      ' DLAMC2,', / ' otherwise supply EMIN explicitly.', / )

c     End of DLAMC2

      END

c***********************************************************************

      DOUBLE PRECISION FUNCTION DLAMC3( A, B )

c  -- LAPACK auxiliary routine (version 2.0) --
c     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
c     Courant Institute, Argonne National Lab, and Rice University
c     October 31, 1992
      implicit   none

c     .. Scalar Arguments ..
      DOUBLE PRECISION   A, B
c     ..

c  Purpose
c  =======

c  DLAMC3  is intended to force  A  and  B  to be stored prior to doing
c  the addition of  A  and  B ,  for use in situations where optimizers
c  might hold one of these in a register.

c  Arguments
c  =========

c  A, B    (input) DOUBLE PRECISION
c          The values A and B.

c =====================================================================

c     .. Executable Statements ..

      DLAMC3 = A + B

      RETURN

c     End of DLAMC3

      END

c***********************************************************************

      SUBROUTINE DLAMC4( EMIN, START, BASE )

c  -- LAPACK auxiliary routine (version 2.0) --
c     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
c     Courant Institute, Argonne National Lab, and Rice University
c     October 31, 1992
      implicit   none

c     .. Scalar Arguments ..
      INTEGER            BASE, EMIN
      DOUBLE PRECISION   START
c     ..

c  Purpose
c  =======

c  DLAMC4 is a service routine for DLAMC2.

c  Arguments
c  =========

c  EMIN    (output) EMIN
c          The minimum exponent before (gradual) underflow, computed by
c          setting A = START and dividing by BASE until the previous A
c          can not be recovered.

c  START   (input) DOUBLE PRECISION
c          The starting point for determining EMIN.

c  BASE    (input) INTEGER
c          The base of the machine.

c =====================================================================

c     .. Local Scalars ..
      INTEGER            I
      DOUBLE PRECISION   A, B1, B2, C1, C2, D1, D2, ONE, RBASE, ZERO
c     ..
c     .. External Functions ..
      DOUBLE PRECISION   DLAMC3
      EXTERNAL           DLAMC3
c     ..
c     .. Executable Statements ..

      A = START
      ONE = 1
      RBASE = ONE / BASE
      ZERO = 0
      EMIN = 1
      B1 = DLAMC3( A*RBASE, ZERO )
      C1 = A
      C2 = A
      D1 = A
      D2 = A
c+    WHILE( ( C1.EQ.A ).AND.( C2.EQ.A ).AND.
c    $       ( D1.EQ.A ).AND.( D2.EQ.A )      )LOOP
   10 CONTINUE
      IF( ( C1.EQ.A ) .AND. ( C2.EQ.A ) .AND. ( D1.EQ.A ) .AND.
     $    ( D2.EQ.A ) ) THEN
         EMIN = EMIN - 1
         A = B1
         B1 = DLAMC3( A / BASE, ZERO )
         C1 = DLAMC3( B1*BASE, ZERO )
         D1 = ZERO
         DO I = 1, BASE
            D1 = D1 + B1
         END DO ! I
         B2 = DLAMC3( A*RBASE, ZERO )
         C2 = DLAMC3( B2 / RBASE, ZERO )
         D2 = ZERO
         DO I = 1, BASE
            D2 = D2 + B2
         END DO ! I
         GO TO 10
      END IF
c+    END WHILE

      RETURN

c     End of DLAMC4

      END

c***********************************************************************

      SUBROUTINE DLAMC5( BETA, P, EMIN, IEEE, EMAX, RMAX )

c  -- LAPACK auxiliary routine (version 2.0) --
c     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
c     Courant Institute, Argonne National Lab, and Rice University
c     October 31, 1992
      implicit   none

c     .. Scalar Arguments ..
      LOGICAL            IEEE
      INTEGER            BETA, EMAX, EMIN, P
      DOUBLE PRECISION   RMAX
c     ..

c  Purpose
c  =======

c  DLAMC5 attempts to compute RMAX, the largest machine floating-point
c  number, without overflow.  It assumes that EMAX + abs(EMIN) sum
c  approximately to a power of 2.  It will fail on machines where this
c  assumption does not hold, for example, the Cyber 205 (EMIN = -28625,
c  EMAX = 28718).  It will also fail if the value supplied for EMIN is
c  too large (i.e. too close to zero), probably with overflow.

c  Arguments
c  =========

c  BETA    (input) INTEGER
c          The base of floating-point arithmetic.

c  P       (input) INTEGER
c          The number of base BETA digits in the mantissa of a
c          floating-point value.

c  EMIN    (input) INTEGER
c          The minimum exponent before (gradual) underflow.

c  IEEE    (input) LOGICAL
c          A logical flag specifying whether or not the arithmetic
c          system is thought to comply with the IEEE standard.

c  EMAX    (output) INTEGER
c          The largest exponent before overflow

c  RMAX    (output) DOUBLE PRECISION
c          The largest machine floating-point number.

c =====================================================================

c     .. Parameters ..
      DOUBLE PRECISION   ZERO        , ONE
      PARAMETER        ( ZERO = 0.0D0, ONE = 1.0D0 )
c     ..
c     .. Local Scalars ..
      INTEGER            EXBITS, EXPSUM, I, LEXP, NBITS, TRY, UEXP
      DOUBLE PRECISION   OLDY, RECBAS, Y, Z
c     ..
c     .. External Functions ..
      DOUBLE PRECISION   DLAMC3
      EXTERNAL           DLAMC3
c     ..
c     .. Intrinsic Functions ..
      INTRINSIC          MOD
c     ..
c     .. Executable Statements ..

c     First compute LEXP and UEXP, two powers of 2 that bound
c     abs(EMIN). We then assume that EMAX + abs(EMIN) will sum
c     approximately to the bound that is closest to abs(EMIN).
c     (EMAX is the exponent of the required number RMAX).

      LEXP = 1
      EXBITS = 1
   10 CONTINUE
      TRY = LEXP*2
      IF( TRY.LE.( -EMIN ) ) THEN
         LEXP   = TRY
         EXBITS = EXBITS + 1
         GO TO 10
      END IF
      IF( LEXP.EQ.-EMIN ) THEN
         UEXP = LEXP
      ELSE
         UEXP = TRY
         EXBITS = EXBITS + 1
      END IF

c     Now -LEXP is less than or equal to EMIN, and -UEXP is greater
c     than or equal to EMIN. EXBITS is the number of bits needed to
c     store the exponent.

      IF( ( UEXP+EMIN ).GT.( -LEXP-EMIN ) ) THEN
         EXPSUM = 2*LEXP
      ELSE
         EXPSUM = 2*UEXP
      END IF

c     EXPSUM is the exponent range, approximately equal to
c     EMAX - EMIN + 1 .

      EMAX  = EXPSUM + EMIN - 1
      NBITS = 1 + EXBITS + P

c     NBITS is the total number of bits needed to store a
c     floating-point number.

      IF( ( MOD( NBITS, 2 ).EQ.1 ) .AND. ( BETA.EQ.2 ) ) THEN

c        Either there are an odd number of bits used to store a
c        floating-point number, which is unlikely, or some bits are
c        not used in the representation of numbers, which is possible,
c        (e.g. Cray machines) or the mantissa has an implicit bit,
c        (e.g. IEEE machines, Dec Vax machines), which is perhaps the
c        most likely. We have to assume the last alternative.
c        If this is true, then we need to reduce EMAX by one because
c        there must be some way of representing zero in an implicit-bit
c        system. On machines like Cray, we are reducing EMAX by one
c        unnecessarily.

         EMAX = EMAX - 1
      END IF

      IF( IEEE ) THEN

c        Assume we are on an IEEE machine which reserves one exponent
c        for infinity and NaN.

         EMAX = EMAX - 1
      END IF

c     Now create RMAX, the largest machine number, which should
c     be equal to (1.0 - BETA**(-P)) * BETA**EMAX .

c     First compute 1.0 - BETA**(-P), being careful that the
c     result is less than 1.0 .

      RECBAS = ONE / BETA
      Z      = BETA - ONE
      Y      = ZERO
      DO I = 1, P
         Z = Z*RECBAS
         IF( Y.LT.ONE )
     $      OLDY = Y
         Y = DLAMC3( Y, Z )
      END DO ! I
      IF( Y.GE.ONE )
     $   Y = OLDY

c     Now multiply by BETA**EMAX to get RMAX.

      DO I = 1, EMAX
         Y = DLAMC3( Y*BETA, ZERO )
      END DO ! I

      RMAX = Y

c     End of DLAMC5

      END
