      SUBROUTINE XERBLA( SRNAME, INFO )

c  -- LAPACK auxiliary routine (version 2.0) --
c     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
c     Courant Institute, Argonne National Lab, and Rice University
c     September 30, 1994
      implicit   none

c     .. Scalar Arguments ..
      CHARACTER*6        SRNAME
      INTEGER            INFO
c     ..

c  Purpose
c  =======

c  XERBLA  is an error handler for the LAPACK routines.
c  It is called by an LAPACK routine if an input parameter has an
c  invalid value.  A message is printed and execution stops.

c  Installers may consider modifying the STOP statement in order to
c  call system-specific exception-handling facilities.

c  Arguments
c  =========

c  SRNAME  (input) CHARACTER*6
c          The name of the routine which called XERBLA.

c  INFO    (input) INTEGER
c          The position of the invalid parameter in the parameter list
c          of the calling routine.

c =====================================================================

c     .. Executable Statements ..

      WRITE( *, FMT = 9999 )SRNAME, INFO

      STOP

 9999 FORMAT( ' ** On entry to ', A6, ' parameter number ', I2, ' had ',
     $      'an illegal value' )

c     End of XERBLA

      END
