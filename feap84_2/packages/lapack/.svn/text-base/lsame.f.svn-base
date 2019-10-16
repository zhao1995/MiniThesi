      LOGICAL          FUNCTION LSAME( CA, CB )

c  -- LAPACK auxiliary routine (version 2.0) --
c     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
c     Courant Institute, Argonne National Lab, and Rice University
c     September 30, 1994
      implicit   none

c     .. Scalar Arguments ..
      CHARACTER          CA, CB
c     ..

c  Purpose
c  =======

c  LSAME returns .TRUE. if CA is the same letter as CB regardless of
c  case.

c  Arguments
c  =========

c  CA      (input) CHARACTER*1
c  CB      (input) CHARACTER*1
c          CA and CB specify the single characters to be compared.

c =====================================================================

c     .. Intrinsic Functions ..
      INTRINSIC          ICHAR
c     ..
c     .. Local Scalars ..
      INTEGER            INTA, INTB, ZCODE
c     ..
c     .. Executable Statements ..

c     Test if the characters are equal

      LSAME = CA.EQ.CB
      IF( LSAME )
     $   RETURN

c     Now test for equivalence if both characters are alphabetic.

      ZCODE = ICHAR( 'Z' )

c     Use 'Z' rather than 'A' so that ASCII can be detected on Prime
c     machines, on which ICHAR returns a value with bit 8 set.
c     ICHAR('A') on Prime machines returns 193 which is the same as
c     ICHAR('A') on an EBCDIC machine.

      INTA = ICHAR( CA )
      INTB = ICHAR( CB )

      IF( ZCODE.EQ.90 .OR. ZCODE.EQ.122 ) THEN

c        ASCII is assumed - ZCODE is the ASCII code of either lower or
c        upper case 'Z'.

         IF( INTA.GE.97 .AND. INTA.LE.122 ) INTA = INTA - 32
         IF( INTB.GE.97 .AND. INTB.LE.122 ) INTB = INTB - 32

      ELSE IF( ZCODE.EQ.233 .OR. ZCODE.EQ.169 ) THEN

c        EBCDIC is assumed - ZCODE is the EBCDIC code of either lower or
c        upper case 'Z'.

         IF( INTA.GE.129 .AND. INTA.LE.137 .OR.
     $       INTA.GE.145 .AND. INTA.LE.153 .OR.
     $       INTA.GE.162 .AND. INTA.LE.169 ) INTA = INTA + 64
         IF( INTB.GE.129 .AND. INTB.LE.137 .OR.
     $       INTB.GE.145 .AND. INTB.LE.153 .OR.
     $       INTB.GE.162 .AND. INTB.LE.169 ) INTB = INTB + 64

      ELSE IF( ZCODE.EQ.218 .OR. ZCODE.EQ.250 ) THEN

c        ASCII is assumed, on Prime machines - ZCODE is the ASCII code
c        plus 128 of either lower or upper case 'Z'.

         IF( INTA.GE.225 .AND. INTA.LE.250 ) INTA = INTA - 32
         IF( INTB.GE.225 .AND. INTB.LE.250 ) INTB = INTB - 32
      END IF
      LSAME = INTA.EQ.INTB

c     End of LSAME

      END
