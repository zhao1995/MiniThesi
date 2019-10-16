      SUBROUTINE DLARF( SIDE, M, N, V, INCV, TAU, C, LDC, WORK )

c  -- LAPACK auxiliary routine (version 2.0) --
c     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
c     Courant Institute, Argonne National Lab, and Rice University
c     February 29, 1992
      implicit   none

c     .. Scalar Arguments ..
      CHARACTER          SIDE
      INTEGER            INCV, LDC, M, N
      DOUBLE PRECISION   TAU
c     ..
c     .. Array Arguments ..
      DOUBLE PRECISION   C( LDC, * ), V( * ), WORK( * )
c     ..

c  Purpose
c  =======

c  DLARF applies a real elementary reflector H to a real m by n matrix
c  C, from either the left or the right. H is represented in the form

c        H = I - tau * v * v'

c  where tau is a real scalar and v is a real vector.

c  If tau = 0, then H is taken to be the unit matrix.

c  Arguments
c  =========

c  SIDE    (input) CHARACTER*1
c          = 'L': form  H * C
c          = 'R': form  C * H

c  M       (input) INTEGER
c          The number of rows of the matrix C.

c  N       (input) INTEGER
c          The number of columns of the matrix C.

c  V       (input) DOUBLE PRECISION array, dimension
c                     (1 + (M-1)*abs(INCV)) if SIDE = 'L'
c                  or (1 + (N-1)*abs(INCV)) if SIDE = 'R'
c          The vector v in the representation of H. V is not used if
c          TAU = 0.

c  INCV    (input) INTEGER
c          The increment between elements of v. INCV <> 0.

c  TAU     (input) DOUBLE PRECISION
c          The value tau in the representation of H.

c  C       (input/output) DOUBLE PRECISION array, dimension (LDC,N)
c          On entry, the m by n matrix C.
c          On exit, C is overwritten by the matrix H * C if SIDE = 'L',
c          or C * H if SIDE = 'R'.

c  LDC     (input) INTEGER
c          The leading dimension of the array C. LDC >= max(1,M).

c  WORK    (workspace) DOUBLE PRECISION array, dimension
c                         (N) if SIDE = 'L'
c                      or (M) if SIDE = 'R'

c  =====================================================================

c     .. Parameters ..
      DOUBLE PRECISION   ONE         , ZERO
      PARAMETER        ( ONE = 1.0D+0, ZERO = 0.0D+0 )
c     ..
c     .. External Subroutines ..
      EXTERNAL           DGEMV, DGER
c     ..
c     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
c     ..
c     .. Executable Statements ..

      IF( LSAME( SIDE, 'L' ) ) THEN

c        Form  H * C

         IF( TAU.NE.ZERO ) THEN

c           w := C' * v

            CALL DGEMV( 'T', M, N, ONE, C, LDC, V, INCV, ZERO,
     $                  WORK, 1 )

c           C := C - v * w'

            CALL DGER( M, N, -TAU, V, INCV, WORK, 1, C, LDC )
         END IF
      ELSE

c        Form  C * H

         IF( TAU.NE.ZERO ) THEN

c           w := C * v

            CALL DGEMV( 'N', M, N, ONE, C, LDC, V, INCV,
     $                  ZERO, WORK, 1 )

c           C := C - w * v'

            CALL DGER( M, N, -TAU, WORK, 1, V, INCV, C, LDC )
         END IF
      END IF

c     End of DLARF

      END
