      SUBROUTINE DSTEQR( COMPZ, N, D, E, Z, LDZ, WORK, INFO )

c  -- LAPACK routine (version 2.0) --
c     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
c     Courant Institute, Argonne National Lab, and Rice University
c     September 30, 1994
      implicit   none

c     .. Scalar Arguments ..
      CHARACTER          COMPZ
      INTEGER            INFO, LDZ, N
c     ..
c     .. Array Arguments ..
      DOUBLE PRECISION   D( * ), E( * ), WORK( * ), Z( LDZ, * )
c     ..

c  Purpose
c  =======

c  DSTEQR computes all eigenvalues and, optionally, eigenvectors of a
c  symmetric tridiagonal matrix using the implicit QL or QR method.
c  The eigenvectors of a full or band symmetric matrix can also be found
c  if DSYTRD or DSPTRD or DSBTRD has been used to reduce this matrix to
c  tridiagonal form.

c  Arguments
c  =========

c  COMPZ   (input) CHARACTER*1
c          = 'N':  Compute eigenvalues only.
c          = 'V':  Compute eigenvalues and eigenvectors of the original
c                  symmetric matrix.  On entry, Z must contain the
c                  orthogonal matrix used to reduce the original matrix
c                  to tridiagonal form.
c          = 'I':  Compute eigenvalues and eigenvectors of the
c                  tridiagonal matrix.  Z is initialized to the identity
c                  matrix.

c  N       (input) INTEGER
c          The order of the matrix.  N >= 0.

c  D       (input/output) DOUBLE PRECISION array, dimension (N)
c          On entry, the diagonal elements of the tridiagonal matrix.
c          On exit, if INFO = 0, the eigenvalues in ascending order.

c  E       (input/output) DOUBLE PRECISION array, dimension (N-1)
c          On entry, the (n-1) subdiagonal elements of the tridiagonal
c          matrix.
c          On exit, E has been destroyed.

c  Z       (input/output) DOUBLE PRECISION array, dimension (LDZ, N)
c          On entry, if  COMPZ = 'V', then Z contains the orthogonal
c          matrix used in the reduction to tridiagonal form.
c          On exit, if INFO = 0, then if  COMPZ = 'V', Z contains the
c          orthonormal eigenvectors of the original symmetric matrix,
c          and if COMPZ = 'I', Z contains the orthonormal eigenvectors
c          of the symmetric tridiagonal matrix.
c          If COMPZ = 'N', then Z is not referenced.

c  LDZ     (input) INTEGER
c          The leading dimension of the array Z.  LDZ >= 1, and if
c          eigenvectors are desired, then  LDZ >= max(1,N).

c  WORK    (workspace) DOUBLE PRECISION array, dimension (max(1,2*N-2))
c          If COMPZ = 'N', then WORK is not referenced.

c  INFO    (output) INTEGER
c          = 0:  successful exit
c          < 0:  if INFO = -i, the i-th argument had an illegal value
c          > 0:  the algorithm has failed to find all the eigenvalues in
c                a total of 30*N iterations; if INFO = i, then i
c                elements of E have not converged to zero; on exit, D
c                and E contain the elements of a symmetric tridiagonal
c                matrix which is orthogonally similar to the original
c                matrix.

c  =====================================================================

c     .. Parameters ..
      DOUBLE PRECISION   ZERO        , ONE        , TWO
      PARAMETER        ( ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0 )
      DOUBLE PRECISION   THREE
      PARAMETER        ( THREE = 3.0D0 )
      INTEGER            MAXIT
      PARAMETER        ( MAXIT = 30 )
c     ..
c     .. Local Scalars ..
      INTEGER            I, ICOMPZ, II, ISCALE, J, JTOT, K, L, L1, LEND,
     $                   LENDM1, LENDP1, LENDSV, LM1, LSV, M, MM, MM1,
     $                   NM1, NMAXIT
      DOUBLE PRECISION   ANORM, B, C, EPS, EPS2, F, G, P, R, RT1, RT2,
     $                   S, SAFMAX, SAFMIN, SSFMAX, SSFMIN, TST
c     ..
c     .. External Functions ..
      LOGICAL            LSAME
      DOUBLE PRECISION   DLAMCH, DLANST, DLAPY2
      EXTERNAL           LSAME, DLAMCH, DLANST, DLAPY2
c     ..
c     .. External Subroutines ..
      EXTERNAL           DLAE2, DLAEV2, DLARTG, DLASCL, DLASET, DLASR,
     $                   DLASRT, DSWAP, XERBLA
c     ..
c     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, SIGN, SQRT
c     ..
c     .. Executable Statements ..

c     Test the input parameters.

      INFO = 0

      IF( LSAME( COMPZ, 'N' ) ) THEN
         ICOMPZ = 0
      ELSE IF( LSAME( COMPZ, 'V' ) ) THEN
         ICOMPZ = 1
      ELSE IF( LSAME( COMPZ, 'I' ) ) THEN
         ICOMPZ = 2
      ELSE
         ICOMPZ = -1
      END IF
      IF( ICOMPZ.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( ( LDZ.LT.1 ) .OR. ( ICOMPZ.GT.0 .AND. LDZ.LT.MAX( 1,
     $         N ) ) ) THEN
         INFO = -6
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DSTEQR', -INFO )
         RETURN
      END IF

c     Quick return if possible

      IF( N.EQ.0 )
     $   RETURN

      IF( N.EQ.1 ) THEN
         IF( ICOMPZ.EQ.2 )
     $      Z( 1, 1 ) = ONE
         RETURN
      END IF

c     Determine the unit roundoff and over/underflow thresholds.

      EPS    = DLAMCH( 'E' )
      EPS2   = EPS**2
      SAFMIN = DLAMCH( 'S' )
      SAFMAX = ONE / SAFMIN
      SSFMAX = SQRT( SAFMAX ) / THREE
      SSFMIN = SQRT( SAFMIN ) / EPS2

c     Compute the eigenvalues and eigenvectors of the tridiagonal
c     matrix.

      IF( ICOMPZ.EQ.2 )
     $   CALL DLASET( 'Full', N, N, ZERO, ONE, Z, LDZ )

      NMAXIT = N*MAXIT
      JTOT = 0

c     Determine where the matrix splits and choose QL or QR iteration
c     for each block, according to whether top or bottom diagonal
c     element is smaller.

      L1 = 1
      NM1 = N - 1

   10 CONTINUE
      IF( L1.GT.N )
     $   GO TO 160
      IF( L1.GT.1 )
     $   E( L1-1 ) = ZERO
      IF( L1.LE.NM1 ) THEN
         DO M = L1, NM1
            TST = ABS( E( M ) )
            IF( TST.EQ.ZERO )
     $         GO TO 30
            IF( TST.LE.( SQRT( ABS( D( M ) ) )*SQRT( ABS( D( M+
     $          1 ) ) ) )*EPS ) THEN
               E( M ) = ZERO
               GO TO 30
            END IF
         END DO ! M
      END IF
      M = N

   30 CONTINUE
      L = L1
      LSV = L
      LEND = M
      LENDSV = LEND
      L1 = M + 1
      IF( LEND.EQ.L )
     $   GO TO 10

c     Scale submatrix in rows and columns L to LEND

      ANORM = DLANST( 'I', LEND-L+1, D( L ), E( L ) )
      ISCALE = 0
      IF( ANORM.EQ.ZERO )
     $   GO TO 10
      IF( ANORM.GT.SSFMAX ) THEN
         ISCALE = 1
         CALL DLASCL( 'G', 0, 0, ANORM, SSFMAX, LEND-L+1, 1, D( L ), N,
     $                INFO )
         CALL DLASCL( 'G', 0, 0, ANORM, SSFMAX, LEND-L, 1, E( L ), N,
     $                INFO )
      ELSE IF( ANORM.LT.SSFMIN ) THEN
         ISCALE = 2
         CALL DLASCL( 'G', 0, 0, ANORM, SSFMIN, LEND-L+1, 1, D( L ), N,
     $                INFO )
         CALL DLASCL( 'G', 0, 0, ANORM, SSFMIN, LEND-L, 1, E( L ), N,
     $                INFO )
      END IF

c     Choose between QL and QR iteration

      IF( ABS( D( LEND ) ).LT.ABS( D( L ) ) ) THEN
         LEND = LSV
         L = LENDSV
      END IF

      IF( LEND.GT.L ) THEN

c        QL Iteration

c        Look for small subdiagonal element.

   40    CONTINUE
         IF( L.NE.LEND ) THEN
            LENDM1 = LEND - 1
            DO M = L, LENDM1
               TST = ABS( E( M ) )**2
               IF( TST.LE.( EPS2*ABS( D( M ) ) )*ABS( D( M+1 ) )+
     $             SAFMIN )GO TO 60
            END DO ! M
         END IF

         M = LEND

   60    CONTINUE
         IF( M.LT.LEND )
     $      E( M ) = ZERO
         P = D( L )
         IF( M.EQ.L )
     $      GO TO 80

c        If remaining matrix is 2-by-2, use DLAE2 or SLAEV2
c        to compute its eigensystem.

         IF( M.EQ.L+1 ) THEN
            IF( ICOMPZ.GT.0 ) THEN
               CALL DLAEV2( D( L ), E( L ), D( L+1 ), RT1, RT2, C, S )
               WORK( L ) = C
               WORK( N-1+L ) = S
               CALL DLASR( 'R', 'V', 'B', N, 2, WORK( L ),
     $                     WORK( N-1+L ), Z( 1, L ), LDZ )
            ELSE
               CALL DLAE2( D( L ), E( L ), D( L+1 ), RT1, RT2 )
            END IF
            D( L ) = RT1
            D( L+1 ) = RT2
            E( L ) = ZERO
            L = L + 2
            IF( L.LE.LEND )
     $         GO TO 40
            GO TO 140
         END IF

         IF( JTOT.EQ.NMAXIT )
     $      GO TO 140
         JTOT = JTOT + 1

c        Form shift.

         G = ( D( L+1 )-P ) / ( TWO*E( L ) )
         R = DLAPY2( G, ONE )
         G = D( M ) - P + ( E( L ) / ( G+SIGN( R, G ) ) )

         S = ONE
         C = ONE
         P = ZERO

c        Inner loop

         MM1 = M - 1
         DO I = MM1, L, -1
            F = S*E( I )
            B = C*E( I )
            CALL DLARTG( G, F, C, S, R )
            IF( I.NE.M-1 )
     $         E( I+1 ) = R
            G = D( I+1 ) - P
            R = ( D( I )-G )*S + TWO*C*B
            P = S*R
            D( I+1 ) = G + P
            G = C*R - B

c           If eigenvectors are desired, then save rotations.

            IF( ICOMPZ.GT.0 ) THEN
               WORK( I ) = C
               WORK( N-1+I ) = -S
            END IF

         END DO ! I

c        If eigenvectors are desired, then apply saved rotations.

         IF( ICOMPZ.GT.0 ) THEN
            MM = M - L + 1
            CALL DLASR( 'R', 'V', 'B', N, MM, WORK( L ), WORK( N-1+L ),
     $                  Z( 1, L ), LDZ )
         END IF

         D( L ) = D( L ) - P
         E( L ) = G
         GO TO 40

c        Eigenvalue found.

   80    CONTINUE
         D( L ) = P

         L = L + 1
         IF( L.LE.LEND )
     $      GO TO 40
         GO TO 140

      ELSE

c        QR Iteration

c        Look for small superdiagonal element.

   90    CONTINUE
         IF( L.NE.LEND ) THEN
            LENDP1 = LEND + 1
            DO M = L, LENDP1, -1
               TST = ABS( E( M-1 ) )**2
               IF( TST.LE.( EPS2*ABS( D( M ) ) )*ABS( D( M-1 ) )+
     $             SAFMIN )GO TO 110
            END DO ! M
         END IF

         M = LEND

  110    CONTINUE
         IF( M.GT.LEND )
     $      E( M-1 ) = ZERO
         P = D( L )
         IF( M.EQ.L )
     $      GO TO 130

c        If remaining matrix is 2-by-2, use DLAE2 or SLAEV2
c        to compute its eigensystem.

         IF( M.EQ.L-1 ) THEN
            IF( ICOMPZ.GT.0 ) THEN
               CALL DLAEV2( D( L-1 ), E( L-1 ), D( L ), RT1, RT2, C, S )
               WORK( M ) = C
               WORK( N-1+M ) = S
               CALL DLASR( 'R', 'V', 'F', N, 2, WORK( M ),
     $                     WORK( N-1+M ), Z( 1, L-1 ), LDZ )
            ELSE
               CALL DLAE2( D( L-1 ), E( L-1 ), D( L ), RT1, RT2 )
            END IF
            D( L-1 ) = RT1
            D( L ) = RT2
            E( L-1 ) = ZERO
            L = L - 2
            IF( L.GE.LEND )
     $         GO TO 90
            GO TO 140
         END IF

         IF( JTOT.EQ.NMAXIT )
     $      GO TO 140
         JTOT = JTOT + 1

c        Form shift.

         G = ( D( L-1 )-P ) / ( TWO*E( L-1 ) )
         R = DLAPY2( G, ONE )
         G = D( M ) - P + ( E( L-1 ) / ( G+SIGN( R, G ) ) )

         S = ONE
         C = ONE
         P = ZERO

c        Inner loop

         LM1 = L - 1
         DO I = M, LM1
            F = S*E( I )
            B = C*E( I )
            CALL DLARTG( G, F, C, S, R )
            IF( I.NE.M )
     $         E( I-1 ) = R
            G = D( I ) - P
            R = ( D( I+1 )-G )*S + TWO*C*B
            P = S*R
            D( I ) = G + P
            G = C*R - B

c           If eigenvectors are desired, then save rotations.

            IF( ICOMPZ.GT.0 ) THEN
               WORK( I ) = C
               WORK( N-1+I ) = S
            END IF

         END DO ! I

c        If eigenvectors are desired, then apply saved rotations.

         IF( ICOMPZ.GT.0 ) THEN
            MM = L - M + 1
            CALL DLASR( 'R', 'V', 'F', N, MM, WORK( M ), WORK( N-1+M ),
     $                  Z( 1, M ), LDZ )
         END IF

         D( L ) = D( L ) - P
         E( LM1 ) = G
         GO TO 90

c        Eigenvalue found.

  130    CONTINUE
         D( L ) = P

         L = L - 1
         IF( L.GE.LEND )
     $      GO TO 90
         GO TO 140

      END IF

c     Undo scaling if necessary

  140 CONTINUE
      IF( ISCALE.EQ.1 ) THEN
         CALL DLASCL( 'G', 0, 0, SSFMAX, ANORM, LENDSV-LSV+1, 1,
     $                D( LSV ), N, INFO )
         CALL DLASCL( 'G', 0, 0, SSFMAX, ANORM, LENDSV-LSV, 1, E( LSV ),
     $                N, INFO )
      ELSE IF( ISCALE.EQ.2 ) THEN
         CALL DLASCL( 'G', 0, 0, SSFMIN, ANORM, LENDSV-LSV+1, 1,
     $                D( LSV ), N, INFO )
         CALL DLASCL( 'G', 0, 0, SSFMIN, ANORM, LENDSV-LSV, 1, E( LSV ),
     $                N, INFO )
      END IF

c     Check for no convergence to an eigenvalue after a total
c     of N*MAXIT iterations.

      IF( JTOT.LT.NMAXIT )
     $   GO TO 10
      DO 150 I = 1, N - 1
         IF( E( I ).NE.ZERO )
     $      INFO = INFO + 1
  150 CONTINUE
      GO TO 190

c     Order eigenvalues and eigenvectors.

  160 CONTINUE
      IF( ICOMPZ.EQ.0 ) THEN

c        Use Quick Sort

         CALL DLASRT( 'I', N, D, INFO )

      ELSE

c        Use Selection Sort to minimize swaps of eigenvectors

         DO II = 2, N
            I = II - 1
            K = I
            P = D( I )
            DO J = II, N
               IF( D( J ).LT.P ) THEN
                  K = J
                  P = D( J )
               END IF
            END DO ! J
            IF( K.NE.I ) THEN
               D( K ) = D( I )
               D( I ) = P
               CALL DSWAP( N, Z( 1, I ), 1, Z( 1, K ), 1 )
            END IF
         END DO ! II
      END IF

  190 CONTINUE

c     End of DSTEQR

      END
