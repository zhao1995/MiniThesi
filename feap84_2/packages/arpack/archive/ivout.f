c-----------------------------------------------------------------------

c  This source file has been adapted from the arpack distribution and is
c  distributed with feap under the following license conditions

c  Rice BSD Software License

c  Permits source and binary redistribution of the software
c  ARPACK and P_ARPACK  for both non-commercial and commercial use.

c  Copyright (c) 2001, Rice University
c  Developed by D.C. Sorensen, R.B. Lehoucq, C. Yang, and K. Maschhoff.
c  All rights reserved.

c  Redistribution and use in source and binary forms, with or without
c  modification, are permitted provided that the following conditions are met:

c  _ Redistributions of source code must retain the above copyright notice,
c    this list of conditions and the following disclaimer.
c  _ Redistributions in binary form must reproduce the above copyright notice,
c    this list of conditions and the following disclaimer in the documentation
c    and/or other materials provided with the distribution.
c  _ If you modify the source for these routines we ask that you change the
c    name of the routine and comment the changes made to the original.
c  _ Written notification is provided to the developers of  intent to use this
c    software.

c  Also, we ask that use of ARPACK is properly cited in any resulting
c  publications or software documentation.

c  _ Neither the name of Rice University (RICE) nor the names of its
c    contributors may be used to endorse or promote products derived from
c    this software without specific prior written permission.

c  THIS SOFTWARE IS PROVIDED BY RICE AND CONTRIBUTORS "AS IS" AND
c  ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
c  TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
c  A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL RICE OR
c  CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
c  SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT
c  NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
c  LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
c  CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,
c  STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
c  ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
c  ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

c  Routine:    IVOUT

c  Purpose:    Integer vector output routine.

c  Usage:      CALL IVOUT (LOUT, N, IX, IDIGIT, IFMT)

c  Arguments
c     N      - Length of array IX. (Input)
c     IX     - Integer array to be printed. (Input)
c     IFMT   - Format to be used in printing array IX. (Input)
c     IDIGIT - Print up to ABS(IDIGIT) decimal digits / number. (Input)
c              If IDIGIT .LT. 0, printing is done with 72 columns.
c              If IDIGIT .GT. 0, printing is done with 132 columns.

c-----------------------------------------------------------------------

      SUBROUTINE IVOUT (LOUT, N, IX, IDIGIT, IFMT)
c     ...
c     ... SPECIFICATIONS FOR ARGUMENTS
      implicit   none
      INTEGER    IX(*), N, IDIGIT, LOUT
      CHARACTER  IFMT*(*)
c     ...
c     ... SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER    I, NDIGIT, K1, K2, LLL
      CHARACTER*80 LINE
*     ...
*     ... SPECIFICATIONS INTRINSICS
      INTRINSIC          MIN
*

      LLL = MIN ( LEN ( IFMT ), 80 )
      DO I = 1, LLL
          LINE(I:I) = '-'
      END DO ! I

      DO I = LLL+1, 80
          LINE(I:I) = ' '
      END DO ! I

      WRITE ( LOUT, 2000 ) IFMT, LINE(1:LLL)
 2000 FORMAT ( /1X, A  /1X, A )

      IF (N .LE. 0) RETURN
      NDIGIT = IDIGIT
      IF (IDIGIT .EQ. 0) NDIGIT = 4

c=======================================================================
c             CODE FOR OUTPUT USING 72 COLUMNS FORMAT
c=======================================================================

      IF (IDIGIT .LT. 0) THEN

      NDIGIT = -IDIGIT
      IF (NDIGIT .LE. 4) THEN
         DO K1 = 1, N, 10
            K2 = MIN0(N,K1+9)
            WRITE(LOUT,1000) K1,K2,(IX(I),I=K1,K2)
         END DO ! K1

      ELSE IF (NDIGIT .LE. 6) THEN
         DO K1 = 1, N, 7
            K2 = MIN0(N,K1+6)
            WRITE(LOUT,1001) K1,K2,(IX(I),I=K1,K2)
         END DO ! K1

      ELSE IF (NDIGIT .LE. 10) THEN
         DO K1 = 1, N, 5
            K2 = MIN0(N,K1+4)
            WRITE(LOUT,1002) K1,K2,(IX(I),I=K1,K2)
         END DO ! K1

      ELSE
         DO K1 = 1, N, 3
            K2 = MIN0(N,K1+2)
            WRITE(LOUT,1003) K1,K2,(IX(I),I=K1,K2)
         END DO ! K1
      END IF

c=======================================================================
c             CODE FOR OUTPUT USING 132 COLUMNS FORMAT
c=======================================================================

      ELSE

      IF (NDIGIT .LE. 4) THEN
         DO K1 = 1, N, 20
            K2 = MIN0(N,K1+19)
            WRITE(LOUT,1000) K1,K2,(IX(I),I=K1,K2)
         END DO ! K1

      ELSE IF (NDIGIT .LE. 6) THEN
         DO K1 = 1, N, 15
            K2 = MIN0(N,K1+14)
            WRITE(LOUT,1001) K1,K2,(IX(I),I=K1,K2)
         END DO ! K1

      ELSE IF (NDIGIT .LE. 10) THEN
         DO K1 = 1, N, 10
            K2 = MIN0(N,K1+9)
            WRITE(LOUT,1002) K1,K2,(IX(I),I=K1,K2)
         END DO ! K1

      ELSE
         DO K1 = 1, N, 7
            K2 = MIN0(N,K1+6)
            WRITE(LOUT,1003) K1,K2,(IX(I),I=K1,K2)
         END DO ! K1
      END IF
      END IF
      WRITE (LOUT,1004)

 1000 FORMAT(1X,I4,' - ',I4,':',20(1X,I5))
 1001 FORMAT(1X,I4,' - ',I4,':',15(1X,I7))
 1002 FORMAT(1X,I4,' - ',I4,':',10(1X,I11))
 1003 FORMAT(1X,I4,' - ',I4,':',7(1X,I15))
 1004 FORMAT(1X,' ')

      END
