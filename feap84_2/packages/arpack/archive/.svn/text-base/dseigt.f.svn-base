c-----------------------------------------------------------------------
c\BeginDoc

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

c\Name: dseigt

c\Description:
c  Compute the eigenvalues of the current symmetric tridiagonal matrix
c  and the corresponding error bounds given the current residual norm.

c\Usage:
c  call dseigt
c     ( RNORM, N, H, LDH, EIG, BOUNDS, WORKL, IERR )

c\Arguments
c  RNORM   Double precision scalar.  (INPUT)
c          RNORM contains the residual norm corresponding to the current
c          symmetric tridiagonal matrix H.

c  N       Integer.  (INPUT)
c          Size of the symmetric tridiagonal matrix H.

c  H       Double precision N by 2 array.  (INPUT)
c          H contains the symmetric tridiagonal matrix with the
c          subdiagonal in the first column starting at H(2,1) and the
c          main diagonal in second column.

c  LDH     Integer.  (INPUT)
c          Leading dimension of H exactly as declared in the calling
c          program.

c  EIG     Double precision array of length N.  (OUTPUT)
c          On output, EIG contains the N eigenvalues of H possibly
c          unsorted.  The BOUNDS arrays are returned in the
c          same sorted order as EIG.

c  BOUNDS  Double precision array of length N.  (OUTPUT)
c          On output, BOUNDS contains the error estimates corresponding
c          to the eigenvalues EIG.  This is equal to RNORM times the
c          last components of the eigenvectors corresponding to the
c          eigenvalues in EIG.

c  WORKL   Double precision work array of length 3*N.  (WORKSPACE)
c          Private (replicated) array on each PE or array allocated on
c          the front end.

c  IERR    Integer.  (OUTPUT)
c          Error exit flag from dstqrb.

c\EndDoc

c-----------------------------------------------------------------------

c\BeginLib

c\Local variables:
c     xxxxxx  real

c\Routines called:
c     dstqrb  ARPACK routine that computes the eigenvalues and the
c             last components of the eigenvectors of a symmetric
c             and tridiagonal matrix.
c     second  ARPACK utility routine for timing.
c     dvout   ARPACK utility routine that prints vectors.
c     dcopy   Level 1 BLAS that copies one vector to another.

c\Author
c     Danny Sorensen               Phuong Vu
c     Richard Lehoucq              CRPC / Rice University
c     Dept. of Computational &     Houston, Texas
c     Applied Mathematics
c     Rice University
c     Houston, Texas

c\Revision history:
c     xx/xx/92: Version ' 2.4'

c\SCCS Information: @(#)
c FILE: seigt.F   SID: 2.4   DATE OF SID: 8/27/96   RELEASE: 2

c\Remarks
c     None

c\EndLib

c-----------------------------------------------------------------------

      subroutine dseigt
     &   ( rnorm, n, h, ldh, eig, bounds, workl, ierr )
      implicit   none

c     %----------------------------------------------------%
c     | Include files for debugging and timing information |
c     %----------------------------------------------------%

      include   'debug.h'
      include   'stat.h'

c     %------------------%
c     | Scalar Arguments |
c     %------------------%

      integer    ierr, ldh, n
      Double precision
     &           rnorm

c     %-----------------%
c     | Array Arguments |
c     %-----------------%

      Double precision
     &           eig(n), bounds(n), h(ldh,2), workl(3*n)

c     %------------%
c     | Parameters |
c     %------------%

      Double precision
     &           zero
      parameter (zero = 0.0D+0)

c     %---------------%
c     | Local Scalars |
c     %---------------%

      integer    k, msglvl
      real       tary(2), etime, t0,t1
      save       t0,t1

c     %----------------------%
c     | External Subroutines |
c     %----------------------%

c     external   dcopy, dstqrb, dvout, second
      external   dcopy, dstqrb, dvout

c     %-----------------------%
c     | Executable Statements |
c     %-----------------------%

c     %-------------------------------%
c     | Initialize timing statistics  |
c     | & message level for debugging |
c     %-------------------------------%

c     call second (t0)
      t0 = etime(tary)
      msglvl = mseigt

      if (msglvl .gt. 0) then
         call dvout (logfil, n, h(1,2), ndigit,
     &              '_seigt: main diagonal of matrix H')
         if (n .gt. 1) then
         call dvout (logfil, n-1, h(2,1), ndigit,
     &              '_seigt: sub diagonal of matrix H')
         end if
      end if

      call dcopy  (n, h(1,2), 1, eig, 1)
      call dcopy  (n-1, h(2,1), 1, workl, 1)
      call dstqrb (n, eig, workl, bounds, workl(n+1), ierr)
      if (ierr .ne. 0) go to 9000
      if (msglvl .gt. 1) then
         call dvout (logfil, n, bounds, ndigit,
     &              '_seigt: last row of the eigenvector matrix for H')
      end if

c     %-----------------------------------------------%
c     | Finally determine the error bounds associated |
c     | with the n Ritz values of H.                  |
c     %-----------------------------------------------%

      do k = 1, n
         bounds(k) = rnorm*abs(bounds(k))
      end do ! k

c     call second (t1)
      t1 = etime(tary)
      tseigt = tseigt + (t1 - t0)

 9000 continue

c     %---------------%
c     | End of dseigt |
c     %---------------%

      end
