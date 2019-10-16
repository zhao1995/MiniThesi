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

c\Name: dsconv

c\Description:
c  Convergence testing for the symmetric Arnoldi eigenvalue routine.

c\Usage:
c  call dsconv
c     ( N, RITZ, BOUNDS, TOL, NCONV )

c\Arguments
c  N       Integer.  (INPUT)
c          Number of Ritz values to check for convergence.

c  RITZ    Double precision array of length N.  (INPUT)
c          The Ritz values to be checked for convergence.

c  BOUNDS  Double precision array of length N.  (INPUT)
c          Ritz estimates associated with the Ritz values in RITZ.

c  TOL     Double precision scalar.  (INPUT)
c          Desired relative accuracy for a Ritz value to be considered
c          "converged".

c  NCONV   Integer scalar.  (OUTPUT)
c          Number of "converged" Ritz values.

c\EndDoc

c-----------------------------------------------------------------------

c\BeginLib

c\Routines called:
c     second  ARPACK utility routine for timing.
c     dlamch  LAPACK routine that determines machine constants.

c\Author
c     Danny Sorensen               Phuong Vu
c     Richard Lehoucq              CRPC / Rice University
c     Dept. of Computational &     Houston, Texas
c     Applied Mathematics
c     Rice University
c     Houston, Texas

c\SCCS Information: @(#)
c FILE: sconv.F   SID: 2.4   DATE OF SID: 4/19/96   RELEASE: 2

c\Remarks
c     1. Starting with version 2.4, this routine no longer uses the
c        Parlett strategy using the gap conditions.

c\EndLib

c-----------------------------------------------------------------------

      subroutine dsconv (n, ritz, bounds, tol, nconv)
      implicit   none

c     %----------------------------------------------------%
c     | Include files for debugging and timing information |
c     %----------------------------------------------------%

      include   'debug.h'
      include   'stat.h'

c     %------------------%
c     | Scalar Arguments |
c     %------------------%

      integer    n, nconv
      Double precision
     &           tol

c     %-----------------%
c     | Array Arguments |
c     %-----------------%

      Double precision
     &           ritz(n), bounds(n)

c     %---------------%
c     | Local Scalars |
c     %---------------%

      integer    i
      Double precision
     &           temp, eps23
      real       tary(2), etime, t0,t1
      save       t0,t1

c     %-------------------%
c     | External routines |
c     %-------------------%

      Double precision
     &           dlamch
      external   dlamch

c     %---------------------%
c     | Intrinsic Functions |
c     %---------------------%

      intrinsic    abs

c     %-----------------------%
c     | Executable Statements |
c     %-----------------------%

c     call second (t0)
      t0 = etime(tary)

      eps23 = dlamch('Epsilon-Machine')
      eps23 = eps23**(2.0D+0 / 3.0D+0)

      nconv  = 0
      do i = 1, n

c        %-----------------------------------------------------%
c        | The i-th Ritz value is considered "converged"       |
c        | when: bounds(i) .le. TOL*max(eps23, abs(ritz(i)))   |
c        %-----------------------------------------------------%

         temp = max( eps23, abs(ritz(i)) )
         if ( bounds(i) .le. tol*temp ) then
            nconv = nconv + 1
         end if

      end do ! i

c     call second (t1)
      t1 = etime(tary)
      tsconv = tsconv + (t1 - t0)

c     %---------------%
c     | End of dsconv |
c     %---------------%

      end
