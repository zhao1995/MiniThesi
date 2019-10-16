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

c\Name: dsgets

c\Description:
c  Given the eigenvalues of the symmetric tridiagonal matrix H,
c  computes the NP shifts AMU that are zeros of the polynomial of
c  degree NP which filters out components of the unwanted eigenvectors
c  corresponding to the AMU's based on some given criteria.

c  NOTE: This is called even in the case of user specified shifts in
c  order to sort the eigenvalues, and error bounds of H for later use.

c\Usage:
c  call dsgets
c     ( ISHIFT, WHICH, KEV, NP, RITZ, BOUNDS, SHIFTS )

c\Arguments
c  ISHIFT  Integer.  (INPUT)
c          Method for selecting the implicit shifts at each iteration.
c          ISHIFT = 0: user specified shifts
c          ISHIFT = 1: exact shift with respect to the matrix H.

c  WHICH   Character*2.  (INPUT)
c          Shift selection criteria.
c          'LM' -> KEV eigenvalues of largest magnitude are retained.
c          'SM' -> KEV eigenvalues of smallest magnitude are retained.
c          'LA' -> KEV eigenvalues of largest value are retained.
c          'SA' -> KEV eigenvalues of smallest value are retained.
c          'BE' -> KEV eigenvalues, half from each end of the spectrum.
c                  If KEV is odd, compute one more from the high end.

c  KEV      Integer.  (INPUT)
c          KEV+NP is the size of the matrix H.

c  NP      Integer.  (INPUT)
c          Number of implicit shifts to be computed.

c  RITZ    Double precision array of length KEV+NP.  (INPUT/OUTPUT)
c          On INPUT, RITZ contains the eigenvalues of H.
c          On OUTPUT, RITZ are sorted so that the unwanted eigenvalues
c          are in the first NP locations and the wanted part is in
c          the last KEV locations.  When exact shifts are selected, the
c          unwanted part corresponds to the shifts to be applied.

c  BOUNDS  Double precision array of length KEV+NP.  (INPUT/OUTPUT)
c          Error bounds corresponding to the ordering in RITZ.

c  SHIFTS  Double precision array of length NP.  (INPUT/OUTPUT)
c          On INPUT:  contains the user specified shifts if ISHIFT = 0.
c          On OUTPUT: contains the shifts sorted into decreasing order
c          of magnitude with respect to the Ritz estimates contained in
c          BOUNDS. If ISHIFT = 0, SHIFTS is not modified on exit.

c\EndDoc

c-----------------------------------------------------------------------

c\BeginLib

c\Local variables:
c     xxxxxx  real

c\Routines called:
c     dsortr  ARPACK utility sorting routine.
c     ivout   ARPACK utility routine that prints integers.
c     second  ARPACK utility routine for timing.
c     dvout   ARPACK utility routine that prints vectors.
c     dcopy   Level 1 BLAS that copies one vector to another.
c     dswap   Level 1 BLAS that swaps the contents of two vectors.

c\Author
c     Danny Sorensen               Phuong Vu
c     Richard Lehoucq              CRPC / Rice University
c     Dept. of Computational &     Houston, Texas
c     Applied Mathematics
c     Rice University
c     Houston, Texas

c\Revision history:
c     xx/xx/93: Version ' 2.1'

c\SCCS Information: @(#)
c FILE: sgets.F   SID: 2.4   DATE OF SID: 4/19/96   RELEASE: 2

c\Remarks

c\EndLib

c-----------------------------------------------------------------------

      subroutine dsgets ( ishift, which, kev, np, ritz, bounds, shifts )
      implicit   none

c     %----------------------------------------------------%
c     | Include files for debugging and timing information |
c     %----------------------------------------------------%

      include   'debug.h'
      include   'stat.h'

c     %------------------%
c     | Scalar Arguments |
c     %------------------%

      character*2 which
      integer    ishift, kev, np

c     %-----------------%
c     | Array Arguments |
c     %-----------------%

      Double precision
     &           bounds(kev+np), ritz(kev+np), shifts(np)

c     %------------%
c     | Parameters |
c     %------------%

      Double precision
     &           one         , zero
      parameter (one = 1.0D+0, zero = 0.0D+0)

c     %---------------%
c     | Local Scalars |
c     %---------------%

      integer    kevd2, msglvl
      real       tary(2), etime, t0,t1

      integer    idum(1)

c     %----------------------%
c     | External Subroutines |
c     %----------------------%

c     external   dswap, dcopy, dsortr, second
      external   dswap, dcopy, dsortr

c     %---------------------%
c     | Intrinsic Functions |
c     %---------------------%

      intrinsic    max, min

c     %-----------------------%
c     | Executable Statements |
c     %-----------------------%

c     %-------------------------------%
c     | Initialize timing statistics  |
c     | & message level for debugging |
c     %-------------------------------%

c     call second (t0)
      t0 = etime(tary)
      msglvl = msgets

      if (which .eq. 'BE') then

c        %-----------------------------------------------------%
c        | Both ends of the spectrum are requested.            |
c        | Sort the eigenvalues into algebraically increasing  |
c        | order first then swap high end of the spectrum next |
c        | to low end in appropriate locations.                |
c        | NOTE: when np < floor(kev/2) be careful not to swap |
c        | overlapping locations.                              |
c        %-----------------------------------------------------%

         call dsortr ('LA', .true., kev+np, ritz, bounds)
         kevd2 = kev / 2
         if ( kev .gt. 1 ) then
            call dswap ( min(kevd2,np), ritz, 1,
     &                   ritz( max(kevd2,np)+1 ), 1)
            call dswap ( min(kevd2,np), bounds, 1,
     &                   bounds( max(kevd2,np)+1 ), 1)
         end if

      else

c        %----------------------------------------------------%
c        | LM, SM, LA, SA case.                               |
c        | Sort the eigenvalues of H into the desired order   |
c        | and apply the resulting order to BOUNDS.           |
c        | The eigenvalues are sorted so that the wanted part |
c        | are always in the last KEV locations.               |
c        %----------------------------------------------------%

         call dsortr (which, .true., kev+np, ritz, bounds)
      end if

      if (ishift .eq. 1 .and. np .gt. 0) then

c        %-------------------------------------------------------%
c        | Sort the unwanted Ritz values used as shifts so that  |
c        | the ones with largest Ritz estimates are first.       |
c        | This will tend to minimize the effects of the         |
c        | forward instability of the iteration when the shifts  |
c        | are applied in subroutine dsapps.                     |
c        %-------------------------------------------------------%

         call dsortr ('SM', .true., np, bounds, ritz)
         call dcopy (np, ritz, 1, shifts, 1)
      end if

c     call second (t1)
      t1 = etime(tary)
      tsgets = tsgets + (t1 - t0)

      if (msglvl .gt. 0) then
         idum(1) = kev
         call ivout (logfil, 1, idum, ndigit, '_sgets: KEV is')
         idum(1) = np
         call ivout (logfil, 1, idum, ndigit, '_sgets: NP is')
         call dvout (logfil, kev+np, ritz, ndigit,
     &        '_sgets: Eigenvalues of current H matrix')
         call dvout (logfil, kev+np, bounds, ndigit,
     &        '_sgets: Associated Ritz estimates')
      end if

c     %---------------%
c     | End of dsgets |
c     %---------------%

      end
