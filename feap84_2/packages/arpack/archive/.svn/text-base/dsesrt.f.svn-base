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

c\Name: dsesrt

c\Description:
c  Sort the array X in the order specified by WHICH and optionally
c  apply the permutation to the columns of the matrix A.

c\Usage:
c  call dsesrt
c     ( WHICH, APPLY, N, X, NA, A, LDA)

c\Arguments
c  WHICH   Character*2.  (Input)
c          'LM' -> X is sorted into increasing order of magnitude.
c          'SM' -> X is sorted into decreasing order of magnitude.
c          'LA' -> X is sorted into increasing order of algebraic.
c          'SA' -> X is sorted into decreasing order of algebraic.

c  APPLY   Logical.  (Input)
c          APPLY = .TRUE.  -> apply the sorted order to A.
c          APPLY = .FALSE. -> do not apply the sorted order to A.

c  N       Integer.  (INPUT)
c          Dimension of the array X.

c  X      Double precision array of length N.  (INPUT/OUTPUT)
c          The array to be sorted.

c  NA      Integer.  (INPUT)
c          Number of rows of the matrix A.

c  A      Double precision array of length NA by N.  (INPUT/OUTPUT)

c  LDA     Integer.  (INPUT)
c          Leading dimension of A.

c\EndDoc

c-----------------------------------------------------------------------

c\BeginLib

c\Routines
c     dswap  Level 1 BLAS that swaps the contents of two vectors.

c\Authors
c     Danny Sorensen               Phuong Vu
c     Richard Lehoucq              CRPC / Rice University
c     Dept. of Computational &     Houston, Texas
c     Applied Mathematics
c     Rice University
c     Houston, Texas

c\Revision history:
c     12/15/93: Version ' 2.1'.
c               Adapted from the sort routine in LANSO and
c               the ARPACK code dsortr

c\SCCS Information: @(#)
c FILE: sesrt.F   SID: 2.3   DATE OF SID: 4/19/96   RELEASE: 2

c\EndLib

c-----------------------------------------------------------------------

      subroutine dsesrt (which, apply, n, x, na, a, lda)
      implicit   none

c     %------------------%
c     | Scalar Arguments |
c     %------------------%

      character*2 which
      logical    apply
      integer    lda, n, na

c     %-----------------%
c     | Array Arguments |
c     %-----------------%

      Double precision
     &           x(0:n-1), a(lda, 0:n-1)

c     %---------------%
c     | Local Scalars |
c     %---------------%

      integer    i, igap, j
      Double precision
     &           temp

c     %----------------------%
c     | External Subroutines |
c     %----------------------%

      external   dswap

c     %-----------------------%
c     | Executable Statements |
c     %-----------------------%

      igap = n / 2

      if (which .eq. 'SA') then

c        X is sorted into decreasing order of algebraic.

   10    continue
         if (igap .eq. 0) go to 9000
         do i = igap, n-1
            j = i-igap
   20       continue

            if (j.lt.0) go to 30

            if (x(j).lt.x(j+igap)) then
               temp = x(j)
               x(j) = x(j+igap)
               x(j+igap) = temp
               if (apply) call dswap( na, a(1, j), 1, a(1,j+igap), 1)
            else
               go to 30
            endif
            j = j-igap
            go to 20
   30       continue
         end do ! i
         igap = igap / 2
         go to 10

      else if (which .eq. 'SM') then

c        X is sorted into decreasing order of magnitude.

   40    continue
         if (igap .eq. 0) go to 9000
         do i = igap, n-1
            j = i-igap
   50       continue

            if (j.lt.0) go to 60

            if (abs(x(j)).lt.abs(x(j+igap))) then
               temp = x(j)
               x(j) = x(j+igap)
               x(j+igap) = temp
               if (apply) call dswap( na, a(1, j), 1, a(1,j+igap), 1)
            else
               go to 60
            endif
            j = j-igap
            go to 50
   60       continue
         end do ! i
         igap = igap / 2
         go to 40

      else if (which .eq. 'LA') then

c        X is sorted into increasing order of algebraic.

   70    continue
         if (igap .eq. 0) go to 9000
         do i = igap, n-1
            j = i-igap
   80       continue

            if (j.lt.0) go to 90

            if (x(j).gt.x(j+igap)) then
               temp = x(j)
               x(j) = x(j+igap)
               x(j+igap) = temp
               if (apply) call dswap( na, a(1, j), 1, a(1,j+igap), 1)
            else
               go to 90
            endif
            j = j-igap
            go to 80
   90       continue
         end do ! i
         igap = igap / 2
         go to 70

      else if (which .eq. 'LM') then

c        X is sorted into increasing order of magnitude.

  100    continue
         if (igap .eq. 0) go to 9000
         do i = igap, n-1
            j = i-igap
  110       continue

            if (j.lt.0) go to 120

            if (abs(x(j)).gt.abs(x(j+igap))) then
               temp = x(j)
               x(j) = x(j+igap)
               x(j+igap) = temp
               if (apply) call dswap( na, a(1, j), 1, a(1,j+igap), 1)
            else
               go to 120
            endif
            j = j-igap
            go to 110
  120       continue
         end do ! i
         igap = igap / 2
         go to 100
      end if

 9000 continue

c     %---------------%
c     | End of dsesrt |
c     %---------------%

      end
