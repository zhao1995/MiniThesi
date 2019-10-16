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

c\Name: dsortr

c\Description:
c  Sort the array X1 in the order specified by WHICH and optionally
c  applies the permutation to the array X2.

c\Usage:
c  call dsortr
c     ( WHICH, APPLY, N, X1, X2 )

c\Arguments
c  WHICH   Character*2.  (Input)
c          'LM' -> X1 is sorted into increasing order of magnitude.
c          'SM' -> X1 is sorted into decreasing order of magnitude.
c          'LA' -> X1 is sorted into increasing order of algebraic.
c          'SA' -> X1 is sorted into decreasing order of algebraic.

c  APPLY   Logical.  (Input)
c          APPLY = .TRUE.  -> apply the sorted order to X2.
c          APPLY = .FALSE. -> do not apply the sorted order to X2.

c  N       Integer.  (INPUT)
c          Size of the arrays.

c  X1      Double precision array of length N.  (INPUT/OUTPUT)
c          The array to be sorted.

c  X2      Double precision array of length N.  (INPUT/OUTPUT)
c          Only referenced if APPLY = .TRUE.

c\EndDoc

c-----------------------------------------------------------------------

c\BeginLib

c\Author
c     Danny Sorensen               Phuong Vu
c     Richard Lehoucq              CRPC / Rice University
c     Dept. of Computational &     Houston, Texas
c     Applied Mathematics
c     Rice University
c     Houston, Texas

c\Revision history:
c     12/16/93: Version ' 2.1'.
c               Adapted from the sort routine in LANSO.

c\SCCS Information: @(#)
c FILE: sortr.F   SID: 2.3   DATE OF SID: 4/19/96   RELEASE: 2

c\EndLib

c-----------------------------------------------------------------------

      subroutine dsortr (which, apply, n, x1, x2)
      implicit   none

c     %------------------%
c     | Scalar Arguments |
c     %------------------%

      character*2 which
      logical    apply
      integer    n

c     %-----------------%
c     | Array Arguments |
c     %-----------------%

      Double precision
     &           x1(0:n-1), x2(0:n-1)

c     %---------------%
c     | Local Scalars |
c     %---------------%

      integer    i, igap, j
      Double precision
     &           temp

c     %-----------------------%
c     | Executable Statements |
c     %-----------------------%

      igap = n / 2

      if (which .eq. 'SA') then

c        X1 is sorted into decreasing order of algebraic.

   10    continue
         if (igap .eq. 0) go to 9000
         do i = igap, n-1
            j = i-igap
   20       continue

            if (j.lt.0) go to 30

            if (x1(j).lt.x1(j+igap)) then
               temp = x1(j)
               x1(j) = x1(j+igap)
               x1(j+igap) = temp
               if (apply) then
                  temp = x2(j)
                  x2(j) = x2(j+igap)
                  x2(j+igap) = temp
               end if
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

c        X1 is sorted into decreasing order of magnitude.

   40    continue
         if (igap .eq. 0) go to 9000
         do i = igap, n-1
            j = i-igap
   50       continue

            if (j.lt.0) go to 60

            if (abs(x1(j)).lt.abs(x1(j+igap))) then
               temp = x1(j)
               x1(j) = x1(j+igap)
               x1(j+igap) = temp
               if (apply) then
                  temp = x2(j)
                  x2(j) = x2(j+igap)
                  x2(j+igap) = temp
               end if
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

c        X1 is sorted into increasing order of algebraic.

   70    continue
         if (igap .eq. 0) go to 9000
         do i = igap, n-1
            j = i-igap
   80       continue

            if (j.lt.0) go to 90

            if (x1(j).gt.x1(j+igap)) then
               temp = x1(j)
               x1(j) = x1(j+igap)
               x1(j+igap) = temp
               if (apply) then
                  temp = x2(j)
                  x2(j) = x2(j+igap)
                  x2(j+igap) = temp
               end if
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

c        X1 is sorted into increasing order of magnitude.

  100    continue
         if (igap .eq. 0) go to 9000
         do i = igap, n-1
            j = i-igap
  110       continue

            if (j.lt.0) go to 120

            if (abs(x1(j)).gt.abs(x1(j+igap))) then
               temp = x1(j)
               x1(j) = x1(j+igap)
               x1(j+igap) = temp
               if (apply) then
                  temp = x2(j)
                  x2(j) = x2(j+igap)
                  x2(j+igap) = temp
               end if
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
c     | End of dsortr |
c     %---------------%

      end
