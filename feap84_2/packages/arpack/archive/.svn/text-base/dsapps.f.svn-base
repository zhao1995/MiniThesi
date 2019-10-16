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

c\Name: dsapps

c\Description:
c  Given the Arnoldi factorization

c     A*V_{k} - V_{k}*H_{k} = r_{k+p}*e_{k+p}^T,

c  apply NP shifts implicitly resulting in

c     A*(V_{k}*Q) - (V_{k}*Q)*(Q^T* H_{k}*Q) = r_{k+p}*e_{k+p}^T * Q

c  where Q is an orthogonal matrix of order KEV+NP. Q is the product of
c  rotations resulting from the NP bulge chasing sweeps.  The updated Arnoldi
c  factorization becomes:

c     A*VNEW_{k} - VNEW_{k}*HNEW_{k} = rnew_{k}*e_{k}^T.

c\Usage:
c  call dsapps
c     ( N, KEV, NP, SHIFT, V, LDV, H, LDH, RESID, Q, LDQ, WORKD )

c\Arguments
c  N       Integer.  (INPUT)
c          Problem size, i.e. dimension of matrix A.

c  KEV     Integer.  (INPUT)
c          INPUT: KEV+NP is the size of the input matrix H.
c          OUTPUT: KEV is the size of the updated matrix HNEW.

c  NP      Integer.  (INPUT)
c          Number of implicit shifts to be applied.

c  SHIFT   Double precision array of length NP.  (INPUT)
c          The shifts to be applied.

c  V       Double precision N by (KEV+NP) array.  (INPUT/OUTPUT)
c          INPUT: V contains the current KEV+NP Arnoldi vectors.
c          OUTPUT: VNEW = V(1:n,1:KEV); the updated Arnoldi vectors
c          are in the first KEV columns of V.

c  LDV     Integer.  (INPUT)
c          Leading dimension of V exactly as declared in the calling
c          program.

c  H       Double precision (KEV+NP) by 2 array.  (INPUT/OUTPUT)
c          INPUT: H contains the symmetric tridiagonal matrix of the
c          Arnoldi factorization with the subdiagonal in the 1st column
c          starting at H(2,1) and the main diagonal in the 2nd column.
c          OUTPUT: H contains the updated tridiagonal matrix in the
c          KEV leading submatrix.

c  LDH     Integer.  (INPUT)
c          Leading dimension of H exactly as declared in the calling
c          program.

c  RESID   Double precision array of length (N).  (INPUT/OUTPUT)
c          INPUT: RESID contains the the residual vector r_{k+p}.
c          OUTPUT: RESID is the updated residual vector rnew_{k}.

c  Q       Double precision KEV+NP by KEV+NP work array.  (WORKSPACE)
c          Work array used to accumulate the rotations during the bulge
c          chase sweep.

c  LDQ     Integer.  (INPUT)
c          Leading dimension of Q exactly as declared in the calling
c          program.

c  WORKD   Double precision work array of length 2*N.  (WORKSPACE)
c          Distributed array used in the application of the accumulated
c          orthogonal matrix Q.

c\EndDoc

c-----------------------------------------------------------------------

c\BeginLib

c\Local variables:
c     xxxxxx  real

c\References:
c  1. D.C. Sorensen, "Implicit Application of Polynomial Filters in
c     a k-Step Arnoldi Method", SIAM J. Matr. Anal. Apps., 13 (1992),
c     pp 357-385.
c  2. R.B. Lehoucq, "Analysis and Implementation of an Implicitly
c     Restarted Arnoldi Iteration", Rice University Technical Report
c     TR95-13, Department of Computational and Applied Mathematics.

c\Routines called:
c     ivout   ARPACK utility routine that prints integers.
c     second  ARPACK utility routine for timing.
c     dvout   ARPACK utility routine that prints vectors.
c     dlamch  LAPACK routine that determines machine constants.
c     dlartg  LAPACK Givens rotation construction routine.
c     dlacpy  LAPACK matrix copy routine.
c     dlaset  LAPACK matrix initialization routine.
c     dgemv   Level 2 BLAS routine for matrix vector multiplication.
c     daxpy   Level 1 BLAS that computes a vector triad.
c     dcopy   Level 1 BLAS that copies one vector to another.
c     dscal   Level 1 BLAS that scales a vector.

c\Author
c     Danny Sorensen               Phuong Vu
c     Richard Lehoucq              CRPC / Rice University
c     Dept. of Computational &     Houston, Texas
c     Applied Mathematics
c     Rice University
c     Houston, Texas

c\Revision history:
c     12/16/93: Version ' 2.4'

c\SCCS Information: @(#)
c FILE: sapps.F   SID: 2.6   DATE OF SID: 3/28/97   RELEASE: 2

c\Remarks
c  1. In this version, each shift is applied to all the subblocks of
c     the tridiagonal matrix H and not just to the submatrix that it
c     comes from. This routine assumes that the subdiagonal elements
c     of H that are stored in h(1:kev+np,1) are nonegative upon input
c     and enforce this condition upon output. This version incorporates
c     deflation. See code for documentation.

c\EndLib

c-----------------------------------------------------------------------

      subroutine dsapps
     &   ( n, kev, np, shift, v, ldv, h, ldh, resid, q, ldq, workd )
      implicit   none

c     %----------------------------------------------------%
c     | Include files for debugging and timing information |
c     %----------------------------------------------------%

      include   'debug.h'
      include   'stat.h'

c     %------------------%
c     | Scalar Arguments |
c     %------------------%

      integer    kev, ldh, ldq, ldv, n, np

c     %-----------------%
c     | Array Arguments |
c     %-----------------%

      Double precision
     &           h(ldh,2), q(ldq,kev+np), resid(n), shift(np),
     &           v(ldv,kev+np), workd(2*n)

c     %------------%
c     | Parameters |
c     %------------%

      Double precision
     &           one, zero
      parameter (one = 1.0D+0, zero = 0.0D+0)

c     %---------------%
c     | Local Scalars |
c     %---------------%

      integer    i, iend, istart, itop, j, jj, kplusp, msglvl
      logical    first
      Double precision
     &           a1, a2, a3, a4, big, c, epsmch, f, g, r, s
      real       tary(2), etime, t0,t1
      save       epsmch, first,  t0,t1

      integer    idum(1)


c     %----------------------%
c     | External Subroutines |
c     %----------------------%

      external   daxpy, dcopy, dscal, dlacpy, dlartg, dlaset, dvout,
     &           ivout, dgemv
c    &           ivout, second, dgemv

c     %--------------------%
c     | External Functions |
c     %--------------------%

      Double precision
     &           dlamch
      external   dlamch

c     %----------------------%
c     | Intrinsics Functions |
c     %----------------------%

      intrinsic  abs

c     %----------------%
c     | Data statments |
c     %----------------%

      data       first / .true. /

c     %-----------------------%
c     | Executable Statements |
c     %-----------------------%

      if (first) then
         epsmch = dlamch('Epsilon-Machine')
         first = .false.
      end if
      itop = 1

c     %-------------------------------%
c     | Initialize timing statistics  |
c     | & message level for debugging |
c     %-------------------------------%

c     call second (t0)
      t0 = etime(tary)
      msglvl = msapps

      kplusp = kev + np

c     %----------------------------------------------%
c     | Initialize Q to the identity matrix of order |
c     | kplusp used to accumulate the rotations.     |
c     %----------------------------------------------%

      call dlaset ('All', kplusp, kplusp, zero, one, q, ldq)

c     %----------------------------------------------%
c     | Quick return if there are no shifts to apply |
c     %----------------------------------------------%

      if (np .eq. 0) go to 9000

c     %----------------------------------------------------------%
c     | Apply the np shifts implicitly. Apply each shift to the  |
c     | whole matrix and not just to the submatrix from which it |
c     | comes.                                                   |
c     %----------------------------------------------------------%

      do jj = 1, np

         istart = itop

c        %----------------------------------------------------------%
c        | Check for splitting and deflation. Currently we consider |
c        | an off-diagonal element h(i+1,1) negligible if           |
c        |         h(i+1,1) .le. epsmch*( |h(i,2)| + |h(i+1,2)| )   |
c        | for i=1:KEV+NP-1.                                        |
c        | If above condition tests true then we set h(i+1,1) = 0.  |
c        | Note that h(1:KEV+NP,1) are assumed to be non negative.  |
c        %----------------------------------------------------------%

   20    continue

c        %------------------------------------------------%
c        | The following loop exits early if we encounter |
c        | a negligible off diagonal element.             |
c        %------------------------------------------------%

         do i = istart, kplusp-1
            big   = abs(h(i,2)) + abs(h(i+1,2))
            if (h(i+1,1) .le. epsmch*big) then
               if (msglvl .gt. 0) then
                  idum(1) = i
                  call ivout (logfil, 1, idum, ndigit,
     &                 '_sapps: deflation at row/column no.')
                  idum(1) = jj
                  call ivout (logfil, 1, idum, ndigit,
     &                 '_sapps: occured before shift number.')
                  call dvout (logfil, 1, h(i+1,1), ndigit,
     &                 '_sapps: the corresponding off diagonal element')
               end if
               h(i+1,1) = zero
               iend = i
               go to 40
            end if
         end do ! i
         iend = kplusp
   40    continue

         if (istart .lt. iend) then

c           %--------------------------------------------------------%
c           | Construct the plane rotation G'(istart,istart+1,theta) |
c           | that attempts to drive h(istart+1,1) to zero.          |
c           %--------------------------------------------------------%

             f = h(istart,2) - shift(jj)
             g = h(istart+1,1)
             call dlartg (f, g, c, s, r)

c            %-------------------------------------------------------%
c            | Apply rotation to the left and right of H;            |
c            | H <- G' * H * G,  where G = G(istart,istart+1,theta). |
c            | This will create a "bulge".                           |
c            %-------------------------------------------------------%

             a1 = c*h(istart,2)   + s*h(istart+1,1)
             a2 = c*h(istart+1,1) + s*h(istart+1,2)
             a4 = c*h(istart+1,2) - s*h(istart+1,1)
             a3 = c*h(istart+1,1) - s*h(istart,2)
             h(istart,2)   = c*a1 + s*a2
             h(istart+1,2) = c*a4 - s*a3
             h(istart+1,1) = c*a3 + s*a4

c            %----------------------------------------------------%
c            | Accumulate the rotation in the matrix Q;  Q <- Q*G |
c            %----------------------------------------------------%

             do j = 1, min(istart+jj,kplusp)
                a1            =   c*q(j,istart) + s*q(j,istart+1)
                q(j,istart+1) = - s*q(j,istart) + c*q(j,istart+1)
                q(j,istart)   = a1
             end do ! j


c            %----------------------------------------------%
c            | The following loop chases the bulge created. |
c            | Note that the previous rotation may also be  |
c            | done within the following loop. But it is    |
c            | kept separate to make the distinction among  |
c            | the bulge chasing sweeps and the first plane |
c            | rotation designed to drive h(istart+1,1) to  |
c            | zero.                                        |
c            %----------------------------------------------%

             do i = istart+1, iend-1

c               %----------------------------------------------%
c               | Construct the plane rotation G'(i,i+1,theta) |
c               | that zeros the i-th bulge that was created   |
c               | by G(i-1,i,theta). g represents the bulge.   |
c               %----------------------------------------------%

                f = h(i,1)
                g = s*h(i+1,1)

c               %----------------------------------%
c               | Final update with G(i-1,i,theta) |
c               %----------------------------------%

                h(i+1,1) = c*h(i+1,1)
                call dlartg (f, g, c, s, r)

c               %-------------------------------------------%
c               | The following ensures that h(1:iend-1,1), |
c               | the first iend-2 off diagonal of elements |
c               | H, remain non negative.                   |
c               %-------------------------------------------%

                if (r .lt. zero) then
                   r = -r
                   c = -c
                   s = -s
                end if

c               %--------------------------------------------%
c               | Apply rotation to the left and right of H; |
c               | H <- G * H * G',  where G = G(i,i+1,theta) |
c               %--------------------------------------------%

                h(i,1) = r

                a1 = c*h(i,2)   + s*h(i+1,1)
                a2 = c*h(i+1,1) + s*h(i+1,2)
                a3 = c*h(i+1,1) - s*h(i,2)
                a4 = c*h(i+1,2) - s*h(i+1,1)

                h(i,2)   = c*a1 + s*a2
                h(i+1,2) = c*a4 - s*a3
                h(i+1,1) = c*a3 + s*a4

c               %----------------------------------------------------%
c               | Accumulate the rotation in the matrix Q;  Q <- Q*G |
c               %----------------------------------------------------%

                do j = 1, min( i+jj, kplusp )
                   a1       =   c*q(j,i) + s*q(j,i+1)
                   q(j,i+1) = - s*q(j,i) + c*q(j,i+1)
                   q(j,i)   = a1
                end do ! j

             end do ! i

         end if

c        %--------------------------%
c        | Update the block pointer |
c        %--------------------------%

         istart = iend + 1

c        %------------------------------------------%
c        | Make sure that h(iend,1) is non-negative |
c        | If not then set h(iend,1) <-- -h(iend,1) |
c        | and negate the last column of Q.         |
c        | We have effectively carried out a        |
c        | similarity on transformation H           |
c        %------------------------------------------%

         if (h(iend,1) .lt. zero) then
             h(iend,1) = -h(iend,1)
             call dscal(kplusp, -one, q(1,iend), 1)
         end if

c        %--------------------------------------------------------%
c        | Apply the same shift to the next block if there is any |
c        %--------------------------------------------------------%

         if (iend .lt. kplusp) go to 20

c        %-----------------------------------------------------%
c        | Check if we can increase the the start of the block |
c        %-----------------------------------------------------%

         do i = itop, kplusp-1
            if (h(i+1,1) .gt. zero) go to 90
            itop  = itop + 1
         end do ! i
   90    continue

c        %-----------------------------------%
c        | Finished applying the jj-th shift |
c        %-----------------------------------%

      end do ! jj

c     %------------------------------------------%
c     | All shifts have been applied. Check for  |
c     | more possible deflation that might occur |
c     | after the last shift is applied.         |
c     %------------------------------------------%

      do i = itop, kplusp-1
         big   = abs(h(i,2)) + abs(h(i+1,2))
         if (h(i+1,1) .le. epsmch*big) then
            if (msglvl .gt. 0) then
               idum(1) = i
               call ivout (logfil, 1, idum, ndigit,
     &              '_sapps: deflation at row/column no.')
               call dvout (logfil, 1, h(i+1,1), ndigit,
     &              '_sapps: the corresponding off diagonal element')
            end if
            h(i+1,1) = zero
         end if
      end do ! i

c     %-------------------------------------------------%
c     | Compute the (kev+1)-st column of (V*Q) and      |
c     | temporarily store the result in WORKD(N+1:2*N). |
c     | This is not necessary if h(kev+1,1) = 0.         |
c     %-------------------------------------------------%

      if ( h(kev+1,1) .gt. zero )
     &   call dgemv ('N', n, kplusp, one, v, ldv,
     &                q(1,kev+1), 1, zero, workd(n+1), 1)

c     %-------------------------------------------------------%
c     | Compute column 1 to kev of (V*Q) in backward order    |
c     | taking advantage that Q is an upper triangular matrix |
c     | with lower bandwidth np.                              |
c     | Place results in v(:,kplusp-kev:kplusp) temporarily.  |
c     %-------------------------------------------------------%

      do i = 1, kev
         call dgemv ('N', n, kplusp-i+1, one, v, ldv,
     &               q(1,kev-i+1), 1, zero, workd, 1)
         call dcopy (n, workd, 1, v(1,kplusp-i+1), 1)
      end do ! i

c     %-------------------------------------------------%
c     |  Move v(:,kplusp-kev+1:kplusp) into v(:,1:kev). |
c     %-------------------------------------------------%

      call dlacpy ('A', n, kev, v(1,np+1), ldv, v, ldv)

c     %--------------------------------------------%
c     | Copy the (kev+1)-st column of (V*Q) in the |
c     | appropriate place if h(kev+1,1) .ne. zero. |
c     %--------------------------------------------%

      if ( h(kev+1,1) .gt. zero )
     &     call dcopy (n, workd(n+1), 1, v(1,kev+1), 1)

c     %-------------------------------------%
c     | Update the residual vector:         |
c     |    r <- sigmak*r + betak*v(:,kev+1) |
c     | where                               |
c     |    sigmak = (e_{kev+p}'*Q)*e_{kev}  |
c     |    betak = e_{kev+1}'*H*e_{kev}     |
c     %-------------------------------------%

      call dscal (n, q(kplusp,kev), resid, 1)
      if (h(kev+1,1) .gt. zero)
     &   call daxpy (n, h(kev+1,1), v(1,kev+1), 1, resid, 1)

      if (msglvl .gt. 1) then
         call dvout (logfil, 1, q(kplusp,kev), ndigit,
     &      '_sapps: sigmak of the updated residual vector')
         call dvout (logfil, 1, h(kev+1,1), ndigit,
     &      '_sapps: betak of the updated residual vector')
         call dvout (logfil, kev, h(1,2), ndigit,
     &      '_sapps: updated main diagonal of H for next iteration')
         if (kev .gt. 1) then
         call dvout (logfil, kev-1, h(2,1), ndigit,
     &      '_sapps: updated sub diagonal of H for next iteration')
         end if
      end if

c     call second (t1)
      t1 = etime(tary)
      tsapps = tsapps + (t1 - t0)

 9000 continue

c     %---------------%
c     | End of dsapps |
c     %---------------%

      end
