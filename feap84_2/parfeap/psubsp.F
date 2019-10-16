c$Id:$
      subroutine psubsp(v,t,gh,d,dp,dtl,p,proceq,
     &                  nf,nv,pgmx,pneq,imas,shift,tol,prt,its)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Revise send/receive for psubsp                   24/11/2006
c       2. Add "include/finclude/petscksp.h"                24/07/2007
c       3. Change 'include/finclude' to 'finclude'          23/01/2009
c       4. Remove common 'pfeapa' (values in 'setups')      05/02/2009
c       5. Update petsc includes to v3.1                    20/07/2010
c       6. Add idum(1) to fix compile error                 05/01/2013
c       7. Limit the number of requested eigenvalues        07/01/2013
c          Add extra argument to scalev
c          Only use numpeq long vectors
c       8. Reduce nf to nmas; check nv/nmas.eq.1 (its =1)   11/04/2013
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Subspace iteration to extract lowest nf eigenpairs
c               Solves:  (A - shift*B)*V = B*V*d for V and d

c      Inputs:
c         nf        - Number of pairs to converge
c         nv        - Size of subspace problem > or = nf
c         pgmx      - Size of gh arrays
c         pneq      - Size of active vectors in V
c         imas      - Switch: =1 for consistent B; =2 for diagonal B
c         shift     - Value of shift
c         tol       - Tolerance to converge pairs
c         prt       - Flag, output iteration arrays if true
c         its       - Maximum number of subspace iterations

c      Scratch:
c         t(vneq)   - Working vector
c         gh(pgmx,4)- Projection of (1) A - shift*B matrix
c                                   (2) B           matrix
c                                   (3,4) - sdatabuf, rdatabuf
c         dp(*)     - Previous iteration values
c         dtl(*)    - Tolerance of eigenvalue iterations
c         p(nv,*)   - Eigenvectors of G*P = H*P*d
c         proceq(*) - Processor mass/equation list

c      Outputs:
c         v(vneq,*) - Eigenvectors
c         d(*)      - Eigenvalues
c-----[--.----+----.----+----.-----------------------------------------]
      implicit   none

#     include   "finclude/petsc.h"
#     include   "pfeapc.h"

#     include   "pfeapb.h"

#     include   "compas.h"
#     include   "evdata.h"
#     include   "fdata.h"
#     include   "gltran.h"
#     include   "iofile.h"
#     include   "sdata.h"
#     include   "setups.h"

#     include   "pointer.h"
#     include   "comblk.h"

      logical    conv,prt
      integer    i,j,k,n, nf,nv,pneq,imas
      integer    its, itlim, it,itt, num,nmas
      integer    req, nn, pgmx
      integer    proceq(*)
      integer    sdatabuf, rdatabuf,idum(1)
      real*4     etime, tary(2)
      real*8     shift,tol, dm,dr, v(vneq,*),t(vneq)
      real*8     gh(pgmx,*),d(*),dp(*),dtl(*),p(nv,*),dpp(4)

      integer    ueige0_msg
      parameter (ueige0_msg=13)

c     Petsc values

      PetscScalar    sol_array(1)
      PetscOffset    i_sol
      PetscInt       starti,endi
      PetscErrorCode ierr

      integer        msg_stat(MPI_STATUS_SIZE)

      save

c     Compute initial iteration vectors


      call pzero(v,nv*vneq)

c     Get diagonal of mass matrix

      if(fl(1)) then
        call MatGetDiagonal(Mmat, xvec, ierr)
      else
        call VecCopy       (Mdiag, xvec, ierr)
      endif
      call VecGetOwnershipRange(xvec, starti, endi ,ierr)
      call VecGetArray         (xvec, sol_array, i_sol, ierr)

c     Count number of nonzero masses


      nn = 0
      do n = 1,numpeq
        if(sol_array(i_sol+n).ne.0.0d0) nn = nn + 1
      end do ! n

c     Sending data asynchronously

      sdatabuf = nn
      do i = 0,ntasks-1
        if(i.ne.rank) then
          call MPI_Isend( sdatabuf, 1, MPI_INTEGER, i,
     &                    ueige0_msg, MPI_COMM_WORLD, req, ierr)
          if(ierr.ne.0) then
            write(iow,*) ' -->> Send Error[',rank,'->',i,']',ierr
            write(  *,*) ' -->> Send Error[',rank,'->',i,']',ierr
          endif

        endif
      end do ! i

c     Receiving data in blocking mode

      nmas = nn
      proceq(1) = 0
      do i = 0, ntasks-1
        if(i.ne.rank) then
          call MPI_Recv( rdatabuf, 1, MPI_INTEGER, i,
     &                   ueige0_msg, MPI_COMM_WORLD, msg_stat, ierr)
          if(ierr.ne.0) then
            write(iow,*) ' -->> Recv Error[',rank,'->',i,']',ierr
            write(  *,*) ' -->> Recv Error[',rank,'->',i,']',ierr
          endif
          nmas        = nmas + rdatabuf
          proceq(i+2) = rdatabuf
        else
          proceq(i+2) = nn
        endif
      end do ! i
      do i = 1,ntasks
        proceq(i+1) = proceq(i) + proceq(i+1)
      end do ! i

c     Should fix the above lines to avoid the barrier
      call MPI_BARRIER(MPI_COMM_WORLD, ierr)

      nf = min(nmas,nf)
      nv = min(nmas,nv)

      nmas = nmas/nv
      i    = proceq(rank+1)
      j    = proceq(rank+1)/nmas + 1
      do n = 1,numpeq
        dm = sol_array(i_sol + n)
        if(dm.ne.0.0d0) then
          v(n,j) = dm
          i      = i + 1
          if(mod(i,nmas).eq.0) j = min(j+1,nv)
        endif
      end do ! n

      call VecRestoreArray(xvec, sol_array, i_sol, ierr)

      idum(1) = 1
      do i = 1,nv
        dp(i)  = 0.0
        dtl(i) = 1.0d0
        call scalev(v(1,i),idum,1,1,numpeq,.false.)
      end do ! i

c     Compute new vectors and project 'a' onto 'g'

      conv = .false.
      itlim = its
      if(nv.eq.nf .or. nmas.eq.1) itlim = 1
      do it = 1,itlim
        itt = it

c       Project 'b' matrix to form 'h' and compute 'z' vectors

        call psprojb(v,t,gh(1,2),pneq,nv)

c       Project 'a' matrix to form 'g'

        call psproja(v,t,gh,nv,pgmx)

c       Solve reduced eigenproblem

        if(imtyp.eq.1) then

c         Eigenproblem: 'g*p = h*p*d'; e.g., vibration modes

          call geig(gh(1,1),gh(1,2),d,p,t,nv,1,prt)

        else

c         Eigenproblem: 'h*p = g*p*1/d'; e.g., structural buckling

          call geig(gh(1,2),gh(1,1),d,p,t,nv,0,prt)

          do n = 1,nv
            if(d(n).ne.0.0d0) d(n) = 1.d0/d(n)
          end do ! n

c         Resort eigenvalues and vectors to ascending order

          num = nv + 1
          do n = 1,nv/2
            dm       = d(n)
            d(n)     = d(num-n)
            d(num-n) = dm
            do i = 1,nv
              dm         = p(i,n)
              p(i,n)     = p(i,num-n)
              p(i,num-n) = dm
            end do ! i
          end do ! n
        endif

c       Compute new iteration vector 'u' from 'z'

        do i = 1,numpeq

c         Move row of 'v' into temporary vector 't'

          do j = 1,nv
            t(j) = v(i,j)
          end do ! j

c         Compute new iiteration vector entries

          do j = 1,nv
            v(i,j) = 0.0d0
            do k = 1,nv
              v(i,j) = v(i,j) + t(k)*p(k,j)
            end do ! k
          end do ! j
        end do ! i

c       Check for convergence

        do n = 1,nv
          if(d(n).ne.0.0d0) dtl(n) = abs((d(n)-dp(n))/d(n))
          dp(n) = d(n)
        end do ! n

c       Compute maximum error for iteration

        conv = .true.
        dm   =  0.0d0
        do n = 1,nf
          if(dtl(n).gt.tol) conv = .false.
          dm = max(dm,dtl(n))
        end do ! n

        if(prt) then
          if(ior.gt.0) write(iow,2006) it,(d(n),n=1,nv)
          if(ior.lt.0) write(  *,2006) it,(d(n),n=1,nv)
          if(itlim.gt.1) then
            if(ior.gt.0) write(iow,2001) it,(dtl(n),n=1,nv)
            if(ior.lt.0) write(  *,2001) it,(dtl(n),n=1,nv)
          endif
        endif

        if(conv) go to 200
        if(rank.eq.0) then
          dm = 1.d0/gtan(1)
          write(*,2000) it,((1.d0/d(n)+shift)*dm,n=1,nf)
          write(*,2001) it,(dtl(n),n=1,nf)
        endif
        if(ior.lt.0) write(*,2005) it, etime(tary)
      end do ! it

c     Scaled vectors mass orthonormal

200   dm = 1.d0/gtan(1)
      do n = 1,nv
        d(n) = (1.0/d(n) + shift)*dm
        dp(n)= sqrt(abs(d(n)))
      end do ! n

c     Output solution values

      write(iow,2000) itt,(d(n),n=1,nv)
      if(itt.gt.1) write(iow,2001) itt,(dtl(n),n=1,nv)
      if(imtyp.eq.1) then
        write(iow,2003) (dp(n),n=1,nv)
        dm = 0.5d0/acos(-1.0d0)
        write(iow,2004) (dp(n)*dm,n=1,nv)
        write(iow,2007)
        dr = 1.d0/dm
        do i = 0,nv-1,4
          do j = 1,min(4,nv-i)
            if(abs(dp(i+j)).gt.tol*10.d0) then
              dpp(j) = dr/dp(i+j)
            else
              dpp(j) = 0.0d0
            endif
          end do ! j
          write(iow,2008) (dpp(j),j=1,min(4,nv-i))
        end do ! i
      endif
      if(ior.lt.0) then
        write(*,2000) itt,(d(n),n=1,nv)
        if(itt.gt.1) write(*,2001) itt,(dtl(n),n=1,nv)
        if(imtyp.eq.1) then
          write(*,2003) (dp(n),n=1,nv)
          write(*,2004) (dp(n)*dm,n=1,nv)
          write(*,2007)
          do i = 0,nv-1,4
            do j = 1,min(4,nv-i)
              if(abs(dp(i+j)).gt.tol*10.d0) then
                dpp(j) = dr/dp(i+j)
              else
                dpp(j) = 0.0d0
              endif
            end do ! j
            write(*,2008) (dpp(j),j=1,min(4,nv-i))
          end do ! i
        endif
      endif

c     Formats

2000  format(/'  SUBSPACE: Current eigenvalues, iteration',i4/
     +        (5x,1p,4d17.8))

2001  format( '  SUBSPACE: Current residuals,   iteration',i4/
     +        (5x,1p,4d17.8))

2003  format(/'  SUBSPACE: Square root of eigenvalues (rad/t)'/
     +        (5x,1p,4d17.8))

2004  format(/'  SUBSPACE: Square root of eigenvalues (Hz.)'/
     +        (5x,1p,4d17.8))

2005  format('  Completed subspace iteration',i4,' Time =',f9.2:,
     +       '  Max. tol =',1p,e12.5)

2006  format(/'  SUBSPACE: Inverse eigenvalues, iteration',i4/
     +        (5x,1p,4d17.8))

2007  format(/'  SUBSPACE: Period in units of time      (T)')

2008  format(5x,1p,4d17.8)

      end
