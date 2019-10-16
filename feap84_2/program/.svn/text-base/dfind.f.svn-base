c$Id:$
      subroutine dfind(d,ri,r0,nupd,g0,g,s,neq,v,w,nbfgs)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Compute BFGS update vectors and search direction

c      Inputs:
c         ri(*)     - Current residual
c         r0(*)     - Old residual
c         g0        - Line-search initial value
c         g         - Line-search final value
c         s         - Line-search step size
c         neq       - Number of equations

c      Outputs:
c         d(*)      - Solution increment vector
c         v         - BFGS update vectors
c         w         - BFGS update vectors
c         nbfgs     - Number BFGS steps
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'compas.h'
      include  'iofile.h'
      include  'ndata.h'
      include  'part0.h'
      include  'pointer.h'
      include  'print.h'
      include  'pscal.h'

      include  'p_int.h'

      logical   upfl
      integer   i,j,ii,nupd,neq,nbfgs
      real*8    s,g0,g,condmx,delgam,dlkdl,coef,fact,stcond, dot
      real*8    d(*),ri(*),r0(*),v(*),w(*)

      save

      data      condmx /1.0d5/

c     Find a new search direction using bfsg method in factored form

c     delgam = delta-(i) : gamma-(i)
c     dlkdl  = delta-(i) : K-(i-1) : delta-(i)

      delgam = s*(g0-g)
      dlkdl  = s*s*g0
      upfl   = delgam.gt.0.0d0 .and. dlkdl.gt.0.0d0

c     If G(0) > G(s) & G(0) > 0

      if (upfl) then

        stcond=sqrt(abs(delgam/dlkdl))

        if(ior.lt.0) write(*,2000) stcond

c       Compute updating vectors v, w and put residual into d

        fact = s/delgam
        coef = 1.d0 + s*stcond
        do i = 1,neq
          v(i)  = ri(i) - coef*r0(i)
          w(i)  = d(i)*fact
          d(i)  = ri(i)
          r0(i) = ri(i)
        end do ! i

c       Check estimate on condition number

        upfl = upfl.and.stcond.lt.condmx
        if ( upfl ) then

c         Save updating factors

          call store(v,w,nupd+1,1)

c         Compute search direction: d

          coef=fact*g
          do i = 1,neq
            d(i) = d(i) + coef*v(i)
          end do ! i
        endif
      else
        do i = 1,neq
          d(i)  = ri(i)
          r0(i) = ri(i)
        end do ! i
      endif

c     Right half of update

      do i = 1,nupd
        ii = nupd - i + 1
        call store(v,w,ii,2)
        coef = dot(w,d,neq)
        do j = 1,neq
          d(j)=d(j)+coef*v(j)
        end do ! j
      end do ! i

c     Resolution

      fp(1)  = na
      fp(2)  = nau
      fp(3)  = nal
      fp(4)  = np(20+npart)
      call psolve(ittyp,d,fp,.false.,.true.,.true.,prnt)

      if(upfl) nupd = nupd + 1

c     Left half of updating

      do i = 1,nupd
        if(i.gt.1) call store(v,w,i,2)
        coef = dot(v,d,neq)
        do j = 1,neq
          d(j) = d(j) + coef*w(j)
        end do ! j
      end do ! i

      nupd = mod(nupd,nbfgs)

c     Format

2000  format(' ---> Stiffness condition no. ',1p,1e12.5)

      end
