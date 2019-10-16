c$Id:$
      subroutine dasol(al,au,ad,b,jp,neqs, neqt, energy, scale)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Add 'setups.h' and check rank for output         24/02/2009
c       2. Separate np(67) and np(307) for sparse solve     01/05/2012
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Solution of algebraic equations stored in profile form

c         Equations '   1  ' to 'neqs' are symmetric.
c         Equations 'neqs+1' to 'neqt' are unsymmetric.

c         Use:
c          a.) All equations are unsymmetric       : neqs = 1 (or 0)
c              N.B.  The top 1 x 1 submatrix is always symmetric.
c                    Both 'al' and 'au' must be provided.

c          b.) All equations are symmetric         : neqs = neqt
c              N.B.  In this case the array 'al' is not used.

c          c.) First 'neqs' equations are symmetric: 1 < neqs < neqt
c              N.B.  Storage of 'al' for unsymmetric equations only.

c      Coefficient matrix must be decomposed into its triangular
c      factors using 'DATRI' before using 'DASOL'.

c      Inputs:
c         al(*)  - Lower triangular factors of A
c         au(*)  - Upper triangular factors of A
c         ad(*)  - Diagonal factors of A
c         b(*)   - Right hand side vector
c         jp(*)  - Pointer array for row/columns of 'al', 'au'.
c         neqs   - Number of symmetric equations
c         neqt   - Number of equations to solve.
c         scale  - Scaling flag

c      Outputs:
c         b(*)   - Solution vector, x.
c         energy - Energy of solution: x*A*x
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'compas.h'
      include  'complx.h'
      include  'iodata.h'
      include  'iofile.h'
      include  'part0.h'
      include  'setups.h'
      include  'ssolve.h'

      include  'pointer.h'
      include  'comblk.h'

      logical   scale
      character fau*12,fal*12
      integer   is, j, jh, neqs, neqt, neq, iunau,iunal, jp(*)
      real*8    dot, bd, energy, al(*),au(*),ad(*),b(*)

      save

c     Initialize energy

      energy = 0.0d0

c     Sparse solver

      if(ittyp.eq.-2) then

c       Scale  right hand side

        if(scale) then
          call pscalb(b,hr(np(234+npart)),1,neqt)
        endif
        call sdrv(neqt,mr(np(47)),mr(np(48)),mr(np(93)+neqt),mr(np(94)),
     &            ad,b,b,mr(np(307)),hr(np(67)),esp,3,is)
        if(scale) then
          call pscalb(b,hr(np(234+npart)),1,neqt)
        endif

c     Blocked Profile Solver

      elseif(ittyp.eq.-1) then

        iunal = ios
        iunau = ios

        fau   = 'Aupper'
        fal   = 'Alower'

        if(neqs.eq.neqt) fal   = fau

        call xdasol( fau,fal,iunau,iunal,b,jp,neqt,energy,
     &               ad,au,al,mr(np(92)),maxbl,scale )

c     Solve in-core

      else

c       Solve complex equations

        if(cplxfl) then
          is = jp(neqt) + 1
          jh = neqt     + 1
          call cdasol(al,al(is),au,au(is),ad,ad(jh),b,b(jh),jp,
     &                neqs,neqt,energy)

c       Solve real equations

        else

c         Find first non-zero entry in right hand side

          do is = 1,neqt
            if(b(is).ne.0.0d0) go to 100
          end do ! is
          if(rank.eq.0) then
            if(ior.gt.0) write(iow,2000)
            if(ior.lt.0) write(*,2000)
          endif
          return

c         Scale  right hand side

100       if(scale) then
            call pscalb(b,hr(np(234+npart)),is,neqt)
          endif

c         Reduce right hand side

c         Do symmetric part

          neq = max(1, neqs)
          do j = is+1,neq
            jh = jp(j) - jp(j-1)
            if(jh.gt.0) then
              b(j) = b(j) - dot(au(jp(j-1)+1),b(j-jh),jh)
            endif
          end do ! j

c         Do unsymmetric part

          do j = max(is,neq)+1,neqt
            jh = jp(j) - jp(j-1)
            if(jh.gt.0) then
              b(j) = b(j) - dot(al(jp(j-1)-jp(neq)+1),b(j-jh),jh)
            endif
          end do ! j

c         Multiply by inverse of diagonal elements

          do j = is,neqt
            bd     = b(j)
            b(j)   = b(j)*ad(j)
            energy = energy + bd*b(j)
          end do ! j

c         Symmetric and unsymmetric backsubstitution

          do j = neqt,2,-1
            jh = jp(j) - jp(j-1)
            if(jh.gt.0) then
              call colred(au(jp(j-1)+1),b(j),jh, b(j-jh))
            endif
          end do ! j

c         Scale  right hand side

          if(scale) then
            call pscalb(b,hr(np(234+npart)),1,neqt)
          endif

        endif

      endif

2000  format(' *WARNING* Zero right-hand-side vector')

      end
