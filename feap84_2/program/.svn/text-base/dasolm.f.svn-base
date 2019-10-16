c$Id:$
      subroutine dasolm(al,au,ad,b,jp,neqr,neqs,neqt,energy,scale)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
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
c         neqr   - Number of equations to reduce for subspace
c         neqs   - Number of symmetric equations
c         neqt   - Number of equations to solve.
c         scale  - Scale flag

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
      include  'pointer.h'
      include  'comblk.h'

      logical   scale, setvar, palloc
      character fau*12,fal*12
      integer   is, j, jr, jh, neq1, neqj, neqr, neqs, neqt, neql
      integer   iunau,iunal, jp(*)
      real*8    dot, bd, energy, al(*),au(*),ad(*),b(*)

      save

c     Initialize energy

      energy = 0.0d0

c     Blocked Profile Solver

      if(compfl) then

        iunal = ios
        iunau = ios

        fau   = 'Aupper'
        fal   = 'Alower'

        if(neqs.eq.neqt) fal   = fau

        call xdasol( fau,fal,iunau,iunal,b,jp,neqt,energy,
     &               ad,au,al,mr(np(92)),maxbl,scale )

c     Solve in-core

      else

c       Solve for complex equations

        if(cplxfl) then
          is = jp(neqt) + 1
          jh = neqt     + 1
          call cdasol(al,al(is),au,au(is),ad,ad(jh),b,b(jh),jp,
     &                neqs,neqt,energy)
          return
        endif

c       Find first non-zero entry in right hand side

        do is = 1,neqt
          if(b(is).ne.0.0d0) go to 100
        end do ! is
        if(ior.gt.0) write(iow,2000)
        if(ior.lt.0) write(*,2000)
        return

c       Scale equations

100     if(scale) then
          call pscalb(b,hr(np(234+npart)),is,neqt)
        endif

c       Reduce right hand side

c       Do symmetric part

        neql = max(1, neqs)
        neq1 = neqr + 1
        do j = is+1,neql
          jr = jp(j-1)
          jh = jp(j) - jr
          if(jh.gt.0) then
            b(j) = b(j) - dot(au(jr+1),b(j-jh),min(jh,jh-j+neq1))
          endif
        end do ! j

c       Do unsymmetric part

        do j = max(is,neql)+1,neqt
          jr = jp(j-1)
          jh = jp(j) - jr
          if(jh.gt.0) then
            jr   = jr   - jp(neql)
            b(j) = b(j) - dot(al(jr+1),b(j-jh),min(jh,jh-j+neq1))
          endif
        end do ! j

c       Multiply by inverse of diagonal elements

        do j = is,min(neqr,neqt)
          bd = b(j)
          b(j) = b(j)*ad(j)
          energy = energy + bd*b(j)
        end do ! j

c       Export neqj by neqj matrix for SVD if neqr < neqt
c       (i.e., there are joints)
c       neqj = total number of joint equations

        if(neqr.lt.neqt) then
          neqj = neqt - neqr

          setvar = palloc(172,'SVDA', neqj*neqj,2)
          setvar = palloc(173,'SVDV', neqj*neqj,2)
          setvar = palloc(174,'SVDW', neqj*neqj,2)
          setvar = palloc(175,'SVDX', neqj*neqj,2)

          call svdsol(ad,au,al,jp,b,hr(np(172)),hr(np(173)),hr(np(174)),
     &                hr(np(175)),neqj,neqr,neqt,energy)
        endif

c       Symmetric and unsymmetric backsubstitution

        do j = neqt,2,-1
          jr = jp(j-1)
          jh = jp(j) - jr
          if(jh.gt.0) then
            call colred(au(jr+1),b(j),min(jh,jh-j+neq1),b(j-jh))
          endif
        end do ! j

c       Scale equations

        if(scale) then
          call pscalb(b,hr(np(234+npart)),1,neqt)
        endif

      endif

2000  format(' *WARNING* Zero right-hand-side vector')

      end
