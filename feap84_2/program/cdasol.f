c$Id:$
      subroutine cdasol(alr,ali,aur,aui,adr,adi,br,bi,jp,neqs,neqt,enr)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Revise computation of energy measure             18/08/2010
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Solution of algebraic equations stored in profile form

c         Equations '   1  ' to 'neqs' are symmetric.
c         Equations 'neqs+1' to 'neqt' are unsymmetric.

c         Use:
c          a.) All equations are unsymmetric       : neqs = 1 (or 0)
c              N.B.  The top 1 x 1 submatrix is always symmetric.
c                    Both 'al' and 'au' must be provided.

c          b.) All equations are symmetric         : neqs = neqt
c              N.B.  In this case array 'al' is not used.

c          c.) First 'neqs' equations are symmetric: 1 < neqs < neqt
c              N.B.  Storage of 'al' for unsymmetric equations only.

c         Solution of complex equations stored in profile form
c         matrix must be decomposed into its triangular factors
c         using CDATRI before using DASOL.

c      Inputs:
c         alr(*) - Real      part of lower array
c         ali(*) - Imaginary part of lower array
c         aur(*) - Real      part of upper array
c         aui(*) - Imaginary part of upper array
c         adr(*) - Real      part of diagonal array
c         adi(*) - Imaginary part of diagonal array
c         br(*)  - Real      part of right hand side
c         bi(*)  - Imaginary part of right hand side
c         jp(*)  - Pointer to end of column/rows in A array
c         neqs   - Number of symmetric equations
c         neqt   - Number of total     equations

c      Outputs:
c         br(*)  - Real      part of solution x
c         bi(*)  - Imaginary part of solution x
c         enr    - Energy product: x*A*x
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'iofile.h'

      integer   is, j, jr, jh, neqs, neqt, neq, jp(*)
      real*8    enr,eni,  bdr, bdi, dot
      real*8    alr(*),aur(*),adr(*),br(*)
      real*8    ali(*),aui(*),adi(*),bi(*)

      save

c     Find first non-zero entry in right hand side

      do is = 1,neqt
        if(br(is).ne.0.0d0 .or. bi(is).ne.0.0d0) go to 100
      end do ! is
      if(ior.gt.0) write(iow,2000)
      if(ior.lt.0) write(*,2000)
      return

c     Reduce right hand side

c     Symmetric part

100   neq = max(1, neqs)
      do j = is+1,neq
        jr = jp(j-1)
        jh = jp(j) - jr
        if(jh.gt.0) then
          br(j) = br(j) - dot(aur(jr+1),br(j-jh),jh)
     &                  + dot(aui(jr+1),bi(j-jh),jh)
          bi(j) = bi(j) - dot(aur(jr+1),bi(j-jh),jh)
     &                  - dot(aui(jr+1),br(j-jh),jh)
        endif
      end do ! j

c     Unsymmetric part

      do j = max(is,neq)+1,neqt
        jr = jp(j-1)
        jh = jp(j) - jr
        if(jh.gt.0) then
          jr    = jr   - jp(neq)
          br(j) = br(j) - dot(alr(jr+1),br(j-jh),jh)
     &                  + dot(ali(jr+1),bi(j-jh),jh)
          bi(j) = bi(j) - dot(alr(jr+1),bi(j-jh),jh)
     &                  - dot(ali(jr+1),br(j-jh),jh)
        endif
      end do ! j

c     Multiply by inverse of diagonal elements & compute energy measure

      eni = 0.0d0
      do j = is,neqt
        bdr = br(j)
        bdi = bi(j)
        br(j) = bdr*adr(j) - bdi*adi(j)
        bi(j) = bdr*adi(j) + bdi*adr(j)
        enr   = enr + bdr*br(j) - bdi*bi(j)
        eni   = eni + bdi*br(j) + bdr*bi(j)
      end do ! j
      enr = sqrt(enr*enr + eni*eni)       ! Final energy measure

c     Backsubstitution

      do j = neqt,2,-1
        jr = jp(j-1)
        jh = jp(j) - jr
        if(jh.gt.0) then
          call cclred(aur(jr+1),aui(jr+1),br(j),bi(j),
     &                jh, br(j-jh),bi(j-jh))
        endif
      end do ! j

c     Format

2000  format(' *WARNING* Zero right-hand-side vector')

      end
