c$Id:$
      subroutine caprod(ad,ac,p,v,jc,ir,neq)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose:  Peform a matrix vector product (A*p) for a symmetric
c                matrix stored in compact form (rows below diagonals,
c                columns above), diagonal stored separtely in 'ad'.

c      Inputs:
c         ad(*)  - Diagonal entries for A array
c         ac(*)  - Off-diagonal entries for symmetric matrix
c         p(*)   - Specified vector for product
c         jc(*)  - Pointer array to locate entries in rows/columns
c         ir(*)  - Location of non-zero entries in A
c         neq    - Number of equations

c      Outputs:
c         v(*)   - Matrix product of A*p.
c-----[--.----+----.----+----.-----------------------------------------]
      implicit none

      integer  neq,ni,nj, jc(*),ir(*)
      real*8   pni,vni, ad(*),ac(*),p(*),v(*)

      save

c     Diagonal part

      do ni = 1,neq
        v(ni) = ad(ni)*p(ni)
      end do ! ni

c     Perform multiplies for non-zero terms

      do ni = 2,neq

        vni = 0.0d0
        pni = p(ni)

        do nj = jc(ni-1)+1,jc(ni)

          vni       = vni       + ac(nj)*p(ir(nj))   ! Lower part
          v(ir(nj)) = v(ir(nj)) + ac(nj)*pni         ! Upper part

        end do ! nj
        v(ni) = v(ni) + vni

      end do ! ni

      end
