c$Id:$
      subroutine pmodify(ld, p, s, ub, nsiz, nst)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Modify RHS for imposed displacements & compress for
c               assembly of element arrays.
c
c      Inputs:
c        ld(*,2) - Equations for rows and columns
c        p(*)    - Unmodified element residual
c        s(*,*)  - Element tangent array
c        nsiz    - Size of current matrix
c        nst     - Dimension of stiffness array

c      Outputs:
c        s(*,*)  - Compressed tangent
c        p(*)    - Compressed residual with non-zero displacement
c-----[--.----+----.----+----.-----------------------------------------]
      implicit   none

      include   'pfeapb.h'
      include   'pointer.h'

      integer    n,m, nsiz, nst
      integer    ld(nst,*)
      real*8     p(nst), s(nst,*), ub(nst)

c     Modify tangent and residual for B.C. in

      if(pfeap_bcin) then
        do n = 1,nsiz
          if(ld(n,5).ne.0) then  ! inactive equation
            do m = 1,nsiz
              s(n,m) = 0.0d0
              s(m,n) = 0.0d0
            end do ! m
            s(n,n) = 1.0d0
            p(n)   = ub(n)       ! non-zero displacement
          endif
        end do ! n
      endif

c     Compress tangent and residual for parallel assembly

      call pcompress(ld,ld(1,7), p,s, nsiz, nst)

      end
