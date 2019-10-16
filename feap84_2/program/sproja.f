c$Id:$
      subroutine sproja(v,t,g,neq,nvc)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Compute subspace projection of 'a' to form 'g'

c      Inputs:
c         v(neq,*) - Set of iteration vectors
c         neq      - Number of equations in A
c         nvc      - Size of projected matrix

c      Scratch:
c         t(neq)   - Working vector

c      Outputs:
c         g(*)     - Projected matrix V_trans * A * V
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'compas.h'
      include  'part0.h'
      include  'ndata.h'
      include  'pointer.h'
      include  'comblk.h'

      include  'p_int.h'

      integer   neq,nvc, i,j,k
      real*8    dot, v(neq,*),t(neq),g(*)

      save

c     Forward reduce eigenvector estimates

      k = 0
      do j = 1,nvc

c       Copy vector 'v' into 'z' and solve equations

        do i = 1,neq
          t(i) = v(i,j)
        end do ! i
        fp(1)  = na
        fp(2)  = nau
        fp(3)  = nal
        fp(4)  = np(20+npart)
        call psolve(ittyp,v(1,j),fp,.false.,.true.,.true.,.false.)

c       Compute projection of stiffness

        do i = 1,j
          k = k + 1
          g(k) = dot(v(1,i),t(1),neq)
        end do ! i
      end do ! j

      end
