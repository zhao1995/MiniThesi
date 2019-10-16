c$Id:$
      subroutine psproja(v,t,gh,nvc,pmax)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Save only numpeq entries in t() for dot product  07/01/2013
c          Eliminate pneq from argument list
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Compute subspace projection of 'a' to form 'g'

c      Inputs:
c         v(neq,*) - Set of iteration vectors
c         nvc      - Size of projected matrix
c         pmax     - Size of G and H arrays

c      Scratch:
c         t(vneq)  - Working vector

c      Outputs:
c         gh(pmax,*) - (1) Projected matrix V_trans * A * V
c                      (2) Projected matrix V_trans * B * V
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include   'pfeapb.h'
      include   'compas.h'
      include   'part0.h'
      include   'ndata.h'
      include   'pointer.h'
      include   'comblk.h'

      include   'p_int.h'

      integer    nvc,pmax
      integer    i, j, k
      real*8     dot,gh(pmax,*),v(vneq,*),t(vneq)

      save

c     Forward reduce eigenvector estimates

      k = 0
      do j = 1,nvc

c       Copy vector 'v' into 't' and solve equations

        do i = 1,numpeq
          t(i) = v(i,j)
        end do ! i
        fp(1) = na
        fp(2) = nau
        fp(3) = nal
        fp(4) = np(20+npart)
        call psolve(ittyp,v(1,j),fp,.false.,.true.,.true.,.false.)

c       Compute projection of stiffness

        do i = 1,j
          k = k + 1
          gh(k,1) = dot(v(1,i),t(1),numpeq)
        end do ! i

      end do ! j


c     Accumulate g and h from all processors

      call pfeapsr(gh,gh(1,3),2*pmax,.true.)

      end
