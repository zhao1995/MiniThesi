c$Id:$
      subroutine elcnt(numel, nen, neix, id, ix, ic, sgn)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Compute array for performing assembly of stiffness
c               matrices into a compacted data structure;

c               Reference: SESM-84/06: Nour-Omid & Taylor.

c      Inputs:
c          numel  -  Number of elements in the mesh.
c          nen    -  Max. number of nodes per element.
c          neix   -  Dimension of ix array.
c          id     -  Nodal equation numbers.
c          ix     -  Element conectivity array.
c          sgn    -  ( 1) for solid elements
c                    (-1) for contact elements

c      Outputs:
c          ic     -  Array first holds element degree of each equation,
c                    then becomes pointer for array that contains set of
c                    elements connected to each equation.
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'comblk.h'
      include  'part0.h'
      include  'pointer.h'
      include  'sdata.h'

      integer   i,j,k,kk,n, numel, nen, neix, sgn
      integer   id(ndf,*), ix(neix, *), ic(*)

      save

c     Count number of elements attached to each node.

      do i = 1, numel
        do j =  1, nen
          n = ix(j,i)
          if(n.gt.0) then
            do k = 1,ndf
              if(npart.eq.ndfp(k)) then
                kk = id(k,n)
                if(kk.gt.0) ic(kk) = ic(kk) + 1
              endif
            end do ! k
          end if ! n > 0
        end do ! j
      end do ! i

c     Element equation treatment

      if(sgn.gt.0) then

c       Finite element equations

        if(npart.eq.ndfp(1) .and. np(211).ne.0) then
          call elcntl(mr(np(32)),ix,mr(np(211)), ic)
        endif

      elseif(sgn.lt.0) then

c       Contact equations

        if(neix.gt.nen) then
          call celcntl(ix, ic)
        endif

      endif

      end
