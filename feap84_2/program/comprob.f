c$Id:$
      subroutine comprob(numnp, nen, nen1, ndf, ix, id,
     &                   ic, ielc, ir, jc, bycol, wdiag, all)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose:  Compute locations for non-zero terms in coefficient
c                matrix.

c      Inputs:
c         numnp      -  Number of nodes in mesh
c         nen        -  Maximum number of nodes on any element
c         nen1       -  Dimension for 'ix' array
c         ndf        -  Number of unknowns at each node.
c         ix(nen1,*) -  List of nodes connected to each element
c         id         -  Active unknowns at each node.
c         ic         -  Pointer for ielc list
c         ielc(*)    -  Holds set of elements connected to each node.
c         bycol      -  Storage by columns if true
c         wdiag      -  Include diagonal if true
c         all        -  All terms in row/col if true

c      Outputs:
c         ir(*)      -  Row number of each nonzero in stiffness matrix.
c         jc(*)      -  End of enteries in ir from a given column.
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'compac.h'
      include  'part0.h'
      include  'pglob1.h'
      include  'pointer.h'
      include  'comblk.h'

      include  'p_int.h'

      logical   bycol, wdiag, all
      integer   i, j, k, ne, nep, neq, nn
      integer   numnp, nen, nen1, ndf, kp, kpo
      integer   ix(nen1,*), id(ndf,*), ic(*), ir(*), ielc(*), jc(*)

      save

c     Set up compressed profile pointers.

      neq = 0
      do i = 1, numnp
        do j = 1,ndf
          if(npart.eq.ndfp(j)) neq = max(neq,id(j,i))
        end do ! j
      end do ! i

c     Check for element equations

      if(npart.eq.ndfp(1) .and. np(211).ne.0) then
        call comprol(mr(np(32)),ix,mr(np(211)), neq)
      endif

c     Add global equations if necessary

      gneq = neq
      if(gceflg) then
        neq = neq + geqnum
      endif

c     Do all equations

      kp  = 0
      nep = 1
      do i = 1, gneq
        ne    = ic(i)
        jc(i) = kp
        kpo   = kp + 1
        do k = nep, ne
          nn = ielc(k)

c         Check element type(>0: FE, <0: contact)

          if(nn.gt.0) then
            call comelm(id,ix(1,nn), ir, ndf,nen,  kpo,kp,i,
     &                  bycol,wdiag,all)
          else
            fp(1) = np(168) - ncen1*(nn + 1)
            call comelm(id,mr(fp(1)),   ir, ndf,ncen, kpo,kp,i,
     &                  bycol,wdiag,all)
          endif

c         End element tests

        end do ! k

c       Global equation additions

        if(gceflg) then
          call comgeq(ir, kp, bycol, all, neq)
        endif

        jc(i) = kp
        nep   = ne + 1
      end do ! i

c     Add global equations storage

      do i = gneq+1,neq
        do k = 1,neq
          ir(kp+k) = k
        end do ! k
        kp    = kp + neq
        jc(i) = kp
      end do ! i

      end
