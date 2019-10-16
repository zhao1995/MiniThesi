c$Id:$
      subroutine pelcon(numel, nen, neix, ix, id, ic, ielc, icneq, sgn)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Replace old structure by do-while loop           15/07/2010
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose:  Compute number of elements connected to each node
c                matrix.

c      Inputs:
c         numel      -  Number of elements in mesh
c         nen        -  Maximum number of nodes on any element
c         neix       -  Dimension for 'ix' array
c         ix(nen1,*) -  List of nodes connected to each element
c         id(ndf,*)  -  Nodal equation numbers
c         icneq      -  Dimension of IELC (= ic(neq))
c         ic(*)      -  Pointer array
c         sgn        -  ( 1) for finite element;
c                       (-1) for contact elements

c      Outputs:
c         ielc(*)    -  Holds set of elements connected to an equation.
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'debugs.h'
      include  'part0.h'
      include  'sdata.h'

      include  'pointer.h'
      include  'comblk.h'

      logical   feflag
      integer   i, j,k,kk, n, icneq, numel, nen, neix, kp, sgn
      integer   neql
      integer   ix(neix,*), id(ndf,*), ic(*), ielc(*)

      save

c     Initialize array

      feflag = sgn.gt.0

      if(feflag) then
        do i = 1,icneq
          ielc(i) = 0
        end do ! i
      endif

c     Find elements connected to each node

      neql = 0
      do i = 1, numel

c       Skip inactive element

        if(feflag .and. ix(neix-1,i).lt.0) go to 120

c       Test for equations

        do j = 1, nen
          n = ix(j,i)
          if(n.gt.0) then
            do k = 1,ndf
              if(npart.eq.ndfp(k)) then
                kk = id(k,n)
                neql = max(neql,kk)
                if(kk.gt.0) then
                  kp = ic(kk)
                  do while( ielc(kp).ne.0 )
                    kp = kp - 1
                  end do ! while
                  if(feflag) then
                    ielc(kp) =  i
                  else
                    ielc(kp) = -i
                  endif
                endif ! kk > 0
              endif ! npart = ndfp(k)
            end do ! k
          endif ! n > 0
        end do ! j

c       Element equation (Lagrange multiplier treatment)

        if(feflag) then
          if(npart.eq.ndfp(1) .and. np(211).ne.0) then
            call pelconl(i,mr(np(32)),ix(1,i),mr(np(211)),
     &                   ic,ielc,neql)
          endif
        else
          if(neix.gt.nen) then
            call pcelconl(i,ix(1,i), ic,ielc,neql)
          endif
        endif

120     continue

      end do ! i

      if(debug .and. ndebug.gt.1) then
        call ioprof(ic,ielc, neql, 'ELEMENT CONNECTIONS',1)
      endif

      end
