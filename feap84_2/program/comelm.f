c$Id:$
      subroutine comelm(id,ix, ir, ndf,nen, kpo,kp,neq,
     &                  bycol,wdiag,all)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose:  Compute equation numbers from elements

c      Inputs:
c         id(ndf,*)  -  Active unknowns at each node.
c         ix(nen1,*) - List of nodes connected to each element
c         ndf        -  Number of unknowns at each node.
c         nen        -  Maximum number of nodes on any element
c         kpo        -  Initial row entry
c         bycol      -  Storage by columns if true
c         wdiag      -  Include diagonal if true
c         all        -  All terms in row/col if true

c      Outputs:
c         ir(*)      -  Row number of each nonzero in stiffness matrix.
c         kp         -  Last entry in ir array
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'iofile.h'
      include  'part0.h'

      logical   addeq, bycol, wdiag, all
      integer   ndf,nen,kpo,kp,neq, i,l,m, kk,neqj
      integer   id(ndf,*),ix(*),ir(*)

      save

      do l = 1,nen
        kk = ix(l)
        if(kk.gt.0) then
          do m = 1, ndf
            if(npart.eq.ndfp(m)) then
              neqj = id(m,kk)

c             Check if equation to be added

              if(all) then                           ! all terms
                addeq   = neqj.gt.0
              elseif(bycol) then                     ! by columns
                if(wdiag) then
                  addeq = neqj.le.neq.and.neqj.gt.0  ! diagonal in
                else
                  addeq = neqj.lt.neq.and.neqj.gt.0  ! diagonal out
                endif
              else                                   ! by rows
                if(wdiag) then
                  addeq = neqj.ge.neq                ! diagonal in
                else
                  addeq = neqj.gt.neq                ! diagonal out
                endif
              endif

c             Add equation to list

              if(addeq) then

c               Check if equation already in list.

                do i = kpo, kp
                  if(ir(i).eq.neqj) go to 200
                end do ! i

c               New equation, add to list

                kp     = kp + 1
                ir(kp) = neqj
200             continue
              endif
            endif
          end do ! m
        endif
      end do ! l

      end
