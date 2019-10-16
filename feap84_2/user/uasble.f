c$Id:$
      subroutine uasble(s,p,ld,ns,afl,bfl,b)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c     Purpose: User assembly for sparse matrix/vector forms

c     Use:     Set ubycol, udiag, uall, ulfl (in compac.h) before call.

c     Inputs:
c       s(ns,ns)  - element matrix
c       p(ns)     - element vector
c       ld(ns)    - local/global active equation numbers
c       ns        - size of arrays
c       afl       - Assemble s(ns,ns) into global storage
c       bfl       - Assemble p(ns)    into global storage

c     Outputs:
c       b(*)      - Assembled RHS vector
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'compac.h'
      include  'compas.h'
      include  'part0.h'
      include  'rdat1.h'

      include  'pointer.h'
      include  'comblk.h'

      logical   afl,bfl
      integer   i, ns
      integer   ld(ns)
      real*8    s(ns,ns),p(ns),b(*)

      save

c     Assemble matrix

      if(afl) then

        call cassem(hr(np(npart)),hr(np(npart)),hr(np(npart)),s,
     &              mr(np(226)),mr(np(225)),mr(np(48)),ld,ns,ulfl,
     &              ubycol,udiag,uall)
      endif

c     Assemble vector and compute norms

      if(bfl) then
        do i = 1,ns
          if(ld(i).gt.0) b(ld(i)) = b(ld(i)) + p(i)
        end do ! i
        if(compre) then
          do i = 1,ns
            rnorm1 = rnorm1 + abs(p(i))
          end do ! i
          rnormn = rnormn + dble(ns)
        endif
      endif

      end
