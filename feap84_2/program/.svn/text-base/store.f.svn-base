c$Id:$
      subroutine store(v,w,nv,itrn)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Store/retrieve BFGS vectors

c      Inputs:
c         v(*),w(*) - BFGS Vectors to save
c         nv        - Number of vector pair to save
c         itrn      - Switch: = 1 to save; = 2 to retrieve

c      Outputs:
c         none      - Saves in hr(*):  (itrn = 1)
c         v(*),w(*) - BFGS Vector   :  (itrn = 2)
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'cdata.h'
      include  'pointer.h'
      include  'comblk.h'

      include  'p_int.h'

      integer   nv, itrn, ns
      real*8    v(*),w(*)

      save

      ns = (neq+neq)*(nv-1)

      if(itrn.eq.1) then
        fp(1) = np(71) + ns
        call pmove(v,hr(fp(1)),neq)
        fp(1) = fp(1) + neq
        call pmove(w,hr(fp(1)),neq)
      elseif(itrn.eq.2) then
        fp(1) = np(71) + ns
        call pmove(hr(fp(1)),v,neq)
        fp(1) = fp(1) + neq
        call pmove(hr(fp(1)),w,neq)
      endif

      end
