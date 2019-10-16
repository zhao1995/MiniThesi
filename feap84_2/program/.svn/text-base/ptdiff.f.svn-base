c$Id:$
      subroutine ptdiff(k1,du,flg)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Remove allocation of TEMP3 and np(113) from call
c          to pndiff                                        07/12/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Driver for numerical tangent differentiations

c      Inputs:
c        k1     - Element number for tangent
c        du     - Increment to use for tangent (replaces default)
c        flg    - Computes analytical tangent if true

c      Outputs:
c        Tangent output to file
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'eltran.h'
      include  'hdatam.h'
      include  'pointer.h'
      include  'sdata.h'
      include  'comblk.h'

      include  'p_int.h'

      logical   flg,setvar,palloc,hfold,h3old
      integer   k1
      real*8    du

      save

      fp(1) = np(33) + (k1 - 1)*nen1

      setvar = palloc(112,'TEMP2',  2*nst*nst,2)

      hfold  = hflgu
      h3old  = h3flgu
      hflgu  = .false.
      h3flgu = .false.

      call pndiff(mr(fp(1)),mr(np(32)),hr(np(41)),hr(np(112)),k1,du,flg)

      hflgu  = hfold
      h3flgu = h3old

      setvar= palloc(112,'TEMP2',0,2)

      end
