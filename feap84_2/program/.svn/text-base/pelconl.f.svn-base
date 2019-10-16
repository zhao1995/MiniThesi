c$Id:$
      subroutine pelconl(i,ie,ix,lagre, ic, ielc,neql)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Correct set of number in kk                      19/09/2011
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose:  Lagrange multipler for finite elements

c      Inputs:
c         i          -
c         ie(nie,*)  -  Element multiplier locator
c         ix(*)      -  List of nodes connected to each element
c         lagre(*)   -  Multipler equation numbers
c         ic(*)      -  Pointer array

c      Outputs:
c         ielc(*)    -  Holds set of elements connected to each node.
c         neql       -  Largest equation number in list
c-----[--.----+----.----+----.-----------------------------------------]
      implicit   none

      include   'cdata.h'
      include   'cdat1.h'
      include   'sdata.h'

      integer    i,ma,n,nlm,kk,kp, neql
      integer    ie(nie,*), ix(*), lagre(*), ic(*), ielc(*)

      ma  = ix(nen1)       ! Material set number
      nlm = ie(nie-8,ma)
      do n = 1,nlm
        kk  = lagre(ix(nen+4)) + n - 1
        neql = max(neql,kk)
        if(kk.gt.0) then
          kp = ic(kk)
100       if(ielc(kp).eq.0) go to 110
          kp = kp - 1
          go to 100
110       ielc(kp) =  i
        endif ! kk > 0
      end do ! n

      end
