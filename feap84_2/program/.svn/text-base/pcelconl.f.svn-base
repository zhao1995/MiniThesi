c$Id:$
      subroutine pcelconl(i,ix, ic, ielc,neql)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose:  Lagrange multipler for contact elements

c      Inputs:
c         i          -
c         ix(*)      -  List of nodes connected to each element
c         ic(*)      -  Pointer array

c      Outputs:
c         ielc(*)    -  Holds set of elements connected to each node.
c         neql       -  Largest equation number in list
c-----[--.----+----.----+----.-----------------------------------------]
      implicit   none

      include   'compac.h'

      integer    i,n,kk,kp, neql
      integer    ix(*),ic(*),ielc(*)

      do n = 0,ix(ncen1)-1
        kk   = ix(ncen+1) + n
        neql = max(neql,kk)
        if(kk.gt.0) then
          kp = ic(kk)
100       if(ielc(kp).eq.0) go to 110
          kp = kp - 1
          go to 100
110       ielc(kp) = -i
        endif ! kk > 0
      end do ! n

      end
