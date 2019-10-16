c$Id:$
      subroutine elmt34(d,ul,xl,ix,tl,s,p,ndf,ndm,nst,isw)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c     Purpose:

c     Inputs:

c     Outputs:
c-----[--.----+----.----+----.-----------------------------------------]
      implicit none

      integer ix(*), ndf , ndm , nst , isw
      real*8  d(*), ul(*), xl(*), tl(*), s(*), p(*)

      if(isw.gt.0) write(*,2000)
2000  format('    Elmt 34: *WARNING* Dummy subroutine called')
      end
