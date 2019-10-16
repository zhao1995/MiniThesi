c$Id:$
      subroutine ptsumpl(x,r)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Sum reaction components for output

c      Inputs:
c         x(ndm,*)  - Nodal coordinates
c         r(ndf,*)  - Nodal reactions

c      Outputs:
c-----[--.----+----.----+----.-----------------------------------------]
      implicit   none

      include   'cdata.h'
      include   'debugs.h'
      include   'sdata.h'
      include   'ptdat9.h'

      integer    n,nt,i,j, nsum
      real*8     x(ndm,numnp),r(ndf,numnp), xx,xt

      save

      do nt = 1,ntplts

c       Reaction sums on coordinates

        i  = max(1,min(itpl(1,nt),ndf))  ! Reaction component
        j  = max(1,min(itpl(2,nt),ndm))  ! Coordinate direction
        xx = tpld(1,nt)                  ! Coordinate value
        xt = tpld(2,nt)                  ! Coordinate tolerance
        if(xt.eq.0.0d0) then
          xt = 0.01d0
        endif
        nsum = 0
        do n = 1,numnp
          if(abs(x(j,n) - xx).le.xt) then
            tpl(nt) = tpl(nt) - r(i,n)
            nsum    = nsum + 1
          endif
        end do ! n
        if(debug) then
          write(*,*) ' REACTION SUMS: NSUM =',nsum
        endif
      end do ! nt

      end
