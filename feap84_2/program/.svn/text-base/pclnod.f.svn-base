c$Id:$
      subroutine pclnod(xx,x,ntyp, ndm,numnp, ncn)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c     Purpose: Determine number of node closest to specified coordinate
c              in array xx(ndm)

c     Inputs:
c       xx(ndm)  - Specified coordinate array
c       x(ndm,*) - Coordinates of nodes
c       ntyp(*)  - Active node indicator:  > 0 active; < 0 inactive
c       ndm      - Space dimension of mesh
c       numnp    - Number of nodes to serch

c     Outputs:
c       ncn      - Closest node
c-----[--.----+----.----+----.-----------------------------------------]
      implicit   none

      include   'xtout.h'

      logical    clflg
      integer    ndm,numnp, ncn, n, ntyp(*)
      real*8     xx(*),x(ndm,numnp), dotx,xmn,tmn

      save

      clflg = .false.
      do n = 1,numnp
        if(ntyp(n).ge.0) then
          tmn = dotx(xx,x(1,n),ndm)
          if(clflg) then
            if(tmn.lt.xmn) then
              xmn  = tmn
              ncn  = n
            elseif(tmn.gt.xmn) then
              xtol = tmn
            endif
          else
            xmn   =  tmn
            xtol  =  tmn
            ncn   =  n
            clflg = .true.
          endif
        endif
      end do ! n

      end
