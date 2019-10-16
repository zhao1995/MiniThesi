c$Id:$
      subroutine pltelm(x,ie,ix,scale,nie,ndm,nen1,n1,n2)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c      1. Add test on 'ma > 0' for plots                    19/07/2013
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Place element numbers on plots of mesh

c      Inputs:
c         x(ndm,*)  - Nodal coordinates of mesh
c         ie(nie,*) - Assembly data for material sets
c         ix(nen1,*)- Element nodal connections
c         scale     - Plot scale factor
c         nie       - Dimension of ie array
c         ndm       - Dimension of x array
c         nen1      - Dimension of ix array
c         n1        - First element number to display
c         n2        - Last element number to display

c      Outputs:
c         none      - Plot outputs to screen/file
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'cdata.h'
      include  'pbody.h'
      include  'pdata4.h'
      include  'plflag.h'
      include  'ppers.h'

      logical   zoom
      integer   nie,ndm,nen1, i,j,n,ii,jj,nn,ma,nx,ny,nz, n1,n2
      integer   ie(nie,*),ix(nen1,*),iplt(30), pstyp,nel
      real*8    scale,dx1, x(ndm,*),xx(3)

      save

c     Write element labels

      dx1 = .005d0/scale
      if(kpers.ne.0) then
        nx = 0
        ny = 0
        nz = 0
      else
        nx = 5
        ny = 1
        nz = 0
      endif
      do n = n1,n2
        ma = ix(nen1,n)
        if(ma.gt.0) then
          pstyp = ie(1,ma)
          if((ix(nen1-1,n).ge.nreg1 .and. ix(nen1-1,n).le.nreg2) .and.
     &       (maplt.eq.0 .or. ma.eq.maplt) .and. pstyp.ne.0) then
            do i = nen,1,-1
              if(ix(i,n).gt.0) then
                nel = i
                exit
              endif
            end do ! i
            call plftyp(pstyp,nel,ie(nie-1,ma))
            call pltord(ix(1,n),ie(nie-1,ma), nn,iplt)
            xx(1) = 0.0d0
            xx(2) = 0.0d0
            xx(3) = 0.0d0
            jj = 0
            nn = max(1,nn-1)
            do i = 1,nn
              j  = iplt(i)
              ii = ix(j,n)
              if(ii.gt.0) then
                jj = jj + 1
                xx(1) = xx(1) + x(1,ii)
                if(ndm.ge.2) xx(2) = xx(2) + x(2,ii)
                if(ndm.ge.3) xx(3) = xx(3) + x(3,ii)
              endif
            end do ! i
            if(jj.gt.0) then
              xx(1) = xx(1)/jj
              xx(2) = xx(2)/jj
              xx(3) = xx(3)/jj
              if(zoom(xx(1),ndm)) then
                call plotl(xx(1)-dx1*nx,xx(2)-dx1*ny,xx(3)-dx1*nz,3)
                if(clip) call plabl(n)
              endif
            endif
          endif
        endif ! ma > 0
      end do ! n

      end
