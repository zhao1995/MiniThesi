c$Id:$
      subroutine plface(ix,ip,x,ndm,nen1,numnp,numel,iln,ct)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Driver routine to display visible faces of elements
c               in perspective views.

c      Inputs:
c         ix(nen1,*)- Nodal connection list
c         ip(8,*)   - Sorted order for symmetry plots
c         x(ndm,*)  - Nodal coordinates for plot
c         nen1      - Dimension of ix array
c         numnp     - Number of nodes
c         numel     - Number of elements/faces
c         iln(*)    - Line type data
c         ct        - Option to plot faces with negative normals

c      Outputs:
c         none      - Plot outputs to screen/file
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'pbody.h'
      include  'pdatay.h'

      integer   n, nsy, ndm, nen1, numnp, numel, nface
      integer   ix(nen1,numel), ip(8,numel), iln(2)
      real*8    x(ndm,numnp),ct

      save

c     Plot faces which are visible

      do nsy = 1,nsym

        lsym  = nsy
        call pltsym(x,ndm,numnp,lsym)
        nface = 0
        do n = 1,numel
          if(ix(nen1-1,n).ge.nreg1 .and. ix(nen1-1,n).le.nreg2 .and.
     &                                   ix(nen1,n).gt.0) then
            ip(nsy,n) = n
            call pfacev(ix(1,n),x,ndm,iln,ct,ip(nsy,n),nface)
          else
            ip(nsy,n) = 0
          endif
        end do ! n
        nfac(nsy) = nface
        call pltsym(x,ndm,numnp,lsym)

      end do ! nsy

c     Pack ip array

      do nsy = 1,nsym
        nface = 0
        do n = 1,numel
          if(ip(nsy,n).gt.0 .and. ix(nen1,n).gt.0) then
            nface = nface + 1
            ip(nsy,nface) = ip(nsy,n)
            if(n.gt.nface) then
              ip(nsy,n) = 0
            endif
          endif
        end do ! n
        if(nface.ne.nfac(nsy)) then
          write(*,*) ' *ERROR* PLFACE:',nface,nfac(nsy)
        endif
      end do ! nsy

c     Set line type to original

      call plline(iln)

      end
