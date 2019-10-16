c$Id:$
      subroutine plot2d(ie,ix,ip,x,xl,nie,ndm,nen,nen1,
     &                  numnp,numel,n1,n2)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Change 'call plot9' to 'call plotel' for general 31/08/2008
c          element plots.
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Two dimensional mesh plot routine

c      Inputs:
c         ie(nie,*) - Assembly data for material sets
c         ix(nen1,*)- Element nodal connections
c         ip(8,*)   - Sorted element order for each quadrant
c         x(ndm,*)  - Nodal coordinates
c         ndm       - Dimension of x  array
c         nen       - Number of nodes on element
c         nen1      - Dimension of ix array
c         numnp     - Number of nodes
c         numel     - Number of elements
c         n1        - Color for plots
c         n2        - Outline indicator

c      Scratch:
c         xl(ndm,*) - Element nodal coordinates

c      Outputs:
c         none      - Plot output to screen/file
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'eldata.h'
      include  'pbody.h'
      include  'pdata2.h'
      include  'pdatay.h'

      integer   nie, ndm, nen, nen1, numnp, numel, n1,n2, isp
      integer   i, j, ii, nsy, nn, nume
      integer   ie(nie,*),ix(nen1,*), ip(8,numel)
      real*8    x(ndm,*),xl(ndm,*)

      save

c     Loop over elements to draw mesh

      do nsy = 1,nsym
        lsym = isym(nsy)
        call pltsym(x,ndm,numnp,lsym)
        nume = nfac(lsym)
        do nn = 1,nume

          n  = ip(lsym,nn)

          if(ix(nen1-1,n).ge.nreg1 .and. ix(nen1-1,n).le.nreg2) then

            ma = ix(nen1,n)

c           Plot correct material number: ma > 0 active material

            if( ma.gt.0 .and. (maplt.eq.0 .or. ma.eq.maplt)) then
              if(n1.eq.0) then
                icolr = ma + 1
              elseif(n1.lt.0) then
                icolr = mod(ix(nen1-1,n),11) + 1
              else
                icolr = mod(n1-1,16) + 1
              endif
              call pppcol(icolr,0)

              do i = 1,nen
                ii = abs(ix(i,n))
                if(ii.gt.0) then
                  nel = i
                  do j = 1,ndm
                    xl(j,i) = x(j,ii)
                  end do ! j
                else
                  do j = 1,ndm
                    xl(j,i) = 0.0d0
                  end do ! j
                endif
              end do ! i

c             Check for a line element

              if(ix(1,n).eq.ix(4,n) .and. ix(2,n).eq.ix(3,n)) nel = 2

              if(n2.ge.0) then
                isp =  1
              else
                isp = -1
              endif
              call plotel(ie(1,ma),ix(1,n),xl,ndm,nel,isp)

            endif
          endif
        end do ! nn

        call pltsym(x,ndm,numnp,lsym)

      end do ! nsy

      end
