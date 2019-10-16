c$Id:$
      subroutine fptplt(xs,ys,tx,nn,nctr)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Subroutine to place text into PostScript file

c      Inputs:
c         xs(2)     - X-location to start placing text
c         ys(2)     - Y-location to start placing text
c         tx        - Character array qith text
c         nn        - Number of characters in text
c         nctr      - Text counter value

c      Outputs:
c         none      - Output written to PostScript file
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'pdatap.h'
      include  'plflag.h'
      include  'plpost.h'
      include  'psdat2.h'
      include  'psdat5.h'
      include  'psdat6.h'

      integer   nn, j,nctr, x,y
      character tx(nn)*1, coord*10, txsiz*2
      real*4    xs(*), ys(*)

      save

c     Close out stroke if necessary

      if(lstrk) then
        call fppsin('s')
      endif
      call fppsdu()

c     Place text on plot in current size

      write(txsiz,'(a1,i1)') 'H',nsizt
      call fppsin(txsiz)
      call fppsdu()

c     Position text

      x = nint(5400.0*xs(2) + 360.0)
      y = nint(5400.0*ys(2) + 360.0)

      xll = min(x,xll)
      yll = min(y,yll)
      xur = max(x,xur) + nint(2.5*float(nn))
      yur = max(y,yur) + 7

      write(coord,'(i4,1x,i4,1x)') x,y

      call fppsin( coord//'m ')
      call fppsdu()

      if (nctr .eq. 1) then
        call fppsin('(')

        do j = 1,nn
          call fppsin(tx(j))
        end do ! j

        call fppsin(') w')
        call fppsdu()
      endif
      call fppsin('(')

      do j = 1,nn
        call fppsin(tx(j))
      end do ! j

      call fppsin(') '//clin//'show')
      call fppsdu()

      oclin = ' '

      lstrk = .false.
      lfill = .false.

      xold = -9980
      yold = -9980

      end
