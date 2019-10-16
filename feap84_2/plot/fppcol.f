c$Id:$
      subroutine fppcol(icol,jcolr,ncolr,isw)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    11/10/2007
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Set hardcopy color information

c      Inputs:
c         icol     - Color parameter
c         ncolr    - Total number of colors available
c         isw      - Switch:

c      Outputs:
c         none     - Set values returned in common blocks
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'plflag.h'
      include  'psdat2.h'

      integer   icol, isw, jj, jcolr, ncolr

      save

c     Line for postscript

      if(isw.eq.0 .or. isw.eq.1) then

c       Reverse line color for black background

        if(blk) then
          clin = 'G1 '
        else
          clin = 'G0 '
        endif
      endif

c     Gray scale for postscript

      if(isw.eq.0 .or. isw.eq.2) then

        if (icol .lt. 0) then
          jj = 1
        elseif (icol .eq. 0) then
          jj = jcolr
        else
          jj = icol
        endif
        if(jj.lt.10) then
          write(cvar,'(a2,i1)') ' G',jj
        else
          write(cvar,'(a2,a1)') ' G',char(87+jj)
        endif

c       Color for postscript

        if(icol.le.0) then
          if(blk) then
            colv = '0 '
          else
            colv = '1 '
          endif
        else
          jj = max( 1, min(ncolr, icol) )
          if(jj.lt.10) then
            write( colv, '(i1,1x)' ) jj
          else
            write( colv, '(a1,1x)' ) char(87+jj)
          endif
        endif
        if(isw.eq.0) then
          if(jj.lt.10) then
            if(.not.blk .and. jj.eq.1) jj = 0
            write(clin,'(a1,i1,a1)') 'h',jj,' '
          else
            write(clin,     '(3a1)') 'h',char(87+jj),' '
          endif
        endif
      endif

      end
