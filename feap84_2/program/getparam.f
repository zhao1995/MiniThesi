c$Id:$
      subroutine getparam(par, value, error, prt)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Get parameters from stored values in vvv

c      Inputs:
c         par      - Character name of parameter

c      Outputs:
c         value    - Current value of parameter 'par'
c         error    - No parameter found
c         prt      - Output parameter name and value if true
c-----[--.----+----.----+----.-----------------------------------------]

      implicit  none

      include  'conval.h'
      include  'iofile.h'

      logical   prt, error
      character par*(*), y*15
      integer   i,j,n, ial,izl,iau,izu,id,iq, i0,i9, le
      real*8    value

      save

      error = .true.

c     Set numeric and upper/lower case locations

      i0  = ichar('0')
      i9  = ichar('9')
      iq  = ichar('=')
      ial = ichar('a')
      izl = ichar('z')
      iau = ichar('A')
      izu = ichar('Z')
      id  = ial - iau

c     Store one or two letters

      le = min(15,len(par))
      y = ' '
      n = 0
      do i = 1,le
        if(par(i:i).ne.' ') then
          n      = n + 1
          y(n:n) = par(i:i)
        endif
      end do ! i

      if(n.gt.2) then
        write(iow,3000) par
        if(ior.lt.0) then
          write(*,3000) par
        endif
        return
      endif

c     Converts all characters to lower case

      do i = 1,n
        j = ichar(y(i:i))
        if(j.ge.iau .and. j.le.izu) then
          y(i:i) = char(j+id)
        endif
      end do ! i

c     Check for blank character or null character = blank line

      if(y(1:1).eq.' '.or.ichar(y(1:1)).eq.0) then
        error = .true.
        return
      endif

c     Locate correct location for the addition

      n = ichar(y(1:1)) - ial + 1
      if(n.le.26) then
        j = ichar(y(2:2))
        if(    j.ge.ial.and. j.le.izl) then  ! Letter
          j = j - ial + 1
        elseif(j.ge.i0 .and. j.le.i9 ) then  ! Numeral
          j = j - i0 + 27
        endif

c       Set value

        value = vvv(n,j)
        error = .false.

        if(prt) then
          write(iow,2000) y(1:1),y(2:2),value
          if(ior.lt.0) then
            write(*,2000) y(1:1),y(2:2),value
          endif
        endif
      else
        write(iow,3001) par
        if(ior.lt.0) then
          write(*,3001) par
        endif
      endif

c     Formats

2000  format(5x,'Constant ',a1,a1,' = ',e15.8)

3000  format('-->ERROR: Parameter ',a,' longer than two characters')
3001  format('-->ERROR: Illegal parameter name: ',a)

      end
