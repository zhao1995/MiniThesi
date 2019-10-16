c$Id:$
      subroutine aparexp(x,xs,v,nex,error)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Change character x(*)*1 to x*(*) also xs         09/01/2012
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Identify parenthetical expressions and evaluate

c      Inputs:
c         x(*)     - String containing expression to evaluate

c      Scratch:
c         xs(*)    - Array used to temporarily store expression
c         v(*)     - Array to hold values

c      Outputs:
c         x(*)     - Expression replaced by upper case letter
c         nex      - Number of upper case letters used
c         error    - Flag, true if error occurs

c      Common returns:
c         www(*)   - Upper case letters with values assigned
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'conval.h'

      logical   error
      character x*(*),xs*(*)
      integer   i,j,k,l, i1,i2,nex
      real*8    val, v(*)

      save

c     Find parenthetical expressions and remove

      do i = 1,125
        if(x(i:i).eq.'(') then
          i1 = i + 1
          do j = i1,125
            if(x(j:j).eq.'(') then
              call errclr('PAREXP')
              call plstop()
            elseif(x(j:j).eq.')') then
              do l = 1,j-i+1
                xs(l:l) = ' '
              end do ! l
              i2 = j - 1
              if(i2.lt.i1) then
                call errclr('PAREXP')
                call plstop()
              else
                k = 0
                do l = i1,i2
                  k = k + 1
                  xs(k:k) = x(l:l)
                  x(l:l)  = ' '
                end do ! l
                x(i2+1:i2+1)  = ' '

c               Evaluate expression in parenthesis

                call aevalex(xs,v,val,k,error)
                if(error) return
                nex = nex + 1
                www(nex) = val

c               Put upper case letter in expression & close up remainder

                x(i:i) = char(nex +64)
                i2 = i2 -i1 + 2
                do l = i1,125
                  x(l:l) = ' '
                  if(l+i2.le.125) then
                    x(l:l) = x(l+i2:l+i2)
                  endif
                end do ! l
              endif
              go to 100
            endif
          end do ! j
100       continue
        endif
      end do ! i

      end
