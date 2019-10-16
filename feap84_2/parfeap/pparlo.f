c$Id:$
      subroutine pparlo(ld,eq,ix,id,ldsize,lnen)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Set local arrays for each element

c      Inputs:
c        eq(ndf,*)  - Equation numbers
c        ix(*)      - Node connection list for current element
c        id(*)      - Boundary code indicators

c      Outputs:
c        ld(nst,i)  - Element equation assembly information
c                     i = 1: local processor global equation numbers ;
c                     i = 2: all   processor global equation numbers ;
c                     i = 3: assembly start column ;
c                     i = 4: assembly stop  column ;
c                     i = 5: active   equation = 0 ;
c                            boundary equation = 1 ;
c                     i = 6: local processor equation numbers ;
c-----[--.----+----.----+----.-----------------------------------------]
      implicit   none

      include   'cdata.h'
      include   'sdata.h'

      include   'pfeapb.h'

      integer    i,j, ni, ldsize, lnen
      integer    ld(ldsize,*), eq(ndf,*), ix(*), id(ndf,*)

c     Set row/column global equation numbers

      ni = 0
      do i = 1,lnen

        if(ix(i).gt.0) then
          do j = 1,ndf                ! Columns
            ld(ni+j,6) = ld(ni+j,1)
            ld(ni+j,2) = eq(j,ix(i))
          end do ! j

          if(ix(i).le.numpn) then     ! Rows
            do j = 1,ndf
              ld(ni+j,1) = eq(j,ix(i))
              ld(ni+j,3) = 1
              ld(ni+j,4) = nst
              if(id(j,ix(i)).eq.0) then
                ld(ni+j,5) = 0           ! Active equation
              else
                ld(ni+j,5) = 1           ! Known solution equation
              endif
            end do ! j
          else
            do j = 1,ndf
              ld(ni+j,1) = 0
              ld(ni+j,3) = 1
              ld(ni+j,4) = 0
              if(id(j,ix(i)).eq.0) then  ! Set equation modify code
                ld(ni+j,5) = 0
              else
                ld(ni+j,5) = 1
              endif
            end do ! j
          endif
        else
          do j = 1,ndf
            ld(ni+j,1) = 0
            ld(ni+j,2) = 0
            ld(ni+j,3) = 1
            ld(ni+j,4) = 0
            ld(ni+j,5) = 1
            ld(ni+j,6) = 0
          end do ! j

        endif ! ix(i) > 0
        ni = ni + ndf
      end do ! i

      end
