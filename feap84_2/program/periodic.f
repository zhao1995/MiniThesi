c$Id:$
      subroutine periodic(prt)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    17/03/2011
c       1. Remove 'ndm' from argument list                  03/05/2013
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Input periodic gradu

c      Inputs:
c         prt     - Print flag
c-----[--.----+----.----+----.-----------------------------------------]
      implicit   none

      include   'elpers.h'
      include   'iofile.h'
      include   'setups.h'

      logical    prt, errck, pinput
      integer    i,j
      real*8     td(9)

c     Read type of input 'strain' or 'gradu'

      do i = 1,3
        errck = pinput(td,3)
        do j = 1,3
          gradu(i,j) = td(j)
        end do ! j
      end do ! i

c     Output

      if(prt) then
        write(iow,2000) (i,i=1,3),(i,(gradu(i,j),j=1,3),i=1,3)
        if(ior.lt.0) then
          write(*,2000) (i,i=1,3),(i,(gradu(i,j),j=1,3),i=1,3)
        endif
      endif

c     Activate periodic case

      perflg = .true.

c     Formats

2000  format(/5x,'P e r i o d i c   G r a d  u   V a l u e s'/
     &      /4x,3(10x,'x-',i1)/(5x,'u-',i1,1p,3e13.5))

      end
