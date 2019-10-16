c$Id:$
      subroutine pcrbrd(x,ndtyp,prt)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Input all coordinates in binary mode

c      Inputs:
c         prt      - Print results if true

c      Outputs:
c         x(ndm,*) - Nodal coordinates
c         ndtyp(*) - Active node identifier
c-----[--.----+----.----+----.-----------------------------------------]
      implicit   none

      include   'cdata.h'
      include   'chdata.h'
      include   'comfil.h'
      include   'sdata.h'
      include   'iodata.h'
      include   'iofile.h'

      logical    prt
      character  fname*120
      integer    loceq, i,j,n
      integer    ndtyp(numnp)
      real*8     x(ndm,numnp)

c     Set filename for read

      loceq = index(xxx,'=')
      do i = loceq+1,120
        if(xxx(i:i).ne.' ') then
          do j = i+1,120
            if(xxx(j:j).eq.' ') then
             fname(1:j-i) = xxx(i:j-1)
             go to 100
            end if
          end do ! j
        endif
      end do ! i
      j = 120
100   write(iow,*) ' Filename = ',fname(1:j-i)
      write(  *,*) ' Filename = ',fname(1:j-i)

c     Open file for binary input

      open(unit = ios,file = fname(1:j-i), form = 'unformatted')
      read (ios)  x
      close(ios)

      do n = 1,numnp
        ndtyp(n) = 0
      end do ! n

c     Output data

      if(prt) then
        call prtitl(prt)
        write(iow,2000)  (i,i=1,ndm)
        do n = 1,numnp
          write(iow,2001) n,(x(i,n),i=1,ndm)
        end do ! n
      endif

c     Formats

2000  format(5x,'Nodal Coordinates'//6x,'Node',3(i6,' Coord':))
2001  format(i10,1p,3e12.4)

      end
