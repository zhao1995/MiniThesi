c$Id:$
      subroutine pelbrd(ix,rben,prt)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Add set of last element number to last_elm       29/01/2012
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Input all element connections in binary mode

c      Inputs:
c         prt        - Print results if ture

c      Outputs:
c         ix(nen1,*) - Element connection data
c         rben(*)    - Rigid body/flexible element indicator
c-----[--.----+----.----+----.-----------------------------------------]
      implicit   none

      include   'cdata.h'
      include   'chdata.h'
      include   'comfil.h'
      include   'sdata.h'
      include   'iodata.h'
      include   'iofile.h'
      include   'p_ptname.h'
      include   'region.h'
      include   'rigid2.h'

      logical    prt
      character  fname*120
      integer    loceq, i,j,n
      integer    ix(nen1,numel),rben(*)

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
      read (ios)  ix
      close(ios)
      last_elm = numel

c     Set region and rigid body indicators

      do n = 1,numel
        rben(n)      = nrigid
        ix(nen1-1,n) = nreg
      end do ! n

c     Output data
      if(prt) then
        call prtitl(prt)
        write(iow,2000)  (i,i=1,nen)
        do n = 1,numel
          write(iow,2001) n,ix(nen1,n),ix(nen1-1,n),(ix(i,n),i=1,nen)
        end do ! n
      endif

c     Formats

2000  format(5x,'E l e m e n t s'//3x,'Elmt Mat Reg',8(i3,' Node':),
     &      (15x,8(i3,' Node':)))
2001  format(i7,2i4,8i8:/(15x,8i8))

      end
