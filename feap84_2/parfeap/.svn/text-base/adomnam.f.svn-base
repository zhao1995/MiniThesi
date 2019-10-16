c$Id:$
      subroutine adomnam(fin,fout,d)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Set file names for outputs of parallel parts

c      Inputs:
c         fin     - Character string of main input file
c         d       - Integer number of domain

c      Outputs:
c         fout    - Character string of output file for domain 'd'
c-----[--.----+----.----+----.-----------------------------------------]
      implicit   none

      integer    d, n
      character  fin*128, fout*128, fext*5

c     Check if underscore '_' exists in name -- strip off

      fout = fin
      do n = 128,1,-1
        if(fin(n:n).eq.'_') then
          fout(n:128) = ' '
          exit
        endif
      end do ! n

c     Add extender for 'd'

      fext(1:5) = ' '
      if(d.lt.10) then
        write(fext,'(a3,i1)') '000',d
      elseif(d.lt.100) then
        write(fext,'(a2,i2)') '00',d
      elseif(d.lt.1000) then
        write(fext,'(a1,i3)') '0',d
      elseif(d.lt.10000) then
        write(fext,'(i4)') d
      endif

      call addext(fout,fext,128,5)

c     Replace '.' by '_' in file name

      n = index(fout,'.')
      if(n.gt.0) then
        fout(n:n) = '_'
      endif

      end
