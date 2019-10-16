c$Id:$
      subroutine setcount(fnam,nprob,i)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Add multiproblem counter to name

c      Inputs:
c         fnam      - Root file name
c         nprob     - Previous problem name
c         i         - First character to append

c      Outputs:
c         nprob     - Current problem name
c         fnam      - Final file name for current problem 'nprob'
c-----[--.----+----.----+----.-----------------------------------------]
      implicit   none

      character  fnam*(*)
      integer    nprob, i

c     Add problem counter to name

      nprob = nprob + 1
      write(fnam(i:i+2),'(a)') '000'
      if(nprob.lt.10) then
        write(fnam(i+2:i+2),'(i1)') nprob
      elseif(nprob.lt.100) then
        write(fnam(i+1:i+2),'(i2)') nprob
      elseif(nprob.lt.1001) then
        write(fnam(  i:i+2),'(i3)') nprob
      else
        write(*,*) 'Exceed 999 limit on number of files'
      endif

      end
