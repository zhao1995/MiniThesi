c$Id:$
      subroutine pconst(prt)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Input parameter expressions:  let = expression

c      Inputs:
c         prt    - Print input values if true

c      Outputs:
c         Values of parameters a-z are stored in array vvv(26)
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'comfil.h'
      include  'iofile.h'
      include  'ioincl.h'
      include  'iodata.h'

      logical   prt, redo, pconset, lopn
      integer   i

      save

c     Input record from file "ior" or keyboard "*"

      if(prt) then
        write(iow,2000)
      endif
      inquire(unit=iwd, opened = lopn)
      if(lopn) write(iwd,'(a)') 'Parameters'

1     redo = .false.
      record = ' '
      if(ior.gt.0) then
        read (ior,1000,err=901,end=902) record
        irecrd(isf) = irecrd(isf) + 1
      else
        write(*,3000)
        call pprint('  -->')
        read (  *,1000,err=901,end=902) record
      endif

      if(lopn) then
        do i = 255,1,-1
          if(record(i:i).ne.' ') go to 100
        end do ! i
        i = 1
100     write(iwd,'(a)') record(1:i)
      endif

      redo = pconset(prt)
      if(redo) go to 1
      return

c     Error on read

901   call  errclr ('PCONST')
      if (ior.lt.0)  goto 1
      return

c     EOF encountered

902   return

c     Formats

 1000 format(a)

 2000 format(/'  C o n s t a n t    V a l u e s'/1x)

 3000 format(' Use "list" to give current values - <CR> to exit'/
     &       ' Input: letter=expression (no blanks)')

      end
