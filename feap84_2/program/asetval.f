c$Id:$
      subroutine asetval(xi,num, val)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Represent character constant by free inputs in string

c      Inputs:
c         xi        - Input string
c         num       - Length of character string
c      Outputs:
c         val       - Value of string
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'codat.h'
      include  'conval.h'
      include  'errchk.h'
      include  'iofile.h'

      logical   errco
      character xt*125,xs*125,xi*(*)
      integer   i, nex, num, nc
      real*8    val, v(25)

      save

      data      nc  / 25 /

c     Read value if no constants have been set

      if(.not.coflg) then
        xs(1:nc) = ' '
        nex      = nc - num
        do i = 1,num
          xs(i+nex:i+nex) = xi(i:i)
        end do ! i
        errco = .false.
        call apinval(xs,val,errco)
        if(errco) go to 60
        return
      endif

c     Find letter number for this parameter

60    xs(1:125) = ' '
      xs(1:num) = xi(1:num)

c     Evaluate expression

1     nex = 0
      errck = .false.
      call aparexp(xs,xt,v,nex,errck)
      if(errck) go to 150

      call apfuncs(xs,v,val,nex,errck)
      if(errck) go to 150
      return

c     An error has been detected in statement respecify

150   if(ior.lt.0) then
        write(*,2001) xi(1:num)
        call pprint('   >')
151     read (*,1000,err=152,end=153) xt
        goto  154
152     call  errclr ('ASETVAL')
        goto  151
153     call  endclr ('ASETVAL',xt)
154     write(iow,2002) xt
        call apcheck(1,xt,errck)
        xs( 1:124)  = xt(1:124)
        xs(125:125) = ' '
        errck = .false.
        go to 1

c     Error on a read

      else
        call  errclr ('SETVAL')
      endif

c     Formats

 1000 format(125a1)

 2001 format(2x,a1,' = ',124a1)

 2002 format('  Correction:>',125a1)

      end
