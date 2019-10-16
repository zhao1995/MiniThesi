c$Id:$
      subroutine pcangle(pl,x,ang,ndm,numnp,numprt,prt,prth)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c     1. Use 'pi' from 'pconstant.h'                        14/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Set coordinate angle values

c      Inputs:

c      Outputs:
c-----[--.----+----.----+----.-----------------------------------------]
      implicit   none

      include   'iofile.h'
      include   'pconstant.h'

      logical    prt,prth
      integer    ndm,numnp,numprt, n,n1,n2,n3
      real*8     con,pl(*),x(ndm,*),ang(*)

      save

      n2 = nint(pl(1))
      n3 = nint(pl(2))
      n1 = min(n2,n3)
      n2 = max(n2,n3)
      n3 = nint(pl(3))
      if(n2.eq.0) then
        n1 = 1
        n2 = numnp
        n3 = 1
      else
        n1 = max(n1,1)
        n3 = max(n3,1)
      end if

      if(prth .and. numprt.le.0) then
        call prtitl(prth)
        write(iow,2000)
        if(ior.lt.0) then
          write(*,2000)
        endif
        numprt = 50
      endif
      con = 180.d0/pi
      do n = n1,n2,n3
        if(max(abs(x(1,n)),abs(x(2,n))).gt.0.0d0) then
          ang(n) = con*atan2(x(2,n),x(1,n))
          if(prt) then
            write(iow,2001) n,ang(n)
            if(ior.lt.0) then
              write(*,2001) n,ang(n)
            endif
            numprt = numprt - 1
          endif
        endif
      end do ! n

c     Formats

2000  format(6x,'Node      Angle')
2001  format(i10,1p,1e12.4)

      end
