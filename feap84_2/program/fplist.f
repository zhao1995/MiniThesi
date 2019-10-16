c$Id:$
      subroutine fplist(ifp,ndf,numnp,prt,prth)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Set type of proportional load to use for each
c               degree of freedom.

c      Inputs:
c         ndf        - Number dof/node
c         numnp      - Number of nodes in mesh
c         prt        - Output generated data if true
c         prth       - Output title/header if true

c      Outputs:
c         ifp(ndf,*) - Proportional load table numbers for dof's
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'iofile.h'

      logical   prt,prth, errck, pinput
      integer   ndf,numnp, i,n, nn
      integer   ifp(ndf,numnp)
      real*8    td(16)

      save

      nn = 0
1     errck = pinput(td,ndf+1)

      n = nint(td(1))
      if(n.gt.0 .and. n.le.numnp) then

c       Print header

        if(mod(nn,50).eq.0) then
          call prtitl(prth)
          write(iow,2000) (i,i=1,ndf)
          if(ior.lt.1.and.prt) write(*,2000) (i,i=1,ndf)
        endif

        do i = 1,ndf
          ifp(i,n) = nint(td(i+1))
        end do ! i
        write(iow,2001) n,(ifp(i,n),i=1,ndf)
        if(ior.lt.0 .and. prt) then
          write(iow,2001) n,(ifp(i,n),i=1,ndf)
        end if
        nn = nn + 1
        go to 1
      end if

2000  format('   N o d a l   P r o p o r t i o n a l   L o a d'
     &  ,'   N u m b e r s'//6x,'Node',9(i3,'-dof'))
2001   format(i10,9i7)

      end
