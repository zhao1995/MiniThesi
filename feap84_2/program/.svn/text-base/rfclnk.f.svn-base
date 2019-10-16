c$Id:$
      subroutine rfclnk(f1,x,rixt,rlink)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Form nodal load vector for current time

c      Inputs:
c         x(ndm,*)     - Nodal coordinates
c         rixt(*)      - List of master nodes
c         rlink(ndf,*) - Link indicators

c      Outputs:
c         f1(ndf,*)    - Total nodal load for t_n+1
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'cdata.h'
      include  'sdata.h'

      integer   ii, nn, master
      integer   rixt(*),rlink(ndf,*)
      real*8    f1(ndf,*),x(ndm,*), xx(3)

      save

      if(ndm.eq.3 .and. ndf.ge.6) then
        do nn = 1,numnp
          master = rixt(nn)

          if(master.gt.0) then
            do ii = 1,3
              xx(ii) = x(ii,master) - x(ii,nn)
            end do ! ii
            if(rlink(4,nn).eq.0) then
              f1(4,master) = f1(4,master) + xx(3)*f1(2,nn)
     &                                    - xx(2)*f1(3,nn)
            endif

            if(rlink(5,nn).eq.0) then
              f1(5,master) = f1(5,master) + xx(1)*f1(3,nn)
     &                                    - xx(3)*f1(1,nn)
            endif

            if(rlink(6,nn).eq.0) then
              f1(6,master) = f1(6,master) + xx(2)*f1(1,nn)
     &                                    - xx(1)*f1(2,nn)
            endif
          endif
        end do ! nn
      else
        write(*,*) ' *WARNING* Master-Slave dimension error'
      endif

      end
