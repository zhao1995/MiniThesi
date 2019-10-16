c$Id:$
      subroutine xsmth(ic,indc,x,ndm,numnp,nsmth)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Smooth

c      Inputs:

c      Outputs:
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      integer   ndm,numnp,nsmth, i,j,m,n,in1,in2, ic(*),indc(*)
      real*8    x(ndm,*),xx(3),div

      save

      do m = 1,nsmth
        in1 = 0
        do n = 1,numnp
          in2 = ic(n)
          if(in2.gt.in1+1) then
            do j = 1,ndm
              xx(j) = x(j,indc(in1+1))
            end do ! j
            do i = in1+2,in2
              do j = 1,ndm
                xx(j) = xx(j) + x(j,indc(i))
              end do ! j
            end do ! i
            div = 1.d0/dble(in2-in1)
            do j = 1,ndm
              x(j,n) = xx(j)*div
            end do ! j
          endif
          in1 = in2
        end do ! n
      end do ! m

      end
