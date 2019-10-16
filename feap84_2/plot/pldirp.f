c$Id:$
      subroutine pdirp(vex,x,numnp,ndm,xm,k3)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Draw vectors for beam and shell directors

c      Inputs:
c         vex(numnp,*) - Directions of directors
c         x(3,*)       - Nodal coordinates in deformed state
c         numnp        - Number of nodes in mesh
c         ndm          - Spatial dimension of mesh
c         xm           - Scale factor for vectors
c         k3           - Color for plot

c      Outputs:
c         none         - Plot outputs to screen/file
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'pdata3.h'
      include  'pointer.h'
      include  'comblk.h'

      include  'p_int.h'

      integer   ndm,ii,m,numnp,k3
      real*8    vex(numnp,*),xm,xx(4,3),x(3,*)

      save

      call pppcol(k3,1)

      fp(1)  =  npty - 1

      do ii=1,numnp

c       Get node location

        if(mr(fp(1)+ii).ge.0) then
          do m = 1,3
           xx(1,m) = x(m,ii)
          end do ! m

          if(ndm.eq.3) then
            if(.not.hide .or. mr(np(66)+ii-1).eq.1) then

c             Get endpoints of principal vectors scaled

              xx(2,1) = xx(1,1) + xm*vex(ii,1)
              xx(2,2) = xx(1,2) + xm*vex(ii,2)
              xx(2,3) = xx(1,3) + xm*vex(ii,3)
              xx(3,1) = xx(1,1) + xm*vex(ii,4)
              xx(3,2) = xx(1,2) + xm*vex(ii,5)
              xx(3,3) = xx(1,3) + xm*vex(ii,6)
              xx(4,1) = xx(1,1) + xm*vex(ii,7)
              xx(4,2) = xx(1,2) + xm*vex(ii,8)
              xx(4,3) = xx(1,3) + xm*vex(ii,9)

c             Plot vectors

              do m=1,3
                call pppcol(m,1)
                call plotl(xx(1,1),xx(1,2),xx(1,3),3)
                call plotl(xx(m+1,1),xx(m+1,2),xx(m+1,3),2)
              end do ! m
            endif

          elseif(ndm.eq.2) then

c           Get endpoints of principal vectors scaled

            xx(2,1) = xx(1,1) + xm*vex(ii,1)
            xx(2,2) = xx(1,2) + xm*vex(ii,2)
            xx(2,3) = xx(1,3)
            xx(3,1) = xx(1,1) + xm*vex(ii,3)
            xx(3,2) = xx(1,2) + xm*vex(ii,4)
            xx(3,3) = xx(1,3)

c           Plot vectors

            do m=1,2
              call plotl(xx(m+1,1),xx(m+1,2),xx(m+1,3),2)
            end do ! m

          endif
        endif
      end do ! ii

      end
