c$Id:$
      subroutine praxp(vex,x,numnp,ndm,xm,k1,k2,k3)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Draw vectors for principal stress axes

c      Inputs:
c         vex(numnp,*) - Directions of principal axes
c         x(3,*)       - Nodal coordinates in deformed state
c         numnp        - Number of nodes in mesh
c         ndm          - Spatial dimension of mesh
c         xm           - Scale factor for vectors
c         k1           - Number of principal stress vector to plot
c                                = 0 for all
c         k2           - Switch: < 0 plot negative components
c                                > 0 plot positive components
c                                = 0 plot all components
c         k3           - Color for plot

c      Outputs:
c         none         - Plot outputs to screen/file
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'pdata3.h'
      include  'pointer.h'
      include  'comblk.h'

      include  'p_int.h'

      integer   ndm,k1,ii,m,numnp,k2,f(3),k3
      real*8    vex(numnp,*),xm,xx(4,3),x(3,*),high,low
      logical   flg(3)

      save

      call pppcol(k3,1)

      flg(1) = .true.
      flg(2) = .true.
      flg(3) = .true.
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
              high    = max(vex(ii,10),vex(ii,11),vex(ii,12))
              low     = min(vex(ii,10),vex(ii,11),vex(ii,12))
              do m = 1,3
                if(vex(ii,9+m).eq.high) f(1) = m
                if(vex(ii,9+m).eq. low) f(3) = m
              end do ! m
              f(2) = 6 - f(1) - f(3)
              if(k2.lt.0) then
                do m=1,3
                  flg(m) = vex(ii,9+f(m)).lt.0.d0
                end do ! m
              elseif(k2.gt.0) then
                do m=1,3
                  flg(m) = vex(ii,9+f(m)).gt.0.d0
                end do ! m
              endif

c             Plot vectors

              do m=1,3
                if( (k1.eq.0 .or. k1.eq.f(m)) .and. flg(m)) then
                  call plotl(xx(1,1),xx(1,2),xx(1,3),3)
                  call plotl(xx(f(m)+1,1),xx(f(m)+1,2),xx(f(m)+1,3),2)
                endif
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
              if(k1.eq.m .or .k1.eq.0 .and.
     &          (k2.le.0 .and. vex(ii,4+m).le.0.0d0  .or.
     &           k2.ge.0 .and. vex(ii,4+m).ge.0.0d0)) then
                call plotl(xx(  1,1),xx(  1,2),xx(  1,3),3)
                call plotl(xx(m+1,1),xx(m+1,2),xx(m+1,3),2)
              endif
            end do ! m

          endif
        endif
      end do ! ii

      end
