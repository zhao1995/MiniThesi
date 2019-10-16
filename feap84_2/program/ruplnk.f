c$Id:$
      subroutine ruplnk(u,ud,x,rixt,rlink,nvec,ndm,ndf,numnp,flg)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Rigid body update and link controls

c      Inputs:

c      Outputs:
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      logical   flg
      integer   nvec,ndm,ndf,numnp, i,j,nn, master

      integer   rixt(*),rlink(ndf,*)
      real*8    u(ndf,numnp,*),ud(ndf,numnp,*),x(ndm,*), xx(3)

      if(ndm.eq.3 .and. ndf.ge.6) then

        do nn = 1,numnp
          master = rixt(nn)
          if(master.gt.0) then
            do i = 1,3
              xx(i) = x(i,master) - x(i,nn)
            end do ! i

c           Solution vector updates

            do j = 1,6

              if(rlink(j,nn).eq.0) then
                do i = 1,nvec
                  u(j,nn,i) = u(j,master,i)
                end do ! i
                if(flg) then
                  do i = 1,5
                    ud(j,nn,i) = ud(j,master,i)
                  end do ! i
                endif
              endif

            end do ! j

            if(rlink(4,nn).eq.0) then
              do i = 1,nvec
                u(2,nn,i) = u(2,nn,i) + xx(3)*u(4,master,i)
                u(3,nn,i) = u(3,nn,i) - xx(2)*u(4,master,i)
              end do ! i
              if(flg) then
                do i = 1,5
                  ud(2,nn,i) = ud(2,nn,i) + xx(3)*ud(4,master,i)
                  ud(3,nn,i) = ud(3,nn,i) - xx(2)*ud(4,master,i)
                end do ! i
              endif
            endif

            if(rlink(5,nn).eq.0) then
              do i = 1,nvec
                u(1,nn,i) = u(1,nn,i) - xx(3)*u(5,master,i)
                u(3,nn,i) = u(3,nn,i) + xx(1)*u(5,master,i)
              end do ! i
              if(flg) then
                do i = 1,5
                  ud(1,nn,i) = ud(1,nn,i) - xx(3)*ud(5,master,i)
                  ud(3,nn,i) = ud(3,nn,i) + xx(1)*ud(5,master,i)
                end do ! i
              endif
            endif

            if(rlink(6,nn).eq.0) then
              do i = 1,nvec
                u(1,nn,i) = u(1,nn,i) + xx(2)*u(6,master,i)
                u(2,nn,i) = u(2,nn,i) - xx(1)*u(6,master,i)
              end do ! i
              if(flg) then
                do i = 1,5
                  ud(1,nn,i) = ud(1,nn,i) + xx(2)*ud(6,master,i)
                  ud(2,nn,i) = ud(2,nn,i) - xx(1)*ud(6,master,i)
                end do ! i
              endif
            endif

          endif

        end do ! nn
      else

        write(*,*) ' *WARNING* Master-Slave not done'

      end if ! ndm & ndf test

      end
