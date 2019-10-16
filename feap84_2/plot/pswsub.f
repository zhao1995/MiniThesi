c$Id:$
      subroutine pswsub(ix, extnd, ip, x, nxd,nxn, swang, nsinc )

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--+---------+---------+---------+---------+---------+---------+-]
c      Purpose: Plot sweep of 2-d mesh through angle 'swang'

c      Inputs:
c        ix(nxd,*) - Nodal connection array
c        extnd(*)  - External node numbers
c        ip(8,*)   - Element numbers in quadrants
c        x(3,*)    - Nodal coordinates
c        nxd       - Dimension of ix array
c        nxn       - Number of nodes on array
c        swang     - Sweep angle
c        nsinc     - Sweep increments

c      Outputs:
c        to screen
c-----[--+---------+---------+---------+---------+---------+---------+-]

      implicit   none

      include   'cdata.h'
      include   'pconstant.h'
      include   'pdatay.h'
      include   'sdata.h'

      integer    nxd,nxn
      integer    ix(nxd,*), extnd(numnp), ip(8,*), nsinc
      integer    i,n,ne, nn, nume, nod1,nod2, nel, nsy
      real*8     x(3,numnp), xx(3), swang, dang, rad, ang, rwang

      save

      rwang = swang*pi/180.0d0
      dang  = rwang/dble(nsinc)

      do nsy = 1,nsym
        lsym = isym(nsy)
        nume = nfac(lsym)
        call pltsym(x,ndm,numnp,lsym)
        do ne = 1,nume
          n = ip(lsym,ne)
          do i = nxn,1,-1
            if(ix(i,n).ne.0) then
              nod2 = ix(i,n)
              nel  = i
              go to 100
            endif
          end do ! i
100       continue
          do i = 1,nel
            if(ix(i,n).gt.0) then
              nod1 = ix(i,n)
              if(extnd(nod1).gt.0) then
                xx(1) = x(1,nod1)
                xx(2) = x(2,nod1)
                xx(3) = x(3,nod1)
                rad   = sqrt(xx(1)*xx(1) + xx(3)*xx(3))
                ang   = atan2(xx(3),xx(1))
                call plotl(xx(1),xx(2),xx(3),3)
                do nn = 1,nsinc
                  xx(1) = rad*cos(ang-dble(nn)*dang)
                  xx(3) = rad*sin(ang-dble(nn)*dang)
                  call plotl(xx(1),xx(2),xx(3),2)
                end do ! nn
                if(extnd(nod2).gt.0) then
                  call plotl(xx(1),xx(2),xx(3),3)
                  rad   = sqrt(x(1,nod2)**2 + x(3,nod2)**2)
                  ang   = atan2(x(3,nod2),x(1,nod2))
                  xx(1) = rad*cos(ang-rwang)
                  xx(3) = rad*sin(ang-rwang)
                  call plotl(xx(1),x(2,nod2),xx(3),2)
                endif
              endif
            endif
            nod2 = nod1
          end do ! i

        end do ! ne

        call pltsym(x,ndm,numnp,lsym)

      end do ! nsy

      end
