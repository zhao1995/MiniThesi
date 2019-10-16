c$Id:$
      subroutine protnd(rnod, prp,prv,theta,nn,xc,v0, x,u,id, ndlist)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Accumulate displacement u(*,*)                   06/04/2009
c       2. Add id to argument and use to set values of u    07/04/2009
c-----[--.----+----.----+----.-----------------------------------------]
c     Purpose: Set displacements for spin condition

c     Input:
c        rnod      - Number of nodes
c        prp       - Spin proportional load
c        prv       - Tranlsational proportional load
c        theta     - Spin velocity
c        nn(3)     - Normal for spin axis
c        xc(3)     - Center of spin axis
c        v0(3)     - Translational velocity
c        x(ndm,*)  - Reference coordinates
c        ndlist(*) - List of nodes to set

c     Output:
c        u(ndf,*)  - Displacement of node
c-----[--.----+----.----+----.-----------------------------------------]
      implicit   none

      include   'cdata.h'
      include   'pload1.h'
      include   'sdata.h'

      include   'pointer.h'
      include   'comblk.h'

      integer    rnod, j,n, nd, ndlist(*), id(ndf,numnp)
      real*8     prp,prv,theta, zz, cs,sn, tol, xsiz
      real*8     nn(3),xc(3),v0(3), x(ndm,numnp), u(ndf,numnp)
      real*8     t(3,3), dx(3)

      data       tol  / 1.d-8 /

c     Determine sizing measure

      xsiz = 0.0d0
      do n = 1,rnod
        nd = ndlist(n)
        do j = 1,3
          xsiz  = max(xsiz,abs(x(j,nd) - xc(j)))
        end do ! j
      end do ! n
      xsiz = tol*xsiz

c     Rotate face

      cs = cos(prp*theta) - 1.0d0
      sn = sin(prp*theta)
      do n = 1,rnod

c       Set global node number

        nd = ndlist(n)

c       Compute local frame direction cosines

        do j = 1,3
          t(j,3) = nn(j)
          dx(j)  = x(j,nd) - xc(j)
        end do ! j
        zz = nn(1)*dx(1) + nn(2)*dx(2) + nn(3)*dx(3)
        do j = 1,3
          dx(j) = dx(j) - zz*nn(j)
        end do
        zz = sqrt(dx(1)**2 + dx(2)**2 + dx(3)**2)
        if(zz.gt.xsiz) then
          do j = 1,3
            t(j,1) = dx(j)/zz
          end do ! j
          t(1,2) = t(2,3)*t(3,1) - t(3,3)*t(2,1)
          t(2,2) = t(3,3)*t(1,1) - t(1,3)*t(3,1)
          t(3,2) = t(1,3)*t(2,1) - t(2,3)*t(1,1)

c         Compute rotation and translation of node

          do j = 1,3
            if(id(j,nd).ne.0) then
              u(j,nd) = u(j,nd) + cs*dx(j) + sn*zz*t(j,2) + prv*v0(j)
            endif
          end do ! j

c       Translation only

        else
          do j = 1,3
            if(id(j,nd).ne.0) then
              u(j,nd) = u(j,nd) + prv*v0(j)
            endif
          end do ! j
        endif

      end do ! n

      end
