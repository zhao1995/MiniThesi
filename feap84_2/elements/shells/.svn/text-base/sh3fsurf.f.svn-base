c$Id:$
      subroutine sh3fsurf(d,dthk,xl,ul,xjw, ndf,ndm,nst, p,s)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c      1. Convert body loads to unit values by dividing by  31/05/2013
c         thickness.
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose:  Nodal force and tangent array for pressure loading

c      INPUT variables
c        d(10)      Value of constant pressure on face
c        dthk       Shell element thickness
c        xl(4,*)    Nodal coordinates
c        ul(4,*)    Nodal displacements
c        xjw(*)     Jacobian times weight
c        ndf        Number of DOF / node
c        ndm        Space dimension
c        nst        Dimension of residual vector

c      OUTPUT variables
c        p(nst)     Contribution to residual
c        s(nst,nst) Contribution to striffness matrix

c      PARAMATER set up:

c         ma   =  0  (for 3-d analysis in reference coordinates)
c              =  1  (for 3-d analysis in current coordinates;
c                     an unsymmetric tangent matrix also computed.)

c         nel  =  4  (for 4-node surface of bi-linear element)

c                        4                3
c                        o----------------o
c                        |       A        |
c                        |       |eta     |
c                        |       |        |
c                        |       +---> xi |
c                        |                |
c                        |                |
c                        |                |
c         Nodes are:     o----------------o
c                        1                2

c       Author:         J.C. Simo

c       Revised:        R.L. Taylor  - - February 1997
c-----[--.----+----.----+----.-----------------------------------------]

      implicit  none

      include  'ddata.h'
      include  'eldata.h'
      include  'eltran.h'

      logical   bflg
      integer   l, lint, ndf, ndm, nst
      integer   ii, jj, i1, j1, i, j
      real*8    dthk, pn, pp, fact

      real*8    sg(3,4), shp(3,4), xl(ndm,*), xu(3,4), xjw(*)
      real*8    d(*), ul(ndf,*), p(ndf,*), s(nst,*), dx(3,2)
      real*8    xii(4), eti(4), bg(3),bl(3)
      real*8    bf(3),bt(3,3),xx(3)

      save

      data      xii / -0.5d0, 0.5d0, 0.5d0,-0.5d0 /
      data      eti / -0.5d0,-0.5d0, 0.5d0, 0.5d0 /
      data      bl  /  3*0.0d0 /

c     Set body loading factors

      call sbodyf(d, bg)

c     Multiply by thickness to get total load on shell surface

      bg(:) = bg(:)*dthk

c     Set pressure loading

      if(int(d(76)).gt.0) then
        bl(3) = d(10)
      else
        bl(3) = d(10)*dm
      endif

c     Compute nodal coordinates in correct reference frame

      if(d(68).eq.0.0d0) then
        fact = 0.d0
      else
        fact = theta(3)
      endif
      do ii = 1,4
        xu(1,ii) = xl(1,ii) + fact*ul(1,ii)
        xu(2,ii) = xl(2,ii) + fact*ul(2,ii)
        xu(3,ii) = xl(3,ii) + fact*ul(3,ii)
      end do ! ii

c     Get quadrature information

      l  = 2
      call int2d (l, lint, sg)

      bflg = d(4).gt.0.0d0 .and. d(65).gt.0.0d0  ! angular velocity

c     First loop over quadrature points

      do l = 1,lint

c       Compute shape functions and geometric factors

        do ii = 1,4
          shp(1,ii) = xii(ii)*(0.5d0+eti(ii)*sg(2,l))
          shp(2,ii) = eti(ii)*(0.5d0+xii(ii)*sg(1,l))
          shp(3,ii) = (0.5d0+xii(ii)*sg(1,l))*(0.5d0 + eti(ii)*sg(2,l))
        end do ! ii

        do ii = 1,3
          dx(ii,1) = shp(1,1)*xu(ii,1) + shp(1,2)*xu(ii,2)
     &             + shp(1,3)*xu(ii,3) + shp(1,4)*xu(ii,4)
          dx(ii,2) = shp(2,1)*xu(ii,1) + shp(2,2)*xu(ii,2)
     &             + shp(2,3)*xu(ii,3) + shp(2,4)*xu(ii,4)
        end do ! ii

c       Angular velocity: d(4) = rho; d(65) = omega

        do i = 1,3
          bf(i) = 0.0d0
        end do ! i
        if(bflg) then
          do ii = 1,3
            xx(ii) = 0.0d0
            do jj = 1,4
              xx(ii) = xx(ii) + shp(3,jj)*xu(ii,jj)
            end do ! jj
          end do ! ii
          call sbodyw(d(4),d(65),xx, bf,bt, .true.)
          do ii = 1,3
            do jj = 1,3
              bt(jj,ii) = bt(jj,ii)*xjw(l)
            end do ! jj
          end do ! ii
        endif

c       Compute nodal loads for pressures & gravity loading

        pn = bl(3)*sg(3,l)
        do ii = 1,4
          pp      = shp(3,ii)*pn
          p(1,ii) = p(1,ii) + pp*(dx(2,1)*dx(3,2) - dx(3,1)*dx(2,2))
     &                      + (bg(1)+bf(1))*shp(3,ii)*xjw(l)
          p(2,ii) = p(2,ii) + pp*(dx(3,1)*dx(1,2) - dx(1,1)*dx(3,2))
     &                      + (bg(2)+bf(2))*shp(3,ii)*xjw(l)
          p(3,ii) = p(3,ii) + pp*(dx(1,1)*dx(2,2) - dx(2,1)*dx(1,2))
     &                      + (bg(3)+bf(3))*shp(3,ii)*xjw(l)
        end do ! ii

c       Compute follower load tangent if necessary

        if(d(68).gt.0.0d0) then
          i1 = 0
          do ii = 1,4
            pp = shp(3,ii)*pn*ctan(1)
            j1 = 0
            do jj = 1,4
              s(i1+1,j1+2) = s(i1+1,j1+2)
     &                     - pp*(shp(1,jj)*dx(3,2) - dx(3,1)*shp(2,jj))
              s(i1+2,j1+3) = s(i1+2,j1+3)
     &                     - pp*(shp(1,jj)*dx(1,2) - dx(1,1)*shp(2,jj))
              s(i1+3,j1+1) = s(i1+3,j1+1)
     &                     - pp*(shp(1,jj)*dx(2,2) - dx(2,1)*shp(2,jj))
              s(i1+1,j1+3) = s(i1+1,j1+3)
     &                     - pp*(shp(2,jj)*dx(2,1) - dx(2,2)*shp(1,jj))
              s(i1+2,j1+1) = s(i1+2,j1+1)
     &                     - pp*(shp(2,jj)*dx(3,1) - dx(3,2)*shp(1,jj))
              s(i1+3,j1+2) = s(i1+3,j1+2)
     &                     - pp*(shp(2,jj)*dx(1,1) - dx(1,2)*shp(1,jj))
              j1 = j1 + ndf
            end do ! jj
            i1 = i1 + ndf
          end do ! ii
        endif

c       Compute rotational body force tangent if necessary

        if(bflg) then
          i1 = 0
          do ii = 1,4
            pp = shp(3,ii)*ctan(1)
            j1 = 0
            do jj = 1,4
              do i = 1,3
                do j = 1,3
                  s(i1+1,j1+2) = s(i1+1,j1+2) + pp*bt(i,j)*shp(3,jj)
                end do ! j
              end do ! i
              j1 = j1 + ndf
            end do ! jj
            i1 = i1 + ndf
          end do ! ii
        endif

      end do ! l

      end
