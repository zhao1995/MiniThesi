c$Id:$
      subroutine sh3fgeo6 ( lint , shp1 , shp2 , shx1 , shx2 ,
     &                      cphm , cpx1 , cpx2 , dir  ,
     &                      sq   , sm   , s    )

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c        Description:    SH3FGEO6 adds additional terms
c                        static geometric tangent stiffness,
c                        in case of a 6 dof formulation.

c        Author:         N. Tarnow and J.C. Simo

c        Date:           February 1992
c-----[--.----+----.----+----.-----------------------------------------]
c        Routine Input:
c        --------------
c        shp1,shp2 ..... Nodal shape function natural derivatives.
c        shx1,shx2 ..... Nodal shape function cartesian derivatives.
c        cpx1,cpx2 ..... Current coordinate global derivatives
c                        at Gauss points.
c        cphm .......... Current coordinate local derivatives
c                        at midside nodes.
c        cdrm .......... Current local directors
c                        at midside nodes.
c        sn ............ Membrane stress.
c        sq ............ Shear stress.
c        sm ............ Bending stress.

c        Routine Output:
c        ---------------
c        s ............. Element stiffness.
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'crotas.h'

      integer   ni,k,l,i,lint
      real*8    shp1(4,4),shp2(4,4),shx1(4,4),shx2(4,4)
      real*8    cphm(3,4),cpx1(3,4),cpx2(3,4),dir(3,9)
      real*8    sq(2,4),sm(3,4),s(24,24),fac(4),cfac(4,3)
      real*8    xm(3,4),xn1(4),xn2(4),xn3(4),xn4(4), fact

      save

c     Evaluate Shape Functions:

      data      xn1  / 0.5d0, 0.0d0, 0.0d0, 0.5d0 /
      data      xn2  / 0.5d0, 0.5d0, 0.0d0, 0.0d0 /
      data      xn3  / 0.0d0, 0.5d0, 0.5d0, 0.0d0 /
      data      xn4  / 0.0d0, 0.0d0, 0.5d0, 0.5d0 /

c     Loop over Gauss Points

      do k = 1,4
        fac(k) = 0.0d0
      end do ! k
      do l = 1,lint
        fac(1) = fac(1) + 2.d0 * shp1(2,l) * sq(1,l)
        fac(2) = fac(2) + 2.d0 * shp1(3,l) * sq(1,l)
        fac(3) = fac(3) + 2.d0 * shp2(4,l) * sq(2,l)
        fac(4) = fac(4) + 2.d0 * shp2(3,l) * sq(2,l)
      end do ! l

      do k = 1,3
        cfac(1,k) = fac(1)*cphm(k,2)
        cfac(2,k) = fac(2)*cphm(k,4)
        cfac(3,k) = fac(3)*cphm(k,1)
        cfac(4,k) = fac(4)*cphm(k,3)
      end do ! k

c     Loop over Nodes

      ni = 0
      do i = 1,4
        if(ndof(i).eq.6)then
          do k = 1,3
            xm(k,i) = cfac(1,k)*xn2(i) + cfac(2,k)*xn4(i)
     &              + cfac(3,k)*xn1(i) + cfac(4,k)*xn3(i)
          end do ! k

c         Loop over Gauss Points

          do l = 1,lint
            do k = 1,3
              xm(k,i) = xm(k,i) + sm(1,l)*cpx1(k,l)*shx1(i,l)
     &                          + sm(2,l)*cpx2(k,l)*shx2(i,l)
     &                          + sm(3,l)*cpx1(k,l)*shx2(i,l)
     &                          + sm(3,l)*cpx2(k,l)*shx1(i,l)
            end do ! k
          end do ! l

c         Compute Projection

          fact = xm(1,i)*dir(1,i) + xm(2,i)*dir(2,i) + xm(3,i)*dir(3,i)
          do k = 1,3
            xm(k,i) = xm(k,i) - fact*dir(k,i)
          end do ! k

c         Add Contributions to Element Stiffness

          do k = 1,3
            do l = 1,3
              s(ni+3+k,ni+3+l) = s(ni+3+k,ni+3+l)
     &                         + 0.5d0*(xm(k,i)*dir(l,i)
     &                                + xm(l,i)*dir(k,i))
            end do ! l
          end do ! k
        endif
        ni = ni + 6
      end do ! i

      end
