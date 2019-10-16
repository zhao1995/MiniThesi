c$Id:$
        subroutine sh3fgeom ( shp1 , shp2 , shx1 , shx2 ,
     &                        cphm , cpx1 , cpx2 , dir  ,
     &                        sn   , sq   , sm   , s           )

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c        Description:    SH3FGEOM is the subroutine which constructs
c                        the static geometric tangent stiffness.

c        Authors:        M.S. Rifai, J.C. Simo, & D.D. Fox.

c        Date:           January 1991.
c-----[--.----+----.----+----.-----------------------------------------]
c        Routine Input:
c        --------------
c        shp1,shp2 ..... Nodal shape function natural derivatives.
c        shx1,shx2 ..... Nodal shape function cartesian derivatives.
c        cpx1,cpx2 ..... Current coordinate global derivatives
c                        at the Gauss points.
c        cphm .......... Current coordinate local derivatives
c                        at the midside nodes.
c                        at the midside nodes.
c        sn ............ Membrane stress.
c        sq ............ Shear stress.
c        sm ............ Bending stress.

c        Routine Output:
c        ---------------
c        s ............. Element stiffness.
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      integer   i  , j  , k  , ni , nj
      integer   ns1  (4), ns2  (4),  ns3  (4), ns4  (4)
      real*8    termb   , termm   , termd   , sq1h    , sq2h
      real*8    sm1i    , sm2i    , sn1i    , sn2i
      real*8    shp1 (4), shp2 (4), shx1 (4), shx2 (4), dir(3,9)
      real*8    cphm (3,4)   , cpx1 (3)     , cpx2 (3)
      real*8    sn   (3)     , sq   (2)     , sm   (3)
      real*8    s    (24,24) , qt   (4,4)   , qd   (4)

c     Scatter Data:

      data      ns1  / 2 , 2 , 4 , 4 /
      data      ns2  / 1 , 3 , 3 , 1 /
      data      ns3  / 1 , 1 , 4 , 4 /
      data      ns4  / 1 , 2 , 2 , 1 /

c     Set Up Shear Terms:

      sq1h    = 0.5d0 * sq(1)
      sq2h    = 0.5d0 * sq(2)

      qt(1,1) = shp1(1)*sq1h + shp2(1)*sq2h
      qt(1,2) = shp1(1)*sq1h
      qt(1,3) = 0.0d0
      qt(1,4) = shp2(1)*sq2h

      qt(2,1) = shp1(2)*sq1h
      qt(2,2) = shp1(2)*sq1h + shp2(2)*sq2h
      qt(2,3) = shp2(2)*sq2h
      qt(2,4) = 0.0d0

      qt(3,1) = 0.0d0
      qt(3,2) = shp2(3)*sq2h
      qt(3,3) = shp1(3)*sq1h + shp2(3)*sq2h
      qt(3,4) = shp1(3)*sq1h

      qt(4,1) = shp2(4)*sq2h
      qt(4,2) = 0.0d0
      qt(4,3) = shp1(4)*sq1h
      qt(4,4) = shp1(4)*sq1h + shp2(4)*sq2h

      do i = 1,4
        qd(i)   = sq(1)*shp1(ns3(i)) * (cphm(1,ns1(i))*dir(1,i)
     &                               +  cphm(2,ns1(i))*dir(2,i)
     &                               +  cphm(3,ns1(i))*dir(3,i))
     &          + sq(2)*shp2(ns4(i)) * (cphm(1,ns2(i))*dir(1,i)
     &                               +  cphm(2,ns2(i))*dir(2,i)
     &                               +  cphm(3,ns2(i))*dir(3,i))
      end do ! i

      ni = 0
      do i = 1,4
        sm1i= sm(1)*shx1(i) + sm(3)*shx2(i)
        sm2i= sm(2)*shx2(i) + sm(3)*shx1(i)

        sn1i= sn(1)*shx1(i) + sn(3)*shx2(i)
        sn2i= sn(2)*shx2(i) + sn(3)*shx1(i)

        nj = 0
        do j = 1,4

c         Bending Term

          termb = sm1i*shx1(j) + sm2i*shx2(j)

c         Membrane Term

          termm = sn1i*shx1(j) + sn2i*shx2(j)

c         Assemble Membrane, Bending and Shear Parts

          do k = 1,3
            s(ni+k,nj+k  ) = s(ni+k,nj+k  ) + termm
            s(ni+k,nj+3+k) = s(ni+k,nj+3+k) + termb + qt(i,j)
            s(ni+3+k,nj+k) = s(ni+3+k,nj+k) + termb + qt(j,i)
          end do ! k
          nj = nj + 6
        end do ! j

c       Diagonal Bending Term

        termd = (sm1i*cpx1(1) + sm2i*cpx2(1))*dir(1,i)
     &        + (sm1i*cpx1(2) + sm2i*cpx2(2))*dir(2,i)
     &        + (sm1i*cpx1(3) + sm2i*cpx2(3))*dir(3,i)

c       Assemble Diagonal Bending and Shear Parts

        do k = 1,3
          s(ni+3+k,ni+3+k) = s(ni+3+k,ni+3+k) + qd(i) - termd
        end do ! k

        ni = ni + 6
      end do ! l

      end
