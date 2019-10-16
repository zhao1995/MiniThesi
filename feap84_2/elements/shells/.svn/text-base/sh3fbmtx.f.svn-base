c$Id:$
        subroutine sh3fbmtx ( shp1 , shp2 , shx1 , shx2 ,
     &                        cphm , cpx1 , cpx2 , cdrm ,
     &                        cdx1 , cdx2 , b           )

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]

c      Description:    SH3FBMTX is the subroutine which computes the
c                      discrete strain-displacement operator (matrix)
c                      for the general shell inextensible element.
c                      The membrane and bending strains are assumed
c                      to be defined in a cartesian reference frame.
c                      The shear strains are computed in the natural
c                      frame using Bathe-Dvorkin interpolations.

c      Authors:        M.S. Rifai, J.C. Simo, & D.D. Fox.
c      Date:           January, 1991.
c      Revised:        R.L. Taylor  - - February 1997
c-----[--.----+----.----+----.-----------------------------------------]
c      Routine Input:
c      --------------
c      shp1,shp2 ..... Nodal shape function natural derivatives.
c      shx1,shx2 ..... Nodal shape function cartesian derivatives.
c      cpx1,cpx2 ..... Current coordinate global derivatives
c                      at the Gauss points.
c      cphm .......... Current coordinate local derivatives
c                      at the midside nodes.
c      cdx1,cdx2 ..... Current director global derivatives
c                      at the Gauss points.
c      cdrm .......... Current local directors
c                      at the midside nodes.
c      xln ........... Localized nodal orthogonal transformation
c                      matrices, xln(3,3,ElementNode).

c      Routine Output:
c      ---------------
c      b ............. Discrete strain-displacement operator.
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      integer   i            , j        ,  nm      , nb
      integer   nsb1 (4)     , nsb2 (4) , nsm1 (4) , nsm2 (4)
      real*8    shp1 (4)     , shp2 (4) , shx1 (4) , shx2 (4)
      real*8    cphm (3,4)   , cpx1 (3) , cpx2 (3)
      real*8    cdrm (3,4)   , cdx1 (3) , cdx2 (3) , b(8,24)

      save

c     Shear Scatter Data:

      data      nsb1 / 2 , 2 , 3 , 3 /
      data      nsb2 / 4 , 3 , 3 , 4 /
      data      nsm1 / 2 , 2 , 4 , 4 /
      data      nsm2 / 1 , 3 , 3 , 1 /

c     [Bmm] Part (Membrane), [Bbb] Part (Bending), and
c     [Bbm] Part (Bending-Membrane Coupling Terms)

c     [Bsm] Part (Shear - Displacement)
c     [Bsb] Part (Shear - Rotation)

      nm = 0
      nb = 3
      do j = 1,4
        do i = 1,3
          b(1,nm+i) = shx1(j)*cpx1(i)
          b(2,nm+i) = shx2(j)*cpx2(i)
          b(3,nm+i) = shx2(j)*cpx1(i) + shx1(j)*cpx2(i)

          b(6,nm+i) = shx1(j)*cdx1(i)
          b(7,nm+i) = shx2(j)*cdx2(i)
          b(8,nm+i) = shx1(j)*cdx2(i) + shx2(j)*cdx1(i)

          b(6,nb+i) = b(1,nm+i)
          b(7,nb+i) = b(2,nm+i)
          b(8,nb+i) = b(3,nm+i)

          b(4,nm+i) = shp1(j) * cdrm(i,nsm1(j))
          b(5,nm+i) = shp2(j) * cdrm(i,nsm2(j))

          b(4,nb+i) = shp1(nsb1(j)) * cphm(i,nsm1(j))
          b(5,nb+i) = shp2(nsb2(j)) * cphm(i,nsm2(j))

        end do ! i
        nm = nm + 6
        nb = nb + 6
      end do ! j

      end
