c$Id:$
      subroutine shpm3d(shpi,xl,ndm,finc)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Use informaiton from qudshp.h                    11/11/2008
c       2. Remove xsj array                                 10/12/2008
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose:  Modified shape function routine for 3D Q1-quad

c      Inputs:

c      Outputs:
c-----[--.----+----.----+----.-----------------------------------------]
      implicit   none

      include   'qudshp.h'

      logical    finc
      integer    ndm, i, j, k, l
      real*8     cc,c1m,c1p,c2m,c2p,c3m,c3p
      real*8     sgl1(8), sgl2(8), sgl3(8), shpi(3,4,9)
      real*8     xl(ndm,*), hgc(3,4),xji(3,3)
      real*8     gamma(4,8),fact(3,4), h1(8),h2(8),h3(8),h4(8)

      save

      data     h1/ 0.125d0, 0.125d0,-0.125d0,-0.125d0,
     &            -0.125d0,-0.125d0, 0.125d0, 0.125d0/
      data     h2/ 0.125d0,-0.125d0,-0.125d0, 0.125d0,
     &            -0.125d0, 0.125d0, 0.125d0,-0.125d0/
      data     h3/ 0.125d0,-0.125d0, 0.125d0,-0.125d0,
     &             0.125d0,-0.125d0, 0.125d0,-0.125d0/
      data     h4/-0.125d0, 0.125d0,-0.125d0, 0.125d0,
     &             0.125d0,-0.125d0, 0.125d0,-0.125d0/

c     COMPUTE ONE-POINT SHAPE FUNCTIONS

      call shp1pt(jac(9),shp3(1,1,9),xji,xl,ndm)

c     COMPUTE STABILIZATION COEFFICIENTS

c     Compute dot products

      do i = 1,3

        fact(i,1) = 0.125d0*(  xl(i,1) + xl(i,2) - xl(i,3) - xl(i,4)
     &                       - xl(i,5) - xl(i,6) + xl(i,7) + xl(i,8))
        fact(i,2) = 0.125d0*(  xl(i,1) - xl(i,2) - xl(i,3) + xl(i,4)
     &                       - xl(i,5) + xl(i,6) + xl(i,7) - xl(i,8))
        fact(i,3) = 0.125d0*(  xl(i,1) - xl(i,2) + xl(i,3) - xl(i,4)
     &                       + xl(i,5) - xl(i,6) + xl(i,7) - xl(i,8))
        fact(i,4) = 0.125d0*(- xl(i,1) + xl(i,2) - xl(i,3) + xl(i,4)
     &                       + xl(i,5) - xl(i,6) + xl(i,7) - xl(i,8))

      end do ! i

c     Compute gamma factors

      do i = 1,8

        gamma(1,i) = h1(i) - fact(1,1)*shp3(1,i,9)
     &                     - fact(2,1)*shp3(2,i,9)
     &                     - fact(3,1)*shp3(3,i,9)

        gamma(2,i) = h2(i) - fact(1,2)*shp3(1,i,9)
     &                     - fact(2,2)*shp3(2,i,9)
     &                     - fact(3,2)*shp3(3,i,9)

        gamma(3,i) = h3(i) - fact(1,3)*shp3(1,i,9)
     &                     - fact(2,3)*shp3(2,i,9)
     &                     - fact(3,3)*shp3(3,i,9)

        gamma(4,i) = h4(i) - fact(1,4)*shp3(1,i,9)
     &                     - fact(2,4)*shp3(2,i,9)
     &                     - fact(3,4)*shp3(3,i,9)

      end do ! i

c     COMPUTE MODIFIED SHAPE FUNTIONS

c     Natural derivatives of hourglass modes

      do l = 1,8

        call jacb3m(sg3(1,l),jac(l),xl,ndm)

        sgl1(l) = sg3(1,l)/jac(l)
        sgl2(l) = sg3(2,l)/jac(l)
        sgl3(l) = sg3(3,l)/jac(l)

c       Modified cartesian derivatives of hourglass modes

        do i = 1,3

          hgc(i,1) = xji(i,2)*sgl3(l)  + xji(i,3)*sgl2(l)
          hgc(i,2) = xji(i,3)*sgl1(l)  + xji(i,1)*sgl3(l)
          hgc(i,3) = xji(i,1)*sgl2(l)  + xji(i,2)*sgl1(l)
          hgc(i,4) = hgc(i,3)*sg3(3,l) + xji(i,3)*sgl1(l)*sg3(2,l)

        end do ! i

c       Modified derivatives of shape function

        do i = 1,3
          do j = 1,8
            shp3(i,j,l) = shp3(i,j,9)

            do k = 1,4
              shp3(i,j,l) = shp3(i,j,l) + gamma(k,j)*hgc(i,k)
            end do ! k

          end do ! j
          shpi(i,4,l) = hgc(i,4)
        end do ! i

c       Shape functions

        c1p = 0.125d0 + 0.125d0*sg3(1,l)
        c1m = 0.125d0 - 0.125d0*sg3(1,l)
        c2p = 1.d0 + sg3(2,l)
        c2m = 1.d0 - sg3(2,l)
        c3p = 1.d0 + sg3(3,l)
        c3m = 1.d0 - sg3(3,l)

        cc         = c1m*c2m
        shp3(4,1,l) = cc*c3m
        shp3(4,5,l) = cc*c3p

        cc         = c1p*c2m
        shp3(4,2,l) = cc*c3m
        shp3(4,6,l) = cc*c3p

        cc         = c1p*c2p
        shp3(4,3,l) = cc*c3m
        shp3(4,7,l) = cc*c3p

        cc         = c1m*c2p
        shp3(4,4,l) = cc*c3m
        shp3(4,8,l) = cc*c3p

      end do ! l

c     COMPUTE ENHANCED SHAPE FUNCTIONS

c     First three enhanced modes

      if (finc) then

        do l = 1,8

          do i = 1,3
            shpi(i,1,l) = xji(i,1)*sgl1(l)
            shpi(i,2,l) = xji(i,2)*sgl2(l)
            shpi(i,3,l) = xji(i,3)*sgl3(l)
          end do ! i

        end do ! l

        do j = 1,4
          do i = 1,3
            shpi(i,j,9) = 0.0d0
          end do ! i
        end do ! j

      endif

      end
