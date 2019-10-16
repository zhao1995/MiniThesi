c$Id:$
        subroutine sh3fpian ( xl  , xjw  , sg  , eh  , eh12 ,
     &                        xnu , f1   , f2  , y1  , y2  , xnt  )

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c      1. Initialize g array for all 6 terms                21/06/2013
c-----[--.----+----.----+----.-----------------------------------------]
c        Description:    SH3FPIAN is a subroutine which computes some
c                        objects needed by Pian-Summihara-type
c                        Hellinger-Reissner mixed element.

c        Authors:        M.S. Rifai, J.C. Simo, & D.D. Fox.

c        Date:           February 1991.
c-----[--.----+----.----+----.-----------------------------------------]
c        Routine Input:
c        --------------
c        xl ............ Initial nodal coordinates
c        xjw ........... Reference jacobian weighted for integration.
c        sg ............ Gauss point coordinates.
c        eh ............ Membrane inverse constitutive coefficient.
c        eh12 .......... bending inverse constitutive coefficient.
c        xnu ........... Poisson's ratio.

c        Routine Output:
c        ---------------
c        f1,f2 ......... Tensor transformations in vector format
c                        at center of element.
c        y1,y2, xnt .... Commonly used parts of [H].
c-----[--.----+----.----+----.-----------------------------------------]

        implicit  none

        integer   i
        real*8    xlcnrm , eta   , xi     , etaxi  , xjweh   , xjwE2
        real*8    eh     , eh12  , xnu    , xnu12
        real*8    xl(3,4), xjw(4), zpc1(3), zpc2(3), xjc(2,2), xlc(3,3)
        real*8    sg(3,4), f1(3) , f2(3)  , y1(4)  , y2(4)   , xnt(13)
c       real*8    dot

c       Compute central Local Derivatives of Coordinates

        do i = 1,3
          zpc1(i) = 0.25d0 * (-xl(i,1)+xl(i,2)+xl(i,3)-xl(i,4))
          zpc2(i) = 0.25d0 * (-xl(i,1)-xl(i,2)+xl(i,3)+xl(i,4))
        end do ! i

c       Compute central Surface Normal

        call vecp ( zpc1(1) , zpc2(1) , xlc(1,3) )

c       xlcnrm   = 1.d0/sqrt ( dot(xlc(1,3),xlc(1,3),3) )
        xlcnrm   = 1.d0/sqrt (xlc(1,3)*xlc(1,3)
     &                      + xlc(2,3)*xlc(2,3)
     &                      + xlc(3,3)*xlc(3,3))
        xlc(1,3) = xlc(1,3) * xlcnrm
        xlc(2,3) = xlc(2,3) * xlcnrm
        xlc(3,3) = xlc(3,3) * xlcnrm

c       Compute Local-Global Jacobian

        call sh3flmda ( xlc(1,3) , xlc  )

c       xjc(1,1) = dot ( xlc(1,1) , zpc1(1) ,3 )
c       xjc(2,1) = dot ( xlc(1,2) , zpc1(1) ,3 )
c       xjc(1,2) = dot ( xlc(1,1) , zpc2(1) ,3 )
c       xjc(2,2) = dot ( xlc(1,2) , zpc2(1) ,3 )

        xjc(1,1) = xlc(1,1) * zpc1(1)
     &           + xlc(2,1) * zpc1(2)
     &           + xlc(3,1) * zpc1(3)

        xjc(2,1) = xlc(1,2) * zpc1(1)
     &           + xlc(2,2) * zpc1(2)
     &           + xlc(3,2) * zpc1(3)

        xjc(1,2) = xlc(1,1) * zpc2(1)
     &           + xlc(2,1) * zpc2(2)
     &           + xlc(3,1) * zpc2(3)

        xjc(2,2) = xlc(1,2) * zpc2(1)
     &           + xlc(2,2) * zpc2(2)
     &           + xlc(3,2) * zpc2(3)

c       Compute central Transformation Tensor in Vector Form

        f1(1) = xjc(1,1) * xjc(1,1)
        f1(2) = xjc(2,1) * xjc(2,1)
        f1(3) = xjc(2,1) * xjc(1,1)
        f2(1) = xjc(1,2) * xjc(1,2)
        f2(2) = xjc(2,2) * xjc(2,2)
        f2(3) = xjc(2,2) * xjc(1,2)

c       Compute Element Integrals of 1,xi,eta and Multiples thereof

        do i = 1 , 13
          xnt(i) = 0.0d0
        end do ! i
        do i = 1 , 4
          xnt(1) = xnt(1) + xjw(i)
          xnt(2) = xnt(2) + xjw(i) * sg(2,i)
          xnt(3) = xnt(3) + xjw(i) * sg(1,i)
          xnt(4) = xnt(4) + xjw(i) * sg(2,i) * sg(1,i)
        end do ! i

        xnt(2) = xnt(2) / xnt(1)
        xnt(3) = xnt(3) / xnt(1)
        xnt(4) = xnt(4) / xnt(1)

        do i = 1 , 4
          eta     = sg(2,i)         - xnt(2)
          xi      = sg(1,i)         - xnt(3)
          etaxi   = sg(2,i)*sg(1,i) - xnt(4)
          xjweh   = xjw(i) / eh
          xjwE2   = xjw(i) / eh12
          xnt(5 ) = xnt(5 ) + xjweh * eta ** 2
          xnt(6 ) = xnt(6 ) + xjweh * xi  ** 2
          xnt(7 ) = xnt(7 ) + xjweh * eta *  xi
          xnt(8 ) = xnt(8 ) + xjweh * xi  *  etaxi
          xnt(9 ) = xnt(9 ) + xjweh * eta *  etaxi
          xnt(10) = xnt(10) + xjweh * etaxi ** 2
          xnt(11) = xnt(11) + xjwE2 * eta ** 2
          xnt(12) = xnt(12) + xjwE2 * xi  ** 2
          xnt(13) = xnt(13) + xjwE2 * eta *  xi
        end do ! i

c       Compute Common Parts of [H]

        xnu12 = 2.d0  * (1.d0  + xnu)
        y1(1) = f1(1) - xnu * f1(2)
        y1(2) = f1(2) - xnu * f1(1)
        y1(3) = f1(3) * xnu12
        y2(1) = f2(1) - xnu * f2(2)
        y2(2) = f2(2) - xnu * f2(1)
        y2(3) = f2(3) * xnu12
        y1(4) = (f1(1)+f1(2)) * xnu
        y2(4) = (f2(1)+f2(2)) * xnu

        end

        subroutine sh3fmbrn ( xjw , sg  , eh  , xnu ,
     &                        f1  , f2  , y1  , y2  , xnt ,
     &                        ndf , nst , ce  , b   , sn  ,
     &                        p   , s   )
c-----[--.----+----.----+----.-----------------------------------------]
c        Description:    SH3FMBRN is a subroutine which computes the
c                        stiffness and residual contributions of the
c                        membrane field, using Pian-Sumihara
c                        stress interpolations for Shell.

c        Authors:        M.S. Rifai, J.C. Simo, & D.D. Fox.

c        Date:           February 1991.
c-----[--.----+----.----+----.-----------------------------------------]
c        Routine Input:
c        --------------
c        xjw ........... Reference jacobian weighted for integration.
c        sg ............ Gauss point coordinates.
c        eh ............ Membrane inverse constitutive coefficient.
c        xnu ........... Poisson's ratio.
c        f1,f2 ......... Tensor transformations in vector format
c                        at center of element.
c        y1,y2, xnt .... Commonly used parts of [H].
c        ndf ........... Number of degrees of freedom per node.
c        nst ........... Number of degrees of freedom per element.
c        ce ............ Current membrane strain measure.
c        b ............. Discrete strain-displacement operator.

c        Routine Output:
c        ---------------
c        sn ............ Membrane stress.
c        s ............. Element stiffness.
c        p ............. Element residual.
c-----[--.----+----.----+----.-----------------------------------------]

        implicit  none

        integer   i      , j      , k        , l       , ndf
        integer   k1     , k2     ,  nst     , n
        integer   ni     , nj
        real*8    eta    , xi     , fac0     , fac1    , fac2
        real*8    fach1  , fach2  , eh       , xnu     , hdet
        real*8    xjw(4) , sg(3,4), f1(3)    , f2(3)
        real*8    xnt(13), ce(3,4), sn(3,4)  , htm(2,2), hin(5,5)
        real*8    ee(5)  , be(5)  , b(8,24,4), g(4,5,3), gh(3,5)
        real*8    y1(4)  , y2(4)  , p(nst)   , s(nst,nst)
c       real*8    dot

        save

c       Zero Matrices

        do i = 1,5
          ee(i) = 0.0d0
          do j = 1,5
            hin(j,i) = 0.0d0
          end do ! j
          do j = 1,4
            g(j,i,1) = 0.0d0
            g(j,i,2) = 0.0d0
            g(j,i,3) = 0.0d0
          end do ! j
        end do ! i

c       Gauss Point Loop

        do i = 1 , 4

c         Common Factors

          eta   = sg(2,i) - xnt(2)
          xi    = sg(1,i) - xnt(3)
          fac0  = xjw(i)
          fac1  = fac0 * eta
          fac2  = fac0 * xi

c         Modified Strains [ee] = [S][ce]

          ee(1) = ee(1) + fac0 * ce(1,i)
          ee(2) = ee(2) + fac0 * ce(2,i)
          ee(3) = ee(3) + fac0 * ce(3,i)

          ee(4) = ee(4) + fac1 * (f1(1) * ce(1,i)
     &                         +  f1(2) * ce(2,i)
     &                         +  f1(3) * ce(3,i))
          ee(5) = ee(5) + fac2 * (f2(1) * ce(1,i)
     &                         +  f2(2) * ce(2,i)
     &                         +  f2(3) * ce(3,i))

c         ee(4) = ee(4) + fac1 * dot ( f1 , ce(1,i) ,3 )
c         ee(5) = ee(5) + fac2 * dot ( f2 , ce(1,i) ,3 )

c         Integrated Strain-Displacement Matrix [G] = [S][B]

          do k = 1 , 4
            do j = 1 , 3
              l = (k-1)*ndf + j
              g(k,1,j) = g(k,1,j) + fac0 * b(1,l,i)
              g(k,2,j) = g(k,2,j) + fac0 * b(2,l,i)
              g(k,3,j) = g(k,3,j) + fac0 * b(3,l,i)

              g(k,4,j) = g(k,4,j) + fac1 * (f1(1)*b(1,l,i)
     &                                   +  f1(2)*b(2,l,i)
     &                                   +  f1(3)*b(3,l,i))

              g(k,5,j) = g(k,5,j) + fac2 * (f2(1)*b(1,l,i)
     &                                   +  f2(2)*b(2,l,i)
     &                                   +  f2(3)*b(3,l,i))

c             g(k,4,j) = g(k,4,j) + fac1 * dot(f1,b(1,l,i),3 )
c             g(k,5,j) = g(k,5,j) + fac2 * dot(f2,b(1,l,i),3 )
            end do ! j
          end do ! k

        end do ! i

c       Compute Second Block of [H-inv]

c       htm(1,1) = xnt(6) * dot ( f1 , y1 ,3 )
c       htm(2,2) = xnt(5) * dot ( f2 , y2 ,3 )
c       htm(1,2) = xnt(7) * dot ( f1 , y2 ,3 )

        htm(1,1) = xnt(6) * (f1(1)*y1(1) + f1(2)*y1(2) + f1(3)*y1(3))
        htm(2,2) = xnt(5) * (f2(1)*y2(1) + f2(2)*y2(2) + f2(3)*y2(3))
        htm(1,2) = xnt(7) * (f1(1)*y2(1) + f1(2)*y2(2) + f1(3)*y2(3))
        hdet     = 1.d0 / (htm(1,1)*htm(2,2) - htm(1,2)**2)

c       Constitutive Relations

c       Compute [H-inv]

        fach1    = eh / (xnt(1) * (1.d0 - xnu * xnu))
        fach2    = eh / (xnt(1) * (1.d0 + xnu) * 2.d0)
        hin(1,1) = fach1
        hin(1,2) = fach1 * xnu
        hin(2,1) = hin(1,2)
        hin(2,2) = hin(1,1)
        hin(3,3) = fach2
        hin(4,4) = htm(2,2) * hdet
        hin(5,5) = htm(1,1) * hdet
        hin(4,5) =-htm(1,2) * hdet
        hin(5,4) = hin(4,5)

c       Compute Stress Variables

        do i = 1,5
          be(i) = hin(i,1)*ee(1)
     1          + hin(i,2)*ee(2)
     2          + hin(i,3)*ee(3)
     3          + hin(i,4)*ee(4)
     4          + hin(i,5)*ee(5)
        end do ! i

c       Gauss Point Loop

        do i = 1 , 4

c         Compute Stresses

          eta     = sg(2,i) - xnt(2)
          xi      = sg(1,i) - xnt(3)
          sn(1,i) = be(1) + eta * f1(1)*be(4) + xi * f2(1)*be(5)
          sn(2,i) = be(2) + eta * f1(2)*be(4) + xi * f2(2)*be(5)
          sn(3,i) = be(3) + eta * f1(3)*be(4) + xi * f2(3)*be(5)

c         Weight Stresses by Jacobian

          sn(1,i) = sn(1,i) * xjw(i)
          sn(2,i) = sn(2,i) * xjw(i)
          sn(3,i) = sn(3,i) * xjw(i)

        end do ! i

c       Node Loop #1

        do i = 1 , 4
          ni = (i-1)*ndf

c         Compute Residual

          do k = 1 , 3
            p(ni+k) = p(ni+k) - g(i,1,k)*be(1)
     &                        - g(i,2,k)*be(2)
     &                        - g(i,3,k)*be(3)
     &                        - g(i,4,k)*be(4)
     &                        - g(i,5,k)*be(5)
          end do ! i

c         Node Loop #2

          do j = 1 , 4
            nj = (j-1)*ndf

c           Compute Stiffness (m-m) Part

            do n = 1,5
              do k1 = 1,3
                gh(k1,n) = g(i,1,k1)*hin(1,n)
     &                   + g(i,2,k1)*hin(2,n)
     &                   + g(i,3,k1)*hin(3,n)
     &                   + g(i,4,k1)*hin(4,n)
     &                   + g(i,5,k1)*hin(5,n)
              end do ! k1
            end do ! n

            do k1 = 1 , 3
              do k2 = 1 , 3
                s(ni+k1,nj+k2) = s(ni+k1,nj+k2)
     &                         + gh(k1,1) * g(j,1,k2)
     &                         + gh(k1,2) * g(j,2,k2)
     &                         + gh(k1,3) * g(j,3,k2)
     &                         + gh(k1,4) * g(j,4,k2)
     &                         + gh(k1,5) * g(j,5,k2)
              end do ! k2
            end do ! k1

          end do ! j

        end do ! i

        end

        subroutine sh3fbend ( xjw  , sg  , eh12, xnu ,
     &                        f1   , f2  , y1  , y2  , xnt ,
     &                        ndf  , nst , cr  , b   , sm  ,
     &                        p    , s   )
c*********************************************************************
c        Description:   SH3FBEND is a subroutine which computes the
c                       stiffness and residual contributions of the
c                       bending field, using Pian-Sumihara
c                       stress interpolations for Shell.

c        Authors:       M.S. Rifai, J.C. Simo, & D.D. Fox.

c        Date:          February 1991.
c-----[--.----+----.----+----.-----------------------------------------]

c        Routine Input:
c        --------------
c        xjw ........... Reference jacobian weighted for integration.
c        sg ............ Gauss point coordinates.
c        eh12 .......... Membrane inverse constitutive coefficient.
c        xnu ........... Poisson's ratio.
c        f1,f2 ......... Tensor transformations in vector format
c                        at center of element.
c        y1,y2, xnt .... Commonly used parts of [H].
c        ndf ........... Number of degrees of freedom per node.
c        nst ........... Number of degrees of freedom per element.
c        cr ............ Current bending strain measure.
c        B ............. Discrete strain-displacement operator.

c        Routine Output:
c        ---------------
c        sm ............ Bending stress.
c        S ............. Element stiffness.
c        P ............. Element residual.
c-----[--.----+----.----+----.-----------------------------------------]
        implicit  none

        integer   i      , j       , k       , l        , ndf
        integer   k1     , k2      , nst     , n
        integer   ni     , nj
        real*8    eta    , xi      , fac0    , fac1     , fac2
        real*8    fach1    , fach2 , hdet    , eh12     , xnu
        real*8    xjw(4) , sg(3,4) , f1(3)   , f2(3)
        real*8    xnt(13), cr(3,4) , sm(3,4) , htm(2,2) , hin(5,5)
        real*8    br(5)  , y1(4)   , y2(4)   , b(8,24,4), g(4,5,6)
        real*8    er(5)  , p(nst)  , s(nst,nst), gh(6,5)
c       real*8    dot

        save

c       Zero Matrices

        do i = 1,5
          er(i) = 0.0d0
          do j = 1,5
            hin(j,i) = 0.0d0
          end do ! j
          do j = 1,4
            g(j,i,1) = 0.0d0
            g(j,i,2) = 0.0d0
            g(j,i,3) = 0.0d0
            g(j,i,4) = 0.0d0
            g(j,i,5) = 0.0d0
            g(j,i,6) = 0.0d0
          end do ! j
        end do ! i

c       Gauss Point Loop

        do i = 1 , 4

c         Common Factors

          eta   = sg(2,i) - xnt(2)
          xi    = sg(1,i) - xnt(3)
          fac0  = xjw(i)
          fac1  = fac0 * eta
          fac2  = fac0 * xi

c         Modified Strains [er] = [S][cr]

          er(1) = er(1) + fac0 * cr(1,i)
          er(2) = er(2) + fac0 * cr(2,i)
          er(3) = er(3) + fac0 * cr(3,i)
          er(4) = er(4) + fac1 * (f1(1)*cr(1,i)
     &                         +  f1(2)*cr(2,i)
     &                         +  f1(3)*cr(3,i))
          er(5) = er(5) + fac2 * (f2(1)*cr(1,i)
     &                         +  f2(2)*cr(2,i)
     &                         +  f2(3)*cr(3,i))

c         er(4) = er(4) + fac1 * dot ( f1 , cr(1,i) ,3 )
c         er(5) = er(5) + fac2 * dot ( f2 , cr(1,i) ,3 )

c         Integrated Strain-Displacement Matrix [G] = [S][B]

          do k = 1 , 4
            do j = 1 , 6
              l = (k-1)*ndf + j
              g(k,1,j) = g(k,1,j) + fac0 * b(6,l,i)
              g(k,2,j) = g(k,2,j) + fac0 * b(7,l,i)
              g(k,3,j) = g(k,3,j) + fac0 * b(8,l,i)

              g(k,4,j) = g(k,4,j) + fac1 * (f1(1)*b(6,l,i)
     &                                   +  f1(2)*b(7,l,i)
     &                                   +  f1(3)*b(8,l,i))
              g(k,5,j) = g(k,5,j) + fac2 * (f2(1)*b(6,l,i)
     &                                   +  f2(2)*b(7,l,i)
     &                                   +  f2(3)*b(8,l,i))

c             g(k,4,j) = g(k,4,j) + fac1 * dot(f1,b(6,l,i),3 )
c             g(k,5,j) = g(k,5,j) + fac2 * dot(f2,b(6,l,i),3 )
            end do ! j
          end do ! k

        end do ! i

c       Compute Second Block of [H-inv]

c       htm(1,1) = xnt(12) * dot ( f1 , y1 ,3 )
c       htm(2,2) = xnt(11) * dot ( f2 , y2 ,3 )
c       htm(1,2) = xnt(13) * dot ( f1 , y2 ,3 )

        htm(1,1) = xnt(12) * (f1(1)*y1(1) + f1(2)*y1(2) + f1(3)*y1(3))
        htm(2,2) = xnt(11) * (f2(1)*y2(1) + f2(2)*y2(2) + f2(3)*y2(3))
        htm(1,2) = xnt(13) * (f1(1)*y2(1) + f1(2)*y2(2) + f1(3)*y2(3))
        hdet = 1.d0 / (htm(1,1)*htm(2,2) - htm(1,2)**2)

c       Constitutive Relations

c       Compute [H-inv]

        fach1    = eh12 / (xnt(1) * (1.d0 - xnu * xnu))
        fach2    = eh12 / (xnt(1) * (1.d0 + xnu) * 2.d0)
        hin(1,1) = fach1
        hin(1,2) = fach1 * xnu
        hin(2,1) = hin(1,2)
        hin(2,2) = hin(1,1)
        hin(3,3) = fach2
        hin(4,4) = htm(2,2) * hdet
        hin(5,5) = htm(1,1) * hdet
        hin(4,5) =-htm(1,2) * hdet
        hin(5,4) = hin(4,5)

c       Compute Stress Variables

        do i = 1,5
          br(i) = hin(i,1)*er(1)
     &          + hin(i,2)*er(2)
     &          + hin(i,3)*er(3)
     &          + hin(i,4)*er(4)
     &          + hin(i,5)*er(5)
        end do ! i

c       Gauss Point Loop

        do i = 1 , 4

c         Compute Stresses

          eta   = sg(2,i) - xnt(2)
          xi    = sg(1,i) - xnt(3)
          sm(1,i) = br(1) + eta * f1(1)*br(4) + xi * f2(1)*br(5)
          sm(2,i) = br(2) + eta * f1(2)*br(4) + xi * f2(2)*br(5)
          sm(3,i) = br(3) + eta * f1(3)*br(4) + xi * f2(3)*br(5)

c         Weight Stresses by Jacobian

          sm(1,i) = sm(1,i) * xjw(i)
          sm(2,i) = sm(2,i) * xjw(i)
          sm(3,i) = sm(3,i) * xjw(i)

        end do ! i

c       Node Loop #1

        do i = 1 , 4
          ni = (i-1)*ndf

c         Compute Residual

          do k = 1 , 6
            p(ni+k) = p(ni+k) - g(i,1,k)*br(1)
     &                        - g(i,2,k)*br(2)
     &                        - g(i,3,k)*br(3)
     &                        - g(i,4,k)*br(4)
     &                        - g(i,5,k)*br(5)
          end do ! k

c         Node Loop #2

          do j = 1 , 4
            nj = (j-1)*ndf

c           Compute Stiffness

            do n = 1,5
              do k1 = 1,6
                gh(k1,n) = g(i,1,k1)*hin(1,n)
     &                   + g(i,2,k1)*hin(2,n)
     &                   + g(i,3,k1)*hin(3,n)
     &                   + g(i,4,k1)*hin(4,n)
     &                   + g(i,5,k1)*hin(5,n)
              end do ! k1
            end do ! n

            do k1 = 1 , 6
              do k2 = 1 , 6
                s(ni+k1,nj+k2) = s(ni+k1,nj+k2)
     &                         + gh(k1,1) * g(j,1,k2)
     &                         + gh(k1,2) * g(j,2,k2)
     &                         + gh(k1,3) * g(j,3,k2)
     &                         + gh(k1,4) * g(j,4,k2)
     &                         + gh(k1,5) * g(j,5,k2)
              end do ! k2
            end do ! k1

          end do ! j

        end do ! i

        end
