c$Id:$
        subroutine sh3fintr ( cphi , dir , cphm , cdrm )

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c        Description:    SH3FINTR is the subroutine which interpolates
c                        the midsurface derivatives and directors
c                        at the midside points.

c        Authors:        M.S. Rifai, J.C. Simo, & D.D. Fox.

c        Date:           January 1991.
c-----[--.----+----.----+----.-----------------------------------------]
c        Routine Input:
c        --------------
c        cphi .......... Current local nodal coordinates.
c        dir ........... Localized nodal director
c                        vectors, dir(3,ElementNode).

c        Routine Output:
c        ---------------
c        cphm .......... Current coordinate local derivatives
c                        at the midside nodes.
c        cdrm .......... Initial local directors
c                        at the midside nodes.
c-----[--.----+----.----+----.-----------------------------------------]

        implicit  none

        integer   i

        real*8    cphi (3,*) , cphm (3,4) , dir(3,9) , cdrm (3,4)

        save

c       Interpolate Midsurface Positions

        do i = 1,3
          cphm(i,1) = ( cphi(i,4) - cphi(i,1) ) * 0.5d0
          cphm(i,2) = ( cphi(i,2) - cphi(i,1) ) * 0.5d0
          cphm(i,3) = ( cphi(i,3) - cphi(i,2) ) * 0.5d0
          cphm(i,4) = ( cphi(i,3) - cphi(i,4) ) * 0.5d0
        end do ! i

c       Interpolate Directors to Midsurface

        do i = 1,3
          cdrm(i,1) = ( dir(i,4) + dir(i,1) ) * 0.5d0
          cdrm(i,2) = ( dir(i,1) + dir(i,2) ) * 0.5d0
          cdrm(i,3) = ( dir(i,2) + dir(i,3) ) * 0.5d0
          cdrm(i,4) = ( dir(i,3) + dir(i,4) ) * 0.5d0
        end do ! i

        end

        subroutine sh3fshap ( xi , shp , shp1 , shp2 )
c-----[--.----+----.----+----.-----------------------------------------]
c        Description:    SH3FSHAP is the subroutine which computes the
c                        element shape functions and local derivatives
c                        for the 4-noded quadrilateral element.
c-----[--.----+----.----+----.-----------------------------------------]
c        Routine Input:
c        --------------
c        xi ............ Local coordinates of the current point.

c        Routine Output:
c        ---------------
c        shp ........... Nodal shape function values.
c        shp1,shp2 ..... Nodal shape function natural derivatives.
c-----[--.----+----.----+----.-----------------------------------------]

        implicit  none

        real*8    xi(2) , shp(4) , shp1(4) , shp2(4)

c       Evaluate (xi) Derivative:

        shp1(2) =  0.25d0 * (1.d0-xi(2))
        shp1(3) =  0.25d0 * (1.d0+xi(2))
        shp1(4) = -shp1(3)
        shp1(1) = -shp1(2)

c       Evaluate (xi) Derivative:

        shp2(3) =  0.25d0 * (1.d0+xi(1))
        shp2(4) =  0.25d0 * (1.d0-xi(1))
        shp2(1) = -shp2(4)
        shp2(2) = -shp2(3)

c       Evaluate Shape Function:

        shp (1) =  shp2(4) * (1.d0-xi(2))
        shp (2) =  shp2(3) * (1.d0-xi(2))
        shp (3) =  shp2(3) * (1.d0+xi(2))
        shp (4) =  shp2(4) * (1.d0+xi(2))

c       Exit

        end

        subroutine sh3fresd ( b      , sn     , sq     , sm ,
     2                        optmix , p                    )
c-----[--.----+----.----+----.-----------------------------------------]
c        Description:    SH3FRESD is the routine which computes the
c                        static residual.
c-----[--.----+----.----+----.-----------------------------------------]
c        Routine Input:
c        --------------
c        B ............. Discrete strain-displacement operator.
c        sn ............ Membrane stress.
c        sq ............ Shear stress.
c        sm ............ Bending stress.
c        optmix ........ Mixed treatment option:
c                        0 => Galerkin (displacement formulation),
c                        1 => Hellinger-Reissner (Pian-Summihara).

c        Routine Output:
c        ---------------
c        P ............. Element residual.
c-----[--.----+----.----+----.-----------------------------------------]

        implicit   none

        logical    optmix
        integer    i
        real*8     b(8,24) , sn(3) , sq(2) , sm(3) , p(24)

        if (.not.optmix) then

c         Calculate Membrane, Shear and Bending Parts

          do i = 1,24
            p(i) = p(i) - b(1,i)*sn(1) - b(2,i)*sn(2) - b(3,i)*sn(3)
     &                  - b(4,i)*sq(1) - b(5,i)*sq(2)
     &                  - b(6,i)*sm(1) - b(7,i)*sm(2) - b(8,i)*sm(3)
          end do ! i

        else

c         Calculate Shear Part

          do i = 1,24
            p(i) = p(i) - b(4,i)*sq(1) - b(5,i)*sq(2)
          end do ! i

        endif


        end

        subroutine sh3fstif ( b      , b1     , c      ,
     &                       optmix , s      )
c-----[--.----+----.----+----.-----------------------------------------]
c        Description:    SH3FSTIF is the routine which computes the
c                        material tangent, i.e., [B-t][C][B].
c-----[--.----+----.----+----.-----------------------------------------]
c        Routine Input:
c        --------------
c        B ............. Discrete strain-displacement operator.
c        C ............. Resultant stress constitutive matrix.
c        OptMix ........ Mixed treatment option:
c                        0 => Galerkin (displacement formulation),
c                        1 => Hellinger-Reissner (Pian-Summihara).

c        Routine Output:
c        ---------------
c        s ............. Element stiffness.
c-----[--.----+----.----+----.-----------------------------------------]

        implicit  none

        logical   optmix
        integer   i       , j

        real*8    b(8,24) , b1(8,24) ,c(8,8) , s(24,24) , temp(24,8)

        if (.not.optmix) then

c         Calculate Membrane [B-transpose] . [C]

          do j = 1,3
            do i = 1,24
              temp(i,j) = b(1,i)*c(1,j) + b(2,i)*c(2,j) + b(3,i)*c(3,j)
            end do ! i
          end do ! j

c         Calculate Bending [B-transpose] . [C]

          do j = 6,8
            do i = 1,24
              temp(i,j) = b(6,i)*c(6,j) + b(7,i)*c(7,j) + b(8,i)*c(8,j)
            end do ! i
          end do ! j
        endif

c       Calculate Shear [B-transpose] . [C]

        do j = 4,5
          do i = 1,24
            temp(i,j) = b(4,i)*c(4,j) + b(5,i)*c(5,j)
          end do ! i
        end do ! j

c       Calculate [B-transpose] . [C] . [B]

        if (.not.optmix) then

          do j = 1,24
            do i = 1,24
              s(i,j) = s(i,j) + temp(i,1)*b1(1,j) + temp(i,2)*b1(2,j)
     &                        + temp(i,3)*b1(3,j) + temp(i,4)*b1(4,j)
     &                        + temp(i,5)*b1(5,j) + temp(i,6)*b1(6,j)
     &                        + temp(i,7)*b1(7,j) + temp(i,8)*b1(8,j)
            end do ! i
          end do ! j

        else

c         Calculate Shear [B-transpose] . [C] . [B]

          do j = 1,24
            do i = 1,24
              s(i,j) = s(i,j) + temp(i,4)*b1(4,j) + temp(i,5)*b1(5,j)
            end do ! i
          end do ! j

        endif

        end

        subroutine sh3fplst ( sn , sq , sm  ,
     &                        dt , st , shp , xjw , lint )
c-----[--.----+----.----+----.-----------------------------------------]
c        Description:    SH3FPLST is the stress projection subroutine
c                        that computes the element contributions to the
c                        global equation system of the least-squares
c                        problem. The `mass' matrix is lumped with
c                        row-sum.
c                        The stresses are arranged as:
c                        1-3 := Membrane resultants.
c                        4-5 := Shear resultants.
c                        6-8 := Bending resultants.

c        Authors:        M.S. Rifai, J.C. Simo, & D.D. Fox.

c        Date:           January, 1991.
c-----[--.----+----.----+----.-----------------------------------------]
c        Routine Input:
c        --------------
c        sn ............ Membrane stress.
c        sq ............ Shear stress.
c        sm ............ Bending stress.
c        shp ........... Nodal shape function values.
c        xjw ........... Reference jacobian weighted for integration.
c        lint........... Number of quadrature points

c        Routine Output:
c        ---------------
c        dt ............ R.H.S. of stress projection equation system.
c        st ............ L.H.S. of stress projection equation system,
c                        lumped by row-sum.
c-----[--.----+----.----+----.-----------------------------------------]

        implicit  none

        include  'cdata.h'
        include  'strnum.h'
        include  'eldata.h'
        include  'iofile.h'

        integer   i , l      , lint
        real*8    xjw(4)     , xjshp

        real*8    sn (3,4)   , sq(2,4)     , sm(3,4) ,
     &            dt (*)     , st(nen,*)   , shp(4,4)

c       Loop over quadrature points

        do l = 1,lint

c         Loop over Nodes:

          do i = 1,4
            xjshp   = xjw(l)*shp(i,l)
            dt(i)   = dt(i) + xjshp

c           Membrane:

            st(i, 1) = st(i, 1) + sn(1,l)*xjshp
            st(i, 2) = st(i, 2) + sn(2,l)*xjshp
            st(i, 4) = st(i, 4) + sn(3,l)*xjshp

c           Shear:

            st(i, 5) = st(i, 5) + sq(1,l)*xjshp
            st(i, 6) = st(i, 6) + sq(2,l)*xjshp

c           Bending:

            st(i, 7) = st(i, 7) + sm(1,l)*xjshp
            st(i, 8) = st(i, 8) + sm(2,l)*xjshp
            st(i,10) = st(i,10) + sm(3,l)*xjshp
          end do ! i

        end do ! l

        iste = 10

        end
