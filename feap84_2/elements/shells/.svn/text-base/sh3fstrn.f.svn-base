c$Id:$
        subroutine sh3fstrn ( shp1 , shp2 , cphm , cpx1 , cpx2 , cdrm ,
     &                        cdx1 , cdx2 , ce   , cx   , cr          )

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c        Description:    SH3FSTRN is the subroutine which computes the
c                        strain measures for the inextensible shell
c                        element. The membrane and bending strains
c                        are defined in a cartesian reference frame.
c                        The shear strains are computed in the natural
c                        frame using Bathe-Dvorkin interpolations.

c        Authors:        M.S. Rifai, J.C. Simo, & D.D. Fox.

c        Date:           January 1991.
c-----[--.----+----.----+----.-----------------------------------------]
c        Routine Input:
c        --------------
c        cpx1,cpx2 ..... Current coordinate global derivatives
c                        at the Gauss points.
c        cphm .......... Current coordinate local derivatives
c                        at the midside nodes.
c        cdx1,cdx2 ..... Current director global derivatives
c                        at the Gauss points.
c        cdrm .......... Current local directors
c                        at the midside nodes.

c        Routine Output:
c        ---------------
c        ce............. Current membrane strain measure.
c        cx............. Current shear strain measure.
c        cr............. Current bending strain measure.
c-----[--.----+----.----+----.-----------------------------------------]

        implicit  none

        real*8    shp1 (4)   , shp2 (4)
        real*8    cphm (3,4) , cpx1 (3) , cpx2 (3)
        real*8    cdrm (3,4) , cdx1 (3) , cdx2 (3)
        real*8    ce   (3)   , cx   (2) , cr   (3)

        save

c       Compute Membrane Strains:

        ce(1) = (cpx1(1)*cpx1(1)+cpx1(2)*cpx1(2)+cpx1(3)*cpx1(3))*0.5d0
        ce(2) = (cpx2(1)*cpx2(1)+cpx2(2)*cpx2(2)+cpx2(3)*cpx2(3))*0.5d0
        ce(3) =  cpx1(1)*cpx2(1)+cpx1(2)*cpx2(2)+cpx1(3)*cpx2(3)

c       Compute Shear Strains:

        cx(1) = 2.d0 * ( shp1(2) * (cphm(1,2)*cdrm(1,2)
     &                           +  cphm(2,2)*cdrm(2,2)
     &                           +  cphm(3,2)*cdrm(3,2))
     &                 + shp1(3) * (cphm(1,4)*cdrm(1,4)
     &                           +  cphm(2,4)*cdrm(2,4)
     &                           +  cphm(3,4)*cdrm(3,4)) )

        cx(2) = 2.d0 * ( shp2(4) * (cphm(1,1)*cdrm(1,1)
     &                           +  cphm(2,1)*cdrm(2,1)
     &                           +  cphm(3,1)*cdrm(3,1))
     &                 + shp2(3) * (cphm(1,3)*cdrm(1,3)
     &                           +  cphm(2,3)*cdrm(2,3)
     &                           +  cphm(3,3)*cdrm(3,3)) )

c       Compute Bending Strains:

        cr(1) = cpx1(1)*cdx1(1) + cpx1(2)*cdx1(2) + cpx1(3)*cdx1(3)
        cr(2) = cpx2(1)*cdx2(1) + cpx2(2)*cdx2(2) + cpx2(3)*cdx2(3)
        cr(3) = cpx1(1)*cdx2(1) + cpx1(2)*cdx2(2) + cpx1(3)*cdx2(3)
     &        + cpx2(1)*cdx1(1) + cpx2(2)*cdx1(2) + cpx2(3)*cdx1(3)

        end

        subroutine sh3fcmtx ( d      , dthk , xjw  ,
     &                        optmix , zph1 , zph2 , c   )
c-----[--.----+----.----+----.-----------------------------------------]
c        Description:    SH3FCMTX is the subroutine which sets-up the
c                        linear resultant constitutive relations.

c        Authors:        M.S. Rifai, J.C. Simo, & D.D. Fox.

c        Date:           January, 1991.

c        Version:        This routine was tested in FEAP version 6.01
c-----[--.----+----.----+----.-----------------------------------------]
c        Routine Input:
c        --------------
c        d(*) .......... Material parameters from inmate
c        dthk .......... Element thickness
c        xjw ........... Jacobian times quadrature weight
c        optmix ........ Mixed treatment option:
c                        0 => Galerkin (displacement formulation),
c                        1 => Hellinger-Reissner (Pian-Summihara).
c        zph1,zph2 ..... Initial coordinate local derivatives.

c        Routine Output:
c        ---------------
c        c ............. Resultant stress constitutive matrix.
c-----[--.----+----.----+----.-----------------------------------------]
        implicit  none

        logical   optmix
        real*8    a11     , a22     , a12     , adet
        real*8    dthk    , dthk3   , xjw
        real*8    d(*)    , zph1(3) , zph2(3) , c(8,8) , ain(2,2)

c       Calculate Metric

        a11 = zph1(1)*zph1(1) + zph1(2)*zph1(2) + zph1(3)*zph1(3)
        a22 = zph2(1)*zph2(1) + zph2(2)*zph2(2) + zph2(3)*zph2(3)
        a12 = zph1(1)*zph2(1) + zph1(2)*zph2(2) + zph1(3)*zph2(3)

c       Calculate Inverse of Metric

        adet     = 1.d0/(a11*a22-a12**2)
        ain(1,1) = a22*adet
        ain(1,2) =-a12*adet
        ain(2,1) = ain(1,2)
        ain(2,2) = a11*adet

        if (.not.optmix) then

c         [Cm] Part (Membrane):

          c(1,1) = d(21)*dthk*xjw
          c(1,2) = d(24)*dthk*xjw
          c(1,3) = 0.d0
          c(2,2) = d(22)*dthk*xjw
          c(2,3) = 0.d0
          c(3,3) = d(27)*dthk*xjw
          c(2,1) = c(1,2)
          c(3,1) = c(1,3)
          c(3,2) = c(2,3)

c         [Cb] Part (Bending):

          dthk3  = dthk**3/12.d0*xjw
          c(6,6) = d(21)*dthk3
          c(6,7) = d(24)*dthk3
          c(6,8) = 0.d0
          c(7,7) = d(22)*dthk3
          c(7,8) = 0.d0
          c(8,8) = d(27)*dthk3
          c(7,6) = c(6,7)
          c(8,6) = c(6,8)
          c(8,7) = c(7,8)
        endif

c       [Cs] Part (Shear):

        dthk3  = d(37)*d(28)*dthk*xjw
        c(4,4) = dthk3 * ain(1,1)
        c(4,5) = dthk3 * ain(1,2)
        c(5,5) = dthk3 * ain(2,2)
        c(5,4) = c(4,5)

        end

        subroutine sh3fstrs ( c      , ce     , cx     , cr ,
     &                        optmix , sn     , sq     , sm )
c-----[--.----+----.----+----.-----------------------------------------]
c        Description:    SH3FSTRS is a subroutine which computes the
c                        stress resultants, given the strain measures
c                        and a linear constitutive matrix.

c        Authors:        M.S. Rifai, J.C. Simo, & D.D. Fox.

c        Date:           January, 1991.

c        Version:        This routine was tested in FEAP version 6.3
c-----[--.----+----.----+----.-----------------------------------------]
c        Routine Input:
c        --------------
c        c ............. Resultant stress constitutive matrix.
c        ce............. Current membrane strain measure.
c        cx............. Current shear strain measure.
c        cr............. Current bending strain measure.
c        optmix ........ Mixed treatment option:
c                        0 => Galerkin (displacement formulation),
c                        1 => Hellinger-Reissner (Pian-Summihara).

c        Routine Output:
c        ---------------
c        sn ............ Membrane stress.
c        sq ............ Shear stress.
c        sm ............ Bending stress.
c-----[--.----+----.----+----.-----------------------------------------]

        implicit  none

        logical   optmix
        real*8    c  (8,8) ,
     1            ce (3)   , cx (2) , cr (3) ,
     2            sn (3)   , sq (2) , sm (3)

        if (.not.optmix) then

c         Membrane Stresses:

          sn(1) = c(1,1)*ce(1) + c(1,2)*ce(2) + c(1,3)*ce(3)
          sn(2) = c(2,1)*ce(1) + c(2,2)*ce(2) + c(2,3)*ce(3)
          sn(3) = c(3,1)*ce(1) + c(3,2)*ce(2) + c(3,3)*ce(3)

c         Calculate Bending Stresses:

          sm(1) = c(6,6)*cr(1) + c(6,7)*cr(2) + c(6,8)*cr(3)
          sm(2) = c(7,6)*cr(1) + c(7,7)*cr(2) + c(7,8)*cr(3)
          sm(3) = c(8,6)*cr(1) + c(8,7)*cr(2) + c(8,8)*cr(3)
        endif

c       Calculate Shear Stresses:

        sq(1) = c(4,4)*cx(1) + c(4,5)*cx(2)
        sq(2) = c(5,4)*cx(1) + c(5,5)*cx(2)

        end
