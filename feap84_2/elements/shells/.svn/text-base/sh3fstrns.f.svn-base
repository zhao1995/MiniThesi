c$Id:$
        subroutine sh3fstrns( shp1  , shp2  ,
     &                        zphm  , zpx1  , zpx2 ,
     &                        zdrm  , zdx1  , zdx2 ,
     &                        uphm  , upx1  , upx2 ,
     &                        udrm  , udx1  , udx2 ,
     &                        ce    , cx    , cr          )

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
c        zpx1,zpx2 ..... Reference coordinate global derivatives
c                        at the Gauss points.
c        zphm .......... Reference coordinate local derivatives
c                        at the midside nodes.
c        zdx1,zdx2 ..... Reference director global derivatives
c                        at the Gauss points.
c        zdrm .......... Reference local directors
c                        at the midside nodes.
c        upx1,upx2 ..... Current coordinate global derivatives
c                        at the Gauss points.
c        uphm .......... Current coordinate local derivatives
c                        at the midside nodes.
c        udx1,udx2 ..... Current director global derivatives
c                        at the Gauss points.
c        udrm .......... Current local directors
c                        at the midside nodes.

c        Routine Output:
c        ---------------
c        ce............. Current membrane strain measure.
c        cx............. Current shear strain measure.
c        cr............. Current bending strain measure.
c-----[--.----+----.----+----.-----------------------------------------]

        implicit  none

        real*8    shp1 (4)   , shp2 (4)
        real*8    zphm (3,4) , zpx1 (3) , zpx2 (3)
        real*8    zdrm (3,4) , zdx1 (3) , zdx2 (3)
        real*8    uphm (3,4) , upx1 (3) , upx2 (3)
        real*8    udrm (3,4) , udx1 (3) , udx2 (3)
        real*8    ce   (3)   , cx   (2) , cr   (3)

        save

c       Compute Membrane Strains:

        ce(1) =  zpx1(1)*upx1(1)+zpx1(2)*upx1(2)+zpx1(3)*upx1(3)
     &        + (upx1(1)*upx1(1)+upx1(2)*upx1(2)+upx1(3)*upx1(3))*0.5d0

        ce(2) =  zpx2(1)*upx2(1)+zpx2(2)*upx2(2)+zpx2(3)*upx2(3)
     &        + (upx2(1)*upx2(1)+upx2(2)*upx2(2)+upx2(3)*upx2(3))*0.5d0

        ce(3) =  zpx1(1)*upx2(1)+zpx1(2)*upx2(2)+zpx1(3)*upx2(3)
     &        +  upx1(1)*zpx2(1)+upx1(2)*zpx2(2)+upx1(3)*zpx2(3)
     &        +  upx1(1)*upx2(1)+upx1(2)*upx2(2)+upx1(3)*upx2(3)

c       Compute Shear Strains:

        cx(1) = 2.d0 * ( shp1(2) * (zphm(1,2)*udrm(1,2)
     &                           +  zphm(2,2)*udrm(2,2)
     &                           +  zphm(3,2)*udrm(3,2)

     &                           +  uphm(1,2)*zdrm(1,2)
     &                           +  uphm(2,2)*zdrm(2,2)
     &                           +  uphm(3,2)*zdrm(3,2)

     &                           +  uphm(1,2)*udrm(1,2)
     &                           +  uphm(2,2)*udrm(2,2)
     &                           +  uphm(3,2)*udrm(3,2))

     &                 + shp1(3) * (zphm(1,4)*udrm(1,4)
     &                           +  zphm(2,4)*udrm(2,4)
     &                           +  zphm(3,4)*udrm(3,4)

     &                           +  uphm(1,4)*zdrm(1,4)
     &                           +  uphm(2,4)*zdrm(2,4)
     &                           +  uphm(3,4)*zdrm(3,4)

     &                           +  uphm(1,4)*udrm(1,4)
     &                           +  uphm(2,4)*udrm(2,4)
     &                           +  uphm(3,4)*udrm(3,4)) )

        cx(2) = 2.d0 * ( shp2(4) * (zphm(1,1)*udrm(1,1)
     &                           +  zphm(2,1)*udrm(2,1)
     &                           +  zphm(3,1)*udrm(3,1)

     &                           +  uphm(1,1)*zdrm(1,1)
     &                           +  uphm(2,1)*zdrm(2,1)
     &                           +  uphm(3,1)*zdrm(3,1)

     &                           +  uphm(1,1)*udrm(1,1)
     &                           +  uphm(2,1)*udrm(2,1)
     &                           +  uphm(3,1)*udrm(3,1))

     &                 + shp2(3) * (zphm(1,3)*udrm(1,3)
     &                           +  zphm(2,3)*udrm(2,3)
     &                           +  zphm(3,3)*udrm(3,3)

     &                           +  uphm(1,3)*zdrm(1,3)
     &                           +  uphm(2,3)*zdrm(2,3)
     &                           +  uphm(3,3)*zdrm(3,3)

     &                           +  uphm(1,3)*udrm(1,3)
     &                           +  uphm(2,3)*udrm(2,3)
     &                           +  uphm(3,3)*udrm(3,3)) )

c       Compute Bending Strains:

        cr(1) = zpx1(1)*udx1(1) + zpx1(2)*udx1(2) + zpx1(3)*udx1(3)
     &        + upx1(1)*zdx1(1) + upx1(2)*zdx1(2) + upx1(3)*zdx1(3)
     &        + upx1(1)*udx1(1) + upx1(2)*udx1(2) + upx1(3)*udx1(3)

        cr(2) = zpx2(1)*udx2(1) + zpx2(2)*udx2(2) + zpx2(3)*udx2(3)
     &        + upx2(1)*zdx2(1) + upx2(2)*zdx2(2) + upx2(3)*zdx2(3)
     &        + upx2(1)*udx2(1) + upx2(2)*udx2(2) + upx2(3)*udx2(3)

        cr(3) = zpx1(1)*udx2(1) + zpx1(2)*udx2(2) + zpx1(3)*udx2(3)
     &        + zpx2(1)*udx1(1) + zpx2(2)*udx1(2) + zpx2(3)*udx1(3)
     &        + upx1(1)*zdx2(1) + upx1(2)*zdx2(2) + upx1(3)*zdx2(3)
     &        + upx2(1)*zdx1(1) + upx2(2)*zdx1(2) + upx2(3)*zdx1(3)
     &        + upx1(1)*udx2(1) + upx1(2)*udx2(2) + upx1(3)*udx2(3)
     &        + upx2(1)*udx1(1) + upx2(2)*udx1(2) + upx2(3)*udx1(3)

        end
