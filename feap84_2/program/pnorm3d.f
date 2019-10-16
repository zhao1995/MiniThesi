c$Id:$
      subroutine pnorm3d(enorm,xl,ndm,nel)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Compute normal vectors (unscaled) for 3-d tetrahedral
c               and brick elements

c      Inputs:
c         xl(ndm,*)  - Element nodal coordinates
c         ndm        - Spatial dimension of mesh
c         nel        - Number of nodes on element

c      Outputs:
c         enorm(3,*) - Unscaled normal vector components
c-----[--.----+----.----+----.-----------------------------------------]
      implicit   none

      integer    ndm,nel, i,nface, faceq(4,6), facet(3,4)
      real*8     v1(3),v2(3), enorm(3,*), xl(ndm,*)

      save

      data       faceq / 1,4,3,2, 1,2,6,5, 2,3,7,6,
     &                   3,4,8,7, 4,1,5,8, 5,6,7,8 /

      data       facet / 1,2,4, 2,3,4, 3,1,4, 1,3,2 /

c     Tetrahedra

      if(nel.eq.4) then

        do nface = 1,4
          do i = 1,3
            v1(i) = xl(i,facet(2,nface)) - xl(i,facet(1,nface))
            v2(i) = xl(i,facet(3,nface)) - xl(i,facet(1,nface))
          end do ! i
          enorm(1,nface) = v1(2)*v2(3) - v1(3)*v2(2)
          enorm(2,nface) = v1(3)*v2(1) - v1(1)*v2(3)
          enorm(3,nface) = v1(1)*v2(2) - v1(2)*v2(1)
        end do ! nface

c     Brick

      elseif(nel.eq.8) then

        do nface = 1,8
          do i = 1,3
            enorm(i,nface) = 0.0d0
          end do ! i
        end do ! nface

        do nface = 1,6
          do i = 1,3
            v1(i) = xl(i,faceq(3,nface)) - xl(i,faceq(1,nface))
            v2(i) = xl(i,faceq(4,nface)) - xl(i,faceq(2,nface))
          end do ! i
          do i = 1,4
            enorm(1,faceq(i,nface)) = enorm(1,faceq(i,nface))
     &                              + v1(2)*v2(3) - v1(3)*v2(2)
            enorm(2,faceq(i,nface)) = enorm(2,faceq(i,nface))
     &                              + v1(3)*v2(1) - v1(1)*v2(3)
            enorm(3,faceq(i,nface)) = enorm(3,faceq(i,nface))
     &                              + v1(1)*v2(2) - v1(2)*v2(1)
          end do ! i
        end do ! nface

      endif

      end
