c$Id:$
        subroutine sh3flmda ( v , xlm )

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c      1. Revise computation of orthogonal axes to be more  31/05/2013
c         consistent over the mesh
c-----[--.----+----.----+----.-----------------------------------------]
c        Description:    SH3FLMDA is subroutine which computes unique
c                        orthogonal transformation matrix which rotates
c                        the 3-rd canonical basis vector, E3 into the
c                        vector v --- without drill.
c                        The singularity of the formula at (or near) the
c                        case: v = -E3 is avoided by computing matrix
c                        which transforms -E3 into v when v.E3 < 0, and
c                        then applying a 180-degree rotation to it.

c        Authors:        M.S. Rifai, J.C. Simo, & D.D. Fox.

c        Date:           January 1991.
c-----[--.----+----.----+----.-----------------------------------------]
c        Routine Input:
c        --------------
c        v ............. Arbitrary vector in R-3.

c        Routine Output:
c        ---------------
c        xlm ........... Orthogonal transformation matrix, which
c                        rotates E3 into v without drill.
c-----[--.----+----.----+----.-----------------------------------------]
        implicit  none

        include  'crotas.h'

        integer   imin , imax , i   , option
        real*8    vmin , vmax , fac , tol , v(3) , xlm(3,2)

        data      tol / 1.0d-10 /

c       Option 1

        option = rotopt

        if(option.eq.1) then

c         Compute location of maximum and minimum entries in v(:)

          vmin = abs(v(1))
          vmax = vmin
          imin = 1
          imax = 1
          do i = 2,3
            if(abs(v(i)).lt.vmin) then
              vmin = abs(v(i))
              imin = i
            endif
            if(abs(v(i)).gt.vmax) then
              vmax = abs(v(i))
              imax = i
            endif
          end do ! i

c         Location of intermediate component

          i = 6 - imin - imax

c         Compute first surface vector: xlm(:1)

          xlm(imax,1) = -v(imin)
          xlm(i   ,1) =  0.0d0
          xlm(imin,1) =  v(imax)

c       Use 1-2 components

        elseif(option.eq.2) then

          if(max(abs(v(1)),abs(v(2))).gt.1.d-03) then
            xlm(1,1) =  v(2)
            xlm(2,1) = -v(1)
            xlm(3,1) =  0.0d0
          else
            xlm(1,1) =  1.0d0
            xlm(2,1) =  0.0d0
            xlm(3,1) =  0.0d0
          endif

c       Original algorithm

        elseif(option.eq.3) then

c         Compute [XLm], when v.E3 > 0

          if (v(3).gt.tol) then
             fac      =   1.d0 / ( 1.d0 + v(3) )
             xlm(1,1) =   v(3) + fac * v(2) * v(2)
             xlm(2,1) =        - fac * v(1) * v(2)
             xlm(3,1) = - v(1)
c            xlm(1,2) =        - fac * v(2) * v(1)
c            xlm(2,2) =   v(3) + fac * v(1) * v(1)
c            xlm(3,2) = - v(2)

c         Compute [XLm], when v.E3 < 0

          else
             fac      =   1.d0 / ( 1.d0 - v(3) )
             xlm(1,1) = - v(3) + fac * v(2) * v(2)
             xlm(2,1) =        - fac * v(1) * v(2)
             xlm(3,1) =   v(1)
c            xlm(1,2) =          fac * v(2) * v(1)
c            xlm(2,2) =   v(3) - fac * v(1) * v(1)
c            xlm(3,2) = - v(2)
          endif

        endif

        vmax = 1.d0/sqrt(xlm(1,1)**2 +xlm(2,1)**2 +xlm(3,1)**2)
        xlm(1,1) = xlm(1,1)*vmax
        xlm(2,1) = xlm(2,1)*vmax
        xlm(3,1) = xlm(3,1)*vmax

c       Compute second surface vector: xlm(:2)

        xlm(1,2) = v(2)*xlm(3,1) - v(3)*xlm(2,1)
        xlm(2,2) = v(3)*xlm(1,1) - v(1)*xlm(3,1)
        xlm(3,2) = v(1)*xlm(2,1) - v(2)*xlm(1,1)

        end
