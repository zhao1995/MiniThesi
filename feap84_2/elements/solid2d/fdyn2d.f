c$Id:$
      subroutine fdyn2d(ul,shp,s,r,is,xcur,cfac,lfac,dmas0,isw)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose:  2-D Dynamics: Computation of residual and tangent.

c      Inputs:
c         d(*)      - Material set parameters
c         ul(ndf,*) - Nodal solution parameters for element
c         xcur      - Current radius
c         cfac      - Consistent factor
c         lfac      - Lumped factor
c         dmas0     - Lumped factor
c         is        - Degree-of-freedom range
c         isw       - Switch to control action

c      Outputs:
c         s(nst,*)  - Element matrix
c         p(nst)    - Element vector
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'cdata.h'
      include  'eldata.h'
      include  'eltran.h'
      include  'pmod2d.h'
      include  'sdata.h'

      integer   isw,is, i,i1,j,j1
      real*8    xcur,cfac,lfac,bdb,dl,dc,di,dmas0,umas0,vmas0
      real*8    ul(ndf,nen,*),shp(3,*),s(nst,*),r(ndf,*)
      real*8    al(3),ac(3),vl(3),vc(3)

c     Compute accelerations

      do i = 1,is
        al(i) = 0.0d0
        vl(i) = 0.0d0
        do j = 1,nel
          al(i) = al(i) + ul(i,j,5)*shp(3,j)
          vl(i) = vl(i) + ul(i,j,4)*shp(3,j)
        end do ! j
        al(i) = al(i)*cfac
        vl(i) = vl(i)*cfac
      end do ! i

c     COMPUTE INERTIA TERMS

      do i = 1,nel

c       Compute inertial effects

        do j = 1,is
          ac(j) = al(j) + lfac*ul(j,i,5)
          vc(j) = vl(j) + lfac*ul(j,i,4)
        end do ! j

c       Element residual

        r(1,i) = r(1,i) - ac(1)*dmas0*shp(3,i)
        r(2,i) = r(2,i) - ac(2)*dmas0*shp(3,i)

c       Torsion residual

        if(stype.eq.8) then
          r(1,i)  = r(1,i) + xcur*vc(3)**2*shp(3,i)*dmas0
          r(3,i)  = r(3,i) - xcur*(xcur*ac(3)
     &                    + 2.d0*vc(1)*vc(3))*shp(3,i)*dmas0
        endif

      end do ! i

c     COMPUTE K (s(nst,nst) = K)

      if(isw.eq.3) then

c       Inertial part.

        dc  = cfac*dmas0
        dl  = lfac*dmas0
        i1  = 0
        do i = 1,nel

c         Lumped contribution

          do j = 1,is
            s(i1+j,i1+j) = s(i1+j,i1+j) + shp(3,i)*dl*ctan(3)
          end do ! j

c         Consistent contribution

          di  = dc*shp(3,i)

          j1  = 0
          do j = 1,nel
            bdb          = di*shp(3,j)
            s(i1+1,j1+1) = s(i1+1,j1+1) + bdb*ctan(3)
            s(i1+2,j1+2) = s(i1+2,j1+2) + bdb*ctan(3)
            if(stype.eq.8) then
              s(i1+3,j1+3) = s(i1+3,j1+3) + xcur*xcur*bdb*ctan(3)

              vmas0 = ctan(2)*2.d0*xcur*bdb

              s(i1+1,j1+3) = s(i1+1,j1+3) - vc(3)*vmas0
              s(i1+3,j1+1) = s(i1+3,j1+1) + vc(3)*vmas0
              s(i1+3,j1+3) = s(i1+3,j1+3) + vc(1)*vmas0

              umas0 = ctan(1)*bdb
              s(i1+1,j1+1) = s(i1+1,j1+1) - vc(3)**2*umas0
              s(i1+3,j1+1) = s(i1+3,j1+1) + 2.d0*(xcur*ac(3)
     &                                    + vc(1)*vc(3))*umas0
            endif
            j1 = j1 + ndf
          end do ! j
          i1 = i1 + ndf
        end do ! i

      endif ! isw = 3

      end
