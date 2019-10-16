c$Id:$
      subroutine convec2d(d,ul,xl,ix,s,r,ndf,ndm,nst,isw)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Set pstyp = -1 for no plots                      02/01/2007
c       2. Set pstyp =  0 for no plots                      31/08/2008
c       3. Revise for higher order elements                 27/03/2009
c       4. Add time option for q_0                          18/10/2010
c-----[--.----+----.----+----.-----------------------------------------]
c     Two dimensional (plane/axisymmetric) thermal surface element

c     N.B. Surface loading for solutions using THERMAL (2d) Element

c-----[--.----+----.----+----.-----------------------------------------]
c     1. Parameters input by this routine when isw = 1 (i.e., mate)

c           H   - Surface parameter
c           T_o - Equilibrium temperature
c           q_o - Normal flux to boundary
c           nn  - Exponent to convection/radiation b.c
c           p(t)- Proportional loading function

c                 flux = q_o*p(t) + H*(T^n - T_o^n)
c                 nn = 1 - convection b.c.
c                 nn = 4 - Stefan-Boltzman radiation b.c.

c           kat - geometry type
c                 1 = plane
c                 2 = axisymmetric
c           N.B. If kat is not = 2 it is set to 1.
c-----[--.----+----.----+----.-----------------------------------------]

c     2. Control parameters

c        This is a two dimensional element which can analyze plane
c        or axisymmetric geometries.  Set control parameters as
c        follows:

c           ndm - set to 2     (x,y or r,z-coords)
c           ndf - set > or = 1 (nodal temperatures)
c           nel - set > or = 2

c....  OUTPUT variables

c        r(1,nel)     Contribution to residual

c....  PARAMATER set up:

c         kat  =  1  (for plane analysis, constant flux)

c              =  2  (for axisymmetric analysis, constant flux)

c         nel  =  2  (for 2-node edge of bi-linear element)

c                                                         ^
c                           ^^^                n (normal) |
c                           |||  + q_n                    |
c                           |||                   <-------+
c         Nodes are:     o----------------o           t (tangent)
c                        2                1
c                        |  interior      |
c                        |     of mesh    |

c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'eldata.h'
      include  'eltran.h'
      include  'iofile.h'
      include  'prld1.h'

      logical   errck, tinput, pcomp
      character tx*15, wlab(2)*6
      integer   ndm,ndf,nst, isw, i,j, i1,j1, kat, l, lint, nn, ix(*)
      real*8    da, xsj, pt, hh, rr, qq, tt, uu, shp1(2,4), sg(2,5)
      real*8    d(*), xl(ndm,*),ul(ndf,*), r(ndf,*), s(nst,*)

      save

      data wlab/'Plane ','Axisym'/

c     Output element type

      if(isw.eq.0) then
        if(ior.lt.0) then
          write(*,*) '   2-d Plane/Axisym Surface Thermal Flux'
        endif

c     Input material properties

      elseif(isw.eq.1) then
    1   nn = 0
        errck = tinput(tx,1,d(11),4)
        if(errck) go to 1

        if    (pcomp(tx,'surf',4)) then
           d(1) = d(11)
           d(2) = d(12)
        elseif(pcomp(tx,'flux',4)) then
           d(3) = d(11)
           d(6) = d(12)
        elseif(pcomp(tx,'expo',4)) then
           d(4) = max(1,nint(d(11)))
        elseif(pcomp(tx,'axis',4)) then
           d(5) = 2.d0
        elseif(pcomp(tx,'plan',4)) then
           d(5) = 1.d0
        elseif(pcomp(tx,'    ',4)) then
           go to 11
        endif
        go to 1

c       Set final parameters

   11   nn  = max(1,nint(d(4)))
        d(4)= nn
        kat = nint(d(5))
        if(kat.ne.2) kat = 1
        d(5) = dble(kat)
        if(ior.lt.0) then
          write(*,2000) d(1),d(2),d(3),nint(d(6)),nn,wlab(kat)
        end if
        write(iow,2000) d(1),d(2),d(3),nint(d(6)),nn,wlab(kat)

        do i = 2,ndf
          ix(i) = 0
        end do ! i

c       Set for no plots

        pstyp = 0

c     Compute conductivity (stiffness) matrix

      elseif(isw.eq.3 .or. isw.eq.6) then

c       Set quadrature

        lint = min(5,nint(d(5)))
        if(lint.eq.0) then
          lint = nel
        endif
        if(nint(d(182)).gt.0) then
          call int1dn(lint,sg)
        else
          call int1d (lint,sg)
        endif

        nn    = nint(d(4))
        kat   = nint(d(5))

c       Set time function for q_o

        if(nint(d(6)).eq.0) then
          pt = 1.0d0
        else
          pt = prldv(nint(d(6)))
        endif

c       Loop over quadrature points

        do l = 1,lint ! {

c         Compute geometric factors

          call shp1d (sg(1,l), xl,  shp1, ndm,nel,xsj)

          rr = 0.0d0
          uu = 0.0d0
          do i = 1,nel ! {
            rr = rr + shp1(2,i)*xl(1,i)
            uu = uu + shp1(2,i)*ul(1,i)
          end do ! i   }

c         Compute surface area ( no radius effect yet)

          da = xsj*sg(2,l)

c         Axisymmetric problems: Modify terms for radius

          if(kat.eq.2) then
            da = da*rr
          endif

c         Thermal properties and loads for flux on face point

          if(nn.eq.1) then
            qq = ( d(3)*pt + d(1)*( uu - d(2) ) )*da
            tt = d(1)*da*ctan(1)
          else
            qq = ( d(3)*pt + d(1)*( uu**nn - d(2)**nn ) )*da
            tt = d(1)*dble(nn)*uu**(nn-1)*da*ctan(1)
          endif

          i1 = 1
          do i = 1,nel ! {
            r(1,i) = r(1,i) - qq*shp1(2,i)

c           Compute tangent matrix for surface convections

            hh = tt*shp1(2,i)
            j1 = 1
            do j = 1,nel ! {
              s(i1,j1) = s(i1,j1) + hh*shp1(2,j)
              j1 = j1 + ndf
            end do ! j   }
            i1 = i1 + ndf
          end do ! i   }

        end do ! l   }

      endif

c     Formats

c-----[--.----+----.----+----.-----------------------------------------]

2000  format(5x,'Two Dimensional Heat Conduction Boundary Element'//
     &      10x,'Surface Parameter ',1p,e12.5/
     &      10x,'Equilibrium Temp. ',1p,e12.5/
     &      10x,'Boundary flux     ',1p,e12.5/
     &      10x,'Proportional load ',i7/
     &      10x,'Temperature exp. n',i7/
     &      10x,a,' Analysis')

      end
