c$Id:$
      subroutine convec1d(d,ul,xl,ix,s,r,ndf,ndm,nst,isw)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Set pstyp = -1 for no plots                      02/01/2007
c       2. Set pstyp =  0 for no plots                      31/08/2008
c       3. Add time option for q_o                          18/10/2010
c-----[--.----+----.----+----.-----------------------------------------]
c     One dimensional (plane/axisymmetric) thermal surface element

c     N.B. Surface loading for solutions using THERMAL (1d) Element

c-----[--.----+----.----+----.-----------------------------------------]
c     1. Parameters input by this routine when isw = 1 (i.e., mate)

c           H   - Surface parameter
c           T_o - Equilibrium temperature
c           q_o - Normal flux to boundary
c           p(t)- Proportional load
c           nn  - Exponent to convection/radiation b.c

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

c           ndm - set to 1     (x or r-coord)
c           ndf - set > or = 1 (nodal temperature)
c           nel - set > or = 1

c      OUTPUT variables

c        r(1,1)     Contribution to residual

c      PARAMATER set up:

c         kat  =  1  (for plane analysis, constant flux)

c              =  2  (for axisymmetric analysis, constant flux)

c         nel  =  1  (for edge of line element)

c                                 ^
c                      ^  normal) |
c                      |  + q_n   |
c                      |  <-------+
c         Node is:     o           t (tangent)
c                      1
c                      |  interior
c                      |     of mesh

c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'eldata.h'
      include  'eltran.h'
      include  'iofile.h'
      include  'prld1.h'

      logical   errck, tinput, pcomp
      character tx*15, wlab(2)*6
      integer   ndm,ndf,nst, isw, i, kat, nn, ix(*)
      real*8    da, hh, rr, qq, pt, uu
      real*8    d(*), xl(ndm,*),ul(ndf,*), r(ndf,*), s(nst,*)

      save

      data wlab/'Plane ','Axisym'/

c     Output element type

      if(isw.eq.0) then
        if(ior.lt.0) then
          write(*,*) '   1-d Plane/Axisym Surface Thermal Flux'
        endif

c     Input material properties

      elseif(isw.eq.1) then

        tx = 'start'
        do while(.not.pcomp(tx,'    ',4))

 1        if(ior.lt.0) then
            write(*,3000)
            call pprint('   >')
          endif
          errck = tinput(tx,1,d(11),4)
          if(errck) go to 1

          if    (pcomp(tx,'surf',4)) then
            d(1) = d(11)
            d(2) = d(12)
          elseif(pcomp(tx,'flux',4)) then
            d(3) = d(11)
            d(6) = d(11)
          elseif(pcomp(tx,'expo',4)) then
            d(4) = d(11)
          elseif(pcomp(tx,'plan',4)) then
            d(5) = 1.d0
          elseif(pcomp(tx,'axis',4)) then
            d(5) = 2.d0
          endif
        end do ! while

c       Set final parameters

        nn  = max(1,nint(d(4)))
        d(4)= nn
        kat = nint(d(5))
        if(kat.ne.2) kat = 1
        d(5) = kat
        if(ior.lt.0) then
          write(*,2000) d(1),d(2),d(3),nint(d(6)),nn,wlab(kat)
        end if
        write(iow,2000) d(1),d(2),d(3),nint(d(6)),nn,wlab(kat)

c       Set ix to eliminate unused dof

        do i = 2,ndf
          ix(i) = 0
        end do ! i

c       Block plotting

        pstyp = 0

c     Compute conductivity (stiffness) matrix

      elseif(isw.eq.3 .or. isw.eq.6) then

        nn    = nint(d(4))
        kat   = nint(d(5))

c       Compute geometric factors

        rr = xl(1,i)
        uu = ul(1,i)

c       Surface area

        if(kat.eq.2) then
          da = rr*d(14)
        else
          da = d(14)
        endif

c       Set time function for q_o

        if(nint(d(6)).eq.0) then
          pt = 1.0d0
        else
          pt = prldv(nint(d(6)))
        endif

c       Thermal properties and loads for flux on face point

        if(nn.eq.1) then
          qq = ( d(3)*pt + d(1)*( uu - d(2) ) )*da
          hh = d(1)*da*ctan(1)
        else
          qq = ( d(3)*pt + d(1)*( uu**nn - d(2)**nn ) )*da
          hh = d(1)*dble(nn)*uu**(nn-1)*da*ctan(1)
        endif

        r(1,1) = r(1,1) - qq

c       Compute tangent matrix for surface convections

        s(1,1) = hh

      endif

c     Formats

c-----[--.----+----.----+----.-----------------------------------------]

2000  format(5x,'One Dimensional Heat Conduction Boundary Element'//
     &      10x,'Surface Parameter ',1p,e12.5/
     &      10x,'Equilibrium Temp. ',1p,e12.5/
     &      10x,'Boundary flux     ',1p,e12.5/
     &      10x,'Proportional load ',i7/
     &      10x,'Temperature exp. n',i7/
     &      10x,a,' Analysis')

3000  format(' Input: SURFace , H, T_0'/
     &       '        FLUX    , q_n'/
     &       '        AXISymmetric or PLANe'/
     &       '   or   EXPOnent, n')

      end
