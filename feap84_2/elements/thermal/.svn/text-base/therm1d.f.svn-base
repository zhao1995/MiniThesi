c$Id:$
      subroutine therm1d(d,ul,xl,ix,s,p,ndf,ndm,nst,isw)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Correct lumped s(j1,j1) (was i1) for isw=5       15/08/2007
c       2. Revise quadrature sets for d(5)                  27/03/2009
c       3. Save xref in 'elcoor.h'                          26/07/2009
c       4. Correct format 2003                              29/03/2011
c-----[--.----+----.----+----.-----------------------------------------]
c     One dimensional (plane/axisymmetric) Linear Thermal Element

c-----[--.----+----.----+----.-----------------------------------------]

c        This is a one dimensional element which can analyze plane
c        or axisymmetric geometries.  Set control parameters as
c        follows:

c           ndm - set to 1     (x or r-coord)
c           ndf - set > or = 1 (nodal temperature)
c           nel - set > or = 2

c              o-----------o ----> x or r
c             1             2

c               Node numbering
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'bdata.h'
      include  'cdata.h'
      include  'eldata.h'
      include  'elplot.h'
      include  'eltran.h'
      include  'fdata.h'
      include  'iofile.h'
      include  'mdata.h'
      include  'part0.h'
      include  'pmod2d.h'
      include  'rdata.h'
      include  'strnum.h'
      include  'comblk.h'

      integer   ndf,ndm,nst,isw, i,j, i1,j1, l,lint, tdof, ix(*)
      real*8    xx, xsj, a1,a3,a4,shj,tdot,cfac,lfac, hh,tinf
      real*8    d(*),ul(ndf,nen,*),xl(ndm,*),s(nst,*),p(*)
      real*8    temp,gradt,flux,dd,shp(2,4),sg(2,5)

      save

c     Set mass factors

      if(d(7).ge.0.0d0 .or. d(183).ne.0.0d0) then
        cfac = d(7)
        lfac = 1.d0 - cfac
      else
        cfac = 0.0d0
        lfac = 0.0d0
      endif

c     Input material properties

      if(isw.eq.1) then

        if(ior.lt.0) write(*,2000)
        write(iow,2000)
        call inmate(d,tdof,0,6)

c       Delete unused parameters

        do i = 2,ndf
          ix(i) = 0
        end do ! i

c       Set to preclude sloping boundary transformations

        ea(1,-iel) = 0
        ea(2,-iel) = 0

c       Set plot sequence

        pstyp = 1
        istv  = max(istv,15)

c     Check of mesh if desired (chec)

      elseif(isw.eq.2) then

        if(nel.eq.3 .or. nel.eq.6 .or. nel.eq.7) then
          call cktris(ix,xl,shp,ndm)
        else
          call ckisop(ix,xl,shp,ndm)
        endif

c     Compute conductivity (stiffness) matrix

      elseif(isw.eq.3 .or. isw.eq.6) then

        lint = min(5,nint(d(5)))
        if(lint.eq.0) then
          lint = nel
        endif
        if(nint(d(182)).gt.0) then
          call int1dn(lint,sg)
        else
          call int1d (lint,sg)
        endif

c       Get global dof for thermal variable

        tdof = max(1,nint(d(19)))
        hh   = d(127)
        tinf = d(128)

        do l = 1,lint

          call shp1d(sg(1,l),xl,shp,ndm,nel,xsj)
          xsj = xsj*sg(2,l)*d(14)

c         Compute flux

          call thfx1d(d,xl,ul,xx,shp, temp,gradt,flux,dd,ndm,ndf,nel)

c         Save data for tplot

          j       = 4*(l-1)
          tt(j+1) = flux
          tt(j+3) = gradt

c         Compute thermal rate

          tdot = 0.0d0
          do j = 1,nel
            tdot = tdot + shp(2,j)*ul(1,j,4)
          end do ! j

          if(stype.eq.3) then
            xsj = xsj*xx
          endif

          j1 = 1
          do j = 1,nel

            a1 = dd*shp(1,j)*xsj
            a3 = d(4)*d(64)*shp(2,j)*xsj
            a4 = d(127)*shp(2,j)*xsj

c           Compute residual

            p(j1) = p(j1) - a1*gradt
     &                    - a3*(cfac*tdot + lfac*ul(1,j,4))
     &                    - a4*(temp - d(128))
     &                    + d(66)*shp(2,j)*xsj*dm
c           Compute tangent

            a1 = a1*ctan(1)
            if(shflg) then
              a3 = a3*ctan(3)
            elseif(ndfo(tdof).gt.0) then
              a3 = a3*ctan(2)
            else
              a3 = 0.0d0
            endif
            a4 = a4*ctan(1) + a3*cfac

c           Lumped rate terms

            s(j1,j1) = s(j1,j1) + a3*lfac

c           Consistent rate and conductivity terms

            i1 = 1
            do i = 1,nel
              s(i1,j1) = s(i1,j1) + a1*shp(1,i) + a4*shp(2,i)
              i1 = i1 + ndf
            end do ! i
            j1 = j1 + ndf
          end do ! j
        end do ! l

c     Output heat flux

      elseif(isw.eq.4) then

        lint = min(5,nint(d(5)))
        if(lint.eq.0) then
          lint = nel
        endif
        if(nint(d(182)).gt.0) then
          call int1dn(lint,sg)
        else
          call int1d (lint,sg)
        endif

        do l=1,lint

          call shp1d(sg(1,l),xl,shp,ndm,nel,xsj)

c         Compute flux and gradients

          call thfx1d(d,xl,ul, xx,shp,temp,gradt,flux,dd, ndm,ndf,nel)

          a4 = -d(127)*(temp - d(128))

          mct = mct - 1
          if(mct.le.0) then
            write(iow,2002) o,head
            if(d(127).gt.0.0d0) then
              write(iow,2004)
            endif
            if(ior.lt.0 .and. pfr) then
              write(*,2002) o,head
              if(d(127).gt.0.0d0) then
                write(*,2004)
              endif
            endif
            mct = 50
          endif
          write(iow,2003) n,ma,xx,flux,gradt
          if(d(127).gt.0.0d0) then
            write(iow,2005) a4
          endif
          if(ior.lt.0 .and. pfr) then
            write(*,2003) n,ma,xx,flux,gradt
            if(d(127).gt.0.0d0) then
              write(*,2005) a4
            endif
          endif

        end do ! l

c     Compute heat capacity (mass) matrix

      elseif(isw.eq.5) then

        lint = min(5,nint(d(5)))
        if(lint.eq.0) then
          lint = nel
        endif
        if(nint(d(182)).gt.0) then
          call int1dn(lint,sg)
        else
          call int1d (lint,sg)
        endif

        do l=1,lint

          call shp1d(sg(1,l),xl,shp,ndm,nel,xsj)
          xsj = xsj*sg(2,l)*d(14)
          if(stype.eq.3) then
            xx = 0.0d0
            do i = 1,nel
              xx = xx + shp(2,i)*xl(1,i)
            end do ! i
            xsj = xsj*xx
          endif
          j1 = 1
          do j = 1,nel
            shj = d(4)*d(64)*shp(2,j)*xsj

c           Lumped capacity (lmas)

            p(j1) = p(j1) + shj
            i1 = 1

c           Consistent (interpolated ) capacity (mass)

            s(j1,j1) = s(j1,j1) + shj*lfac
            do i = 1,nel
              s(i1,j1) = s(i1,j1) + shj*shp(2,i)*cfac
              i1 = i1 + ndf
            end do ! i
            j1 = j1 + ndf
          end do ! j
        end do ! l

c     Compute surface flux loading (not implemented)

c     elseif(isw.eq.7) then

c     Compute nodal heat flux for output/plots

      elseif(isw.eq.8) then

        call thcn1d(d,xl,ul,shp,p,s,p(nen+1),ndf,ndm,nel)

c     Compute error data for heat flux

      elseif(isw.eq.11) then

        call ther1d(d,xl,ul,shp,s,ndf,ndm,nel,nen)

c     External node check

      elseif(isw.eq.26) then

      endif

c     Formats

2000  format(5x,'F o u r i e r   H e a t   C o n d u c t i o n')

2002  format(a1,20a4//5x,'Element Flux'//'  Elmt Mat 1-Coord  2-Coord'
     &            ,'      1-Flux      2-Flux      1-Grad      2-Grad')
2003  format(i8,i4,1p,4e12.3)

2004  format(28x,' Surf. Conv.')
2005  format(28x,1p,1e12.3)

      end

      subroutine thcn1d(d,xl,ul,shp,dt,st,ser,ndf,ndm,nel)

      implicit  none

      include  'iodata.h'
      include  'cdata.h'
      include  'prstrs.h'
      include  'strnum.h'

      integer   ndf,ndm,nel, j,l,lint
      real*8    xx,xsj,xg,d(*), sg(2,3)
      real*8    dt(*),st(nen,*),ser(*),xl(ndm,*),shp(2,*)
      real*8    temp,gradt,flux,dd,ul(ndf,*)

      save

c     Lumped projection routine

      lint = min(5,nint(d(5)))
      if(lint.eq.0) then
        lint = nel
      endif
      if(nint(d(182)).gt.0) then
        call int1dn(lint,sg)
      else
        call int1d (lint,sg)
      endif

      do l=1,lint

        call shp1d(sg(1,l),xl,shp,ndm,nel,xsj)
        xsj = xsj*sg(2,l)*d(14)

        call thfx1d(d,xl,ul, xx,shp,temp,gradt,flux,dd, ndm,ndf,nel)

        temp = -d(127)*(temp - d(128))

c       Compute lumped projection and assemble stress integrals

        do j = 1,nel
          xg      = xsj*shp(2,j)
          dt(j)   = dt(j) + xg
          st(j,13) = st(j,13) + flux*xg
          st(j,15) = st(j,15) + temp*xg
          ser(j)  = ser(j)  + erav*xg
        end do ! j
      end do ! l

      iste = 15

      end

      subroutine thfx1d(d,xl,ul, xx,shp, temp,gradt,flux,dd,
     &                  ndm,ndf,nel)

c     Compute thermal gradient and flux

      implicit  none

      include  'elcoor.h'

      integer   ndm,ndf,nel, i
      real*8    d(*),xl(ndm,*),ul(ndf,*), shp(2,*)
      real*8    xx, temp,gradt,flux,dd

      save

      temp  = 0.0d0
      gradt = 0.0d0
      xx    = 0.0d0
      do i = 1,nel
        gradt = gradt + shp(1,i)*ul(1,i)
        temp  = temp  + shp(2,i)*ul(1,i)
        xx    = xx    + shp(2,i)*xl(1,i)
      end do ! i

c     Compute thermal flux

      dd   = d(61)

      flux = -dd*gradt

c     Save coordinates

      xref(1) = xx
      do i = 2,3
        xref(i) = 0.0d0
      end do ! i

      end

      subroutine ther1d(d,xl,ul,shp,st,ndf,ndm,nel,nen)

      implicit  none

      include  'adapt1.h'
      include  'adapt2.h'
      include  'errind.h'

      integer   ndf,ndm,nel,nen, i,ii
      real*8    g,xx,xsj,detd, st(nen,*),xl(ndm,*),shp(2,*)
      real*8    d(*),gradt,flux,dd, temp,gradp,fluxp,sg
      real*8    ul(ndf,*),ss(3)

      save

      data      ss/-1.d0, 1.d0, 0.d0/

c     Simple routine

      vfem   = 0.d0
      vproj  = 0.d0
      verror = 0.d0
      vener  = 0.d0
      venere = 0.d0
      heta   = 0.0d0
      g      = 1.d0/sqrt(3.0d0)
      do ii = 1,2
        sg = ss(ii)*g
        call shp1d(sg,xl,shp,ndm,nel,xsj)
        call thfx1d(d,xl,ul, xx,shp,temp,gradt,flux,dd, ndm,ndf,nel)
        fluxp = 0.0d0
        do i = 1,nel
          fluxp = fluxp + shp(2,i)*st(i,7)
        end do ! i

c       Compute integral of stress squares for error indicator use

        detd  =  1.d0/dd
        gradp = -detd*fluxp

        heta = heta + xsj
        vfem    = vfem   + flux*flux*xsj
        vproj   = vproj  + fluxp*fluxp*xsj
        verror  = verror + ((fluxp-flux)**2)*xsj
        vener   = vener  + flux*gradt*xsj
        venere  = venere + (fluxp-flux)*(gradp-gradt)*xsj
      end do ! ii
      arsq  =  arsq  + heta
      efem  =  efem  + vfem
      eproj =  eproj + vproj
      eerror=  eerror+ verror
      eener =  eener + vener
      eenere=  eenere+ venere
      areai = heta

c     Check for triangle

      heta  =  d(50)*sqrt(heta)

      end
