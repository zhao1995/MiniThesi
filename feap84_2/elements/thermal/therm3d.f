c$Id:$
      subroutine therm3d(d,ul,xl,ix,s,p,ndf,ndm,nst,isw)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Add interpolated capacity for isw = 5            15/08/2007
c       2. Change 'ord' to 'nel' on 'tetshp' calls          05/11/2007
c       3. Remove 'nel' from call to 'quadr3d'              23/01/2009
c       4. Save xref in 'elcoor.h'                          26/07/2009
c       5. Add history variable treatment and RVE model     19/11/2009
c       6. Add isw.eq.14 for history initialization         28/11/2009
c       7. Add 'pfeapr.h', compute el_vol, etc. use rhoc    30/11/2009
c       8. Output flux at each quadrature point             08/05/2012
c       9. Uset 'oelmt.h' instead of 'pfeapr.h'             09/04/2012
c      10. Remove unused sg0 and sv0 arrays                 11/05/2012
c-----[--.----+----.----+----.-----------------------------------------]
c     Three dimensional Linear Thermal Element

c-----[--.----+----.----+----.-----------------------------------------]

c        This is a three dimensional element which can analyze
c        general geometries.  Set control parameters as
c        follows:

c           ndm - set to 3      (x,y or r,z-coords)
c           ndf - set > or =  1 (nodal temperatures)
c           nel - set > or =  8 for  8-node linear brick
c                     > or = 27 for 27-node quadratic brick
c                     > or = 64 for 64-node cubic brick
c                     > or =  4 for  4-node linear tetrahedron
c                     > or = 10 for 10-node qudratic tetrahedron

c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'bdata.h'
      include  'cdata.h'
      include  'eldata.h'
      include  'elplot.h'
      include  'eltran.h'
      include  'fdata.h'
      include  'hdata.h'
      include  'iofile.h'
      include  'mdata.h'
      include  'oelmt.h'
      include  'part0.h'
      include  'pmod2d.h'
      include  'qudshp.h'
      include  'rdata.h'
      include  'strnum.h'
      include  'comblk.h'

      include  'setups.h'

      integer   ndf,ndm,nst,isw, i,j, i1,j1, l, nn,nhv, tdof, ix(*)
      real*8    a0,a1,a2,a3,shj,tdot,lfac,cfac, rhoc
      real*8    d(*),ul(ndf,nen,*),xl(ndm,*),s(nst,*),p(*)
      real*8    xx(3),gradt(3),flux(3),dd(3,3)

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

        pstyp = 3

        istv  = max(istv,15)

c     Check of mesh if desired (chec)

      elseif(isw.eq.2) then

        call ckbrk8(n,ix,xl,ndm,nel,shp3)

c     Compute conductivity (stiffness) matrix

      elseif(isw.eq.3 .or. isw.eq.6) then

        call quadr3d(d,.true.)

c       Check for thermal dof

        tdof = max(1,nint(d(19)))

        nhv = nint(d(15))
        nn  = 0
        do l = 1,lint

          call interp3d(l, xl, ndm,nel)

c         Set the local element constants for multi-scale

          el_vol = jac(l)
          el_rho = jac(l)*d(4)
          el_c   = jac(l)*d(64)

c         Compute flux

          call thfx3d(d,xl,ul, xx,shp3(1,1,l),gradt,rhoc,
     &                hr(nh1+nn),hr(nh2+nn),nhv,flux,dd,ndm,ndf,nel,isw)

c         Save data for tplot

          j       = 6*(l-1)
          tt(j+1) = flux(1)
          tt(j+2) = flux(2)
          tt(j+3) = flux(3)
          tt(j+4) = gradt(1)
          tt(j+5) = gradt(2)
          tt(j+6) = gradt(3)

c         Compute thermal rate

          tdot = 0.0d0
          do j = 1,nel
            tdot = tdot + shp3(4,j,l)*ul(1,j,4)
          end do ! j

          j1 = 1
          do j = 1,nel

            a1 = (dd(1,1)*shp3(1,j,l) + dd(1,2)*shp3(2,j,l)
     &         +  dd(1,3)*shp3(3,j,l))*jac(l)
            a2 = (dd(2,1)*shp3(1,j,l) + dd(2,2)*shp3(2,j,l)
     &         +  dd(2,3)*shp3(3,j,l))*jac(l)
            a3 = (dd(3,1)*shp3(1,j,l) + dd(3,2)*shp3(2,j,l)
     &         +  dd(3,3)*shp3(3,j,l))*jac(l)

            a0 = rhoc*shp3(4,j,l)*jac(l)

c           Compute residual

            p(j1) = p(j1) - a1*gradt(1) - a2*gradt(2) - a3*gradt(3)
     &                    - a0*(cfac*tdot + lfac*ul(1,j,4))
     &                    + d(66)*shp3(4,j,l)*jac(l)*dm

c           Compute tangent

            if(shflg) then
              a0 = a0*ctan(3)
            elseif(ndfo(tdof).gt.0) then
              a0 = a0*ctan(2)
            else
              a0 = 0.0d0
            endif
            a1 = a1*ctan(1)
            a2 = a2*ctan(1)
            a3 = a3*ctan(1)

c           Lumped rate terms

            s(j1,j1) = s(j1,j1) + a0*lfac

            i1 = 1
            do i = 1,nel
              s(i1,j1) = s(i1,j1) + a1*shp3(1,i,l) + a2*shp3(2,i,l)
     &                            + a3*shp3(3,i,l) + a0*shp3(4,i,l)*cfac
              i1 = i1 + ndf
            end do ! i
            j1 = j1 + ndf
          end do ! j

          nn = nn + nhv

        end do ! l

c     Output heat flux

      elseif(isw.eq.4) then

        call quadr3d(d,.true.)

c       Check for thermal dof

        tdof = max(1,nint(d(19)))

        nhv = nint(d(15))
        nn  = 0
        do l = 1,lint

          call interp3d(l, xl, ndm,nel)

c         Compute flux

          call thfx3d(d,xl,ul, xx,shp3(1,1,l),gradt,rhoc,
     &                hr(nh1+nn),hr(nh2+nn),nhv,flux,dd,ndm,ndf,nel,isw)
          nn = nn + nhv

          mct = mct - 1
          if(mct.le.0) then
            write(iow,2002) o,head
            if(ior.lt.0 .and. pfr) then
              write(*,2002) o,head
            endif
            mct = 50
          endif
          write(iow,2003) n,ma,xx,flux,gradt
          if(ior.lt.0 .and. pfr) then
            write(*,2003) n,ma,xx,flux,gradt
          endif

        end do ! l

c     Compute heat capacity (mass) matrix

      elseif(isw.eq.5) then

        call quadr3d(d,.false.)

        do l = 1,lint
          call interp3d(l, xl, ndm,nel)
          j1 = 1
          do j = 1,nel
            shj = d(4)*d(64)*shp3(4,j,l)*jac(l)

c           Lumped capacity (lmas)

            p(j1) = p(j1) + shj
            i1 = 1

c           Consistent (interpolated) capacity (cmas)

            s(j1,j1) = s(j1,j1) + shj*lfac
            do i = 1,nel
              s(i1,j1) = s(i1,j1) + shj*shp3(4,i,l)*cfac
              i1 = i1 + ndf
            end do ! i
            j1 = j1 + ndf
          end do ! j
        end do ! l

c     Compute surface flux loading (not implemented)

c     elseif(isw.eq.7) then

c     Compute nodal heat flux for output/plots

      elseif(isw.eq.8) then

        call thcn3d(d,xl,ul,p,s,p(nen+1),ndf,ndm,nel,isw)

c     Compute error data for heat flux

      elseif(isw.eq.11) then

        call ther3d(d,xl,ul,s,ndf,ndm,nel,nen,isw)

c     Initialize history data for thermal material

      elseif(isw.eq.14) then

        call quadr3d(d,.true.)

c       Check for thermal dof

        tdof = max(1,nint(d(19)))

        nhv = nint(d(15))
        nn  = 0
        do l = 1,lint

          call interp3d(l, xl, ndm,nel)

c         Compute flux

          call thfx3d(d,xl,ul, xx,shp3(1,1,l),gradt,rhoc,
     &                hr(nh1+nn),hr(nh2+nn),nhv,flux,dd,ndm,ndf,nel,isw)
          nn = nn + nhv

        end do ! l

c     External node determination

      elseif(isw.eq.26) then

        call pcorner3d()

      endif

c     Formats

2000  format(5x,'F o u r i e r   H e a t   C o n d u c t i o n')

2002  format(a1,20a4//5x,'Element Fluxes'//
     &      ' Elmt Matl  1-coord  2-coord  3-coord',
     &      '      1-flux      2-flux      3-flux'/
     &  37x,'      1-grad      2-grad      3-grad')

2003  format(2i5,0p,3f9.3,1p,3e12.3/37x,1p,3e12.3)

      end

      subroutine thcn3d(d,xl,ul,dt,st,ser,ndf,ndm,nel,isw)

      implicit  none

      include  'cdata.h'
      include  'hdata.h'
      include  'iodata.h'
      include  'prstrs.h'
      include  'qudshp.h'
      include  'strnum.h'
      include  'comblk.h'

      integer   ndf,ndm,nel,isw, j,l, nn,nhv
      real*8    xg, rhoc
      real*8    d(*),dt(*),st(nen,*),ser(*),xl(ndm,*)
      real*8    xx(3),gradt(3),flux(3),dd(3,3),ul(ndf,*)

      save

c     Lumped projection routine

      call quadr3d(d,.false.)

      nhv = nint(d(15))
      nn  = 0
      do l = 1,lint
        call interp3d(l, xl, ndm,nel)

        call thfx3d(d,xl,ul, xx,shp3(1,1,l),gradt,rhoc,
     &              hr(nh1+nn),hr(nh2+nn),nhv,flux,dd,ndm,ndf,nel,isw)

c       Compute lumped projection and assemble stress integrals

        do j = 1,nel
          xg      = jac(l)*shp3(4,j,l)
          dt(j)   = dt(j) + xg
          st(j,13) = st(j,13) + flux(1)*xg
          st(j,14) = st(j,14) + flux(2)*xg
          st(j,15) = st(j,15) + flux(3)*xg
          ser(j)  = ser(j)  + erav*xg
        end do ! j
        nn = nn + nhv
      end do ! l

      iste = 15

      end

      subroutine thfx3d(d,xl,ul, xx,shp,gradt,rhoc,hn,hn1,nh,flux,dd,
     &                  ndm,ndf,nel,isw)

c     Compute thermal gradient and flux

      implicit  none

      include  'elcoor.h'

      integer   nh, ndm,ndf,nel,isw, i
      real*8    ta, rhoc
      real*8    d(*),xl(ndm,*),ul(ndf,*), shp(4,*), hn(*),hn1(*)
      real*8    xx(3),gradt(3),flux(3),dd(3,3)

      save

c     Compute temperature and thermal gradient

      ta       = 0.0d0
      gradt(1) = 0.0d0
      gradt(2) = 0.0d0
      gradt(3) = 0.0d0
      xx(1)    = 0.0d0
      xx(2)    = 0.0d0
      xx(3)    = 0.0d0
      do i = 1,nel
        ta       = ta       + shp(4,i)*ul(1,i)
        gradt(1) = gradt(1) + shp(1,i)*ul(1,i)
        gradt(2) = gradt(2) + shp(2,i)*ul(1,i)
        gradt(3) = gradt(3) + shp(3,i)*ul(1,i)
        xx(1)    = xx(1)    + shp(4,i)*xl(1,i)
        xx(2)    = xx(2)    + shp(4,i)*xl(2,i)
        xx(3)    = xx(3)    + shp(4,i)*xl(3,i)
      end do ! i

c     Compute thermal flux and conductivity

      call modltd(d, ta,gradt, hn,hn1,nh, dd,flux,rhoc, isw)

c     Save coordinates

      do i = 1,3
        xref(i) = xx(i)
      end do ! i

      end

      subroutine ther3d(d,xl,ul,st,ndf,ndm,nel,nen,isw)

      implicit  none

      include  'adapt1.h'
      include  'adapt2.h'
      include  'errind.h'
      include  'hdata.h'
      include  'qudshp.h'
      include  'comblk.h'

      integer   ndf,ndm,nel,nen,isw, i,j,l, nn,nhv
      real*8    st(nen,*),xl(ndm,*)
      real*8    d(*),gradt(3),flux(3),dd(3,3),ul(ndf,*), rhoc
      real*8    xx(3),gradp(3),fluxp(3)

      save

      vfem   = 0.d0
      vproj  = 0.d0
      verror = 0.d0
      vener  = 0.d0
      venere = 0.d0
      heta   = 0.0d0
      call quadr3d(d,.false.)

      nhv = nint(d(15))
      nn  = 0
      do l = 1,lint
        call interp3d(l, xl, ndm,nel)

        call thfx3d(d,xl,ul, xx,shp3(1,1,l),gradt,rhoc,
     &              hr(nh1+nn),hr(nh2+nn),nhv,flux,dd,ndm,ndf,nel,isw)
        do i = 1,3
          fluxp(i) = 0.0d0
        end do ! i
        do i = 1,nel
          do j = 1,3
            fluxp(j) = fluxp(j) + shp3(4,i,l)*st(i,j+6)
          end do ! j
        end do ! i

c       Compute integral of stress squares for error indicator use

        call invert(dd,3,3)

        gradp(1) = -dd(1,1)*fluxp(1) - dd(1,2)*fluxp(2)
        gradp(2) = -dd(2,1)*fluxp(1) - dd(2,2)*fluxp(2)

        heta = heta + jac(l)
        do i = 1,3
         vfem    = vfem   + flux(i)*flux(i)*jac(l)
         vproj   = vproj  + fluxp(i)*fluxp(i)*jac(l)
         verror  = verror + ((fluxp(i)-flux(i))**2)*jac(l)
         vener   = vener  + flux(i)*gradt(i)*jac(l)
         venere  = venere + (fluxp(i)-flux(i))
     &                    * (gradp(i)-gradt(i))*jac(l)
        end do ! i

        nn = nn + nhv

      end do ! l
      arsq  =  arsq  + heta
      efem  =  efem  + vfem
      eproj =  eproj + vproj
      eerror=  eerror+ verror
      eener =  eener + vener
      eenere=  eenere+ venere
      areai = heta

      heta  =  d(50)*sqrt(heta)

      end
