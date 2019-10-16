c$Id:$
      subroutine therm2d(d,ul,xl,ix,s,p,ndf,ndm,nst,isw)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Correct lumped s(j1,j1) (was i1) for isw=5       15/08/2007
c       2. Remove 'nel' from call to 'quadr2d'              23/01/2009
c       3. Correct arguments on shp2 for call to thfx2d     13/02/2009
c       4. Save xref in 'elcoor.h'                          26/07/2009
c       5. Add history variable treatment and RVE model     20/11/2009
c       6. Add isw.eq.14 for history initialization         28/11/2009
c       7. Add 'pfeapr.h', compute el_vol, etc. use rhoc    30/11/2009
c       8. Output only 2 flux and gradt not all             06/09/2010
c       9. Uset 'oelmt.h' instead of 'pfeapr.h'             09/04/2012
c-----[--.----+----.----+----.-----------------------------------------]
c     Two dimensional (plane/axisymmetric) Linear Thermal Element

c     N.B.  Surface flux loading may be specified using: convec2d.f

c-----[--.----+----.----+----.-----------------------------------------]
c        This is a two dimensional element which can analyze plane
c        or axisymmetric geometries.  Set control parameters as
c        follows:

c           ndm - set to 2     (x,y or r,z-coords)
c           ndf - set > or = 1 (nodal temperatures)
c           nel - set > or = 4

c                    A eta
c             4      |      3
c              o-----|-----o
c              |     |     |
c              |     |     |
c              |     +-----|----> xi
c              |           |
c              |           |
c              o-----------o
c             1             2

c               Node numbering
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'bdata.h'
      include  'cdata.h'
      include  'complx.h'
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

      integer   ndf,ndm,nst,isw, i,j, i1,j1, l, nn,nhv, tdof, ix(*)
      real*8    xx,yy, a1,a2,a3,a4,shj,tdot,cfac,lfac, hh,tinf
      real*8    d(*),ul(ndf,nen,*),xl(ndm,*),s(nst,*),p(*)
      real*8    temp,rhoc, gradt(3),flux(3),dd(3,3)

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

        pstyp = 2
        istv  = max(istv,15)

c     Check of mesh if desired (chec)

      elseif(isw.eq.2) then

        if(nel.eq.3 .or. nel.eq.6 .or. nel.eq.7) then
          call cktris(ix,xl,shp2,ndm)
        else
          call ckisop(ix,xl,shp2,ndm)
        endif

c     Compute conductivity (stiffness) matrix

      elseif(isw.eq.3 .or. isw.eq.6) then

        call quadr2d(d,.true.)

c       Set the local element constants for multi-scale

        el_vol = jac(l)
        el_rho = jac(l)*d(4)
        el_c   = jac(l)*d(64)

c       Get global dof for thermal variable

        tdof = max(1,nint(d(19)))
        hh   = d(127)
        tinf = d(128)

        nhv  = nint(d(15))
        nn   = 0
        do l = 1,lint

          call interp2d(l, xl,ix, ndm,nel, .false.)
          jac(l) = jac(l)*d(14)

c         Compute flux

          call thfx2d(d,xl,ul,xx,yy,shp2(1,1,l), temp,gradt,rhoc,
     &                hr(nh1+nn),hr(nh2+nn),nhv,flux,dd,ndm,ndf,nel,isw)

c         Save data for tplot

          j       = 4*(l-1)
          tt(j+1) = flux(1)
          tt(j+2) = flux(2)
          tt(j+3) = gradt(1)
          tt(j+4) = gradt(2)

c         Compute thermal rate

          tdot = 0.0d0
          do j = 1,nel
            tdot = tdot + shp2(3,j,l)*ul(1,j,4)
          end do ! j

          if(stype.eq.3) then
            jac(l) = jac(l)*xx
          endif

          j1 = 1
          do j = 1,nel

            a1 = (dd(1,1)*shp2(1,j,l) + dd(1,2)*shp2(2,j,l))*jac(l)
            a2 = (dd(2,1)*shp2(1,j,l) + dd(2,2)*shp2(2,j,l))*jac(l)
            a3 = rhoc*shp2(3,j,l)*jac(l)
            a4 = d(127)*shp2(3,j,l)*jac(l)

c           Compute residual

            p(j1) = p(j1) - a1*gradt(1) - a2*gradt(2)
     &                    - a3*(cfac*tdot + lfac*ul(1,j,4))
     &                    - a4*(temp - d(128))
     &                    + d(66)*shp2(3,j,l)*jac(l)*dm

c           Compute tangent

            a1 = a1*ctan(1)
            a2 = a2*ctan(1)
            if(shflg) then
              if(cplxfl) then
                a3 = 0.0d0
              else
                a3 = a3*ctan(3)
              endif
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
              s(i1,j1) = s(i1,j1) + a1*shp2(1,i,l) + a2*shp2(2,i,l)
     &                            + a4*shp2(3,i,l)
              i1 = i1 + ndf
            end do ! i
            j1 = j1 + ndf
          end do ! j

c         Complex solution part

          if(cplxfl) then
            call ctherm2d(d,xl,ul, shp2(1,1,l),jac(l),cfac,lfac,
     &                    s(1,nst+1),p(nst+1),nn,nhv,isw)
          endif

          nn = nn + nhv
        end do ! l

c     Output heat flux

      elseif(isw.eq.4) then

        call quadr2d(d,.true.)

        nhv  = nint(d(15))
        nn   = 0
        do l=1,lint

          call interp2d(l, xl,ix, ndm,nel, .false.)

c         Compute flux and gradients

          call thfx2d(d,xl,ul, xx,yy,shp2(1,1,l),temp,gradt,rhoc,
     &                hr(nh1+nn),hr(nh2+nn),nhv,flux,dd,ndm,ndf,nel,isw)

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
          write(iow,2003) n,ma,xx,yy,(flux(i),i=1,2),(gradt(i),i=1,2)
          if(d(127).gt.0.0d0) then
            write(iow,2005) a4
          endif
          if(ior.lt.0 .and. pfr) then
            write(*,2003) n,ma,xx,yy,(flux(i),i=1,2),(gradt(i),i=1,2)
            if(d(127).gt.0.0d0) then
              write(*,2005) a4
            endif
          endif

          nn = nn + nhv

        end do ! l

c     Compute heat capacity (mass) matrix

      elseif(isw.eq.5) then

        call quadr2d(d,.false.)

        do l = 1,lint

          call interp2d(l, xl,ix, ndm,nel, .false.)
          jac(l) = jac(l)*d(14)

          if(stype.eq.3) then
            xx = 0.0d0
              do i = 1,nel
              xx = xx + shp2(3,i,l)*xl(1,i)
            end do ! i
            jac(l) = jac(l)*xx
          endif
          j1 = 1
          do j = 1,nel
            shj = d(4)*d(64)*shp2(3,j,l)*jac(l)

c           Lumped capacity (lmas)

            p(j1) = p(j1) + shj

c           Consistent (interpolated ) capacity (mass)

            s(j1,j1) = s(j1,j1) + shj*lfac
            i1 = 1
            do i = 1,nel
              s(i1,j1) = s(i1,j1) + shj*shp2(3,i,l)*cfac
              i1 = i1 + ndf
            end do ! i
            j1 = j1 + ndf
          end do ! j
        end do ! l

c     Compute surface flux loading (not implemented)

c     elseif(isw.eq.7) then

c     Compute nodal heat flux for output/plots

      elseif(isw.eq.8) then

        call thcn2d(ix,d,xl,ul,p,s,p(nen+1),ndf,ndm,nel,isw)

c     Compute error data for heat flux

      elseif(isw.eq.11) then

        call ther2d(ix,d,xl,ul,shp2,s,ndf,ndm,nel,nen)

c     Initialize history variables

      elseif(isw.eq.14) then

        call quadr2d(d,.true.)

c       Get global dof for thermal variable

        tdof = max(1,nint(d(19)))
        hh   = d(127)
        tinf = d(128)

        nhv  = nint(d(15))
        nn   = 0
        do l = 1,lint

          call interp2d(l, xl,ix, ndm,nel, .false.)

c         Compute flux

          call thfx2d(d,xl,ul,xx,yy,shp2(1,1,l), temp,gradt,rhoc,
     &                hr(nh1+nn),hr(nh2+nn),nhv,flux,dd,ndm,ndf,nel,isw)
          nn = nn + nhv
        end do ! l

c     External node determination

      elseif(isw.eq.26) then

        call pcorner2d()

      endif

c     Formats

2000  format(5x,'F o u r i e r   H e a t   C o n d u c t i o n')

2002  format(a1,20a4//5x,'Element Flux'//'  Elmt Mat 1-Coord  2-Coord'
     &            ,'      1-Flux      2-Flux      1-Grad      2-Grad')
2003  format(i6,i4,0p,2f9.3,1p,4e12.3)

2004  format(28x,' Surf. Conv.')
2005  format(28x,1p,1e12.3)

      end

      subroutine thcn2d(ix,d,xl,ul,dt,st,ser,ndf,ndm,nel,isw)

      implicit  none

      include  'iodata.h'
      include  'cdata.h'
      include  'hdata.h'
      include  'prstrs.h'
      include  'qudshp.h'
      include  'strnum.h'
      include  'comblk.h'

      integer   ndf,ndm,nel, j,l, isw, nn,nhv,ne8, ix(*), ixl(16)
      real*8    xx,yy,xg,rhoc, d(*)
      real*8    dt(*),st(nen,*),ser(*),xl(ndm,*)
      real*8    temp,gradt(3),flux(3),dd(3,3),ul(ndf,*)

      save

      data      ixl/ 1,2,3,4,5,6,7,8,9, 7*0 /

c     Lumped projection routine

      call quadr2d(d,.false.)
      if(quad .and. nel.eq.8) then
        call meanx(ix,xl,ndm)
        ne8 = 9
      else
        ne8 = nel
        do j = 1,nel
          ixl(j) = ix(j)
        end do
      endif

      nhv  = nint(d(15))
      nn   = 0
      do l = 1,lint
        call interp2d(l, xl,ixl, ndm,ne8, .false.)
        jac(l) = jac(l)*d(14)

        call thfx2d(d,xl,ul, xx,yy,shp2(1,1,l),temp,gradt,rhoc,
     &              hr(nh1+nn),hr(nh2+nn),nhv,flux,dd,ndm,ndf,nel,isw)

        temp = -d(127)*(temp - d(128))

c       Compute lumped projection and assemble stress integrals

        do j = 1,nel
          xg      = jac(l)*shp2(3,j,l)
          dt(j)   = dt(j)   + xg
          st(j,13) = st(j,13) + flux(1)*xg
          st(j,14) = st(j,14) + flux(2)*xg
          st(j,15) = st(j,15) + temp   *xg
          ser(j)  = ser(j)  + erav   *xg
        end do ! j
        nn = nn + nhv
      end do ! l

      iste = 15

      end

      subroutine thfx2d(d,xl,ul, xx,yy,shp, temp,gradt,rhoc, hn,hn1,nhv,
     &                  flux,dd,ndm,ndf,nel,isw)

c     Compute thermal gradient and flux

      implicit  none

      include  'elcoor.h'

      integer   ndm,ndf,nel,nhv,isw, i
      real*8    d(*),xl(ndm,*),ul(ndf,*), shp(3,*), hn(*),hn1(*)
      real*8    xx,yy, temp,rhoc, gradt(3),flux(3),dd(3,3)

      save

      temp     = 0.0d0
      gradt(1) = 0.0d0
      gradt(2) = 0.0d0
      gradt(3) = 0.0d0
      xx       = 0.0d0
      yy       = 0.0d0
      do i = 1,nel
        gradt(1) = gradt(1) + shp(1,i)*ul(1,i)
        gradt(2) = gradt(2) + shp(2,i)*ul(1,i)
        temp     = temp     + shp(3,i)*ul(1,i)
        xx       = xx       + shp(3,i)*xl(1,i)
        yy       = yy       + shp(3,i)*xl(2,i)
      end do ! i

c     Compute thermal flux and conductivity

      call modltd(d, temp,gradt, hn,hn1, nhv, dd,flux,rhoc, isw)

c     Save coordinates

      xref(1) = xx
      xref(2) = yy
      xref(3) = 0.0d0

      end

      subroutine ther2d(ix,d,xl,ul,shp,st,ndf,ndm,nel,nen)

      implicit  none

      include  'adapt1.h'
      include  'adapt2.h'
      include  'errind.h'
      include  'hdata.h'
      include  'comblk.h'

      integer   ndf,ndm,nel,nen, i,j,ii, nn,nhv,ix(*)
      real*8    g,xx,yy,xsj,detd,rhoc, st(nen,*),xl(ndm,*),shp(3,*)
      real*8    d(*),gradt(3),flux(3),dd(3,3),ul(ndf,*),ss(9),tt(9)
      real*8    temp,gradp(2),fluxp(2),ddp(2,2),sg(2)

      save

      data      ss/-1.d0, 1.d0,1.d0,-1.d0, 0.d0,1.d0,0.d0,-1.d0,0.d0/
      data      tt/-1.d0,-1.d0,1.d0, 1.d0,-1.d0,0.d0,1.d0, 0.d0,0.d0/

c     Simple routine

      vfem   = 0.d0
      vproj  = 0.d0
      verror = 0.d0
      vener  = 0.d0
      venere = 0.d0
      heta   = 0.0d0
      g      = 1.d0/sqrt(3.0d0)
      nhv    = nint(d(15))
      nn     = 0
      do ii = 1,4
        sg(1) = ss(ii)*g
        sg(2) = tt(ii)*g
        call shp2d(sg,xl,shp,xsj,ndm,nel,ix,.false.)
        call thfx2d(d,xl,ul, xx,yy,shp,temp,gradt,rhoc,
     &              hr(nh1+nn),hr(nh2+nn),nhv,flux,dd, ndm,ndf,nel,11)
        do i = 1,2
          fluxp(i) = 0.0d0
        end do ! i
        do i = 1,nel
          do j = 1,2
            fluxp(j) = fluxp(j) + shp(3,i)*st(i,j+6)
          end do ! j
        end do ! i

c       Compute integral of stress squares for error indicator use

        detd     =  1.d0/(dd(1,1)*dd(2,2) - dd(1,2)*dd(2,1))
        ddp(1,1) =  dd(2,2)*detd
        ddp(1,2) = -dd(1,2)*detd
        ddp(2,1) = -dd(2,1)*detd
        ddp(2,2) =  dd(1,1)*detd
        gradp(1) = -ddp(1,1)*fluxp(1) - ddp(1,2)*fluxp(2)
        gradp(2) = -ddp(2,1)*fluxp(1) - ddp(2,2)*fluxp(2)

        heta = heta + xsj
        do i = 1,2
         vfem    = vfem   + flux(i)*flux(i)*xsj
         vproj   = vproj  + fluxp(i)*fluxp(i)*xsj
         verror  = verror + ((fluxp(i)-flux(i))**2)*xsj
         vener   = vener  + flux(i)*gradt(i)*xsj
         venere  = venere + (fluxp(i)-flux(i))*(gradp(i)-gradt(i))*xsj
        end do ! i

        nn = nn + nhv

      end do ! ii
      arsq  =  arsq  + heta
      efem  =  efem  + vfem
      eproj =  eproj + vproj
      eerror=  eerror+ verror
      eener =  eener + vener
      eenere=  eenere+ venere
      areai = heta

c     Check for triangle

      if(nel.lt.4 .or. ix(1).eq.ix(2) .or. ix(2).eq.ix(3)
     1            .or. ix(3).eq.ix(4) .or. ix(4).eq.ix(1) ) then
        heta = heta*2.d0
      endif
      heta  =  d(50)*sqrt(heta)

      end

      subroutine ctherm2d(d,xl,ul, shp,xsj,cfac,lfac, s,p,nn,nhv,isw)

      implicit   none

      include   'cdata.h'
      include   'eldata.h'
      include   'hdata.h'
      include   'pmod2d.h'
      include   'rdata.h'
      include   'sdata.h'
      include   'comblk.h'

      real*8     d(*),xl(ndm,*),ul(ndf,nen,*), shp(3,*)
      real*8     s(nst,*),p(*), xsj, cfac,lfac

      integer    i,j, i1,j1, nn,nhv, isw
      real*8     xx,yy,tdot,temp, a1,a2,a3,a4, omega,rhoc
      real*8     gradt(3),flux(3),dd(3,3)

c     Shift frequency

      omega = sqrt(abs(shift))

c     Compute flux

      call thfx2d(d,xl,ul(1,1,8),xx,yy,shp, temp,gradt,rhoc,
     &            hr(nh1+nn),hr(nh2+nn),nhv,flux,dd,ndm,ndf,nel,isw)

c     Compute thermal rate

      tdot = 0.0d0
      do j = 1,nel
        tdot = tdot + shp(3,j)*ul(1,j,11)
      end do ! j

      if(stype.eq.3) then
        xsj = xsj*xx
      endif

      j1 = 1
      do j = 1,nel

        a1 = (dd(1,1)*shp(1,j) + dd(1,2)*shp(2,j))*xsj
        a2 = (dd(2,1)*shp(1,j) + dd(2,2)*shp(2,j))*xsj
        a3 = rhoc*shp(3,j)*xsj*omega
        a4 = d(127)*shp(3,j)*xsj

c       Compute residual

        p(j1) = p(j1) - a1*gradt(1) - a2*gradt(2)
     &                - a3*(cfac*tdot + lfac*ul(1,j,4))
     &                - a4*(temp - d(128))
     &                + d(66)*shp(3,j)*xsj*dm

c       Compute tangent

        a4 = a3*cfac

c       Lumped rate terms

        s(j1,j1) = s(j1,j1) + a3*lfac

c       Consistent rate and conductivity terms

        i1 = 1
        do i = 1,nel
          s(i1,j1) = s(i1,j1) + a4*shp(3,i)
          i1 = i1 + ndf
        end do ! i
        j1 = j1 + ndf
      end do ! j

      end
