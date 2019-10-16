c$Id:$
      subroutine sld2d2(d,ul,xl,ix,tl,s,r,ndf,ndm,nst,isw, ther)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Use constants from 'pconstant.h'                 14/11/2006
c       2. Revise form of augmented statement               26/01/2007
c       3. Add line feed to format 2003                     07/03/2007
c       4. Modify augmentation for higher order elements;   14/03/2007
c          Use d(185) for augmenting factor
c       5. Add radial body loading                          18/05/2007
c       6. Add direct call to quadr2d and interp2d          11/11/2008
c       7. Remove 'nel' from call to 'quadr2d'              23/01/2009
c       8. Increase order of element and quadrature to 6    17/02/2009
c       9. Add principal strains to stress/strain outputs   02/03/2009
c      10. Revise augment: nhi -> fp(1), augfa -> augf      14/04/2009
c      11. Dimension eps(9,3,64) and also other arrays      02/06/2009
c      12. Move computation of psil to material models      26/07/2009
c      13. For augment set npm by call to quadr2d           15/08/2009
c      14. Remove restriction on plane stress               29/11/2010
c      15. Compute v_avg and sig_33 for multiscale use      03/12/2010
c      16. Add 'l' to modlsd call                           05/01/2012
c      17. Add average of density for multiscale            10/05/2012
c      18. Introduce jvol(64) for computations so stress    03/10/2012
c          projected on areas of elements for axisymmetry.
c      19. Add eps on call to slcn2d                        01/01/2013
c      20. Pass strains to stcn2z for z-zhu projections     01/01/2014
c-----[--.----+----.----+----.-----------------------------------------]
c     Purpose: Two-dimensional mixed u-p-theta small deformation element

c      Inputs:
c         d(*)  - Element parameters
c         ul(ndf,*) - Current nodal solution parameters
c         xl(ndm,*) - Nodal coordinates
c         ix(*)     - Global nodal connections
c         tl(*)     - Nodal temp vector
c         ndf       - Degree of freedoms/node
c         ndm       - Mesh coordinate dimension
c         nst       - Element array dimension
c         isw       - Solution option switch

c      Outputs:
c         s(nst,*)  - Element array
c         r(ndf,*)  - Element vector
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'augdat.h'
      include  'bdata.h'
      include  'cdata.h'
      include  'elbody.h'
      include  'eldata.h'
      include  'elengy.h'
      include  'elplot.h'
      include  'eltran.h'
      include  'fdata.h'
      include  'hdata.h'
      include  'iofile.h'
      include  'oelmt.h'
      include  'pconstant.h'
      include  'part0.h'
      include  'pmod2d.h'
      include  'ptdat6.h'
      include  'p_int.h'
      include  'qudshp.h'
      include  'rdata.h'
      include  'comblk.h'

      logical   ther
      integer   ndf,ndm,nst,isw,i,i1,j,jj,j1,l
      integer   ncp,ncs,nhv,istrt,nn, ix(*)
      real*8    alam,ha, epp,dv,dl,d1,type,bt(3), aj0,aj1
      real*8    dsigtr, mpress, dmass, dmshp, dtheta, cfac,lfac, fac
      real*8    d(*),       ul(ndf,nen,*),  xl(ndm,*),   tl(*), s(nst,*)
      real*8    r(ndf,*),   xx(2,64)
      real*8    bbd(3,7),   aa(6,6,5,64),   dd(7,7)
      real*8    sigm(9),    sigl(16,64),    bpra(3)  ,   bbar(2,64,64)
      real*8    al(3),      ac(3),          vl(3)    ,   x0(2)
      real*8    phi(6,64),  theta(3,64),    hh(6,6)
      real*8    press(64),  pbar(64),       hsig(6),     eps(9,3,64)
      real*8    irad(64),   ta(64),         epsd(6),     epsv(64)
      real*8    dther(7),   body(3),        mom(3,64),   rr(6)
      real*8    bdy(3),     xr(2),          epsm(9)  ,   jvol(64)
      real*8    epsp(6,64)

      save

c     TEMPORARY SET OF TEMPERATURE

      data    ta      / 64*0.0d0 /
      data    alam,ha /  2*0.0d0 /

c     Data inputs: Set controls and errors

      if( isw.eq.1 ) then

c       Adjust history variable storage for element types

        nhv = nint(d(15))
        if(nen.le.4) then
          nh1 = nh1 + 2
        elseif(nen.le.6) then
          nh1 = nh1 + 2 - nhv*2
        elseif(nen.le.7) then
          nh1 = nh1 + 6 - nhv*2
        elseif(nen.le.9) then
          nh1 = nh1 + 6
        else
          nh1 = nh1 + 12
        endif

c       Write error if plane stress requested

c       if(stype.eq.1) then
c         write(  *,4001)
c         write(ilg,4001)
c         write(iow,4001)
c         call plstop()
c       endif

c     Augmented Lagrangian update for nested iteration

      elseif(isw.eq.10) then

        call quadr2d(d,.true.)

        d1    = augf*d(185)
        fp(1) = nh2 - 1
        fp(2) = fp(1) + npm
        do i = 1,npm
          hr(fp(1)+i) = hr(fp(1)+i) + d1*hr(fp(2)+i)
        end do ! i

c     Compute tangent stiffness and residual force vector

      elseif(isw.eq. 3 .or. isw.eq. 4 .or. isw.eq. 6 .or.
     &       isw.eq. 8 .or. isw.eq.13 .or. isw.eq.14 .or.
     &       isw.eq.16 .or. isw.eq.25) then

c       Integration order set to static

        if((d(7).ge.0.0 .or. d(183).ne.0.0d0) .and.
     &             (ndfo(1).gt.0 .or. shflg)) then
          cfac = d(7)
          lfac = 1.d0 - cfac
          if(nel.eq.3) then
            call masst3(stype,cfac,d(4),xl,ul(1,1,5),r,s)
            cfac = 0.0d0
          endif
        else
          cfac = 0.0d0
          lfac = 0.0d0
        endif

c       Initialize momenta

        if(isw.eq.13) then
          do j = 1,nel
            do i = 1,3
              mom(i,j) = 0.0d0
            end do ! i
          end do ! j
        endif

c       Proportional body forces

        call sbodyf(d, body)
        do i = 1,3
          bdy(i) = body(i)
        end do ! i

        estore = 0.0d0

c       Set element quadrature order

        call quadr2d(d,.true.)

c       Set number of history terms / quadradure point

        type   = max(0,stype - 2)
        nhv    = nint(d(15))
        istrt  = nint(d(84))

c       Center estimate

        if(npm.gt.1) then
          do i = 1,2
            x0(i) = 0.0d0
            do j = 1,nel
              x0(i) = x0(i) + xl(i,j)
            end do ! j
            x0(i) = x0(i)/dble(nel)
          end do ! i
        endif

c       MECHANICAL ELEMENT

        do l = 1,lint

c         Shape functions and derivatives

          call interp2d(l, xl,ix, ndm,nel, .false.)

c         Mixed volume effect and temperature projection

          do i = 1,3
            theta(i,l) = 0.0d0
          end do ! i

          ta(l) = -d(9)

c         Compute coordinates

          xx(1,l) = 0.0d0
          xx(2,l) = 0.0d0
          do i = 1,nel
            xx(1,l) = xx(1,l) + shp2(3,i,l)*xl(1,i)
            xx(2,l) = xx(2,l) + shp2(3,i,l)*xl(2,i)
          end do ! i

c         Compute volumetric strain from displacements

          jvol(l) = jac(l)
          if(stype.lt.3) then
            irad(l) = 0.0d0
            ncp     = 2
            ncs     = 4
          else
            jvol(l) = jvol(l)*xx(1,l)
            irad(l) = 1.d0/xx(1,l)
            ncp     = 2
            ncs     = 4
            if(stype.eq.8) then
              ncp   = 3
              ncs   = 6
            endif
          endif

c         Compute div u (volumetric strain) from displacements

          do i = 1,nel
            fac        = shp2(1,i,l) + shp2(3,i,l) * irad(l)
            theta(1,l) = theta(1,l)  + fac         * ul(1,i,1)
     &                               + shp2(2,i,l) * ul(2,i,1)
            theta(2,l) = theta(2,l)  + fac         * ul(1,i,2)
     &                               + shp2(2,i,l) * ul(2,i,2)
            theta(3,l) = theta(3,l)  + fac         * ul(1,i,3)
     &                               + shp2(2,i,l) * ul(2,i,3)
            ta(l)      = ta(l)       + shp2(3,i,l) * tl(i)
          end do ! i

c         Set the pressure functions

          phi(1,l) = 1.d0
          if(npm.gt.1) then
            phi(2,l) = xx(1,l) - x0(1)
            phi(3,l) = xx(2,l) - x0(2)
            phi(4,l) = phi(2,l)**2
            phi(5,l) = phi(2,l)*phi(3,l)
            phi(6,l) = phi(3,l)**2
          endif
        end do ! l

c       Mixed model for volumetric and temperature response

        call bbar2s(phi,shp2,jvol,lint,nel,npm,hh,irad,theta,bbar)

c       Initialize storage for augmented function

        if(npm.gt.1) then
          do i = 1,npm
            rr(i) = 0.0d0
          end do ! i
        endif

c       Compute strains and stresses at quadrature points

        nn = 2*npm
        do l = 1,lint
          call strn2m(shp2(1,1,l),xl,ul,theta(1,l),irad(l),
     &                ndm,ndf,nel,nen,eps(1,1,l))

          epsv(l) = theta(1,l)

c         Set augmented multipler

          alam = hr(nh2)
          do i = 1,npm-1
            alam = alam + phi(i+1,l)*hr(nh2+i)
          end do ! i

          call modlsd(l,d,ta(l),eps(1,1,l),hr(nn+nh1),hr(nn+nh2),nhv,
     &                istrt,aa(1,1,1,l),sigl(1,l),alam,ha,isw)
          if(npm.gt.1) then
            do i = 1,npm
              rr(i) = rr(i) + phi(i,l)*ha*jvol(l)
            end do ! i
          endif

c         Volumetric stress

          pbar(l) = one3*(sigl(1,l) + sigl(2,l) + sigl(3,l))

c         Compute energy

          if(isw.eq.13) then
            aj0 = d(14)*jvol(l)
            aj1 = d( 4)*aj0
            call sengy(ul,shp2(1,1,l), aj0,aj1, lfac,cfac, ncp,3, mom)
          endif

          nn = nn + nhv
        end do ! l

c       Accumulate energy

        if(isw.eq.13) then
          do j = 1,nel
            do i = 1,ncp
              epl(i) = epl(i) + mom(i,j)
            end do ! i
            epl(4) = epl(4) + xl(1,j)*mom(3,j)
            epl(5) = epl(5) - xl(2,j)*mom(3,j)
            epl(6) = epl(6) + xl(1,j)*mom(2,j) - xl(2,j)*mom(1,j)
          end do ! j
          return
        endif

c       Integrate constant pressure over volume, set augmented values

        if(npm.eq.1) then

          hr(nh2+1) = ha
          mpress    = 0.0d0
          do l = 1,lint
            mpress  = mpress  + pbar(l)*jvol(l)
          end do ! l

c         Divide pressure by volume

          press(1) = mpress * hh(1,1)
          do l = 2,lint
            press(l) = press(1)
          end do ! l

c       Higher order element pressures

        else

          fp(1) = nh2 + npm  - 1
          do i = 1,npm
            hr(fp(1)+i) = 0.0d0
            do j = 1,npm
              hr(fp(1)+i) = hr(fp(1)+i) + hh(i,j)*rr(j)
            end do ! j
          end do ! i

          do i = 1,npm
            sigm(i) = 0.0d0
          end do ! i
          do l = 1,lint
            mpress  = pbar(l)*jvol(l)
            sigm(1) = sigm(1) + mpress
            do i = 2,npm
              sigm(i) = sigm(i) + mpress*phi(i,l)
            end do ! i
          end do ! l

c         Divide pressure by reference volume

          do i = 1,npm
            hsig(i) = 0.0d0
            do j = 1,npm
              hsig(i) = hsig(i) + hh(i,j)*sigm(j)
            end do ! j
          end do ! i

          do l = 1,lint
            press(l) = hsig(1)
            do i = 2,npm
              press(l) = press(l) + hsig(i)*phi(i,l)
            end do ! i
          end do ! l

        endif

c       Compute mixed stress and multiply by volume element

        do l = 1,lint
          dsigtr    =  press(l)  - pbar(l)
          sigl(1,l) =  sigl(1,l) + dsigtr
          sigl(2,l) =  sigl(2,l) + dsigtr
          sigl(3,l) =  sigl(3,l) + dsigtr
        end do ! l

c       Tangent and residual computations

        if(isw.eq.3 .or. isw.eq.6 .or. isw.eq.14) then

c         Compute mixed pressure

          if(isw.eq.3 .or. isw.eq.6) then

            do l = 1,lint

c             Store time history plot data for element

              i = 6*(l-1)
              do j = 1,6
                tt(j+i) = sigl(j,l)
                sigm(j) = sigl(j,l)*jvol(l)
              end do ! j

c             Compute parameters for multiscale plane strain

              v_avg  = v_avg + jvol(l)
              v_rho  = v_rho + jvol(l)*d(4)
              sig_33 = sig_33 + sigm(3)

c             Compute acceleration

              dmass = d(4)*jvol(l)
              do i = 1,ncp
                al(i) = 0.0d0
                do j = 1,nel
                  al(i) = al(i) + shp2(3,j,l)*ul(i,j,5)
                end do ! j
                al(i) = al(i)*cfac
              end do ! i

c             Rayleigh damping

              if(d(77).ne.0.0d0) then
                do i = 1,ncp
                  vl(i) = 0.0d0
                  do j = 1,nel
                    vl(i) = vl(i) + shp2(3,j,l)*ul(i,j,4)
                  end do ! j
                  vl(i)   = vl(i)*cfac
                end do ! i

c               Compute mass damping residual

                do i = 1,nel
                  fac    = shp2(3,i,l)*jvol(l)*d(77)*d(4)
                  do j = 1,ncp
                  r(j,i) = r(j,i) - (vl(j) + lfac*ul(j,i,4))*fac
                  end do ! j
                end do ! i
              endif

              if(d(78).ne.0.0d0) then
                do i = 1,ncs
                  epsd(i) = 0.0d0
                end do ! i
                do j = 1,nel
                  epsd(1) = epsd(1) + shp2(1,j,l)*ul(1,j,4)
                  epsd(2) = epsd(2) + shp2(2,j,l)*ul(2,j,4)
                  epsd(3) = epsd(3) + shp2(3,j,l)*ul(1,j,4)
                  epsd(4) = epsd(4) + shp2(2,j,l)*ul(1,j,4)
     &                              + shp2(1,j,l)*ul(2,j,4)
                end do ! j
                epsd(3) = epsd(3)*irad(l)
                if(stype.eq.8) then
                  do j = 1,nel
                    epsd(5) = epsd(5) +  shp2(2,j,l)*ul(3,j,4)
                    epsd(6) = epsd(6) + (shp2(1,j,l)
     &                                -  shp2(3,j,l)*irad(l))*ul(3,j,4)
                  end do ! j
                endif
                fac = one3*(theta(2,l) - epsd(1) - epsd(2) - epsd(3))
                epsd(1) = (epsd(1) + fac)*d(78)*jvol(l)
                epsd(2) = (epsd(2) + fac)*d(78)*jvol(l)
                epsd(3) = (epsd(3) + fac)*d(78)*jvol(l)
                do j = 4,ncs
                  epsd(j) =  epsd(j)*d(78)*jvol(l)
                end do ! j

                do i = 1,ncs
                  do j = 1,ncs
                    sigm(j) = sigm(j) + aa(j,i,5,l)*epsd(i)
                  end do ! j
                end do ! i
              endif

              if(stype.lt.3 .and. nint(d(69)).eq.5) then
                do i = 1,2
                  xr(i) = 0.0d0
                  do j = 1,nel
                    xr(i) = xr(i) + shp2(3,j,l)*xl(i,j)
                  end do ! j
                end do ! i
                fac = sqrt(xr(1)**2 + xr(2)**2)
                if(fac.gt.0.0d0) then
                  bdy(1) = (xr(1)*body(1) - xr(2)*body(2))/fac
                  bdy(2) = (xr(2)*body(1) + xr(1)*body(2))/fac
                endif
              endif

c             Compute residual

              do j = 1,nel

                do i = 1,ncp
                  ac(i)  = d(4)*(al(i) + lfac*ul(i,j,5))
                end do ! i

                r(1,j) = r(1,j) + (bdy(1) - ac(1))*shp2(3,j,l)*jvol(l)
     &                          - shp2(1,j,l)*sigm(1)
     &                          - shp2(2,j,l)*sigm(4)
     &                          - shp2(3,j,l)*sigm(3)*irad(l)
                r(2,j) = r(2,j) + (bdy(2) - ac(2))*shp2(3,j,l)*jvol(l)
     &                          - shp2(1,j,l)*sigm(4)
     &                          - shp2(2,j,l)*sigm(2)
              end do ! j
              if(stype.eq.8) then
                do j = 1,nel
                  r(3,j) = r(3,j) + (bdy(3) - ac(3))*shp2(3,j,l)*jvol(l)
     &                   -  shp2(2,j,l)*sigm(5)
     &                   - (shp2(1,j,l)-shp2(3,j,l)*irad(l))*sigm(6)
                end do ! j
              endif

c             Compute mixed tangent stiffness matrix

              if(isw.eq.3) then

c               Multiply tangent moduli by volume element

                call dmatmx( aa(1,1,1,l), dd )
                d1 = jvol(l)*ctan(1)
                do i = 1,7
                  do j = 1,7
                    dd(i,j) = dd(i,j)*d1
                  end do ! j
c                 Thermo-mechanical coupling
                  dther(i) = aa(i,1,2,l)*d1
                end do ! i

                if(d(78).ne.0.0d0) then
                  d1 = jvol(l)*d(78)*ctan(2)
                  do i = 1,6
                    do j = 1,6
                      dd(i,j) = dd(i,j) + aa(i,j,5,l)*d1
                    end do ! j
                  end do ! i
                endif

c               Mass factors

                dv = (ctan(3) + d(77)*ctan(2))*jvol(l)*d(4)*cfac
                dl = (ctan(3) + d(77)*ctan(2))*jvol(l)*d(4)*lfac

c               Compute row terms

                i1 = 0
                do i = 1,nel

c                 Compute bmat-t * dd * dvol

                  do jj = 1,7

                    bbd(1,jj) =  shp2(1,i,l)*dd(1,jj)
     &                        +  shp2(2,i,l)*dd(4,jj)
     &                        +  shp2(3,i,l)*dd(3,jj)*irad(l)
     &                        +  bbar(1,i,l)*dd(7,jj)

                    bbd(2,jj) =  shp2(2,i,l)*dd(2,jj)
     &                        +  shp2(1,i,l)*dd(4,jj)
     &                        +  bbar(2,i,l)*dd(7,jj)

                    bbd(3,jj) =  shp2(2,i,l)*dd(5,jj)
     &                        + (shp2(1,i,l) - shp2(3,i,l)
     &                        * irad(l))*dd(6,jj)
                  end do ! jj
c                          _
c                 Compute: B_trans * D * alpha * j * w

                  if(ther) then
                    fac   = (dther(1) + dther(2) + dther(3))*one3
                    bt(1) =  shp2(1,i,l)*dther(1) + shp2(2,i,l)*dther(4)
     &                    +  shp2(3,i,l)*dther(3) * irad(l)
     &                    + (bbar(1,i,l) - shp2(1,i,l)
     &                                   - shp2(3,i,l)*irad(l))*fac

                    bt(2) =  shp2(2,i,l)*dther(2) + shp2(1,i,l)*dther(4)
     &                    + (bbar(2,i,l) - shp2(2,i,l))*fac

                    bt(3) =  shp2(2,i,l)*dther(5)
     &                    + (shp2(1,i,l)-shp2(3,i,l)*irad(l))*dther(6)
                  endif

c                 Compute tangent stiffness

                  fac = shp2(3,i,l)*dl
                  do jj = 1,ncp
                    s(i1+jj,i1+jj) = s(i1+jj,i1+jj) + fac
                  end do ! jj

                  dmshp = shp2(3,i,l)*dv

                  j1 = 0
                  do j = 1,nel

c                   Inertial tangent

                    do jj = 1,ncp
                      s(i1+jj,j1+jj) = s(i1+jj,j1+jj)+dmshp*shp2(3,j,l)
                    end do ! jj

c                   Compute mechanics part of tangent stiffness

                    do jj = 1,ncp
                      s(i1+jj,j1+1) = s(i1+jj,j1+1)
     &                              + bbd(jj,1)*shp2(1,j,l)
     &                              + bbd(jj,4)*shp2(2,j,l)
     &                              + bbd(jj,3)*shp2(3,j,l)*irad(l)
     &                              + bbd(jj,7)*bbar(1,j,l)

                      s(i1+jj,j1+2) = s(i1+jj,j1+2)
     &                              + bbd(jj,2)*shp2(2,j,l)
     &                              + bbd(jj,4)*shp2(1,j,l)
     &                              + bbd(jj,7)*bbar(2,j,l)
                    end do ! jj

c                   Torsion terms

                    if(stype.eq.8) then
                      do jj = 1,3
                        s(i1+jj,j1+3) = s(i1+jj,j1+3)
     &                                + bbd(jj,5)*shp2(2,j,l)
     &                                + bbd(jj,6)*(shp2(1,j,l)
     &                                           - shp2(3,j,l)*irad(l))
                      end do ! jj
                    endif

c                   Thermo-mechanical coupling matrix

                    if(ther) then
                      do jj = 1,ncp
                      s(i1+jj,j1+3)  = s(i1+jj,j1+3)+bt(jj)*shp2(3,j,l)
                      end do ! jj
                    endif

                    j1 = j1 + ndf
                  end do ! j
                  i1 = i1 + ndf
                end do ! i
              endif ! isw = 3
            end do ! l

          endif ! isw = 3 or 6

c         Multiply by thickness if not unity

          if((isw.eq.3 .or. isw.eq.6) .and. d(14).ne.1.d0) then

            do j = 1,nst
              do i = 1,nst
                s(i,j) = s(i,j)*d(14)
              end do ! i
            end do ! j
            do j = 1,nel
              do i = 1,ndf
                r(i,j) = r(i,j)*d(14)
              end do ! i
            end do ! j

          endif

c       Output stresses.

        elseif(isw.eq.4 .or. isw.eq.8) then

          do i = 1,9
            sigm(i) = 0.0d0
          end do ! i
          do i = 1,3
            bpra(i) = 0.0d0
          end do ! i
          epp    = 0.0d0
          dtheta = 0.0d0

c         Output stresses

          if (isw .eq. 4) then

            do l = 1,lint
              do i = 1,3
                sigm(i  ) = sigl(i  ,l)
                sigm(i+3) = sigl(i+3,l)
                epsm(i  ) = eps (i  ,1,l)
                epsm(i+3) = eps (i+3,1,l)*0.5d0
              end do ! i

              mct = mct - 4
              if(stype.eq.8) then
                call pstr3d(sigm,sigm(7))
                call pstr3d(epsm,epsm(7))
                if(mct.le.0) then
                  write(iow,2001) o,head
                  if(ior.lt.0) write(*,2001) o,head
                  mct = 50
                endif
                write(iow,2002) n,ma,xx(1,l),xx(2,l),
     &                         (sigm(i),i=7,9),(epsm(i),i=7,9),
     &                         (sigm(i),i=1,6),(eps(i,1,l),i=1,6)
                if(ior.lt.0) then
                  write(*,2002) n,ma,xx(1,l),xx(2,l),
     &                         (sigm(i),i=7,9),(epsm(i),i=7,9),
     &                         (sigm(i),i=1,6),(eps(i,1,l),i=1,6)
                endif
              else
                call pstr2d(sigm,sigm(7))
                call pstr2d(epsm,epsm(7))
                if(mct.le.0) then
                  write(iow,2003) o,head
                  if(ior.lt.0) write(*,2003) o,head
                  mct = 50
                endif
                write(iow,2004) n,ma,xx(1,l),xx(2,l),
     &                         (sigm(i),i=7,9),(epsm(i),i=7,9),
     &                         (sigm(i),i=1,4),(eps(i,1,l),i=1,4)
                if(ior.lt.0) then
                  write(*,2004) n,ma,xx(1,l),xx(2,l),
     &                         (sigm(i),i=7,9),(epsm(i),i=7,9),
     &                         (sigm(i),i=1,4),(eps(i,1,l),i=1,4)
                endif
              endif
            end do ! l

c         Project stresses onto nodes

          else

c           Store strains for plots

            do l = 1,lint
              do j = 1,6
                epsp(j,l) = eps(j,1,l)
              end do ! j
            end do ! l

            call slcn2d(ix,sigl,epsp,r,s,r(nen+1,1),nel,16)
          endif

c       Compute J-integrals and material forces

        elseif(isw.eq.16) then

          call pjint2d(d,ul,tl,shp2,jvol,epsv,sigl,r,ndf,ndm,lint,16)

c       Compute Z-Z projections

        elseif(isw.eq.25) then

          call stcn2z(xl,sigl,epsp,shp2,jac,lint,ndm,nel,16)

        endif ! isw = 4 or 8

      endif ! isw = 3 or 4 or 6 or 8 or 14 or 16

c     Formats

2001  format(a1,20a4//5x,'Element Stresses and Strains'//
     &      5x,'Elmt Mat    1-coord    2-coord'/
     &   15x,' 1-stress   2-stress      Angle',
     &    2x,' 1-strain   2-strain      Angle'/
     &    2x,'11-stress  22-stress  33-stress  12-stress',
     &    2x,'23-stress  31-stress'/'11-strain  22-strain  33-strain',
     &    2x,'12-strain  23-strain  31-strain'/39(' -'))
2002  format(i9,i4,1p,2e11.3/13x,2(1p,2e11.3,0p,f11.2)/(13x,1p,6e11.3))

2003  format(a1,20a4//5x,'Element Stresses and Strains'//
     &      5x,'Elmt Mat    1-coord    2-coord'/
     &   15x,' 1-stress   2-stress      Angle',
     &    2x,' 1-strain   2-strain      Angle'/
     &   15x,'11-stress  22-stress  33-stress  12-stress'/
     &   15x,'11-strain  22-strain  33-strain  12-strain'/39(' -'))
2004  format(i9,i4,1p,2e11.3/13x,2(1p,2e11.3,0p,f11.2)/(13x,1p,4e11.3))

c4001  format(' *ERROR* Plane stress not implemented for MIXEd type.')

      end
