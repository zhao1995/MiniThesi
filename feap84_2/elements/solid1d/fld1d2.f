c$Id:$
      subroutine fld1d2(d,ul,xl,s,r,ndf,ndm,nst,isw)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Use constants from 'pconstant.h'                 14/11/2006
c       2. Use d(185) for augmenting                        14/03/2007
c       3. Set quadrature for d(5).eq.0                     26/03/2009
c       4. Use 'quadr1d' and 'interp1d' for solution        01/04/2009
c       5. Revise augment: nhi -> point                     14/04/2009
c       6. Move set of 'npm' to 'quadr1d'                   15/08/2009
c       7. Add flag to call list                            03/03/2010
c       8. Add epsl on call to slcn1d                       01/01/2013
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Finite deformation mixed model: u-p-theta

c      Inputs:
c         d(*)  - Element parameters
c         ul(ndf,*) - Current nodal solution parameters
c         xl(ndm,*) - Nodal coordinates
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
      include  'elcoor.h'
      include  'eldata.h'
      include  'elengy.h'
      include  'elplot.h'
      include  'eltran.h'
      include  'hdata.h'
      include  'iofile.h'
      include  'pconstant.h'
      include  'part0.h'
      include  'pmod2d.h'
      include  'p_point.h'
      include  'qudshp.h'
      include  'rdata.h'
      include  'comblk.h'

      logical   dynflg
      integer   ndf,ndm,nst, isw, i,i1, j,jj,j1, k, l, nhv,nn, istrt
      real*8    augfp,d1,epp,thlog, xxm,x0, xlamd, ha, bd3, body(3)
      real*8    ta,qfact,dsigtr,dpress,mpress,dmass,cmshp,lmshp,dtheta
      real*8    d(*),        ul(ndf,nen,*), xl(ndm,*),   s(nst,*)
      real*8    df(9,8),     fi(9,2,8),     finv(9,8),   detf(2,8)
      real*8    xr(8),       xu(8),         ru(8),       r(ndf,*)
      real*8    bbd(7),      bei(6),        ad(7,7,5,8), dd(7,7)
      real*8    shpr(8,8),   shpbar(8,8),   dvl0(8),     dvol(8)
      real*8    sigm(9),     sigl(16,8),    bpra(3),     rr(3)
      real*8    acc,         theta(2,8),    hh(3,3),     hsig(3)
      real*8    press(8),    weng(8),       ur(8),       phi(3,8)
      real*8    epsm(3),     epsl(3,8),     egreen(6),   ealmansi(6)

      save

c     TEMPORARY SET OF TEMPERATURE

      data    ta    / 0.0d0 /

c     Augmented Lagrangian update for nested iteration

      if(isw.eq.10) then

        d1      = augfp*d(185)
        hr(nh2) = hr(nh2) + d1*hr(nh2+1)

c     Compute tangent stiffness and residual force vector

      elseif(isw.eq. 3 .or. isw.eq. 4 .or. isw.eq. 6 .or.
     &       isw.eq. 8 .or. isw.eq.14 .or. isw.eq.16 .or.
     &       isw.eq.25) then

        augfp  = augf
        estore = 0.0d0

c       Compute current geometry

        do j = 1,nel
          xu(j) = xl(1,j) + ul(1,j,1)
        end do ! j

c       Center coordinate

        x0 = 0.5d0*(xl(1,1) + xl(1,2))

c       Set quadrature order

        call quadr1d(d)

c       Get shape functions and derivatives in geometry at time t_n+1

        do l = 1,lint
          call interp1d(l, xl, ndm,nel,.false.)

c         Compute coordinates at gauss points

          xr(l) = 0.0d0
          ur(l) = 0.0d0
          do i = 1,nel
            xr(l) = xr(l) + xl(1,i)  *shp1(2,i,l)
            ur(l) = ur(l) + ul(1,i,1)*shp1(2,i,l)
          end do ! i

          phi(1,l) = 1.0d0
          phi(2,l) = xr(l) - x0
          phi(3,l) = phi(2,l)**2

        end do ! l
        xref(2) = 0.0d0
        xref(3) = 0.0d0
        xcur(2) = 0.0d0
        xcur(3) = 0.0d0

c       Set number of history terms / quadradure point

        nhv   = nint(d(15))
        istrt = nint(d(84))

c       MECHANICAL ELEMENT

        if(isw.eq.3 .or. isw.eq.6 .or. isw.eq.14) then

c         Compute f, finv, df and det(fei) at conf t-n+1

          call kine1m(shp1,xl,ul,fi,finv,df,detf,ndm,ndf,nel,nen,lint)

c         compute volume at current state

          do l = 1,lint
            dvol(l) = jac(l)*detf(1,l)
          end do ! l

c         Mixed model for volumetric response

          call bbar1m(phi,shp1,shpr,dvol,detf,lint,nel,npm,
     &                hh,theta,shpbar)

c         Compute mixed model deformation gradient

          call fbar1m(fi,detf,theta,lint)

c         Initialicze storage for augmented function

          if(npm.gt.1) then
            do i = 1,npm
              rr(i) = 0.0d0
            end do ! i
          endif

c         Compute Cauchy stresses and spatial tangent tensor at t-n+1

          dynflg = ctan(3).ne.0.0d0
          nn     = 2*npm
          do l = 1,lint

c           Compute coordinates in reference and current configuration

            xref(1) = xr(l)
            xcur(1) = xr(l) + ur(l)

c           Set augmented multiplier
            xlamd = hr(nh2)
            do i = 1,npm-1
              xlamd = xlamd + phi(i+1,l)*hr(nh2+i)
            end do ! i

            call modlfd(l,d,fi(1,1,l),finv(1,l),df(1,l),theta(1,l),ta,
     &                  hr(nn+nh1),hr(nn+nh2),nhv,istrt,ad(1,1,1,l),
     &                  sigl(1,l),bei,xlamd,ha,.true.,isw)

c           Save augmented function
            if(npm.gt.1) then
              do i = 1,npm
                rr(i) = rr(i) + phi(i,l)*ha*jac(l)
              end do ! i
            endif
            nn = nn + nhv
          end do ! l

          if(isw.eq.14) return

c         Set augmented function

          if(npm.eq.1) then
            hr(nh2+1) = ha
          else
            point = nh2 + npm - 1
            do i = 1,npm
              hr(point+i) = 0.0d0
              do j = 1,npm
                hr(point+i) = hr(point+i) + hh(i,j)*rr(j)
              end do ! j
            end do ! i
          endif

c         Compute mixed pressure

          if(isw.eq.3 .or. isw.eq.6) then

c           Set body load levels

            call sbodyf(d, body)

            if(npm.eq.1) then

              press(1) = 0.0d0
              do l = 1,lint

c               Modify volume element and integrate pressure
c               over reference volume

                press(1) = press(1) + one3*(sigl(1,l) + sigl(2,l)
     &                              + sigl(3,l))*jac(l)
                jac(l)  = jac(l) * theta(1,l)

              end do ! l

c             Divide pressure by reference volume

              press(1) = press(1) * hh(1,1)
              do l = 2,lint
                press(l) = press(1)
              end do ! l

            else

              do i = 1,npm
                sigm(i) = 0.0d0
              end do ! i

              do l = 1,lint

c               Modify volume element and integrate pressure
c               over reference volume

                mpress   = one3*(sigl(1,l) + sigl(2,l)
     &                         + sigl(3,l))*jac(l)
                sigm(1) = sigm(1) + mpress
                do k = 2,npm
                  sigm(k) = sigm(k) + mpress*phi(k,l)
                end do ! k
                dvol(l) = jac(l) * theta(1,l)

              end do ! l

c             Divide pressure by reference volume

              do i = 1,npm
                hsig(i) = 0.0d0
                do j = 1,npm
                  hsig(i) = hsig(i) + hh(i,j)*sigm(j)
                end do ! j
              end do ! i
              do l = 1,lint
                press(l) = hsig(1)
                do j = 2,npm
                  press(l) = press(l) + hsig(j)*phi(j,l)
                end do ! j
              end do ! l

            endif

            do l = 1,lint

c             Compute mixed stress and multiply by volume element

              dsigtr  = press(l)*detf(1,l)/theta(1,l)
     &                - (sigl(1,l)+sigl(2,l)+sigl(3,l))*one3
              sigm(1) =  sigl(1,l) + dsigtr
              sigm(2) =  sigl(2,l) + dsigtr
              sigm(3) =  sigl(3,l) + dsigtr

c             Store time history plot data for element

              i = 6*(l-1)
              do j = 1,3
                tt(j+i) = sigm(j)
                sigm(j) = sigm(j)*dvol(l)
              end do ! j

c             Compute acceleration

              if((d(7).ge.0.0 .or. d(183).ne.0.0d0) .and.
     &                   (ndfo(1).gt.0 .or. shflg)) then
                dmass = d(4)*dvl0(l)
              else
                dmass = 0.0d0
              endif

              cmshp = dmass * d(7)    ! Consistent
              lmshp = dmass - cmshp   ! Lumped
              acc = 0.0d0
              do j = 1,nel
                acc = acc + shp1(2,j,l)*ul(1,j,5)
              end do ! j
              acc = acc*cmshp - body(1)*jac(l)

c             Compute residual

              do j = 1,nel
                ru(j)   = shp1(1,j,l)*sigm(1)
                r(1,j)  = r(1,j) - ru(j) - shpr(j,l)*sigm(3)
     &                           - shp1(2,j,l)*(acc - lmshp*ul(1,j,5))
              end do ! j

c             Compute mixed tangent stiffness matrix

              if(isw.eq.3) then

c               Part 1: Geometric tangent matrix

                if(gflag) then
                  i1 = 1
                  do i = 1,nel
                    bd3 = shpr(i,l)*sigm(3)*ctan(1)
                    j1 = 1
                    do j = 1,nel
                      s(i1,j1) = s(i1,j1) + shp1(1,i,l)*ru(j)*ctan(1)
     &                                    + shpr(j,l)*bd3
                      j1 = j1 + ndf
                    end do ! j
                    i1 = i1 + ndf
                  end do ! i
                endif ! gflag

c               Part 2: Material tangent matrix

c               Modify tangent moduli for stress factors

                mpress = press(l)*detf(1,l)/theta(1,l)
                dpress = one3*(sigl(1,l) + sigl(2,l) + sigl(3,l))

                call dmatdx(ad(1,1,1,l),sigl(1,l),dpress,mpress)

c               Multiply tangent moduli by volume element

                d1 = dvol(l)*ctan(1)
                do i = 1,7
                  do j = 1,7
                    dd(i,j) = ad(i,j,1,l)*d1
                  end do ! j
                end do ! i

c               Compute row terms

                dmass = ctan(3)*dmass
                i1    = 1
                do i = 1,nel

c                 Compute bmat-t * dd * jac

                  do jj = 1,7
                    bbd(jj) = shp1(1,i,l)*dd(1,jj)
     &                      + shpr(  i,l)*dd(3,jj)
     &                      + shpbar(i,l)*dd(7,jj)
                  end do ! jj

                  lmshp = shp1(2,i,l)*dmass
                  cmshp = lmshp*d(7)
                  s(i1,i1) = s(i1,i1) + lmshp - cmshp

                  j1 = 1
                  do j = 1,nel
                    s(i1,j1) = s(i1,j1) + cmshp*shp1(2,j,l)
     &                                  + bbd(1)*shp1(1,j,l)
     &                                  + bbd(3)*shpr(  j,l)
     &                                  + bbd(7)*shpbar(j,l)
                    j1 = j1 + ndf
                  end do ! j
                  i1 = i1 + ndf
                end do ! i
              endif ! isw = 3
            end do ! l

          endif ! isw .eq. 3 or 6

c       Output stresses.

        elseif(isw.eq.4.or.isw.eq.8.or.isw.eq.16.or.isw.eq.25) then

          do i = 1,9
            sigm(i) = 0.0d0
          end do ! i
          do i = 1,3
            bpra(i) = 0.0d0
            epsm(i) = 0.0d0
          end do ! i
          xxm    = 0.0d0
          epp    = 0.0d0
          dtheta = 0.0d0
          qfact  = 1.d0/dble(lint)

c         Compute f, finv, df and det(fei) at conf t-n+1

          call kine1m(shp1,xl,ul,fi,finv,df,detf,ndm,ndf,nel,nen,lint)

          call bbar1m(phi,shp1,shpr,dvol,detf,lint,nel,npm,
     &                hh,theta,shpbar)

          call fbar1m(fi,detf,theta,lint)

c         Second loop over Gauss points

          nn  = 2*npm
          do l = 1,lint

c           Compute coordinates in reference and current configuration

            xref(1) = xr(l)
            xcur(1) = xr(l) + ur(l)

c           Set augmented function

            xlamd = hr(nh2)
            do i = 1,npm-1
              xlamd = xlamd + phi(i+1,l)*hr(nh2+i)
            end do ! i

c           Compute Cauchy stresses and spatial tangent tensor at t-n+1

            call modlfd(l,d,fi(1,1,l),finv(1,l),df(1,l),theta(1,l),ta,
     &                  hr(nn+nh1),hr(nn+nh2),nhv,istrt,ad,sigl(1,l),
     &                  bei,hr(nh2),hr(nh2+1),.true.,isw)
            weng(l) = estore

c           Compute Green-Lagrange strains and Almansi strains

            call fstrain(fi(1,1,l),finv(1,l), egreen, ealmansi)

            epsl(1:3,l) = ealmansi(1:3)

            call pstr3d(bei, bpr)

c           Average stresses and stretches for printing

            xxm = xxm + qfact *xu(l)
            do i = 1,3
              bpra(i) = bpra(i) + 0.5d0*qfact*log(bpr(i))
            end do ! i
            do i = 1,3
              sigm(i) = sigm(i) + qfact*sigl(i,l)
              epsm(i) = epsm(i) + qfact*epsl(i,l)
            end do ! i
            epp      = epp      + qfact*sigl( 9,l)
            dtheta   = dtheta   + qfact*theta(1,l)
            nn = nn + nhv
          end do ! l

c         Output stresses

          if (isw .eq. 4) then

            mct = mct - 2
            if(mct.le.0) then
              write(iow,2001) o,head
              if(ior.lt.0) write(*,2001) o,head
              mct = 50
            endif

c           Compute potential damage variable

            thlog = log(abs(dtheta))
            write(iow,2002) n,ma,xxm,(sigm(i),i=1,3),epsm,bpra,thlog,epp
            if(ior.lt.0) then
              write(*,2002) n,ma,xxm,(sigm(i),i=1,3),epsm,bpra,thlog,epp
            endif

c         Project stresses onto nodes

          elseif(isw.eq.8) then

            call slcn1d(sigl,epsl,shp1,jac,r,s,r(nen+1,1),lint,nel,16)

c         Compute fracture indices

          elseif(isw.eq.16) then

            call pfrac1f(fi,theta,sigl,weng, shp1,jac, r,
     &                   lint,ndf,ndm,3)

c         Compute Z-Z projections

          elseif(isw.eq.25) then

            call stcn1z(xl,sigl,shp1,jac,lint,ndm,nel,16)

          endif
        endif ! isw = 4 or 8 or 16 or 25

      endif ! isw

c     Formats

2001  format(a1,20a4//5x,'Element Stresses'//'    Elmt Mat',
     &   '     1-Coord   11-Stress   22-Stress   33-Stress'/12x,
     &   '     Almansi   11-Strain   22-Strain   33-Strain'/12x,
     &   '   log(lam1)   log(lam2)   log(lam3)     log-J     eff-ep')

2002  format(i8,i4,1p,4e12.4/24x,1p3e12.4/12x,1p5e12.4/1x)

      end
