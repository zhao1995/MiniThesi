c$Id:$
      subroutine sld1d2(d,ul,xl,tl,s,r,ndf,ndm,nst,isw, ther)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Use constants from 'pconstant.h'                 14/11/2006
c       2. Use d(185) for augmenting                        14/03/2007
c       3. Set quadrature for d(5).eq.0                     26/03/2009
c       4. Use 'quard1d' and 'interp1d' for solution        01/04/2009
c       5. Revise augment: augfa > augf                     14/04/2009
c       6. Move set of 'npm' to 'quadr1d'                   15/08/2009
c       7. Add flag to call list                            03/03/2010
c       8. Add 'l' to modlsd call                           05/01/2012
c       9. Add epsl on call to slcn1d                       01/01/2013
c-----[--.----+----.----+----.-----------------------------------------]
c     Purpose: Two-dimensional mixed u-p-theta small deformation element

c      Inputs:
c         d(*)  - Element parameters
c         ul(ndf,*) - Current nodal solution parameters
c         xl(ndm,*) - Nodal coordinates
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
      include  'pconstant.h'
      include  'part0.h'
      include  'pmod2d.h'
      include  'qudshp.h'
      include  'rdata.h'
      include  'comblk.h'

      logical   ther
      integer   ndf,ndm,nst,isw,i,i1,j,jj,j1,l,nhi,nhv,istrt,nn
      real*8    alam,ha, epp, dv,dl,d1, type, bt1
      real*8    dsigtr, mpress, dmass, dmshp, dtheta, cfac,lfac, fac
      real*8    al, ac, vl, x0
      real*8    d(*),      ul(ndf,nen,*), xl(ndm,*),    tl(*), s(nst,*)
      real*8    r(ndf,*),  xx(10),        aa(6,6,5,10), dd(7,7)
      real*8    sigm(9),   sigl(16,10),   bpra(3),      bbar(8,8)
      real*8    phi(3,10), theta(3,10),   hh(3,3),      bbd(7)
      real*8    press(10), pbar(10),      hsig(6),      eps(6,10)
      real*8    irad(10),  ta(10),        epsd(3),      epsv(10)
      real*8    dther(7),  bf(3),         epsl(3,10)

      save

      data    nhi   / 2 /

c     TEMPORARY SET OF TEMPERATURE

      data    ta    / 10*0.0d0 /

c     Augmented Lagrangian update for nested iteration

      if(isw.eq.10) then

        hr(nh2) = hr(nh2) + augf*d(185)*hr(nh2+1)

c     Compute tangent stiffness and residual force vector

      elseif(isw.eq. 3 .or. isw.eq. 4 .or. isw.eq. 6 .or.
     &       isw.eq. 8 .or. isw.eq.14 .or. isw.eq.16 .or.
     &       isw.eq.25) then

c       Proportional body forces

        call sbodyf(d, bf)

        estore = 0.0d0

c       Set element quadrature order

        call quadr1d(d)

c       Set number of history terms / quadradure point

        type   = max(0,stype - 2)
        nhv    = nint(d(15))
        istrt  = nint(d(84))

c       Center radius estimate

        x0 = 0.5d0*(xl(1,1) + xl(1,2))

c       MECHANICAL ELEMENT

        do l = 1,lint

c         Shape functions and derivatives

          call interp1d(l, xl, ndm,nel,.false.)

c         Mixed volume effect and temperature projection

          do i = 1,3
            theta(i,l) = 0.0d0
          end do ! i

          ta(l)      = -d(9)

c         Compute radial coordinate

          xx(l) = 0.0d0
          do i = 1,nel
            xx(l) = xx(l) + shp1(2,i,l)*xl(1,i)
          end do ! i

c         Compute volumetric strain from displacements

          if(stype.lt.3) then
            irad(l) = 0.0d0
          else
            jac(l) = jac(l)*xx(l)
            irad(l) = 1.d0/xx(l)
          endif

          do i = 1,nel
            fac        = shp1(1,i,l) + shp1(2,i,l) * irad(l)
            theta(1,l) = theta(1,l) + fac        * ul(1,i,1)
            theta(2,l) = theta(2,l) + fac        * ul(1,i,2)
            theta(3,l) = theta(3,l) + fac        * ul(1,i,3)
            ta(l)      = ta(l)      + shp1(2,i,l) * tl(i)
          end do ! i

c         Set the pressure functions

          phi(1,l) = 1.d0
          phi(2,l) = xx(l) - x0
          phi(3,l) = phi(2,l)**2
        end do ! l

c       Mixed model for volumetric and temperature response

        call bbar1s(phi,shp1,jac,lint,nel,npm,hh,irad,theta,bbar)

c       Compute strains and stresses at quadrature points

        nn = nhi
        do l = 1,lint
          call strn1m(shp1(1,1,l),xl,ul,theta(1,l),irad(l),
     &                ndm,ndf,nel,nen,eps(1,l))

          epsv(l) = theta(1,l)

          call modlsd(l,d,ta(l),eps(1,l),hr(nn+nh1),hr(nn+nh2),nhv,
     &                istrt,aa(1,1,1,l),sigl(1,l),alam,ha,isw)

c         Volumetric stress

          pbar(l) = one3*(sigl(1,l) + sigl(2,l) + sigl(3,l))

          nn = nn + nhv
        end do ! l

c       Integrate pressure over volume

        if(nel.eq.2) then

          mpress = 0.0d0
          do l = 1,lint
            mpress  = mpress  + pbar(l)*jac(l)
          end do ! l

c         Divide pressure by volume

          press(1) = mpress * hh(1,1)
          do l = 2,lint
            press(l) = press(1)
          end do ! l

c       Quadratic and cubic element

        else

          do i = 1,npm
            sigm(i) = 0.0d0
          end do ! i
          do l = 1,lint
            mpress  = pbar(l)*jac(l)
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

c           Integration order set to static

            if((d(7).ge.0.0 .or. d(183).ne.0.0d0) .and.
     &                 (ndfo(1).gt.0 .or. shflg)) then
              cfac = d(7)
              lfac = 1.d0 - cfac
            else
              cfac = 0.0d0
              lfac = 0.0d0
            endif

            do l = 1,lint

c             Store time history plot data for element

              i = 6*(l-1)
              do j = 1,3
                tt(j+i) = sigl(j,l)
                sigm(j) = sigl(j,l)*jac(l)
              end do ! j

c             Compute acceleration

              dmass = d(4)*jac(l)
              al = 0.0d0
              do j = 1,nel
                al = al + shp1(2,j,l)*ul(1,j,5)
              end do ! j
              al = al*cfac

c             Rayleigh damping

              if(d(77).ne.0.0d0) then
                vl = 0.0d0
                do j = 1,nel
                  vl = vl + shp1(2,j,l)*ul(1,j,4)
                end do ! j
                vl = cfac*vl

c               Compute mass damping residual

                do i = 1,nel
                  fac    = shp1(2,i,l)*jac(l)*d(77)*d(4)
                  r(1,i) = r(1,i) - (vl + lfac*ul(1,i,4))*fac
                end do ! i
              endif

              if(d(78).ne.0.0d0) then
                do i = 1,3
                  epsd(i) = 0.0d0
                end do ! i
                do j = 1,nel
                  epsd(1) = epsd(1) + shp1(1,j,l)*ul(1,j,4)
                end do ! j
                fac = one3*(theta(2,l) - epsd(1) - epsd(2) - epsd(3))
                epsd(1) = (epsd(1) + fac)*d(78)*jac(l)
                epsd(2) = (epsd(2) + fac)*d(78)*jac(l)
                epsd(3) = (epsd(3) + fac)*d(78)*jac(l)

                do i = 1,3
                  do j = 1,3
                    sigm(j) = sigm(j) + aa(j,i,5,l)*epsd(i)
                  end do ! j
                end do ! i
              endif

c             Compute residual

              do j = 1,nel
                ac     = d(4)*(al + lfac*ul(1,j,5))
                r(1,j) = r(1,j) + (bf(1) - ac)*shp1(2,j,l)*jac(l)
     &                          - shp1(1,j,l)*sigm(1)
     &                          - shp1(2,j,l)*sigm(3)*irad(l)
              end do ! j

c             Compute mixed tangent stiffness matrix

              if(isw.eq.3) then

c               Multiply tangent moduli by volume element

                call dmatmx( aa(1,1,1,l), dd )
                d1 = jac(l)*ctan(1)
                do i = 1,7
                  do j = 1,7
                    dd(i,j) = dd(i,j)*d1
                  end do ! j
c                 Thermo-mechanical coupling
                  dther(i) = aa(i,1,2,l)*d1
                end do ! i

                if(d(78).ne.0.0d0) then
                  d1 = jac(l)*d(78)*ctan(2)
                  do i = 1,6
                    do j = 1,6
                      dd(i,j) = dd(i,j) + aa(i,j,5,l)*d1
                    end do ! j
                  end do ! i
                endif

c               Mass factors

                dv = (ctan(3) + d(77)*ctan(2))*jac(l)*d(4)*cfac
                dl = (ctan(3) + d(77)*ctan(2))*jac(l)*d(4)*lfac

c               Compute row terms

                i1 = 1
                do i = 1,nel

c                 Compute bmat-t * dd * dvol

                  do jj = 1,7

                    bbd(jj) =  shp1(1,i,l)*dd(1,jj)
     &                      +  shp1(2,i,l)*dd(3,jj)*irad(l)
     &                      +   bbar(i,l)*dd(7,jj)
                  end do ! jj
c                          _
c                 Compute: B_trans * D * alpha * j * w

                  if(ther) then
                    fac = (dther(1) + dther(2) + dther(3))*one3
                    bt1 =  shp1(1,i,l)*dther(1)
     &                  +  shp1(2,i,l)*dther(3) * irad(l)
     &                  + (bbar(i,l) - shp1(1,i,l)
     &                               - shp1(2,i,l)*irad(l))*fac
                  endif

c                 Compute tangent stiffness

                  fac = shp1(2,i,l)*dl
                  s(i1,i1) = s(i1,i1) + fac

                  dmshp = shp1(2,i,l)*dv

                  j1 = 1
                  do j = 1,nel

c                   Inertial tangent

                    s(i1,j1) = s(i1,j1) + dmshp*shp1(2,j,l)

c                   Compute mechanics part of tangent stiffness

                    s(i1,j1) = s(i1,j1) + bbd(1)*shp1(1,j,l)
     &                                  + bbd(3)*shp1(2,j,l)*irad(l)
     &                                  + bbd(7)*bbar(j,l)
c                   Thermo-mechanical coupling matrix

                    if(ther) then
                      s(i1,j1+1)  = s(i1,j1+1) + bt1*shp1(2,j,l)
                    endif

                    j1 = j1 + ndf
                  end do ! j
                  i1 = i1 + ndf
                end do ! i
              endif ! isw = 3
            end do ! l

          endif ! isw = 3 or 6

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
                sigm(i) = sigl(i,l)
              end do ! i

              mct = mct - 2
              if(mct.le.0) then
                write(iow,2001) o,head
                if(ior.lt.0) write(*,2001) o,head
                mct = 50
              endif

c             Compute potential damage variable

              write(iow,2002)  n,ma,(sigm(i),i=1,3),
     &                        xx(l),(eps(i,l),i=1,3)
              if(ior.lt.0 .and. pfr) then
                write(*,2002)  n,ma,(sigm(i),i=1,3),
     &                        xx(l),(eps(i,l),i=1,3)
              endif
            end do ! l

c         Project stresses onto nodes

          else
            do l = 1,lint
              epsl(1:3,l) = eps(1:3,l)
            end do ! l
            call slcn1d(sigl,epsl,shp1,jac,r,s,r(nen+1,1),lint,nel,16)
          endif

c       Compute J-integrals and material forces

        elseif(isw.eq.16) then

          call pjint1d(d,ul,tl,shp1,jac,epsv,sigl,r,ndf,ndm,lint)

c       Compute Z-Z projections

        elseif(isw.eq.25) then

          call stcn1z(xl,sigl,shp1,jac,lint,ndm,nel,16)

        endif ! isw = 4 or 8

      endif ! isw = 3 or 4 or 6 or 8 or 14 or 16

c     Formats

2001  format(a1,20a4//5x,'Element Stresses'//'    Elmt Mat',
     &   '   11-stress   22-stress   33-stress'/
     &   '     1-coord   11-strain   22-strain   33-strain')
2002  format(i8,i4,1p,3e12.3/1p,4e12.3/1x)

      end
