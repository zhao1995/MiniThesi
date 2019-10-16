c$Id:$
      subroutine sld2d1(d,ul,xl,ix,tl,s,p,ndf,ndm,nst,isw, ther)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Increase quadrature order for 3-node to l=3      21/03/2007
c       2. Add direct call to quadr2d and interp2d          11/11/2008
c       3. Remove 'nel' from call to 'quadr2d'              23/01/2009
c       4. Increase order of element and quadrature to 6    17/02/2009
c       5. Add principal strains to stress/strain outputs   02/03/2009
c       6. Correct thermal coupling term location 'tdf'     30/04/2009
c          Correct index on vv.
c       7. Allow up to 64 node elements                     02/06/2009
c       8. Add geometric stiffness                          24/06/2009
c       9. Move computation of psil to material models      26/07/2009
c      10. Add 'ndf' to call of geomnd                      14/08/2010
c      11. Compute v_avg and sig_33 for multiscale use      03/12/2010
c      12. Add 'l' to modlsd call                           05/01/2012
c      13. Add average of density for multiscale            10/05/2012
c      14. Add eps on call to slcn2d                        01/01/2013
c      15. Pass strains to stcn2z for z-zhu projections     01/01/2014
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Plane and axisymmetric linear elastic element routine

c      Inputs:
c         d(*)  - Element parameters
c         ul(ndf,*) - Current nodal solution parameters
c         xl(ndm,*) - Nodal coordinates
c         ix(*)     - Global nodal connections
c         ndf       - Degree of freedoms/node
c         ndm       - Mesh coordinate dimension
c         nst       - Element array dimension
c         isw       - Solution option switch

c      Outputs:
c         s(nst,*)  - Element array
c         p(ndf,*)  - Element vector
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'bdata.h'
      include  'cdata.h'
      include  'complx.h'
      include  'elbody.h'
      include  'eldata.h'
      include  'elplot.h'
      include  'eltran.h'
      include  'fdata.h'
      include  'hdata.h'
      include  'iofile.h'
      include  'oelmt.h'
      include  'part0.h'
      include  'pmod2d.h'
      include  'ptdat6.h'
      include  'qudshp.h'
      include  'rdata.h'
      include  'strnum.h'
      include  'comblk.h'

      logical   ther
      integer   ndf,ndm,nst,isw,j,k,l,j1,k1,nhv,nn,istrt,ncp,ncs
      integer   ix(*), tdof, tdf
      real*8    alam,ha, aj1,aj2,aj3,aj0,xx,yy,lfac,cfac,sfac,jac0,dv,ta
      real*8    bd(3,6), vv(3), bt1,bt2,bt3, sh3
      real*8    d(*),ul(ndf,nen,*),xl(ndm,*),tl(*),s(nst,*),p(ndf,*)
      real*8    sig(16,64),eps(9,3),epsp(6,64),epsv(64)
      real*8    dd(6,6,5),shpr(64), mom(3,64),epp(6)
      real*8    peps(3)

      save

      data      eps     / 27*0.0d0 /
      data      alam,ha /  2*0.0d0 /

c     Compute stress-divergence vector (p) and stiffness matrix (s)

      nhv   = nint(d(15))
      istrt = nint(d(84))

c     Set number of stress components

      if(stype.eq.8) then
        ncp = 3
        ncs = 6
      else
        ncp = 2
        ncs = 4
      endif

c     Data inputs

      if( isw.eq.1 ) then

      elseif(isw.eq. 3 .or. isw.eq. 5 .or. isw.eq. 6 .or.
     &       isw.eq.13 .or. isw.eq.14) then

c       Integration order set to static

        if((isw.eq. 3 .or. isw.eq. 6 .or. isw.eq.13) .and.
     &     (d(7).ge.0.0 .or. d(183).ne.0.0d0)        .and.
     &     (ndfo(1).gt.0 .or. shflg)               ) then
          cfac = d(7)
          lfac = 1.d0 - cfac
          if(nel.eq.3) then
            call masst3(stype,cfac,d(4),xl,ul(1,1,5),p,s)
            cfac = 0.0d0
          endif
        else
          cfac = 0.0d0
          lfac = 0.0d0
        endif

c       Initialize momenta

        if(isw.eq.13) then
          do j = 1,nel
            do k = 1,3
              mom(k,j) = 0.0d0
            end do ! k
          end do ! j
        endif

c       Compute Gauss quadrature points and weights

        call quadr2d(d,.true.)

c       Zero shpr matrix

        do j = 1,16
          shpr(j) = 0.0d0
        end do ! j

c       Numerical integration loop

        if(stype.eq.1)  then
          nn = lint
        else
          nn = 0
        endif

        do l = 1,lint

c         Get shape functions for point

          call interp2d(l, xl,ix, ndm,nel, .false.)

c         Compute stresses and strains

          call strn2d(d,xl,ul,tl,shp2(1,1,l),ndf,ndm,nel,
     &                xx,yy,ta,eps)

c         Plane stress modification

          if(stype.eq.1) then
            eps(3,1) = hr(nh2-1+l)
            eps(3,2) = hr(nh1-1+l)
          endif

          call modlsd(l,d,ta,eps,hr(nh1+nn),hr(nh2+nn),nhv,istrt,
     &                dd,sig(1,l),alam,ha,isw)

c         Store plane stress strain

          if(stype.eq.1) then
            hr(nh2-1+l) = eps(3,1)
          endif

c         Multiply jacobian by radius for axisymmetry

          if(stype.eq.3 .or. stype.eq.8) then
            jac0   = jac(l)
            jac(l) = jac(l)*xx
            do j = 1,nel
              shpr(j) = shp2(3,j,l)/xx
            end do ! j
          else
            jac0 = 0.0d0
          end if

          if(isw.eq.3 .or. isw.eq.6) then

c           Store time history plot data for element

            k = 10*(l-1)
            do j = 1,6
              tt(j+k) = sig(j,l)
            end do ! j
            k = k + 6
            do j = 1,4
              tt(j+k) = eps(j,1)
            end do ! j

c           Set element parameters for multiscale

            v_avg  = v_avg + jac(l)
            v_rho  = v_rho + jac(l)*d(4)
            sig_33 = sig_33 + jac(l)*sig(3,l)

c           Rayleigh damping effects

            dv = d(4)*(ctan(3) + d(77)*ctan(2))*jac(l)

            if(d(78).ne.0.0d0) then
              call rays2d(d,shp2(1,1,l),sig(1,l),dd(1,1,5),ul(1,1,4),xl,
     &                    ndf,ndm,nel)
              sfac = d(78)*ctan(2)
            else
              sfac = 0.0d0
            endif

c           Compute gravity, thermal, inertia, and stress contributions

            call resid2d(cfac,lfac,jac(l),jac0,xx,yy,shp2(1,1,l),
     &                   sig(1,l),d,ul(1,1,4),ul(1,1,5),p,ndf)

c           Thermal coupling property

            if(ther .and. d(129).ne.0.0d0) then
              tdof = nint(d(19))

              if(tdof.gt.0) then

c               Velocity

                do k = 1,ncp
                  vv(k) = 0.0d0
                  do j = 1,nel
                    vv(k) = vv(k) + shp2(3,j,l)*ul(k,j,4)
                  end do ! j
                end do ! k

                do j = 1,nel

c                 Compute B_trans * D * alpha * j * w

                  aj1  = shp2(1,j,l)*jac(l)*d(129)
                  aj2  = shp2(2,j,l)*jac(l)*d(129)
                  aj3  = shp2(3,j,l)*jac0*d(129)
                  bt1  = aj1*dd(1,1,2) + aj3*dd(3,1,2) * aj2*dd(4,1,2)
                  bt2  = aj2*dd(2,1,2) + aj1*dd(4,1,2)
                  p(tdof,j) = p(tdof,j) - bt1*vv(1) - bt2*vv(2)
                  if(stype.eq.8) then
                    bt3       = aj2*dd(5,1,2) + (aj1-aj3)*dd(6,1,2)
                    p(tdof,j) = p(tdof,j) - bt3*vv(3)
                  endif

                end do ! j
              endif

            endif

c           Tangent stiffness computation

            if(isw.eq.3) then

c             Modify tangent for stiffness rayleigh damping

              do j = 1,ncs
                do k = 1,ncs
                  dd(k,j,1) = dd(k,j,1)*ctan(1) + dd(k,j,5)*sfac
                end do ! k
c               Thermo-mechanical coupling
c               dd(k,1,2) = dd(k,1,2)*ctan(1)
              end do ! j

              j1 = 1
              do j = 1,nel

                aj1 = shp2(1,j,l)*jac(l)
                aj2 = shp2(2,j,l)*jac(l)
                aj3 = shp2(3,j,l)*jac0

c               Compute B_trans * D * j * w

                do k = 1,6
                  bd(1,k) = aj1*dd(1,k,1) + aj2*dd(4,k,1)
     &                    + aj3*dd(3,k,1)
                  bd(2,k) = aj2*dd(2,k,1) + aj1*dd(4,k,1)
                  bd(3,k) = aj2*dd(5,k,1) + (aj1-aj3)*dd(6,k,1)
                end do ! k

c               Compute B_trans * D * alpha * j * w

                bt1  = aj1*dd(1,1,2) + aj3*dd(3,1,2) * aj2*dd(4,1,2)
                bt2  = aj2*dd(2,1,2) + aj1*dd(4,1,2)
                bt3  = aj2*dd(5,1,2) + (aj1-aj3)*dd(6,1,2)

c               Compute lumped mass matrix

                aj0          = shp2(3,j,l)*dv
                s(j1  ,j1  ) = s(j1  ,j1  ) + aj0*lfac
                s(j1+1,j1+1) = s(j1+1,j1+1) + aj0*lfac

c               Loop over columns (symmetry noted)

                k1 = 1
                do k = 1,nel
                  s(j1  ,k1  ) = s(j1  ,k1  ) + bd(1,1)*shp2(1,k,l)
     &                                        + bd(1,4)*shp2(2,k,l)
     &                                        + bd(1,3)*shpr(k)
     &                                        + cfac*aj0*shp2(3,k,l)

                  s(j1  ,k1+1) = s(j1  ,k1+1) + bd(1,2)*shp2(2,k,l)
     &                                        + bd(1,4)*shp2(1,k,l)

                  s(j1+1,k1  ) = s(j1+1,k1  ) + bd(2,1)*shp2(1,k,l)
     &                                        + bd(2,4)*shp2(2,k,l)
     &                                        + bd(2,3)*shpr(k)

                  s(j1+1,k1+1) = s(j1+1,k1+1) + bd(2,2)*shp2(2,k,l)
     &                                        + bd(2,4)*shp2(1,k,l)
     &                                        + cfac*aj0*shp2(3,k,l)

                  if(stype.eq.8) then
                    s(j1  ,k1+2) = s(j1  ,k1+2) + bd(1,5)*shp2(2,k,l)
     &                           + bd(1,6)*(shp2(1,k,l) - shpr(k))
                    s(j1+1,k1+2) = s(j1+1,k1+2) + bd(2,5)*shp2(2,k,l)
     &                           + bd(2,6)*(shp2(1,k,l) - shpr(k))
                    s(j1+2,k1+2) = s(j1+2,k1+2) + bd(3,5)*shp2(2,k,l)
     &                           + bd(3,6)*(shp2(1,k,l) - shpr(k))
     &                                          + cfac*aj0*shp2(3,k,l)
                    s(j1+2,k1  ) = s(j1+2,k1  ) + bd(3,1)*shp2(1,k,l)
     &                                          + bd(3,4)*shp2(2,k,l)
     &                                          + bd(3,3)*shpr(k)
                    s(j1+2,k1+1) = s(j1+2,k1+1) + bd(3,2)*shp2(2,k,l)
     &                                          + bd(3,4)*shp2(1,k,l)
                  endif

c                 Add thermo-mechanical coupling term

                  if(ther) then
                    tdf          = tdof - 1
                    sh3          = shp2(3,k,l)*ctan(1)
                    s(j1  ,k1+tdf) = s(j1  ,k1+tdf) + bt1*sh3
                    s(j1+1,k1+tdf) = s(j1+1,k1+tdf) + bt2*sh3
                    if(stype.eq.8) then
                      s(j1+2,k1+tdf) = s(j1+2,k1+tdf) + bt3*sh3
                    endif

c                   Rate terms

                    if(d(129).ne.0.0d0) then
                      sh3            = shp2(3,k,l)*d(129)*ctan(2)
                      s(j1+tdf,k1  ) = s(j1+tdf,k1  ) + bt1*sh3
                      s(j1+tdf,k1+1) = s(j1+tdf,k1+1) + bt2*sh3
                      if(stype.eq.8) then
                        s(j1+tdf,k1+2) = s(j1+tdf,k1+2) + bt3*sh3
                      endif
                    endif

                  endif

                  k1 = k1 + ndf
                end do ! k

                j1 = j1 + ndf
              end do ! j
            end if
            if(cplxfl) then
              call cst2d1(shp2(1,1,l),shpr,jac(l),jac0,xx,eps,ul(1,1,8),
     &                    dd,dd(1,1,3),s(1,nst+1),p,p(1,nen+1))
            endif

c         Compute geometric stiffness

          elseif(isw.eq.5) then

            call geomnd(shp2(1,1,l),shpr,jac(l), sig(1,l),s,
     &                  ndm,ndf,nel,nst)

c         Compute energy

          elseif(isw.eq.13) then

            if(d(14).ne.1.d0) then
              aj0 = d(14)*jac(l)
            else
              aj0 = jac(l)
            endif
            aj1 = d(4)*aj0
            call sengy(ul,shp2(1,1,l), aj0,aj1, lfac,cfac, ncp,3, mom)

          end if
          nn = nn + nhv
        end do ! l

c       Accumulate momenta

        if(isw.eq.13) then
          do j = 1,nel
            do k = 1,ncp
              epl(k) = epl(k) + mom(k,j)
            end do ! i
            epl(4) = epl(4) + xl(1,j)*mom(3,j)
            epl(5) = epl(5) - xl(2,j)*mom(3,j)
            epl(6) = epl(6) + xl(1,j)*mom(2,j) - xl(2,j)*mom(1,j)
          end do ! j
        endif

c       Multiply by thickness if not unity

        if((isw.eq.3 .or. isw.eq.6) .and. d(14).ne.1.d0) then

          do j = 1,nst
            do k = 1,nst
              s(k,j) = s(k,j)*d(14)
            end do ! k
          end do ! j
          do j = 1,nel
            do k = 1,ndf
              p(k,j) = p(k,j)*d(14)
            end do ! k
          end do ! j

        endif

c     Output of element quantities

      elseif(isw.eq.4 .or. isw.eq.8 .or. isw.eq.16 .or. isw.eq.25) then

        call quadr2d(d,.true.)

c       Numerical integration loop

        if(stype.eq.1)  then
          nn = lint
        else
          nn = 0
        endif

c       Compute element stresses, strains, and forces

        do l = 1,lint

c         Compute element shape functions

          call interp2d(l, xl,ix, ndm,nel, .false.)

c         Compute strains and coordinates

          call strn2d(d,xl,ul,tl,shp2(1,1,l),ndf,ndm,nel,
     &                xx,yy,ta,eps)

c         Plane stress modification

          if(stype.eq.1) then
            eps(3,1) = hr(nh2-1+l)
            eps(3,2) = hr(nh1-1+l)
          endif

          call modlsd(l,d,ta,eps,hr(nh1+nn),hr(nh2+nn),nhv,istrt,
     &                dd,sig(1,l),alam,ha,isw)

c         Compute volumetric strain

          epsv(l) = eps(1,1) + eps(2,1) + eps(3,1)
          do j = 1,3
            epp(j  ) = eps(j  ,1)
            epp(j+3) = eps(j+3,1)*0.5d0
          end do ! j

          if(isw.eq.4) then

            mct = mct - 4
            if(stype.eq.8) then
              call pstr3d(sig(1,l),sig(7,l))
              call pstr3d(epp,peps)
              if(mct.le.0) then
                write(iow,2001) o,head
                if(ior.lt.0) write(*,2001) o,head
                mct = 50
              endif
              write(iow,2002) n,ma,xx,yy,(sig(j,l),j=7,9),peps,
     &                       (sig(j,l),j=1,6),(eps(j,1),j=1,6)
              if(ior.lt.0) then
                write(*,2002) n,ma,xx,yy,(sig(j,l),j=7,9),peps,
     &                       (sig(j,l),j=1,6),(eps(j,1),j=1,6)
              endif
            else
              call pstr2d(sig(1,l),sig(7,l))
              call pstr2d(epp,peps)
              if(mct.le.0) then
                write(iow,2003) o,head
                if(ior.lt.0) write(*,2003) o,head
                mct = 50
              endif
              write(iow,2004) n,ma,xx,yy,(sig(j,l),j=7,9),peps,
     &                       (sig(j,l),j=1,4),(eps(j,1),j=1,4)
              if(ior.lt.0) then
                write(*,2004) n,ma,xx,yy,(sig(j,l),j=7,9),peps,
     &                       (sig(j,l),j=1,4),(eps(j,1),j=1,4)
              endif
            endif

c         Store strains for plots

          elseif(isw.eq.8) then

            do j = 1,6
              epsp(j,l) = eps(j,1)
            end do ! j

          endif
          nn = nn + nhv
        end do ! l

c       Compute nodal stress values

        if(isw.eq.8) then

          call slcn2d(ix,sig,epsp,p,s,p(nen+1,1),nel,16)

c       Compute J-integrals and material forces

        elseif(isw.eq.16) then

          call pjint2d(d,ul,tl,shp2,jac,epsv,sig,p,ndf,ndm,lint,16)

c       Compute Z-Z projections

        elseif(isw.eq.25) then

          call stcn2z(xl,sig,epsp,shp2,jac,lint,ndm,nel,16)

        endif

c     Compute nodal stress error values

      elseif(isw.eq.11) then

        call ster2d(ix,d,xl,ul,tl,s,ndf,ndm,nel,nen)

      endif

c     Formats for input-output

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

      end

      subroutine cst2d1(shp,shpr,jac,jac0,xx,epsr,ul,dr,di,
     &                  si,pr,pi)

      implicit   none
      include   'cdata.h'
      include   'eldata.h'
      include   'eltran.h'
      include   'pmod2d.h'
      include   'sdata.h'

      integer    j,j1,k,k1
      real*8     shp(3,*),shpr(*),epsr(*),epsi(6),ul(ndf,*)
      real*8     dr(6,6),di(6,6), sigr(6),sigi(6),si(nst,*),pr(*),pi(*)
      real*8     jac,jac0,xx, aj1,aj2,aj3
      real*8     bd11,bd21,bd12,bd22,bd13,bd23,bd14,bd24

c     Compute imaginary strains

      do j = 1,6
        epsi(j) = 0.0d0
      end do ! j
      do j = 1,nel
        epsi(1) = epsi(1) + shp(1,j)*ul(1,j)
        epsi(2) = epsi(2) + shp(2,j)*ul(2,j)
        epsi(3) = epsi(3) + shp(3,j)*ul(1,j)
        epsi(4) = epsi(4) + shp(1,j)*ul(2,j) + shp(2,j)*ul(1,j)
      end do ! j

c     Multiply jacobian by radius for axisymmetry

      if(stype.eq.3 .or. stype.eq.8) then
        epsi(3) = epsi(3)/xx
      else
        epsi(3) = 0.0d0
      end if

c     compute imaginary stress

      do j = 1,4
        sigr(j) = 0.0d0
        sigi(j) = 0.0d0
        do k = 1,4
          sigr(j) = sigr(j) + di(j,k)*epsi(k)
          sigi(j) = sigi(j) + dr(j,k)*epsi(k) + di(j,k)*epsr(k)
        end do ! k
      end do ! j

      j1 = 1
      do j = 1,nel

        aj1 = shp(1,j)*jac
        aj2 = shp(2,j)*jac
        aj3 = shp(3,j)*jac0

c       Compute B_trans * D * j * w

        bd11 = aj1*di(1,1) + aj3*di(3,1) + aj2*di(4,1)
        bd12 = aj1*di(1,2) + aj3*di(3,2) + aj2*di(4,2)
        bd13 = aj1*di(1,3) + aj3*di(3,3) + aj2*di(4,3)
        bd14 = aj1*di(1,4) + aj3*di(3,4) + aj2*di(4,4)

        bd21 = aj2*di(2,1) + aj1*di(4,1)
        bd22 = aj2*di(2,2) + aj1*di(4,2)
        bd23 = aj2*di(2,3) + aj1*di(4,3)
        bd24 = aj2*di(2,4) + aj1*di(4,4)

c       Loop over columns (symmetry noted)

        k1 = 1
        do k = 1,nel
          si(j1  ,k1  ) = si(j1  ,k1  ) + (bd11*shp(1,k)
     &                                  +  bd14*shp(2,k)
     &                                  +  bd13*shpr(k))*ctan(1)

          si(j1  ,k1+1) = si(j1  ,k1+1) + (bd12*shp(2,k)
     &                                  +  bd14*shp(1,k))*ctan(1)

          si(j1+1,k1  ) = si(j1+1,k1  ) + (bd21*shp(1,k)
     &                                  +  bd24*shp(2,k)
     &                                  +  bd23*shpr(k))*ctan(1)

          si(j1+1,k1+1) = si(j1+1,k1+1) + (bd22*shp(2,k)
     &                                  +  bd24*shp(1,k))*ctan(1)

          k1 = k1 + ndf
        end do ! k

c       Residual for added real part

        pr(j1  ) = pr(j1  ) + aj1*sigr(1) + aj2*sigr(4) + aj3*sigr(3)
        pr(j1+1) = pr(j1+1) + aj2*sigr(2) + aj1*sigr(4)

c       Residual for imaginary part

        pi(j1  ) = pi(j1  ) - aj1*sigi(1) - aj2*sigi(4) - aj3*sigi(3)
        pi(j1+1) = pi(j1+1) - aj2*sigi(2) - aj1*sigi(4)

        j1 = j1 + ndf
      end do ! j

      end
