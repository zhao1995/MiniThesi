c$Id:$
      subroutine fld3d1tm(d,ul,xl,s,r,ndf,ndm,nst,isw)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    13/01/2010
c       1. Add isw to call of iner3d                        27/01/2011
c       2. Add epsp(6,125) to pass strains to slcn3d        15/08/2013
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: 3-D finite deformation thermo-mech displacement element
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'bdata.h'
      include  'cdata.h'
      include  'cdamag.h'
      include  'defgrd.h'
      include  'elcoor.h'
      include  'eldata.h'
      include  'elengy.h'
      include  'elplot.h'
      include  'eltran.h'
      include  'hdata.h'
      include  'iofile.h'
      include  'ndata.h'
      include  'part0.h'
      include  'plast3f.h'
      include  'pmod2d.h'
      include  'prstrs.h'
      include  'ptdat6.h'
      include  'qudshp.h'
      include  'rdata.h'
      include  'tdata.h'
      include  'comblk.h'

      logical   bflg
      integer   ndf,ndm,nst,isw, i,ii,i1, j,jj,j1,  k
      integer   l, nhv,nn,istrt
      real*8    bdb, rho,rhol,rhom,pote, xipa,wena,xlamd,ha
      real*8    qfact,epp,thlog, qbody
      real*8    d(*),ul(ndf,nen,*),xl(ndm,*),s(nst,*),r(ndf,*)
      real*8    xx(3),bf(3),bt(3,3)
      real*8    be(6), body(3), bep(3,125),bepa(3)
      real*8    aa(7,7,5),dd(7,7,125)
      real*8    sigv(13),sigl(10,125), sigp(10,125), epsp(6,125)
      real*8    fluxv(4),fluxl(3,125), fluxp(3,125),kt(3,3)
      real*8    dtherm(3,3,125), bbc(6)
      real*8    bbd(3,6),bba(3),bbt(4), r1(3,125), tbi(3)
      real*8    dvol(125),weng(125),xr(3,125),ur(3,125), tht(125)

      save

c     MATERIAL DATA

      if(isw.eq.1) then

        xlamd = 0.0d0
        ha    = 0.0d0

c     N.B. Augmenting not allowed in displacement model

      elseif(isw.eq.10) then


c     COMPUTE TANGENT STIFFNESS/RESIDUAL FORCE AND OUTPUTS

      else

c       Get quadrature information

        call quadr3d(d,.true.)

c       Compute shape functions

        do l = 1,lint
          call interp3d(l, xl, ndm,nel)
        end do ! l

c       Compute coordinates

        do l = 1,lint
          do i = 1,3
            xr(i,l) = 0.0d0
            ur(i,l) = 0.0d0
            do j = 1,nel
              xr(i,l) = xr(i,l) + shp3(4,j,l)*xl(i,j)
              ur(i,l) = ur(i,l) + shp3(4,j,l)*ul(i,j,1)
            end do ! j
            ur(i,l) = ur(i,l) + xr(i,l)
          end do ! i
        end do ! l

c       Set history type and size

        nhv   = nint(d(15))
        istrt = nint(d(84))

c       Initialize history variables

        if(isw.eq.14) then

c         Initialize deformation gradients

          call pfinit(f,df,finv,detf, lint)

        else

c         Compute deformation gradient, inverse and determinant.

          call kine3d(shp3,ul,f,finv,df,detf,ndf,nel,nen,lint)

          do l = 1,lint

            dvol(l) = jac(l)*detf(1,l)

c           Compute spatial gradient of temperature

            call tgrad3df(shp3(1,1,l),ul,gradt(1,l),tg(l),ndf,nel)

          end do ! l

        endif

c       Compute Cauchy stresses and spatial tangent tensor at t-n+1

        nn  = 0
        do l = 1,lint

          estore = 0.0d0
          do i = 1,3
            xref(i) = xr(i,l)
            xcur(i) = xr(i,l) + ur(i,l)
          end do ! i
          call modltm(d,f(1,1,l),finv(1,l),df(1,l),detf(1,l),
     &                gradt(1,l),tg(l),hr(nh1+nn),hr(nh2+nn),nhv,istrt,
     &                sigv,fluxv, aa,kt, be,xlamd,ha,.false., isw)
          weng(l) = estore

c         Compute principal stretches

          call pstr3d(be,bep(1,l))

c         Multiply tangent moduli and stresses by volume element.
c         Store time history plot data for element

          k = 9*(l-1)
          do i = 1,6
            tt(i+k)   = sigv(i)
            sigp(i,l) = sigv(i)
            sigl(i,l) = sigv(i)*dvol(l)
            do j = 1,6
              dd(j,i,l) = aa(j,i,1)*dvol(l)*ctan(1)
            end do ! j
            dd(i,7,l) = aa(i,7,1)*dvol(l)*ctan(1)
            dd(7,i,l) = aa(7,i,1)*dvol(l)*ctan(1)
          end do ! i
          sigp( 9,l) = sigv( 9)
          sigp(10,l) = sigv(10)
          k = k + 6
          do i = 1,3
            tt(i+k)    = fluxv(i)
            fluxp(i,l) = fluxv(i)
            fluxl(i,l) = fluxv(i)*dvol(l)
            do j = 1,3
              dtherm(j,i,l) = kt(j,i)*dvol(l)*ctan(1)
            end do ! j
          end do ! i

          dd(7,7,l)   = aa(7,7,1)*dvol(l)*ctan(1)
          tht(l)      = fluxv(4) *dvol(l)

          nn = nn + nhv

        end do ! l

c       ENERGY COMPUTATIONS

        if(isw.eq.13) then

          rhom = 0.0d0
          rhol = 0.0d0
          pote = 0.0d0
          do l = 1,lint

            rho = d(4)*jac(l)
            do i = 1,3
              tbi(i) = 0.0d0
            end do ! i
            do j = 1,nel
              rhol = rhol + (ul(1,j,4)**2 +  ul(2,j,4)**2
     &                    +  ul(3,j,4)**2)*shp3(4,j,l)*rho
              do i = 1,3
                tbi(i) = tbi(i) + shp3(4,j,l)*ul(i,j,4)
              end do ! i
            end do ! j
            rhom = rhom + (tbi(1)**2 + tbi(2)**2 + tbi(3)**2)*rho
            pote = pote +  weng(l)*jac(l)
          end do ! l

c         Accumulate energies

          epl(7) = epl(7) + (rhol + d(7)*(rhom - rhol))*0.5d0
          epl(8) = epl(8) +  pote

c       STIFFNESS AND RESIDUAL

        elseif(isw.eq.3 .or. isw.eq.6) then

c         Compute body forces values

          call sbodyf(d, body)
          bflg = d(4).gt.0.0d0 .and. d(65).gt.0.0d0  ! angular velocity

c         Compute inertia effects: shflg = .true. for eigen shifts

          if(ctan(3).ne.0.0 .and. (ndfo(1).gt.0 .or. shflg)) then
            call iner3d(d,xl,ul(1,1,4),ul(1,1,5),s,r,
     &                  nel,ndf,ndm,nst, isw)
          endif ! ctan(3) test

c         Thermal transients

          if(ctan(2).ne.0.0 .and. (ndfo(4).gt.0 .or. shflg)) then
            call thtrans3d(d,xl,ul(1,1,4),s,r, nel,ndf,ndm,nst)
          endif ! ctan(2) test

          do l = 1,lint

c           Angular velocity: d(4) = rho; d(65) = omega

            do i = 1,3
              bf(i) = 0.0d0
            end do ! i
            if(bflg) then
              call sbodyw(d(4),d(65),ur(1,l), bf,bt, .true.)
              do ii = 1,3
                do jj = 1,3
                  bt(jj,ii) = bt(jj,ii)*jac(l)*ctan(1)
                end do ! jj
              end do ! ii
            endif

c           COMPUTE STRESS DIVERGENCE TERM

c           Compatible internal force.

            qbody = d(66)*jac(l)   ! Heat source term

            do i = 1,nel
              r1(1,i) = shp3(1,i,l)*sigl(1,l)
     &                + shp3(2,i,l)*sigl(4,l)
     &                + shp3(3,i,l)*sigl(6,l)
              r1(2,i) = shp3(1,i,l)*sigl(4,l)
     &                + shp3(2,i,l)*sigl(2,l)
     &                + shp3(3,i,l)*sigl(5,l)
              r1(3,i) = shp3(1,i,l)*sigl(6,l)
     &                + shp3(2,i,l)*sigl(5,l)
     &                + shp3(3,i,l)*sigl(3,l)

c             Mechanical residual

              do j = 1,3
                r(j,i)  = r(j,i) + shp3(4,i,l)*(body(j) + bf(j))*jac(l)
     &                  - r1(j,i)
              end do ! j

c             Thermal residual: Include heating term times abs temp.

              r(4,i) = r(4,i) + shp3(4,i,l)*(qbody - tht(l))
     &                        + shp3(1,i,l)*fluxl(1,l)
     &                        + shp3(2,i,l)*fluxl(2,l)
     &                        + shp3(3,i,l)*fluxl(3,l)
            end do ! i

c           COMPUTE K11

            if(isw.eq.3) then

              i1 = 0
              do i = 1,nel

c               PART 1. - geometric tangent

                j1 = 0
                if(gflag) then
                  do j = 1,nel

c                   Accumulate geometric factor with consistent mass

                    bdb = (r1(1,i)*shp3(1,j,l)
     &                  +  r1(2,i)*shp3(2,j,l)
     &                  +  r1(3,i)*shp3(3,j,l))*ctan(1)

                    do jj = 1,3
                      s(i1+jj,j1+jj) = s(i1+jj,j1+jj) + bdb
                    end do ! jj
                    j1 = j1 + ndf
                  end do ! j
                endif

c               Angular velocity tangent

                if(bflg) then
                  do jj = 1,3
                    do ii = 1,3
                      bdb = shp3(4,i,l)*bt(ii,jj)
                      j1  = 0
                      do j = 1,nel
                        s(i1+ii,j1+jj) = s(i1+ii,j1+jj)
     &                                 + bdb*shp3(4,j,l)
                        j1             = j1 + ndf
                      end do ! j
                    end do ! ii
                  end do ! jj
                endif

c               PART 2. - tangent modulus part (based upon aa-array)

                do jj = 1,6
                  bbd(1,jj) = shp3(1,i,l)*dd(1,jj,l)
     &                      + shp3(2,i,l)*dd(4,jj,l)
     &                      + shp3(3,i,l)*dd(6,jj,l)
                  bbd(2,jj) = shp3(1,i,l)*dd(4,jj,l)
     &                      + shp3(2,i,l)*dd(2,jj,l)
     &                      + shp3(3,i,l)*dd(5,jj,l)
                  bbd(3,jj) = shp3(1,i,l)*dd(6,jj,l)
     &                      + shp3(2,i,l)*dd(5,jj,l)
     &                      + shp3(3,i,l)*dd(3,jj,l)

c                 Thermal coupling with displacements

                  bbc(jj) = shp3(4,i,l)*dd(7,jj,l)

                end do ! jj

c               Thermal conductivity terms

                do jj = 1,3
                  bbt(jj) = shp3(1,i,l)*dtherm(1,jj,l)
     &                    + shp3(2,i,l)*dtherm(2,jj,l)
     &                    + shp3(3,i,l)*dtherm(3,jj,l)
                end do ! jj

c               Thermal heating term

                bbt(4)  = shp3(4,i,l)*dd(7,7,l)

c               Thermal coupling terms

                bba(1) = shp3(1,i,l)*dd(1,7,l)
     &                 + shp3(2,i,l)*dd(4,7,l)
     &                 + shp3(3,i,l)*dd(6,7,l)
                bba(2) = shp3(1,i,l)*dd(4,7,l)
     &                 + shp3(2,i,l)*dd(2,7,l)
     &                 + shp3(3,i,l)*dd(5,7,l)
                bba(3) = shp3(1,i,l)*dd(6,7,l)
     &                 + shp3(2,i,l)*dd(5,7,l)
     &                 + shp3(3,i,l)*dd(3,7,l)

c               Compute tangent stiffness

                j1 = 0
                do j = 1,nel

c                 Mechanical part

                  s(i1+1,j1+1) = s(i1+1,j1+1)
     &                          + bbd(1,1)*shp3(1,j,l)
     &                          + bbd(1,4)*shp3(2,j,l)
     &                          + bbd(1,6)*shp3(3,j,l)
                  s(i1+1,j1+2) = s(i1+1,j1+2)
     &                          + bbd(1,4)*shp3(1,j,l)
     &                          + bbd(1,2)*shp3(2,j,l)
     &                          + bbd(1,5)*shp3(3,j,l)
                  s(i1+1,j1+3) = s(i1+1,j1+3)
     &                          + bbd(1,6)*shp3(1,j,l)
     &                          + bbd(1,5)*shp3(2,j,l)
     &                          + bbd(1,3)*shp3(3,j,l)
                  s(i1+2,j1+1) = s(i1+2,j1+1)
     &                          + bbd(2,1)*shp3(1,j,l)
     &                          + bbd(2,4)*shp3(2,j,l)
     &                          + bbd(2,6)*shp3(3,j,l)
                  s(i1+2,j1+2) = s(i1+2,j1+2)
     &                          + bbd(2,4)*shp3(1,j,l)
     &                          + bbd(2,2)*shp3(2,j,l)
     &                          + bbd(2,5)*shp3(3,j,l)
                  s(i1+2,j1+3) = s(i1+2,j1+3)
     &                          + bbd(2,6)*shp3(1,j,l)
     &                          + bbd(2,5)*shp3(2,j,l)
     &                          + bbd(2,3)*shp3(3,j,l)
                  s(i1+3,j1+1) = s(i1+3,j1+1)
     &                          + bbd(3,1)*shp3(1,j,l)
     &                          + bbd(3,4)*shp3(2,j,l)
     &                          + bbd(3,6)*shp3(3,j,l)
                  s(i1+3,j1+2) = s(i1+3,j1+2)
     &                          + bbd(3,4)*shp3(1,j,l)
     &                          + bbd(3,2)*shp3(2,j,l)
     &                          + bbd(3,5)*shp3(3,j,l)
                  s(i1+3,j1+3) = s(i1+3,j1+3)
     &                          + bbd(3,6)*shp3(1,j,l)
     &                          + bbd(3,5)*shp3(2,j,l)
     &                          + bbd(3,3)*shp3(3,j,l)

c                 Coupling part: Momentum

                  s(i1+1,j1+4) = s(i1+1,j1+4) + bba(1)*shp3(4,j,l)
                  s(i1+2,j1+4) = s(i1+2,j1+4) + bba(2)*shp3(4,j,l)
                  s(i1+3,j1+4) = s(i1+3,j1+4) + bba(3)*shp3(4,j,l)

c                 Coupling part: Thermal

                  s(i1+4,j1+1) = s(i1+4,j1+1) + bbc(1)*shp3(1,j,l)
     &                                        + bbc(4)*shp3(2,j,l)
     &                                        + bbc(6)*shp3(3,j,l)

                  s(i1+4,j1+2) = s(i1+4,j1+2) + bbc(4)*shp3(1,j,l)
     &                                        + bbc(2)*shp3(2,j,l)
     &                                        + bbc(5)*shp3(3,j,l)

                  s(i1+4,j1+3) = s(i1+4,j1+3) + bbc(6)*shp3(1,j,l)
     &                                        + bbc(5)*shp3(2,j,l)
     &                                        + bbc(3)*shp3(3,j,l)

c                 Thermal part

                  s(i1+4,j1+4) = s(i1+4,j1+4)
     &                         + bbt(1)*shp3(1,j,l)
     &                         + bbt(2)*shp3(2,j,l)
     &                         + bbt(3)*shp3(3,j,l)
     &                         + bbt(4)*shp3(4,j,l)
                  j1 = j1 + ndf
                end do ! j
                i1 = i1 + ndf
              end  do ! i

            endif

          end do ! l

c       OUTPUT STRESSES

        elseif(isw.eq.4 .or. isw.eq.8 .or. isw.eq.16) then

          xipa = 0.0d0
          wena = 0.0d0
          do i = 1,3
            xx(i) = 0.0d0
          end do ! i
          do i = 1,13
            sigv(i) = 0.0d0
          end do ! i
          do i = 1,6
            ebig(i) = 0.0d0
            esml(i) = 0.0d0
          end do ! i
          do i = 1,3
            bepa(i) = 0.0d0
          end do ! i
          qfact = 1.d0/dble(lint)
          epp   = 0.0d0

          do l = 1,lint

            do j = 1,3
              do i=1,nel
                xx(j)   = xx(j) + qfact*shp3(4,i,l)*xl(j,i)
              end do ! i
            end do ! j

c           Compute Green and Almansi strains

            call fstrain(f,finv, egreen,ealmansi)

            do i = 1,6
              epsp(i,l) = ealmansi(i)
            end do ! i

c           Move stresses and jacobian for printing

            xipa = xipa + qfact*xipr
            wena = wena + qfact*wengy
            do i = 1,6
              sigv(i) = sigv(i) + qfact*sigp(i,l)
              ebig(i) = ebig(i) + qfact*egreen(i)
              esml(i) = esml(i) + qfact*ealmansi(i)
            end do ! i
            epp      = epp      + qfact*sigp( 9,l)
            sigv(10) = sigv(10) + qfact*sigp(10,l)
            do i = 1,3
              bepa(i) = bepa(i) + 0.5d0*qfact*log(bep(i,l))
            end do ! i

c           Store arrays for plotting.

            sigp(10,l) = 0.0d0
            sigp(9,l)  = 0.0d0

          end do ! l

c         OUTPUT STRESSES

          if (isw .eq. 4) then

            call pstr3d(sigv,sigv(7))

            mct = mct - 2
            if(mct.le.0) then
              write(iow,2001) o,head
              if(ior.lt.0) write(*,2001) o,head
              mct = 50
            endif

            thlog = bepa(1) + bepa(2) + bepa(3)

            write(iow,2002) n,ma,(sigv(ii),ii=1,9),bepa,
     &                      xx,thlog,epp,sigv(10),esml
            if(ior.lt.0) then
              write(*,2002) n,ma,(sigv(ii),ii=1,9),bepa,
     &                      xx,thlog,epp,sigv(10),esml
            end if

c         PROJECT STRESSES ONTO THE NODES FOR PLOTTING

          elseif(isw.eq.8) then

c           Compute current geometry

            call slcn3d(sigp,epsp, r,s, nel)

c         COMPUTE FRACTURE INDICES

          elseif(isw.eq.16) then

            call pfrac3f(f,detf,sigp,weng, dvol, r, ndf,ndm,3)
          endif
        endif
      endif

c     FORMAT STATEMENTS

2001  format(a1,20a4//5x,'Element Stresses'//
     &   '  Elmt  Matl  11-stress  22-stress  33-stress',
     &   '  12-stress  23-stress  13-stress'/4x,'Cauchy  ',
     &   '   1-stress   2-stress   3-stress',
     &   '  log(lam1)  log(lam2)  log(lam3)'/12x,
     &    '   1-coord    2-coord    3-coord',
     &   '       log-J     eff-ep      Yield'/4x,'Almansi '
     &   '  11-Strain  22-Strain  33-Strain  12-Strain',
     &   '  23-Strain  31-Strain')

2002  format(/2i6,1p,6e11.3/(12x,1p,6e11.3:))

      end
