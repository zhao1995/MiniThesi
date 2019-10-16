c$Id:$
      subroutine fld3d1(d,ul,xl,s,r,ndf,ndm,nst,isw)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Remove augmenting statement                      26/01/2007
c       2. Correct computation of 4-node tet inertia;       29/03/2007
c          Change cmom(3) to ac(3); add cfac, lfac definitions
c       3. Add 14/15 node tet computation                   06/07/2007
c       4. Change 'ord' to 'nel' on 'tetshp' calls          05/11/2007
c       5. Revise computation of inertial effects           21/12/2007
c       6. Revise quadrature order for tetrahedra           28/12/2007
c          Modify form for computing intertia
c       7. Add 'qfact' to average stresses correctly        04/03/2008
c       8. Set 'df' to zero for isw.eq.14                   10/10/2008
c       9. Remove 'nel' from call to 'quadr3d'              23/01/2009
c      10. Add 'defgrd' for deformation gradient values     27/04/2009
c      11. Add prints for Almansi strains                   19/10/2009
c      12. Remove xu array                                  26/01/2009
c      13. Modify isw options to restrict computations      11/05/2010
c      14. Move inertial computation after stiffness        01/01/2011
c          Add isw to iner3d to control tangent computation
c      15. Add 'l' to modlfd call                           05/01/2012
c      16. Add average of density for multiscale            10/05/2012
c      17. Add eps to slcn3d call                           01/01/2013
c      18. Use pengy3d to compute momenta and energies      06/08/2013
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: 3-D finite deformation displacement element
c      Remark: This a completely standard mechanical element
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'bdata.h'
      include  'cdata.h'
      include  'cdamag.h'
      include  'debugs.h'
      include  'defgrd.h'
      include  'elcoor.h'
      include  'eldata.h'
      include  'elengy.h'
      include  'elplot.h'
      include  'eltran.h'
      include  'hdata.h'
      include  'iofile.h'
      include  'ndata.h'
      include  'oelmt.h'
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
      real*8    bdb, xipa,wena,xlamd,ha, ta
      real*8    qfact,epp,thlog
      real*8    d(*),ul(ndf,nen,*),xl(ndm,*),s(nst,*),r(ndf,*)
      real*8    xx(3),bf(3),bt(3,3)
      real*8    be(6), body(3), bep(3,125),bepa(3)
      real*8    aa(6,6,5),dd(6,6,125)
      real*8    sigv(13),sigl(10,125), sigp(10,125), epsp(6,125)
      real*8    bbd(3,6),r1(3,125)
      real*8    dvol(125),weng(125),xr(3,125),ur(3,125)

      save

c     MATERIAL DATA

      if(isw.eq.1) then

        xlamd = 0.0d0
        ha    = 0.0d0
        ta    = 0.0d0

c     N.B. Augmenting not allowed in displacement model

      elseif(isw.eq.10) then


c     COMPUTE TANGENT STIFFNESS/RESIDUAL FORCE AND OUTPUTS

      elseif(isw.eq.3 .or. isw.eq.4  .or. isw.eq.6  .or.
     &       isw.eq.8 .or. isw.eq.13 .or. isw.eq.14 .or.
     &       isw.eq.16) then

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
          call modlfd(l,d,f(1,1,l),finv(1,l),df(1,l),detf(1,l),ta,
     &                hr(nh1+nn),hr(nh2+nn),nhv,istrt,aa,sigv,be,
     &                xlamd,ha,.false.,isw)
          weng(l) = estore

c         Compute principal stretches

          call pstr3d(be,bep(1,l))

c         Multiply tangent moduli and stresses by volume element.
c         Store time history plot data for element

          k = 15*(l-1)
          do i = 1,9
            tt(i+6+k) = f(i,3,l) ! Displacement gradient
          end do ! i
          do i = 1,6
            tt(i+k)   = sigv(i)
            sigp(i,l) = sigv(i)
            sigl(i,l) = sigv(i)*dvol(l)
            do j = 1,6
              dd(j,i,l) = aa(j,i,1)*dvol(l)*ctan(1)
            end do ! j
          end do ! i
          sigp( 9,l) = sigv( 9)
          sigp(10,l) = sigv(10)
          nn = nn + nhv

        end do ! l

c       ENERGY COMPUTATIONS

        if(isw.eq.13) then

          call pengy3d(ndm,ndf,nel, weng, d,ul,xl)

c       STIFFNESS AND RESIDUAL

        elseif(isw.eq.3 .or. isw.eq.6) then

c         Compute body forces values

          call sbodyf(d, body)
          bflg = d(4).gt.0.0d0 .and. d(65).gt.0.0d0  ! angular velocity

          do l = 1,lint

c           Compute parameters for multiscale plane strain

            v_avg = v_avg + jac(l)
            v_rho = v_rho + jac(l)*d(4)

c           Angular velocity: d(4) = rho; d(65) = omega

            do i = 1,3
              bf(i) = 0.0d0
            end do ! i
            if(bflg) then
              call sbodyw(d(4),d(65),ur(1,l), bf,bt, .true.)
              do ii = 1,3
                do jj = 1,3
                  bt(jj,ii) = bt(jj,ii)*jac(l)
                end do ! jj
              end do ! ii
            endif

c           COMPUTE STRESS DIVERGENCE TERM

c           Compatible internal force.

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

              do j = 1,3
                r(j,i)  = r(j,i) + shp3(4,i,l)*(body(j) + bf(j))*jac(l)
     &                  - r1(j,i)
              end do ! j
            end do ! i

c           COMPUTE K11

            if(isw.eq.3) then

              i1 = 0
              do i = 1,nel

c               PART 1. - geometric tangenta

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
                end do ! jj

c               Compute tangent stiffness

                j1 = 0
                do j  = 1,nel
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
                  j1 = j1 + ndf
                end do ! j
                i1 = i1 + ndf
              end  do ! i

            endif

          end do ! l

c         Compute inertia effects: shflg = .true. for eigen shifts

          if(ctan(3).ne.0.0 .and. (ndfo(1).gt.0 .or. shflg)) then
            call iner3d(d,xl,ul(1,1,4),ul(1,1,5),s,r,
     &                  nel,ndf,ndm,nst, isw)
          endif ! ctan(3) test

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

            mct = mct - 2
            if(mct.le.0) then
              write(iow,2001) o,head
              if(ior.lt.0) write(*,2001) o,head
              mct = 50
            endif

c           Output quadrature point stresses

            if(qoutfl) then
              do l = 1,lint
                do j = 1,3
                  xx(j) = 0.0d0
                  do i=1,nel
                    xx(j)   = xx(j) + shp3(4,i,l)*xl(j,i)
                  end do ! i
                end do ! j

                call pstr3d(sigp(1,l),sigv(7))

                do i = 1,3
                  bepa(i) = 0.5d0*log(bep(i,l))
                end do ! i
                thlog = bepa(1) + bepa(2) + bepa(3)

                write(iow,2002) n,ma,(sigp(ii,l),ii=1,6),
     &                         (sigv(ii),ii=7,9),bepa,
     &                          xx,thlog,sigp(9,l),sigp(10,l),
     &                         (epsp(ii,l),ii=1,6)
                if(ior.lt.0) then
                  write(*,2002) n,ma,(sigp(ii,l),ii=1,6),
     &                         (sigv(ii),ii=7,9),bepa,
     &                          xx,thlog,sigp(9,l),sigp(10,l),
     &                         (epsp(ii,l),ii=1,6)
                end if
              end do ! l

c           Output averaged stresses

            else

              call pstr3d(sigv,sigv(7))

              thlog = bepa(1) + bepa(2) + bepa(3)

              write(iow,2002) n,ma,(sigv(ii),ii=1,9),bepa,
     &                        xx,thlog,epp,sigv(10),esml
              if(ior.lt.0) then
                write(*,2002) n,ma,(sigv(ii),ii=1,9),bepa,
     &                        xx,thlog,epp,sigv(10),esml
              end if
           endif

c         PROJECT STRESSES ONTO THE NODES FOR PLOTTING

          elseif(isw.eq.8) then

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
