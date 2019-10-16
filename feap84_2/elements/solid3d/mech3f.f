c$Id:$
      subroutine mech3f(d,ul,xl,s,p,ndf,ndm,nst,isw, finc)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Add 'defgrd.h' for deformation gradient values   27/04/2009
c       2. Add prints for Almansi strains                   19/10/2009
c       3. Set arguments to hr(*) to proper ints            17/04/2011
c       4. Add 'l' to modlfd call                           05/01/2012
c       5. Add average of density for multiscale            10/05/2012
c       6. Add eps to slcn3d call                           01/01/2013
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: 3-D finite deformation enhanced strain element
c               F. ARMERO & JCS & RLT 12/91 (Mod. 8/92)
c               Local iteration on enhanced terms - reduced data base

c      Remark: This a completely standard mechanical element. The only
c      ------  differences are:
c              a) Output transfer variables; in this case,
c                 deformation gradients in array `f(9,2,lint)'
c                 f(*,1,*) = F_n+1; f(*,2,*) = F_n

c              b) Reading material properties (isw=1) and
c                 checking for errors (isw=2) is done outside.
c                 However, these functions are never called
c                 and need not be modified.
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

      logical   finc, noconv, bflg
      integer   ndf,ndm,nst,isw, i,ii,i1,i2, nenit, j,jj,j1,j2,  k
      integer   l,lcentr, nn,nhv,ndfi,ninc,nhi,istrt
      real*8    bdb,fact, rho,rhol,rhom,pote, tftau, tempi
      real*8    xx1,xx2,xx3, xipa,wena, temp, tb4, tol1,tolu ,ta
      real*8    d(*),ul(ndf,nen,*),xl(ndm,*),s(nst,*),p(ndf,*)
      real*8    dtf(6),tf(6),tfg(3,3,9), epsp(6,9)
      real*8    f0(3,3,2),shp0(3,8,9),shpi(3,4,9),cmom(3),be(6)
      real*8    aa(6,6,5),dd(6,6,9),sigv(13),sigl(10,9),sigp(10,9)
      real*8    gg(12,24),hh(12,12),bb(12),dui(12), bf(3),bt(3,3)
      real*8    bbd(3,6),r1(3,8), r2(3,3,9), r0(3,8), tbi(3)
      real*8    dvol(9),dvol0(9),weng(9),xr(3,9),ur(3,9), body(3)

      save

c     TEMPORARY TEMPERATURE

      data      ta   / 0.0d0 /
      data      ndfi / 3 /, ninc / 12 /, lcentr / 9 /

c     Go to correct process

      go to(1,1,3,3,1,3,1,3,1,10,1,1,3,3), isw

   1  return

c     COMPUTE TANGENT STIFFNESS AND RESIDUAL FORCE VECTOR

c     Set up parameters

   3  if(finc) then
        l    = -9
        nhi  = 14
      else
        l    = 2
        nhi  = 2
      endif
      call int3d(l,lint,sg3)

c     Compute coordinates

      do l = 1,lint
        call shp3d(sg3(1,l),jac(l),shp3(1,1,l),xl,ndm,nel)
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

c     Degrees of freedom with enhanced modes.

      nhv    = nint(d(15))
      istrt  = nint(d(84))

c     Initialize history variables

      if(isw.eq.14) then

c       Initialize deformation gradient values

        call pfinit(f,df,finv,detf, lint)

        nn =  nhi
        do l = 1,lint

c         Compute Cauchy stresses and spatial tangent tensor at t-n+1

          estore = 0.0d0
          do i = 1,3
            xref(i) = xr(i,l)
            xcur(i) = xr(i,l) + ur(i,l)
          end do ! i
          call modlfd(l,d,f(1,1,l),finv(1,l),df(1,l),detf(1,l),ta,
     &                hr(nh1+nn),hr(nh2+nn),nhv,istrt,aa,sigv,be,
     &                hr(nh2),hr(nh2+1),.false.,isw)
          nn = nn + nhv
        end do ! l
        return
      endif

c     LOCAL ITERATION ON ENHANCED TERMS

c     DO WHILE

      tolu = 1.d-03*tol*rnmax/numel
      noconv = .true.
      nenit  = 0

      do while (noconv)

        do i = 1,ninc
          bb(i) = 0.0d0
          do j = 1,ninc
            hh(j,i) = 0.0d0
          end do ! j
        end do ! i

c       COMPUTE SHAPE FUNCTIONS & DEFORMATION GRADIENT

        if(finc) then
          call shpm3d(shpi,xl,ndm,finc)
        else
          do l = 1,lint
            call shp3d(sg3(1,l),jac(l),shp3(1,1,l),xl,ndm,nel)
          end do ! l
        endif

c       Compute deformation gradient at xi=0

        do i = 1,3
          do j = 1,3
            f0(i,j,1) = 0.0d0
            f0(i,j,2) = 0.0d0
            do k = 1,nel
              f0(i,j,1) = f0(i,j,1) + ul(i,k,1)*shp3(j,k,9)
              f0(i,j,2) = f0(i,j,2) - ul(i,k,2)*shp3(j,k,9)
            end do ! k
            f0(i,j,2) = f0(i,j,2) + f0(i,j,1)
          end do ! j
          f0(i,i,1) = f0(i,i,1) + 1.0d0
          f0(i,i,2) = f0(i,i,2) + 1.0d0
        end do ! i

        nn = nhi
        do l = 1,lint

c         Compute deformation gradient, inverse and determinant.

          call kine3f(shp3(1,1,l),shpi(1,1,l),ul,hr(nh2+2),hr(nh1+2),
     &                f(1,1,l),f0,finv(1,l),df(1,l),detf(1,l),
     &                ndf,ndfi,nel,nen,fact,finc)

c         Modify displacement mode shape functions

          do i = 1,3
            do j = 1,8
              shp0(i,j,l) = shp3(i,j,9)
              shp3(i,j,l) = shp3(i,j,l) + fact * shp3(i,j,9)
            end do ! j
          end do ! i

c         Push forward enhanced shape functions (EXCEPT LAST)

          do i = 1,3
             call pushv3f (finv(1,l),shpi(1,i,l))
          end do ! i

c         Push forward standard (modified) shape functions

          do i = 1,nel
           call pushv3f (finv(1,l),shp3(1,i,l))
          end do ! i

        end do ! l

c       TRANSFER FOR ENERGY COMPUTATION

        if(isw.eq.13) go to 13

c       LOOP OVER GAUSS POINTS

        do l = 1,lint

          dvol0(l) = jac(l)*sg3(4,l)
          dvol(l)  = dvol0(l)*detf(1,l)

c         Compute Cauchy stresses and spatial tangent tensor at t-n+1

          estore = 0.0d0
          do i = 1,3
            xref(i) = xr(i,l)
            xcur(i) = xr(i,l) + ur(i,l)
          end do ! i
          call modlfd(l,d,f(1,1,l),finv(1,l),df(1,l),detf(1,l),ta,
     &                hr(nh1+nn),hr(nh2+nn),nhv,istrt,aa,sigv,be,
     &                hr(nh2),hr(nh2+1),.false.,isw)
          weng(l) = estore

c         Multiply tangent moduli and stresses by volume element.

c         Store time history plot data for element

          k = 6*(l-1)
          do i = 1,6
            tt(i+k)   = sigv(i)
            sigp(i,l) = sigv(i)
            sigv(i)   = sigv(i)*dvol(l)
            sigl(i,l) = sigv(i)
            do j = 1,6
              aa(i,j,1)   = aa(i,j,1)*dvol(l)*ctan(1)
              dd(i,j,l) = aa(i,j,1)
            end do ! j
          end do ! i

c         Enhanced mode only

          if(finc) then

c           Enhanced internal force.

            i2 = 0
            do i = 1,3
              r2(1,i,l) = shpi(1,i,l)*sigv(1)
     &                  + shpi(2,i,l)*sigv(4)
     &                  + shpi(3,i,l)*sigv(6)
              r2(2,i,l) = shpi(1,i,l)*sigv(4)
     &                  + shpi(2,i,l)*sigv(2)
     &                  + shpi(3,i,l)*sigv(5)
              r2(3,i,l) = shpi(1,i,l)*sigv(6)
     &                  + shpi(2,i,l)*sigv(5)
     &                  + shpi(3,i,l)*sigv(3)

              bb(i2+1) = bb(i2+1) - r2(1,i,l)
              bb(i2+2) = bb(i2+2) - r2(2,i,l)
              bb(i2+3) = bb(i2+3) - r2(3,i,l)
              i2 = i2 + ndfi
            end do ! i

c           Last enhanced mode internal force

            do i = 1,3
              tfg(i,1,l) = f0(i,1,1)*finv(1,l)
     &                   + f0(i,2,1)*finv(2,l)
     &                   + f0(i,3,1)*finv(3,l)
              tfg(i,2,l) = f0(i,1,1)*finv(4,l)
     &                   + f0(i,2,1)*finv(5,l)
     &                   + f0(i,3,1)*finv(6,l)
              tfg(i,3,l) = f0(i,1,1)*finv(7,l)
     &                   + f0(i,2,1)*finv(8,l)
     &                   + f0(i,3,1)*finv(9,l)
            end do ! i

            tf(1) = tfg(1,1,l)
            tf(2) = tfg(2,2,l)
            tf(3) = tfg(3,3,l)
            tf(4) = tfg(1,2,l) + tfg(2,1,l)
            tf(5) = tfg(2,3,l) + tfg(3,2,l)
            tf(6) = tfg(3,1,l) + tfg(1,3,l)

            tftau = tf(1)*sigv(1) + tf(2)*sigv(2) + tf(3)*sigv(3)
     &            + tf(4)*sigv(4) + tf(5)*sigv(5) + tf(6)*sigv(6)

            bb(10) = bb(10) - tftau * shpi(1,4,l)
            bb(11) = bb(11) - tftau * shpi(2,4,l)
            bb(12) = bb(12) - tftau * shpi(3,4,l)

c           COMPUTE K22 (hh(9,9) = K22)

c           PART 1. - geometric part.

            if(gflag) then
              i2 = 0
              do i = 1,3
                j2 = 0
                do j = 1,i
                  bdb     = r2(1,i,l)*shpi(1,j,l)
     &                    + r2(2,i,l)*shpi(2,j,l)
     &                    + r2(3,i,l)*shpi(3,j,l)
                  do jj = 1,3
                    hh(i2+jj,j2+jj) = hh(i2+jj,j2+jj) + bdb
                  end do ! jj
                  j2 = j2 + ndfi
                end do ! j

                i2 = i2 + ndfi
              end do ! i
            endif

c           PART 2. - tangent modulus part (based upon aa-array)

            i2 = 0
            do i  = 1,3

c             Compute bmat-t * aa * dvol

              do jj = 1,6
                bbd(1,jj) = shpi(1,i,l)*aa(1,jj,1)
     &                    + shpi(2,i,l)*aa(4,jj,1)
     &                    + shpi(3,i,l)*aa(6,jj,1)
                bbd(2,jj) = shpi(1,i,l)*aa(4,jj,1)
     &                    + shpi(2,i,l)*aa(2,jj,1)
     &                    + shpi(3,i,l)*aa(5,jj,1)
                bbd(3,jj) = shpi(1,i,l)*aa(6,jj,1)
     &                    + shpi(2,i,l)*aa(5,jj,1)
     &                    + shpi(3,i,l)*aa(3,jj,1)
              end do ! jj

c             Compute tangent stiffness

              j2 = 0
              do j  = 1,i
                do ii = 1,3
                  hh(i2+ii,j2+1) = hh(i2+ii,j2+1)
     &                           + bbd(ii,1)*shpi(1,j,l)
     &                           + bbd(ii,4)*shpi(2,j,l)
     &                           + bbd(ii,6)*shpi(3,j,l)
                  hh(i2+ii,j2+2) = hh(i2+ii,j2+2)
     &                           + bbd(ii,4)*shpi(1,j,l)
     &                           + bbd(ii,2)*shpi(2,j,l)
     &                           + bbd(ii,5)*shpi(3,j,l)
                  hh(i2+ii,j2+3) = hh(i2+ii,j2+3)
     &                           + bbd(ii,6)*shpi(1,j,l)
     &                           + bbd(ii,5)*shpi(2,j,l)
     &                           + bbd(ii,3)*shpi(3,j,l)
                end do ! ii
                j2 = j2 + ndfi
              end do ! j

              i2 = i2 + ndfi
            end do ! i

c           LAST ENHANCED MODE

            do i = 1,6
              dtf(i) = 0.0d0
              do j = 1,6
                dtf(i) = dtf(i) + tf(j)*aa(j,i,1)
              end do ! j
            end do ! i

c           Form K32 (Material and Geometric)

            j2 = 0
            do j = 1,3

              tbi(1) = shpi(1,j,l)*dtf(1)
     &               + shpi(2,j,l)*dtf(4)
     &               + shpi(3,j,l)*dtf(6)
              tbi(2) = shpi(2,j,l)*dtf(2)
     &               + shpi(1,j,l)*dtf(4)
     &               + shpi(3,j,l)*dtf(5)
              tbi(3) = shpi(3,j,l)*dtf(3)
     &               + shpi(2,j,l)*dtf(5)
     &               + shpi(1,j,l)*dtf(6)

              do jj = 1,3

                bdb = tfg(jj,1,l)*r2(1,j,l)
     &              + tfg(jj,2,l)*r2(2,j,l)
     &              + tfg(jj,3,l)*r2(3,j,l) + tbi(jj)

                hh(10,j2+jj) = hh(10,j2+jj) + shpi(1,4,l)*bdb
                hh(11,j2+jj) = hh(11,j2+jj) + shpi(2,4,l)*bdb
                hh(12,j2+jj) = hh(12,j2+jj) + shpi(3,4,l)*bdb

              end do ! jj

              j2 = j2 + ndfi

            end do ! j

c           Form K33 (Material and Geometric)

            tftau = 0.0d0
            do i = 1,6
              tftau = tftau + tf(i)*dtf(i)
            end do ! i

            do i = 1,6
              dtf(i) = 0.0d0
            end do ! i

            do i = 1,3
              dtf(1) = dtf(1) + tfg(i,1,l)*tfg(i,1,l)
              dtf(2) = dtf(2) + tfg(i,2,l)*tfg(i,2,l)
              dtf(3) = dtf(3) + tfg(i,3,l)*tfg(i,3,l)
              dtf(4) = dtf(4) + tfg(i,1,l)*tfg(i,2,l)
              dtf(5) = dtf(5) + tfg(i,2,l)*tfg(i,3,l)
              dtf(6) = dtf(6) + tfg(i,3,l)*tfg(i,1,l)
            end do ! i

            bdb = dtf(1)*sigv(1) + dtf(2)*sigv(2) + dtf(3)*sigv(3)
     &          +(dtf(4)*sigv(4) + dtf(5)*sigv(5) + dtf(6)*sigv(6))*2.d0
     &          + tftau

            do jj = 1,3
              temp  = shpi(jj,4,l) * bdb
              do ii = jj,3
                hh(9+ii,9+jj) = hh(9+ii,9+jj) + shpi(ii,4,l) * temp
              end do ! ii
            end do ! jj
          endif

          nn = nn + nhv
        end do ! l

c       Compute symmetric part of hh

        if(finc) then

          do j = 1,ninc
            do i = 1,j
              hh(i,j) = hh(j,i)
            end do ! i
          end do ! j

c         Invert diagonal array

          call invert(hh,ninc,ninc)

c         Compute K22-inverse * b

          do i = 1,ninc
            dui(i) = 0.d0
            do j = 1,ninc
              dui(i) = dui(i) + hh(i,j)*bb(j)
            end do ! j
          end do ! i

c         Check convergence

          tol1  = 0.0d0
          do i = 1,ninc
            tol1 = tol1 + bb(i)*dui(i)
          end do ! i

          if(abs(tol1).le.tolu .and. nenit.ge.1) then
            noconv = .false.
          endif

          nenit = nenit + 1
          if(nenit.ge.3 .or. tolu.eq.0.0d0) then
            noconv = .false.
          endif

c         Update enhanced parameters

          if(noconv) then
            do i = 2,ninc + 1
              hr(nh2+i) = hr(nh2+i) + dui(i-1)
            end do ! i
          endif

c       Displacement solution

        else
          noconv = .false.
        endif

      end do ! while

      if(isw.eq.4 .or. isw.eq.8 .or. isw.eq.16) go to 4

c     LOOP OVER GAUSS POINTS

      do i = 1,24
        do j = 1,ninc
          gg(j,i) = 0.0d0
        end do ! j
      end do ! i

c     Set body force

      call sbodyf(d, body)
      bflg = d(4).gt.0.0d0 .and. d(65).gt.0.0d0

      do l = 1,lint

c       Angular velocity: d(4) = rho; d(65) = omega

        do i = 1,3
          bf(i) = 0.0d0
        end do ! i
        if(bflg) then
          call sbodyw(d(4),d(65),ur(1,l), bf,bt, .true.)
        endif

        dvol0(l) = jac(l)*sg3(4,l)
        dvol(l)  = dvol0(l)*detf(1,l)

c       Multiply tangent moduli and stresses by volume element.

        do i = 1,6
          sigv(i) = sigl(i,l)
          do j = 1,6
            aa(i,j,1) = dd(i,j,l)
          end do ! j
        end do ! i

c       Compute change in momentum

        if(d(7).ge.0.0d0) then
          rho  = d(4)*dvol0(l)
          rhom = rho * d(7)
          rhol = rho - rhom
        else
          rhom = 0.0d0
          rhol = 0.0d0
        endif
        do i = 1,3
          cmom(i) = 0.0d0
          do j = 1,nel
            cmom(i) = cmom(i) + shp3(4,j,l)*ul(i,j,5)
          end do ! j
          cmom(i) = rhom*cmom(i)
        end do ! i

c       COMPUTE STRESS DIVERGENCE TERM

c       Compatible internal force.

        do i = 1,nel
          r1(1,i) = shp3(1,i,l)*sigv(1)
     &            + shp3(2,i,l)*sigv(4)
     &            + shp3(3,i,l)*sigv(6)
          r1(2,i) = shp3(1,i,l)*sigv(4)
     &            + shp3(2,i,l)*sigv(2)
     &            + shp3(3,i,l)*sigv(5)
          r1(3,i) = shp3(1,i,l)*sigv(6)
     &            + shp3(2,i,l)*sigv(5)
     &            + shp3(3,i,l)*sigv(3)

          do j = 1,3
            p(j,i)  = p(j,i) + shp3(4,i,l)*(dvol0(l)*(body(j) + bf(j))
     &              - (cmom(j)+rhol*ul(j,i,5))) - r1(j,i)
          end do ! j
        end do ! i


c       COMPUTE K11 (s(nst,nst) = K11)

        if(isw.eq.3) then

c         Compute parameters for multiscale plane strain

          v_avg = v_avg + dvol0(l)
          v_rho = v_rho + dvol0(l)*d(4)

c         PART 1. - geometric and inertial parts.

          i1 = 0
          do i = 1,nel

c           Inertia Part

            if(d(7).ge.0.0 .and. (ndfo(1).gt.0 .or. shflg)) then
              rhom  = rho*ctan(3)
            else
              rhom = 0.0d0
            endif

            tb4   = shp3(4,i,l)*rhom

c           Lumped mass factor

            rhol  = tb4*(1.d0 - d(7))
            do jj = 1,3
              s(i1+jj,i1+jj) = s(i1+jj,i1+jj) + rhol
            end do ! jj

c           Consistent mass factor

            tb4 = tb4*d(7)

            j1 = 0
            if(gflag) then
              do j = 1,i

c               Accumulate geometric factor with consistent mass

                bdb = (r1(1,i)*shp3(1,j,l) + r1(2,i)*shp3(2,j,l)
     &              +  r1(3,i)*shp3(3,j,l))*ctan(1) + tb4*shp3(4,j,l)

                do jj = 1,3
                  s(i1+jj,j1+jj) = s(i1+jj,j1+jj) + bdb
                end do ! jj
                j1 = j1 + ndf
              end do ! j
            else
              do j = 1,i
                do jj = 1,3
                  s(i1+jj,j1+jj) = s(i1+jj,j1+jj) + tb4*shp3(4,j,l)
                end do ! jj
                j1 = j1 + ndf
              end do ! j
            endif

            i1 = i1 + ndf
          end do ! i

c         PART 2. - tangent modulus part (based upon aa-array)

          i1 = 0
          do i  = 1,nel

c           Angular velocity tangent

            if(bflg) then
              do jj = 1,3
                do ii = 1,3
                  bdb = shp3(4,i,l)*bt(ii,jj)
                  j1  = 0
                  do j = 1,i
                    s(i1+ii,j1+jj) = s(i1+ii,j1+jj) + bdb*shp3(4,j,l)
                    j1             = j1 + ndf
                  end do ! j
                end do ! ii
              end do ! jj
            endif

c           Compute bmat-t * aa * dvol

            do jj = 1,6
              bbd(1,jj) = shp3(1,i,l)*aa(1,jj,1)
     &                  + shp3(2,i,l)*aa(4,jj,1)
     &                  + shp3(3,i,l)*aa(6,jj,1)
              bbd(2,jj) = shp3(1,i,l)*aa(4,jj,1)
     &                  + shp3(2,i,l)*aa(2,jj,1)
     &                  + shp3(3,i,l)*aa(5,jj,1)
              bbd(3,jj) = shp3(1,i,l)*aa(6,jj,1)
     &                  + shp3(2,i,l)*aa(5,jj,1)
     &                  + shp3(3,i,l)*aa(3,jj,1)
            end do ! i

c           Compute tangent stiffness

            j1 = 0
            do j  = 1,i
              do ii = 1,3
                s(i1+ii,j1+1) = s(i1+ii,j1+1) + bbd(ii,1)*shp3(1,j,l)
     &                                        + bbd(ii,4)*shp3(2,j,l)
     &                                        + bbd(ii,6)*shp3(3,j,l)
                s(i1+ii,j1+2) = s(i1+ii,j1+2) + bbd(ii,4)*shp3(1,j,l)
     &                                        + bbd(ii,2)*shp3(2,j,l)
     &                                        + bbd(ii,5)*shp3(3,j,l)
                s(i1+ii,j1+3) = s(i1+ii,j1+3) + bbd(ii,6)*shp3(1,j,l)
     &                                        + bbd(ii,5)*shp3(2,j,l)
     &                                        + bbd(ii,3)*shp3(3,j,l)
              end do ! ii
              j1 = j1 + ndf
            end do ! j

            i1 = i1 + ndf
          end do ! i

c         Enhanced mode calculation

c         COMPUTE K21 ( gg(9,24) = K21)

c         PART 1. - geometric part.

          if(finc .and. l.ne.lcentr) then

            i2 = 0
            do i = 1,3

              j2 = 0
              do j = 1,nel
                bdb = (r2(1,i,l)*shp3(1,j,l)
     &              +  r2(2,i,l)*shp3(2,j,l)
     &              +  r2(3,i,l)*shp3(3,j,l))*ctan(1)
                do jj = 1,3
                  gg(i2+jj,j2+jj) = gg(i2+jj,j2+jj) + bdb
                end do ! jj
                j2 = j2 + ndfi
              end do ! j

              i2 = i2 + ndfi
            end do ! i

c           PART 2. - tangent modulus part (based upon aa-array)

            i2 = 0
            do i  = 1,3

c             Compute bmat-t * aa * dvol

              do jj = 1,6
                bbd(1,jj) = shpi(1,i,l)*aa(1,jj,1)
     &                    + shpi(2,i,l)*aa(4,jj,1)
     &                    + shpi(3,i,l)*aa(6,jj,1)
                bbd(2,jj) = shpi(1,i,l)*aa(4,jj,1)
     &                    + shpi(2,i,l)*aa(2,jj,1)
     &                    + shpi(3,i,l)*aa(5,jj,1)
                bbd(3,jj) = shpi(1,i,l)*aa(6,jj,1)
     &                    + shpi(2,i,l)*aa(5,jj,1)
     &                    + shpi(3,i,l)*aa(3,jj,1)
              end do ! jj

c             Compute tangent stiffness

              j2 = 0
              do j  = 1,nel
                do ii = 1,3
                  gg(i2+ii,j2+1) = gg(i2+ii,j2+1)
     &                           + bbd(ii,1)*shp3(1,j,l)
     &                           + bbd(ii,4)*shp3(2,j,l)
     &                           + bbd(ii,6)*shp3(3,j,l)
                  gg(i2+ii,j2+2) = gg(i2+ii,j2+2)
     &                           + bbd(ii,4)*shp3(1,j,l)
     &                           + bbd(ii,2)*shp3(2,j,l)
     &                           + bbd(ii,5)*shp3(3,j,l)
                  gg(i2+ii,j2+3) = gg(i2+ii,j2+3)
     &                           + bbd(ii,6)*shp3(1,j,l)
     &                           + bbd(ii,5)*shp3(2,j,l)
     &                           + bbd(ii,3)*shp3(3,j,l)
                end do ! ii
                j2 = j2 + ndfi
              end do ! j

              i2 = i2 + ndfi
            end do ! i

c           LAST ENHANCED MODE

c           Form K31 (Material and Geometric)

            tf(1) = tfg(1,1,l)
            tf(2) = tfg(2,2,l)
            tf(3) = tfg(3,3,l)
            tf(4) = tfg(1,2,l) + tfg(2,1,l)
            tf(5) = tfg(2,3,l) + tfg(3,2,l)
            tf(6) = tfg(3,1,l) + tfg(1,3,l)

            do i = 1,6
              dtf(i) = 0.0d0
              do j = 1,6
                dtf(i) = dtf(i) + tf(j)*aa(j,i,1)
              end do ! j
            end do ! i

            do i = 1,nel
              call pushv3f (finv(1,l),shp0(1,i,l))
            end do ! i

            j2 = 0
            do j = 1,8
              r0(1,j) = shp0(1,j,l)*sigv(1) + shp0(2,j,l)*sigv(4)
     &                + shp0(3,j,l)*sigv(6) + shp3(1,j,l)*dtf(1)
     &                + shp3(2,j,l)*dtf(4)  + shp3(3,j,l)*dtf(6)

              r0(2,j) = shp0(1,j,l)*sigv(4) + shp0(2,j,l)*sigv(2)
     &                + shp0(3,j,l)*sigv(5) + shp3(2,j,l)*dtf(2)
     &                + shp3(1,j,l)*dtf(4)  + shp3(3,j,l)*dtf(5)

              r0(3,j) = shp0(1,j,l)*sigv(6) + shp0(2,j,l)*sigv(5)
     &                + shp0(3,j,l)*sigv(3) + shp3(3,j,l)*dtf(3)
     &                + shp3(2,j,l)*dtf(5)  + shp3(1,j,l)*dtf(6)

              do jj = 1,3

                bdb = (tfg(jj,1,l)*r1(1,j) + tfg(jj,2,l)*r1(2,j)
     &              +  tfg(jj,3,l)*r1(3,j) +  r0(jj,j))*ctan(1)

                gg(10,j2+jj) = gg(10,j2+jj) + shpi(1,4,l)*bdb
                gg(11,j2+jj) = gg(11,j2+jj) + shpi(2,4,l)*bdb
                gg(12,j2+jj) = gg(12,j2+jj) + shpi(3,4,l)*bdb

              end do ! jj

              j2 = j2 + ndfi
            end do ! j

          endif

        endif

c       END OF LOOP OVER GAUSS POINTS

      end do ! l

c     Condense internal modes

      if(finc) call stcon3f(hh,gg,dui,ndf,ndfi,nel,nst,ninc,s,p)

c     Compute upper part by symmetry

      do j = 1,nst
        do i = 1,j
          s(i,j) = s(j,i)
        end do ! i
      end do ! j
      return

c     OUTPUT STRESSES

   4  xx1  = 0.d0
      xx2  = 0.d0
      xx3  = 0.d0
      xipa = 0.0d0
      wena = 0.0d0
      do i = 1,13
        sigv(i) = 0.0d0
      end do ! i
      do i = 1,6
        ebig(i) = 0.0d0
        esml(i) = 0.0d0
      end do ! i

      do l = 1,lint

        do i=1,nel
          tempi = 0.125d0*shp3(4,i,l)
          xx1   = xx1 + tempi*xl(1,i)
          xx2   = xx2 + tempi*xl(2,i)
          xx3   = xx3 + tempi*xl(3,i)
        end do ! i

c       Compute Green & Almansi strains

        call fstrain(f(1,1,l),finv(1,l), egreen,ealmansi)

        do i = 1,6
          epsp(i,l) = ealmansi(i)
        end do ! i

c       Move stresses and jacobian for printing

        xipa = xipa + 0.125d0*xipr
        wena = wena + 0.125d0*wengy
        do i = 1,6
          sigv(i) = sigv(i) + 0.125d0*sigp(i,l)
          ebig(i) = ebig(i) + 0.125d0*egreen(i)
          esml(i) = esml(i) + 0.125d0*ealmansi(i)
        end do ! i

c       Store arrays for plotting.

        sigp(10,l) = 0.0d0
        sigp(9,l)  = 0.0d0

      end do ! l

c     OUTPUT STRESSES

      if (isw .eq. 4) then


        mct = mct - 2
        if(mct.le.0) then
          write(iow,2001) o,head
          if(ior.lt.0) then
            write(*,2001) o,head
          endif
          mct = 50
        endif

c       Output results at quadrature point

        if(qoutfl) then
          xx1 = 0.0d0
          xx2 = 0.0d0
          xx3 = 0.0d0
          do i=1,nel
            xx1   = xx1 + shp3(4,i,l)**xl(1,i)
            xx2   = xx2 + shp3(4,i,l)**xl(2,i)
            xx3   = xx3 + shp3(4,i,l)**xl(3,i)
          end do ! i
          call pstr3d(sigp(1,l),sigv(7))
          write(iow,2002) n,ma,(sigp(ii,l),ii=1,6),(sigv(ii),ii=7,9),
     &                    xx1,xx2,xx3,(epsp(ii,l),ii=1,6)
          if(ior.lt.0) then
            write(*,2002) n,ma,(sigp(ii,l),ii=1,6),(sigv(ii),ii=7,9),
     &                    xx1,xx2,xx3,(epsp(ii,l),ii=1,6)
          end if

c       Output averaged results only

        else
          call pstr3d(sigv,sigv(7))
          write(iow,2002) n,ma,(sigv(ii),ii=1,9),xx1,xx2,xx3,esml
          if(ior.lt.0) then
            write(*,2002) n,ma,(sigv(ii),ii=1,9),xx1,xx2,xx3,esml
          end if
        endif

c     PROJECT STRESSES ONTO THE NODES FOR PLOTTING

      elseif(isw.eq.8) then

        call slcn3d(sigp,epsp, p,s, nel)

c     Compute fracture indices

      elseif(isw.eq.16) then

        call pfrac3f(f,detf,sigp,weng, dvol, p, ndf,ndm,3)

      endif
      return

c     AUGMENTED LAGRANGIAN UPDATE FOR NESTED ITERATION

  10  if(hr(nh2).ne.0.0d0) then
        hr(nh2+1) = hr(nh2+1) + d(1)*log(abs(hr(nh2)))
      endif
      return

c     ENERGY COMPUTATIONS
  13  rhom = 0.0d0
      rhol = 0.0d0
      pote = 0.0d0
      nn = nhi
      do l = 1,lint

        dvol0(l) = jac(l)*sg3(4,l)
        rho      = d(4)*dvol0(l)

c       Compute Cauchy stresses and spatial tangent tensor at t-n+1

        estore = -1.0d0
        do i = 1,3
          xref(i) = xr(i,l)
          xcur(i) = xr(i,l) + ur(i,l)
        end do ! i
        call modlfd(l,d,f(1,1,l),finv(1,l),df(1,l),detf(1,l),ta,
     &              hr(nh1+nn),hr(nh2+nn),nhv,istrt,aa,sigv,be,
     &              hr(nh2),hr(nh2+1),.false.,isw)

        do i = 1,3
          tbi(i) = 0.0d0
        end do ! i
        do j = 1,nel
          rhol = rhol + (ul(1,j,4)**2 + ul(2,j,4)**2 + ul(3,j,4)**2)*
     &                  shp3(4,j,l)*rho
          do i = 1,3
            tbi(i) = tbi(i) + shp3(4,j,l)*ul(i,j,4)
          end do ! i
        end do ! j
        rhom = rhom + (tbi(1)**2 + tbi(2)**2 + tbi(3)**2)*rho
        pote = pote +  estore*dvol0(l)
        nn   = nn + nhv
      end do ! l

c     Accumulate energies

      epl(7) = epl(7) + (rhol + d(7)*(rhom - rhol))*0.5d0
      epl(8) = epl(8) +  pote

      return

c     FORMAT STATEMENTS

2001  format(a1,20a4//5x,'element stresses'//'  elmt  matl',
     &   '  11-stress  22-stress  33-stress  12-stress',
     &   '  23-stress  13-stress'/12x,
     &   '   1-stress   2-stress   3-stress',
     &   '    1-coord    2-coord    3-coord '/12x,
     &   '  11-Almansi 22-Almansi 33-Almansi 12-Almansi',
     &    ' 23-Almansi 31-Almansi')

2002  format(/2i6,1p,6e11.3/(12x,1p,6e11.3:))

      end
