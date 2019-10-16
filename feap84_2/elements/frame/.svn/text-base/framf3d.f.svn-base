c$Id:$
      subroutine framf3d(d,ul,xl,ix,s,p,ndf,ndm,nst,isw)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Move set of rot0/rot1 from loop near line 775.   04/16/2007
c       2. Multiply lumped mass by ctan(3) for transients.  30/09/2008
c          Correct indexing on lumper residual and mass.
c       3. Remove call to pltln -- not required for plots   14/10/2008
c       4. Dimension tmp(4) in strefb for quaternion stores 06/09/2010
c-----[--.----+----.----+----.-----------------------------------------]
c     Geometrically nonlinear elastic 3D beam element
c     (Biot stress/strain formulation)
c     Programmed by:                               Date        Version
c       A. Ibrahimbegovic & M. Al Mikdad           1996            1.0
c     Modified by:
c       E. Kasper & R. Taylor                      1998            1.1

c     Pure displ./rotations formulation 2 node elem.
c                              (reduced integration)
c     Rotation parameters = rotation vector
c     Arbitrary position in space
c     Decomposition of rotation field

c     Reference: A. Ibrahimbegovic and M. Al Mikdad, IJNME, [1998],
c                Vol. 41, pp. 781-814
c-----[--.----+----.----+----.-----------------------------------------]
      implicit none

      include 'bm3cn.h'
      include 'bdata.h'
      include 'cdata.h'
      include 'crotas.h'
      include 'eldata.h'
      include 'eltran.h'
      include 'erotas.h'
      include 'ddata.h'
      include 'iofile.h'
      include 'hdata.h'
      include 'part0.h'
      include 'pview.h'
      include 'rdata.h'
      include 'rdat1.h'
      include 'pointer.h'
      include 'comblk.h'

      integer  l,lint,nhv,im,j,jj,jm,k,kk,li,lj,lk
      integer  ix(*),ndf,ndm,nst,isw
      real*8   xjac,xma,xmr,xsj,xsjm,xsjc,xsjk
      real*8   xx,yy,zz, shp11,shp12,shp21,shp22
      real*8   a1(3),v1(3),tm(3),ula(3),urota(3),utaro(3)
      real*8   dya(3,3),tmat(3,3),shp(2,3),sg(2,3),w1(3)
      real*8   xmat(3,3),matm(3,3),matc(3,3),matk(3,3),xlam(3,3)
      real*8   xm(3),xn(3),ugrad(3),omg(3),th(3),thd(3),fi(6,2)
      real*8   dbn(3,3),bn(3,3),dbm(3,3),bm(3,3),sgn(6,6),sgm(6,6)
      real*8   cm(3,3),cn(3,3),cnm(3,3),tmpn(3,3),te(3,3),sa(3,3)
      real*8   d(*),ul(ndf,nen,*),xl(ndm,*),s(nst,*),p(*)

      save

c     Input material properties

      if(isw.eq.1) then

        lint = 0

c       Set history storage - no incompatible modes

        nhv  = 3 + 4
        nh1  = nh1 + nhv

c       Set rotational update type

        rotyp = -3

c     Check element for errors in input data

      elseif(isw.eq.2) then

c     Compute tangent stiffness and Residual

      elseif(isw.eq.3 .or. isw.eq.6) then

        if(nint(d(182)).gt.0) then
          lint = nel
          call int1dn(lint, sg)
        else
          lint = nel - 1
          call int1d(lint, sg)
        endif

c       Add body force to residual

        call fbody3d(d,xl, p, ndm,ndf, isw)

c       Stiffness and residual computation

        nhv = nint(d(15))
        do l = 1,lint

          call bmshp (sg(1,l),xl,shp,ndm,nel,xjac)
          xsj = sg(2,l)*xjac
          call strefb(hr(nh1),hr(nh2),hr(nh1+3),hr(nh2+3),ul,d,
     &                shp,ugrad,omg,th,thd,xl,te,xm,xn,cn,cm,cnm,xlam,
     &                ndm,ndf,nen, hr(nh1+7),hr(nh2+7),nhv,isw)

          if (isw.eq.3) then
            call geomfb(xn,xm,ugrad,omg,th,thd,te,sa,sgm,sgn)

            xsjk = xsj * ctan(1)
            do jm = 1,6
              do im = 1,6
                sgn(im,jm) = sgn(im,jm)*xsjk
                sgm(im,jm) = sgm(im,jm)*xsjk
              end do ! im
            end do ! jm
            do jm = 1,3
              do im = 1,3
                cn(im,jm)  = cn(im,jm)*xsjk
                cm(im,jm)  = cm(im,jm)*xsjk
                cnm(im,jm) = cnm(im,jm)*xsjk
              end do ! im
            end do ! jm
          endif

          do im = 1,3
            w1(im) = d(10+im)*xsj
            xn(im) = xn(im)*xsj
            xm(im) = xm(im)*xsj
          end do ! im

c         Form residual & K matrix

c         loop over rows

          jj = 0
          do j = 1,nel

c           compute gravity load and residual

            call bmbnfb(bm,bn,omg,ugrad,te,sa,shp(2,j),shp(1,j))
            do im = 1,3
              p(jj+im) = p(jj+im) + w1(im)*shp(2,j)
            end do ! im
            do im = 1,3
              p(jj+im) = p(jj+im) - shp(1,j)*xn(im)
              p(jj+im+3) = p(jj+im+3)
     &                   - bm(1,im)*xm(1)
     &                   - bm(2,im)*xm(2)
     &                   - bm(3,im)*xm(3)
     &                   - bn(1,im)*xn(1)
     &                   - bn(2,im)*xn(2)
     &                   - bn(3,im)*xn(3)
            end do ! im
            if (isw.eq.3) then
              do lj = 1,3
                do li = 1,3
                  dbn(li,lj) = 0.0d0
                  dbm(li,lj) = 0.0d0
                  do lk = 1,3
                    dbn(li,lj) = dbn(li,lj) +  cn(li,lk)*bn(lk,lj)
     &                                      + cnm(li,lk)*bm(lk,lj)
                    dbm(li,lj) = dbm(li,lj) +  cm(li,lk)*bm(lk,lj)
     &                                      + cnm(lk,li)*bn(lk,lj)
                  end do ! lk
                end do ! li
              end do ! lk

c             Loop over columns

              kk = 0
              do k = 1,nel

c               Geometric stiffness

                shp12 = shp(1,j)*shp(2,k)
                shp21 = shp(2,j)*shp(1,k)
                shp22 = shp(2,j)*shp(2,k)
                do jm = 1,3
                  do im = 1,3
                    s(  jj+im,3+kk+jm) = s(  jj+im,3+kk+jm) +
     &                                   sgn(  im,3+jm)*shp12
                    s(3+jj+im,  kk+jm) = s(3+jj+im,  kk+jm) +
     &                                   sgn(3+im,  jm)*shp21
                    s(3+jj+im,3+kk+jm) = s(3+jj+im,3+kk+jm) +
     &                                   sgn(3+im,3+jm)*shp22

                    s(3+jj+im,3+kk+jm) = s(3+jj+im,3+kk+jm)
     &                                 + sgm(  im,  jm)*shp22
     &                                 + sgm(3+im,  jm)*shp12
     &                                 + sgm(  im,3+jm)*shp21
                  end do ! im
                end do ! jm

c               Material part of stiffness

                call bmbnfb(bm,bn,omg,ugrad,te,sa,shp(2,k),shp(1,k))
                do lj = 1,3
                  do li = 1,3
                    tmpn(li,lj) = 0.0d0
                    do lk = 1,3
                      tmpn(li,lj) = tmpn(li,lj) +  cn(li,lk)*bn(lk,lj)
     &                                          + cnm(li,lk)*bm(lk,lj)
                    end do ! lk
                  end do ! li
                end do ! lk
                shp11 = shp(1,j)*shp(1,k)
                do jm = 1,3
                  do im = 1,3
                    s(  jj+im,  kk+jm) = s(  jj+im,  kk+jm)
     &                                 + shp11*cn(im,jm)
                    s(  jj+im,3+kk+jm) = s(  jj+im,3+kk+jm)
     &                                 + shp(1,j)*tmpn(im,jm)
                    s(3+jj+im,  kk+jm) = s(3+jj+im,  kk+jm)
     &                                 + shp(1,k)*dbn(jm,im)
                    s(3+jj+im,3+kk+jm) = s(3+jj+im,3+kk+jm)
     &                                 + bn(1,jm)*dbn(1,im)
     &                                 + bn(2,jm)*dbn(2,im)
     &                                 + bn(3,jm)*dbn(3,im)

                    s(3+jj+im,3+kk+jm) = s(3+jj+im,3+kk+jm)
     &                                 + bm(1,jm)*dbm(1,im)
     &                                 + bm(2,jm)*dbm(2,im)
     &                                 + bm(3,jm)*dbm(3,im)
                  end do ! im
                end do ! jm
                kk = kk + ndf
              end do ! k
            endif
            jj = jj + ndf
          end do ! j
        end do ! l

c       Compute rotational body force residual and tangent

        if(d(4).gt.0.0d0 .and. d(65).gt.0.0d0) then
          call fbody3w(d(4)*d(32),d(65),xl,ul, p,s, isw.eq.3)
        endif

c       Modify for inertial effects (consistent mass)

        if (ctan(3).ne.0.0 .and. (ndfo(1).gt.0 .or. shflg)) then

c         Lumped  mass

          if (d(7).eq.0.0) then

            xsj = sqrt((xl(1,2) - xl(1,1))**2
     &               + (xl(2,2) - xl(2,1))**2
     &               + (xl(3,2) - xl(3,1))**2)
            xmr = 0.5d0 * d(4) * d(32) * xsj
            xma = xmr * ctan(3)

            jj = 0
            do j = 1,nel
              do im = 1,3
                p(jj+im)       = p(jj+im) - ul(im,j,5) * xmr
                s(jj+im,jj+im) = s(jj+im,jj+im)        + xma
              end do ! im
              jj = jj + ndf
            end do ! j

c         Consistent mass

          elseif(frotas) then

            lint = nel
            if(nint(d(182)).gt.0) then
              call int1dn(lint, sg)
            else
              call int1d(lint, sg)
            endif
            do l = 1,lint

c             Compute shape functions

              call bmshp(sg(1,l),xl,shp,ndm,nel,xjac)
              xsj  = sg(2,l)*xjac

              call extrfb(shp,tm,th,v1,a1)

              call accefb(ul,d,shp,xl,utaro,tm,th,v1,a1,dya,tmat,
     &                    ula,urota, ndm,ndf,nen)

              xmr = d(4) * d(32) * xsj
              if (isw.eq.3) then

                call massfb(v1,th,dya,tmat,utaro,matk,matc,matm)

                xsjk = xsj * ctan(1)
                xsjc = xsj * ctan(2)
                xsjm = xsj * ctan(3)

                xma  = xmr * ctan(3)

                do jm = 1,3
                  do im = 1,3
                  xmat(im,jm) = xsjk * matk(im,jm)
     &                        + xsjc * matc(im,jm)
     &                        + xsjm * matm(im,jm)
                  end do ! im
                end do ! jm
              endif

              do im = 1,3
                ula(im)   = ula(im)   * xmr
                urota(im) = urota(im) * xsj
              end do ! im

              jj = 0
              do j = 1,nel
                do im = 1,3
                  p(jj+im  ) = p(jj+im  ) - shp(2,j)*ula(im  )
                  p(jj+im+3) = p(jj+im+3) - shp(2,j)*urota(im)
                end do ! im

c               Loop over columns

                if (isw.eq.3) then

                  kk = 0
                  do k = 1,nel

c                   Add consistent mass matrix

                    shp22 = shp(2,j)*shp(2,k)
                    do jm = 1,3
                      s(jj+jm  ,kk+jm  ) = s(jj+jm  ,kk+jm  )
     &                                   + shp22 * xma
                      do im = 1,3
                        s(jj+im+3,kk+jm+3) = s(jj+im+3,kk+jm+3)
     &                                     + shp22 * xmat(im,jm)
                      end do ! im
                    end do ! jm
                    kk = kk + ndf
                  end do ! k
                endif
                jj = jj + ndf
              end do ! j
            end do ! l
          else
            write(iow,3000)
            call plstop()
          endif
        endif

c     Output stresses

      elseif(isw.eq.4) then

        if(nint(d(182)).gt.0) then
          lint = nel
          call int1dn(lint, sg)
        else
          lint = nel - 1
          call int1d(lint, sg)
        endif

c       Compute element stresses, strains

        nhv = nint(d(15))
        do l = 1,lint
          call bmshp (sg(1,l),xl,shp,ndm,nel,xjac)

c         Compute coordinates

          xx = 0.0d0
          yy = 0.0d0
          zz = 0.0d0
          do j = 1,nel
            xx = xx + shp(2,j)*xl(1,j)
            yy = yy + shp(2,j)*xl(2,j)
            zz = zz + shp(2,j)*xl(3,j)
          end do ! j

c         Compute Cauchy stresses

          call strefb(hr(nh1),hr(nh2),hr(nh1+3),hr(nh2+3),ul,d,
     &                shp,ugrad,omg,th,thd,xl,te,xm,xn,cn,cm,cnm,xlam,
     &                ndm,ndf,nen, hr(nh1+7),hr(nh2+7),nhv,isw)

c         Output stresses and strains

          mct = mct - 2
          if(mct.le.0) then
            write(iow,2001) o,head
            if(ior.lt.0) then
              write(*,2001) o,head
            endif
            mct = 50
          endif
          write(iow,2002) n,ma,(xm(li),li=1,3),
     &                    xx,yy,zz,(xn(li),li=1,3)
          if(ior.lt.0) then
            write(*,2002) n,ma,(xm(li),li=1,3),
     &                    xx,yy,zz,(xn(li),li=1,3)
          endif
        end do ! l

c     Compute constant mass matrix for initial acceleration

      elseif(isw.eq.5) then

        lint = nel
        if(nint(d(182)).gt.0) then
          call int1dn(lint, sg)
        else
          call int1d(lint, sg)
        endif

        do l = 1,lint

c         Compute shape functions

          call bmshp(sg(1,l),xl,shp,ndm,nel,xjac)
          xsj  = sg(2,l)*xjac

          call extrfb(shp,tm,th,v1,a1)

          call accefb(ul,d,shp,xl,utaro,tm,th,v1,a1,dya,tmat,
     &                ula,urota, ndm,ndf,nen)

          xma         = xsj  * d(4) * d(32)

          do jm = 1,3
            do im = 1,3
              xmat(im,jm) = xsj  * dya(im,jm)
            end do ! im
          end do ! jm

c         For each node j compute db = rho*shape*dv

          jj = 0
          do j = 1,nel

c           Compute a lumped mass
            kk    = 0
            shp22 = xma*shp(2,j)
            do k = 1,3
              p(kk+k) = p(kk+k) + shp22
            end do ! k

            do k = j,nel

c             Compute a consistant mass (symmetry noted)
              shp22 = shp(2,j)*shp(2,k)
              do jm = 1,3
                s(jj+jm  ,kk+jm  ) = s(jj+jm  ,kk+jm  )
     &                             + shp22 * xma
                do im = 1,3
                  s(jj+im+3,kk+jm+3) = s(jj+im+3,kk+jm+3)
     &                               + shp22 * xmat(im,jm)
                end do ! im
              end do ! jm
              kk = kk + ndf
            end do ! k
            jj = jj + ndf
          end do ! j
        end do ! l

c     Compute surface tractions

      elseif(isw.eq.7) then

c     Compute nodal stress values

      elseif(isw.eq.8) then

        do j = 1,2
          do im = 1,6
            fi(im,j) = 0.0d0
          end do ! im
        end do ! j

        if(nint(d(182)).gt.0) then
          lint = nel
          call int1dn(lint, sg)
        else
          lint = nel - 1
          call int1d(lint, sg)
        endif

c       Residual computation

        nhv = nint(d(15))
        do l = 1,lint

          call bmshp (sg(1,l),xl,shp,ndm,nel,xjac)
          xsj = sg(2,l)*xjac
          call strefb(hr(nh1),hr(nh2),hr(nh1+3),hr(nh2+3),ul,d,
     &                shp,ugrad,omg,th,thd,xl,te,xm,xn,cn,cm,cnm,xlam,
     &                ndm,ndf,nen, hr(nh1+7),hr(nh2+7),nhv,isw)

          do im = 1,3
            xn(im) = xn(im)*xsj
            xm(im) = xm(im)*xsj
          end do ! im

c         Form residual

          do j = 1,2

            call bmbnfb(bm,bn,omg,ugrad,te,sa,shp(2,j),shp(1,j))
            do im = 1,3
              fi(im,j) = fi(im,j) + shp(1,j)*xn(im)
              fi(im+3,j) = fi(im+3,j)
     &                   + bm(1,im)*xm(1) + bm(2,im)*xm(2)
     &                   + bm(3,im)*xm(3) + bn(1,im)*xn(1)
     &                   + bn(2,im)*xn(2) + bn(3,im)*xn(3)
            end do ! im
          end do ! j
        end do ! l

c       Rotate axes to current element frame

        do j = 1,2
          do im = 1,3
            sgn(im  ,j) = xlam(1,im)*fi(1,j)
     &                  + xlam(2,im)*fi(2,j)
     &                  + xlam(3,im)*fi(3,j)
            sgn(im+3,j) = xlam(1,im)*fi(4,j)
     &                  + xlam(2,im)*fi(5,j)
     &                  + xlam(3,im)*fi(6,j)
          end do ! im
        end do ! j
        do im = 1,6
          sgn(im,1) = -sgn(im,1)
        end do ! im
        call frcn3d(sgn,p,s)

c     Compute nodal error values

      elseif(isw.eq.9) then

c     Augment

      elseif(isw.eq.10) then

c     Initialize history parameters

      elseif(isw.eq.14) then

        if (nel.ne.2) stop '2-node elem. only'

c       Quaternion initialization

        hr(nh1+6) = 1.d0
        hr(nh2+6) = 1.d0
        if(nint(d(182)).gt.0) then
          lint = nel
          call int1dn(lint, sg)
        else
          lint = nel - 1
          call int1d(lint, sg)
        endif

c       Initialize constitutive history parameters

        nhv = nint(d(15))
        do l = 1,lint

          call bmshp (sg(1,l),xl,shp,ndm,nel,xjac)
          call strefb(hr(nh1),hr(nh2),hr(nh1+3),hr(nh2+3),ul,d,
     &                shp,ugrad,omg,th,thd,xl,te,xm,xn,cn,cm,cnm,xlam,
     &                ndm,ndf,nen, hr(nh1+7),hr(nh2+7),nhv,isw)

        end do ! ll

c     Project stress quantities

      elseif(isw.eq.20) then

        do im = 1,nel
          mxl(im)   = ix(im)
          xll(1,im) = xl(1,im)
          xll(2,im) = xl(2,im)
          xll(3,im) = xl(3,im)
        end do ! im
        if(cs.gt.0.0d0) then
          do im = 1,nel
            xll(1,im) = xll(1,im) + cs*ul(1,im,1)
            xll(2,im) = xll(2,im) + cs*ul(2,im,1)
            xll(3,im) = xll(3,im) + cs*ul(3,im,1)
          end do ! im
        endif
        mel    = nel
        mxl(4) = ma

        if(nint(d(182)).gt.0) then
          lint = nel
          call int1dn(lint, sg)
        else
          lint = nel - 1
          call int1d(lint, sg)
        endif

        nhv = nint(d(15))
        do l = 1,lint
          call bmshp (sg(1,l),xl,shp,ndm,nel,xjac)
          call strefb(hr(nh1),hr(nh2),hr(nh1+3),hr(nh2+3),ul,d,
     &                shp,ugrad,omg,th,thd,xl,te,xm,xn,cn,cm,cnm,xlam,
     &                ndm,ndf,nen, hr(nh1+7),hr(nh2+7),nhv,isw)
        end do ! ll

      endif

c     Formats for input-output

2001  format(a1,20a4//5x,'Element Stresses'//'  Element Material',9x,
     1   '   Mom. - 1    Mom. - 2    Mom. - 3  '/
     2   '  1-Coord  2-Coord  3-Coord',
     3   '   Nor. - 1    Nor. - 2    Nor. - 3 ')

2002  format(2i9,9x,3e12.3/3f9.3,3e12.3/1x)

3000  format('  *ERROR* Need special rotational update - ROTA // 7. ')

      end

      subroutine strefb(df0,df1,qla0,qla1,ul,d,shp,ugrad,omg,
     &                  th,thd,xl,te,xm,xn,cn,cm,cnm,xlam,
     &                  ndm,ndf,nen,hn,h1,nh,isw)

c     Compute strain and stress

      implicit none

      include 'bm3cn.h'
      include 'eldata.h'
      include 'elplot.h'
      include 'iofile.h'

      integer  ndm,ndf,nen,i,im,j,jm,k,s2i(6),isw,nh
      real*8   tmps, dxsm, tol, theta
      real*8   ul(ndf,nen,6),d(*),shp(2,3)
      real*8   ugrad(3),xlam(3,3),th(3),thd(3),cmod(6)
      real*8   eps(3),xm(3),xn(3),ds0(3),xl(ndm,*)
      real*8   cn(3,3),cm(3,3),cnm(3,3),te(3,3),xlam0(3,3),tmp(4)
      real*8   df0(3),df1(3),qla0(4),qla1(4),xlam1(3,3),omg(3)
      real*8   hn(*),h1(*),strain(6),stress(6),mhook(6,6)
      real*8   ahook(6,6),epsb(3),df1b(3),xnb(3),xmb(3)

      data     s2i /2,3,1,5,6,4/,  tol / 1.d-6 /

c     Length, displacement gradient, rotation vector and derivative

      do k = 1,3
        ds0(k)   = shp(1,1) * xl(k  ,1)   + shp(1,2) * xl(k  ,2)
        ugrad(k) = shp(1,1) * ul(k  ,1,1) + shp(1,2) * ul(k  ,2,1)
        th(k)    = shp(2,1) * ul(k+3,1,2) + shp(2,2) * ul(k+3,2,2)
        thd(k)   = shp(1,1) * ul(k+3,1,2) + shp(1,2) * ul(k+3,2,2)
      end do ! k

      if(nel.eq.3) then
        do k = 1,3
          ds0(k)   = ds0(k)   + shp(1,3) * xl(k  ,3)
          ugrad(k) = ugrad(k) + shp(1,3) * ul(k  ,3,1)
          th(k)    = th(k)    + shp(2,3) * ul(k+3,3,2)
          thd(k)   = thd(k)   + shp(1,3) * ul(k+3,3,2)
        end do ! k
      endif

c     Axial deformation gradient

      do k = 1,3
        ugrad(k) = ugrad(k) +  ds0(k)
      end do ! k

c     Incremental curvatures

      call rottrq(th,te)
      do i = 1,3
        omg(i) = te(i,1)*thd(1) + te(i,2)*thd(2) + te(i,3)*thd(3)
      end do ! i

c     Set reference configuration triad

      if    (nint(d(96)).eq.1) then ! REFE,NODE option
        do j = 1,3
          tmp(j) = d(96+j) - xl(j,1)
        end do ! j
      elseif(nint(d(96)).eq.2) then ! REFE,VECT option
        do j = 1,3
          tmp(j) = d(96+j)
        end do ! j
      elseif(nint(d(96)).eq.3) then ! REFE,POLAr option
        df1(1) = 0.5d0*(xl(1,1) + xl(1,2))
        df1(2) = 0.5d0*(xl(2,1) + xl(2,2))
        theta  = atan2(df1(2),df1(1))
        tmp(1) = cos(theta)
        tmp(2) = sin(theta)
        tmp(3) = 0.0d0
      else
        tmp(1) = 0.0d0
        tmp(2) = 0.0d0
        tmp(3) = 1.0d0
      endif
      tmps = tmp(1)*tmp(1) + tmp(2)*tmp(2) + tmp(3)*tmp(3)

c     Check that axis of beam not parallel to reference vector

      if(isw.eq.14) then

        df1(1) = ds0(2)*tmp(3) - ds0(3)*tmp(2)
        df1(2) = ds0(3)*tmp(1) - ds0(1)*tmp(3)
        df1(3) = ds0(1)*tmp(2) - ds0(2)*tmp(1)

        dxsm = ds0(1)**2 + ds0(2)**2 + ds0(3)**2
        if(tmps.eq.0.0d0) then
          write(iow,3001) n
          call plstop()
        endif
        if(dxsm.eq.0.0d0) then
          write(iow,3002) n
          call plstop()
        endif
        if(df1(1)**2 + df1(2)**2 + df1(3)**2 .lt. tol*(dxsm*tmps)) then
          write(iow,3003) n
          call plstop()
        endif
      endif

      if (tmps.eq.0.0d0) then

        dxsm = abs(ds0(1))
        j    = 1
        do i = 2,3
          if(abs(ds0(i)).lt.dxsm) then
            dxsm = abs(ds0(i))
            j    = i
          endif
        end do ! i
        tmp(j) = 1.0d0

      endif

c     Initial curv.

      do i = 1,3
        xlam0(i,1) = ds0(i)
      end do ! i

c     Compute orthonormal initial transformation triad

      call vecp(tmp,xlam0(1,1), xlam0(1,2))
      tmps = 1.d0/sqrt(xlam0(1,2)**2 + xlam0(2,2)**2 + xlam0(3,2)**2)
      do i = 1,3
        xlam0(i,2) = xlam0(i,2)*tmps
      end do ! i

      call vecp(xlam0(1,1),xlam0(1,2), xlam0(1,3))
      tmps = 1.d0/sqrt(xlam0(1,3)**2 + xlam0(2,3)**2 + xlam0(3,3)**2)
      do i = 1,3
        xlam0(i,3) = xlam0(i,3)*tmps
      end do ! i

c     Update rotation using exponential map

c     Form th unit quarternion & th * lamda quaternion product

      call rotqua(th,tmp)
      call quamul(tmp,qla0,qla1)
      call quanrm(qla1)

c     Form orthog. matrix

      call quamat(qla1,xlam1)
      do j = 1,3
        do i = 1,3
          xlam(i,j) = xlam1(i,1)*xlam0(1,j)
     &              + xlam1(i,2)*xlam0(2,j)
     &              + xlam1(i,3)*xlam0(3,j)
        end do ! i
      end do ! j

c     Spatial strain measures

      do j = 1,3
        do i = 1,3
          rot0(i,j) = xlam0(i,j)
          rot1(i,j) = xlam(i,j)
        end do ! i
      end do ! j

      call lamrot(th,xlam1)
      do j = 1,3
        eps(j) =  ugrad(j) - xlam (j,1)
        df1(j) =  omg(j)   + xlam1(j,1)*df0(1)
     &                     + xlam1(j,2)*df0(2)
     &                     + xlam1(j,3)*df0(3)
      end do ! j

c     Spatial tangent moduli

      if(nint(d(100)).eq.0) then

c       Compute material moduli

        cmod(1) = d(1)*d(32)                         ! EA  - ea
        cmod(2) = d(37)*0.5d0*d(1)/(1.d0+d(2))*d(32) ! GA2 - ga2
        cmod(3) = d(38)*0.5d0*d(1)/(1.d0+d(2))*d(32) ! GA3 - ga3
        cmod(4) = 0.5d0*d(1)/(1.d0+d(2))*d(36)       ! GJ  - gi
        cmod(5) = d(1)*d(33)                         ! EI2 - ei2
        cmod(6) = d(1)*d(34)                         ! EI3 - ei3

        do jm = 1,3
          do im = 1,3
            cm(im,jm) =  cmod(4)*xlam(im,1)*xlam(jm,1)
     &                +  cmod(5)*xlam(im,2)*xlam(jm,2)
     &                +  cmod(6)*xlam(im,3)*xlam(jm,3)
            cn(im,jm) =  cmod(1)*xlam(im,1)*xlam(jm,1)
     &                +  cmod(2)*xlam(im,2)*xlam(jm,2)
     &                +  cmod(3)*xlam(im,3)*xlam(jm,3)
            cnm(im,jm) = 0.d0
          end do ! im
        end do ! jm

c       Spatial stress components

        do k = 1,3
          xn(k) = cn(k,1)*eps(1) + cn(k,2)*eps(2) + cn(k,3)*eps(3)
          xm(k) = cm(k,1)*df1(1) + cm(k,2)*df1(2) + cm(k,3)*df1(3)
        end do ! k

c       Mapping between Simo to Ibrahimbegovic elements
c         SIMO                            IBRAHIMBEGOVIC

c         stress(1) = V1                  xn(1) = N  =  int sig dA
c         stress(2) = V2                  xn(2) = V1
c         stress(3) = N  =  int sig dA    xn(3) = V2
c         stress(4) = M1 =  int sig z dA  xm(1) = T
c         stress(5) = M2 = -int sig y dA  xm(2) = M1 =  int sig z dA
c         stress(6) = T                   xm(3) = M2 = -int sig y dA

c         mhook(1,1) =  k_x GA            cn(1,1) =  E int dA
c         mhook(2,2) =  k_y GA            cn(1,5) =  E int z dA
c         mhook(3,3) =  E int dA          cn(1,6) = -E int y dA
c         mhook(3,4) =  E int z dA        cn(2,2) =  GA_2
c         mhook(3,5) = -E int y dA        cn(3,3) =  GA_2
c         mhook(4,4) =  E int z^2 dA      cm(1,1) =  GJ
c         mhook(4,5) = -E int z y dA      cm(2,2) =  E int z^2 dA
c         mhook(5,5) =  E int y^2 dA      cm(2,3) = -E int z y dA
c         mhook(6,6) =  GJ                cm(3,3) =  E int y^2 dA

c     Call integrated model

      else

c       Rotate strains to material frame

        do i = 1, 3
          epsb(i) = xlam(1,i)*eps(1)+xlam(2,i)*eps(2)+xlam(3,i)*eps(3)
          df1b(i) = xlam(1,i)*df1(1)+xlam(2,i)*df1(2)+xlam(3,i)*df1(3)
        end do !i

c       Map components

        strain(3) = epsb(1)
        strain(1) = epsb(2)
        strain(2) = epsb(3)

        strain(6) = df1b(1)
        strain(4) = df1b(2)
        strain(5) = df1b(3)

c       Material frame constitution for resultants

        call bm3res(d,hn,h1,nh,strain, stress,mhook, isw)

c       Remap components

        xnb(1) = stress(3)
        xnb(2) = stress(1)
        xnb(3) = stress(2)

        xmb(1) = stress(6)
        xmb(2) = stress(4)
        xmb(3) = stress(5)

c       Transform to spatial components

        do i = 1, 3
          xn(i) = xlam(i,1)*xnb(1)+xlam(i,2)*xnb(2)+xlam(i,3)*xnb(3)
          xm(i) = xlam(i,1)*xmb(1)+xlam(i,2)*xmb(2)+xlam(i,3)*xmb(3)
        end do !i

        do i = 1, 6
          do j = 1, 6
            ahook(s2i(i),s2i(j)) = mhook(i,j)
          end do !j
        end do !i

c       Transform coupled material moduli to spatial frame

        do jm = 1,3
          do im = 1,3
            cm(im,jm) =  ahook(4,4)*xlam(im,1)*xlam(jm,1)
     &       + (ahook(5,5)*xlam(im,2)+ahook(6,5)*xlam(im,3))*xlam(jm,2)
     &       + (ahook(6,5)*xlam(im,2)+ahook(6,6)*xlam(im,3))*xlam(jm,3)

            cn(im,jm) =  ahook(1,1)*xlam(im,1)*xlam(jm,1)
     &                +  ahook(2,2)*xlam(im,2)*xlam(jm,2)
     &                +  ahook(3,3)*xlam(im,3)*xlam(jm,3)
            cnm(im,jm) = xlam(im,1)*(ahook(1,5)*xlam(jm,2)
     &                             + ahook(1,6)*xlam(jm,3))
          end do ! im
        end do ! jm

      endif

c     Compute local measures for plots

      do jm = 1,3
        tt(jm  ) = xlam(1,jm)*xn(1)
     &           + xlam(2,jm)*xn(2)
     &           + xlam(3,jm)*xn(3)
        tt(jm+3) = xlam(1,jm)*xm(1)
     &           + xlam(2,jm)*xm(2)
     &           + xlam(3,jm)*xm(3)
      end do ! jm

c     Formats

3001  format(' *ERROR* Reference vector has zero length in element',i8)
3002  format(' *ERROR* Frame element has zero length in element',i8)
3003  format(' *ERROR* Frame element parallel to reference vector in',
     &       ' element',i8)

      end

      subroutine geomfb(xn,xm,ugrad,omg,th,thd,te,sa,sgm,sgn)

c     Form constitutive parts for geometric stiffness

      implicit none

      include 'pmod2d.h'

      integer  i,j,l
      real*8   theta,ugrad(3),omg(3),te(3,3),sa(3,3),th(3),thd(3)
      real*8   xm(3),xn(3),sgm(6,6),sgn(6,6),tmpm(3,3),tmp(3,3)

      do j = 1,6
        do i = 1,6
          sgm(i,j) = 0.0d0
          sgn(i,j) = 0.0d0
        end do ! i
      end do ! j
      if(.not.gflag) return
      theta = th(1)*th(1) + th(2)*th(2) + th(3)*th(3)

      do j = 1,3
        do i = 1,3
          tmpm(i,j) = te(1,i)*te(1,j)
     &              + te(2,i)*te(2,j)
     &              + te(3,i)*te(3,j)
        end do ! i
      end do ! j

c     Form xn-related parts

      do l = 1,3
        sgn(1,l+3) = -xn(2)*te(3,l) + xn(3)*te(2,l)
        sgn(2,l+3) = -xn(3)*te(1,l) + xn(1)*te(3,l)
        sgn(3,l+3) = -xn(1)*te(2,l) + xn(2)*te(1,l)
      end do ! l

      do i = 1,3
        do j = 4,6
          sgn(j,i) = sgn(i,j)
        end do ! i
      end do ! j

      call setnfb(th,tmp,te,xn,ugrad,tmpm)
      do j = 1,3
        do i = 1,3
          sgn(i+3,j+3) = tmp(i,j)
        end do ! i
      end do ! j

      call skgmfb(sgn)

c     Form xm-related parts

      call setnfb(th,tmp,te,xm,omg,tmpm)
      do j = 1,3
        do i = 1,3
          sgm(i,j) = tmp(i,j)
        end do ! i
      end do ! j

      call forsfb(thd,th,sa)

      do l = 1,3
        tmp(1,l) =  xm(2)*te(3,l) - xm(3)*te(2,l)
        tmp(2,l) =  xm(3)*te(1,l) - xm(1)*te(3,l)
        tmp(3,l) =  xm(1)*te(2,l) - xm(2)*te(1,l)
      end do ! l

      do j = 1,3
        do i = 1,3
          sgm(i,j+3) = te(1,i)*tmp(1,j)
     &               + te(2,i)*tmp(2,j)
     &               + te(3,i)*tmp(3,j)
          sgm(i,j)   = sgm(i,j)
     &               - sa(1,i)*tmp(1,j)
     &               - sa(2,i)*tmp(2,j)
     &               - sa(3,i)*tmp(3,j)
     &               - sa(1,j)*tmp(1,i)
     &               - sa(2,j)*tmp(2,i)
     &               - sa(3,j)*tmp(3,i)
        end do ! i
      end do ! j

      call setcfb(th,xm,tmp)
      do j = 1,3
        do i = 1,3
          sgm(i,j+3) = sgm(i,j+3) + tmp(i,j)
          sgm(j+3,i) = sgm(i,j+3)
        end do ! i
      end do ! j

      call setafb(th,thd,xm,tmp)
      do j = 1,3
        do i = 1,3
          sgm(i,j)  = sgm(i,j) + tmp(i,j)
        end do ! i
      end do ! j
      call skgmfb(sgm)

      end

      subroutine setafb(th,thd,xm,tmp)

c     Form additional geom. stiffness part - symmetric part

      implicit none

      integer  i,j
      real*8   a1,a2,a3, b1,b2,b3, ca,cb,c1,c2,c3,c4,c5,c6,c7,c8
      real*8   th(3),thd(3),xm(3),tmp(3,3),a(3)

c     Compute a-vector by cross product of thd with x

      a(1) = thd(2)*xm(3) - thd(3)*xm(2)
      a(2) = thd(3)*xm(1) - thd(1)*xm(3)
      a(3) = thd(1)*xm(2) - thd(2)*xm(1)

c     Compute c1 to c8

      call coeffc(th, c1, c2, c3, c4, c5)

      c6 = th (1)*xm(1) + th (2)*xm(2) + th (3)*xm(3)
      c7 = thd(1)*xm(1) + thd(2)*xm(2) + thd(3)*xm(3)
      c8 = thd(1)*th(1) + thd(2)*th(2) + thd(3)*th(3)
      ca = c1*c7 + c3*c6*c8
     &   + ( th(1)*a(1) +  th(2)*a(2)  + th(3)*a(3))*c2

c     Compute a1, a2, a3

      call coeffa(th, a1, a2, a3)
      cb = a1*c7 + a2*(th(1)*a(1) + th(2)*a(2) + th(3)*a(3))
     &   + a3*c8*c6

      do j = 1,3
        b1 = c5*xm(j) + c3*c6*th(j)
        b2 = cb*th(j) + c3*(c6*thd(j) + c8*xm(j)) + c2*a(j)
        b3 = c5*thd(j) + c3*c8*th(j)
        do i = 1,3
          tmp(i,j) = b1*thd(i) + b2*th(i) + b3*xm(i) + c2*th(j)*a(i)
        end do ! i
        tmp(j,j) = tmp(j,j) + ca
      end do ! j

      end

      subroutine setcfb(th,xm,tmp)

c     Form additional geometric stiffness part

      implicit none

      integer  i,j
      real*8   c0,c1,c2,c3,c4,c5,c6,th(3),xm(3),tmp(3,3),a(3)

      a(1) = xm(2)*th(3) - xm(3)*th(2)
      a(2) = xm(3)*th(1) - xm(1)*th(3)
      a(3) = xm(1)*th(2) - xm(2)*th(1)

      call coeffc(th, c1, c2, c3, c4, c5)

      c6 = th(1)*xm(1) + th(2)*xm(2) + th(3)*xm(3)
      c3 = c3*c6
      c0 = c5*c6

      tmp(1,1) =  c0
      tmp(1,2) =  xm(3)*c4
      tmp(1,3) = -xm(2)*c4

      tmp(2,1) = -tmp(1,2)
      tmp(2,2) =  c0
      tmp(2,3) =  xm(1)*c4

      tmp(3,1) = -tmp(1,3)
      tmp(3,2) = -tmp(2,3)
      tmp(3,3) =  c0

      do j = 1,3
        c0 = c1*xm(j) + c2*a(j) + c3*th(j)
        do i = 1,3
          tmp(i,j) = tmp(i,j) + c0 * th(i)
     &                        + c5 * xm(i) * th(j)
        end do ! i
      end do ! j

      end

      subroutine setnfb(th,tmp,te,x,ugrad,tmpm)

c     Form matrix elaborate part of sgm or sgn (beautiful observation)

      implicit none

      integer  i,j
      real*8   c0,c1,c2,c3,c4,c5,c6,tmps
      real*8   th(3),tmp(3,3),at(3),ugrad(3),te(3,3),x(3)
      real*8   a(3),tmpm(3,3),tmpc(3),tmpr(3)

      a(1) = x(2)*ugrad(3) - x(3)*ugrad(2)
      a(2) = x(3)*ugrad(1) - x(1)*ugrad(3)
      a(3) = x(1)*ugrad(2) - x(2)*ugrad(1)

      call coeffc(th, c1, c2, c3, c4, c5)

      c6 = th(1)*a(1) + th(2)*a(2) + th(3)*a(3)
      c3 = c3*c6
      c0 = c5*c6

      at(1) = th(2)*a(3) - th(3)*a(2)
      at(2) = th(3)*a(1) - th(1)*a(3)
      at(3) = th(1)*a(2) - th(2)*a(1)

      tmp(1,1) =  c0
      tmp(1,2) = -a(3)*c4
      tmp(1,3) =  a(2)*c4

      tmp(2,1) = -tmp(1,2)
      tmp(2,2) =  c0
      tmp(2,3) = -a(1)*c4

      tmp(3,1) = -tmp(1,3)
      tmp(3,2) = -tmp(2,3)
      tmp(3,3) =  c0

      do i = 1,3
        c0 = c1*a(i) - c2*at(i) + c3*th(i)
        do j = 1,3
          tmp(i,j) = tmp(i,j) + c5 * th(i) * a(j)
     &                        + c0 * th(j)
        end do ! i
      end do ! j

      do i = 1,3
        tmpc(i) = te(1,i)*x(1)
     &          + te(2,i)*x(2)
     &          + te(3,i)*x(3)
        tmpr(i) = te(1,i)*ugrad(1)
     &          + te(2,i)*ugrad(2)
     &          + te(3,i)*ugrad(3)
      end do ! i
      tmps = x(1)*ugrad(1)
     &     + x(2)*ugrad(2)
     &     + x(3)*ugrad(3)
      do j = 1,3
        do i = 1,3
          tmp(i,j) = tmp(i,j) + tmpc(i)*tmpr(j) - tmps * tmpm(i,j)
        end do ! i
      end do ! j

      end

      subroutine forsfb(a,th,tmp)

c     Form matrix s-alpha

      implicit none

      integer  i,j
      real*8   c0,c1,c2,c3,c4,c5,c6,a(3),th(3),tmp(3,3),at(3)

      call coeffc(th, c1, c2, c3, c4, c5)

      c6       =  th(1)*a(1) + th(2)*a(2) + th(3)*a(3)
      c3       =  c3*c6

      c0       =  c5*c6
      tmp(1,1) =  c0
      tmp(1,2) =  a(3)*c4
      tmp(1,3) = -a(2)*c4

      tmp(2,1) = -tmp(1,2)
      tmp(2,2) =  c0
      tmp(2,3) =  a(1)*c4

      tmp(3,1) = -tmp(1,3)
      tmp(3,2) = -tmp(2,3)
      tmp(3,3) =  c0

      at(1)    =  th(2)*a(3) - th(3)*a(2)
      at(2)    =  th(3)*a(1) - th(1)*a(3)
      at(3)    =  th(1)*a(2) - th(2)*a(1)

      do i = 1,3
        c0 = c1*a(i) + c2*at(i) +c3*th(i)
        do j = 1,3
          tmp(i,j) = tmp(i,j) + c0*th(j)
     &                        + c5*th(i)*a(j)
        end do ! i
      end do ! j

      end

      subroutine skgmfb(a)

c     symmetrize sgm & sgn - for round-off elimination

      implicit none

      integer  i,j
      real*8   a(6,6)

      do j = 1,6
        do i = 1,j
          a(i,j) = (a(i,j) + a(j,i)) * 0.5d0
        end do ! i
      end do ! j
      do j = 1,6
        do i = 1,j
          a(j,i) = a(i,j)
        end do ! i
      end do ! j

      end

      subroutine bmbnfb(bm,bn,omg,ugrad,te,sa,s,sd)

c     Form Bm-alpha & Bn-alpha matrices

      implicit none

      integer  i,j,l,k
      real*8   s,sd
      real*8   bm(3,3),bn(3,3),omg(3),ugrad(3),te(3,3),sa(3,3),tmp(3,3)

      do l = 1,3
        tmp(1,l) =  omg(2)*te(3,l) - omg(3)*te(2,l)
        tmp(2,l) =  omg(3)*te(1,l) - omg(1)*te(3,l)
        tmp(3,l) =  omg(1)*te(2,l) - omg(2)*te(1,l)
      end do ! l
      do j = 1,3
        do i = 1,3
          bm(i,j) = (sa(i,j) + tmp(i,j))*s + te(i,j)*sd
        end do ! i
      end do ! j

      do k = 1,3
        bn(1,k) = (ugrad(2)*te(3,k) - ugrad(3)*te(2,k)) * s
        bn(2,k) = (ugrad(3)*te(1,k) - ugrad(1)*te(3,k)) * s
        bn(3,k) = (ugrad(1)*te(2,k) - ugrad(2)*te(1,k)) * s
      end do ! k

      end

      subroutine rottrq(rot, t)

c     Compute transformation matrix R between rotation vector and
c     quaternion increments.

c     Added by ai.

      implicit none

      real*8 rotnrm, rotnr2, fac1, fac2, fac3, rot(3), t(3,3)

c     Compute norm of rot

      rotnr2 = rot(1)*rot(1) + rot(2)*rot(2) + rot(3)*rot(3)

c     Check for asymptotic solution if ||rot|| = 0.0:

      if (rotnr2 .lt. 1.d-08) then

        fac1  = 1.0d0     - rotnr2*(1.d0/6.d0
     &                    - rotnr2*(1.d0/120.d0
     &                    - rotnr2/5040.d0))
        fac2  = 0.5d0     - rotnr2*(1.d0/24.d0
     &                    - rotnr2*(1.d0/720.d0
     &                    - rotnr2/40320.d0))
        fac3  = 1.d0/6.d0 - rotnr2/(1.d0/120.d0
     &                    - rotnr2/(1.d0/5040.d0
     &                    - rotnr2/362880.d0))
      else
        rotnrm = sqrt ( rotnr2 )
        fac1   = sin(rotnrm)         /rotnrm
        fac2   = (1.d0 - cos(rotnrm))/rotnr2
        fac3   = (1.d0 - fac1       )/rotnr2
      endif

c     Assemble transformation matrix

      t(1,1)  =         fac1 + fac3*rot(1)*rot(1)
      t(1,2)  = -rot(3)*fac2 + fac3*rot(1)*rot(2)
      t(1,3)  =  rot(2)*fac2 + fac3*rot(1)*rot(3)
      t(2,1)  =  rot(3)*fac2 + fac3*rot(2)*rot(1)
      t(2,2)  =         fac1 + fac3*rot(2)*rot(2)
      t(2,3)  = -rot(1)*fac2 + fac3*rot(2)*rot(3)
      t(3,1)  = -rot(2)*fac2 + fac3*rot(3)*rot(1)
      t(3,2)  =  rot(1)*fac2 + fac3*rot(3)*rot(2)
      t(3,3)  =         fac1 + fac3*rot(3)*rot(3)

      end

      subroutine coeffc(rot, fac1, fac2, fac3, fac4, fac5)

c     Compute coefficients c0 c1 c2 c3 c4 c5 a1 a2 a3
c     Added by ma:                                    February 24 1996

      implicit none

      real*8   rotnrm, rotnr2, rotnr3, rotnr4, rotnr5, rot(3)
      real*8   fac0, fac1, fac2, fac3, fac4, fac5

c     Compute norm of rot

      rotnr2 = rot(1)*rot(1) + rot(2)*rot(2) + rot(3)*rot(3)

      if (rotnr2 .lt .1.d-08) then

         fac0  = 1.0d0 - rotnr2*(1.d0/6.d0
     &                 - rotnr2*(1.d0/120.d0
     &                 - rotnr2/5040.d0))

         fac1  = - 1.d0/3.d0 + rotnr2*(1.d0/30.d0
     &                       - rotnr2*(1.d0/840.d0
     &                       - rotnr2/45360.d0))

         fac2  = - 1.d0/12.d0 + rotnr2*(1.d0/180.d0
     &                        - rotnr2*(1.d0/6720.d0
     &                        - rotnr2/453600.d0))

         fac3  = - 1.d0/60.d0 + rotnr2*(1.d0/1260.d0
     &                        - rotnr2*(1.d0/60480
     &                        - rotnr2/4989600.d0))

         fac4  = 0.5d0 - rotnr2*(1.d0/24.d0
     &                 - rotnr2*(1.d0/720.d0
     &                 - rotnr2/40320.d0))

         fac5  = 1.d0/6.d0 - rotnr2/(1.d0/120.d0
     &                     - rotnr2/(1.d0/5040.d0
     &                     - rotnr2/362880.d0))

      else

        rotnrm = sqrt ( rotnr2 )
        rotnr3 = rotnrm * rotnr2
        rotnr4 = rotnrm * rotnr3
        rotnr5 = rotnrm * rotnr4

        fac0=sin(rotnrm)                       /rotnrm

        fac1=(rotnrm*cos(rotnrm) - sin(rotnrm))/rotnr3

        fac2=(rotnrm*sin(rotnrm)
     &       + 2.0d0*cos(rotnrm) - 2.0d0      )/rotnr4

        fac3=( 3.0d0*sin(rotnrm) - 2.0d0*rotnrm
     &      - rotnrm*cos(rotnrm)              )/rotnr5

        fac4=(1.d0 - cos(rotnrm)              )/rotnr2

        fac5=(1.d0 - fac0                     )/rotnr2

      endif

      end

      subroutine coeffa(rot, a1, a2, a3)

c     Compute coefficients a1 a2 a3

c     Added by ma.                                    February 26 1996

      implicit none

      real*8   rotnrm, rotnr2, rotnr3, rotnr4, rotnr5, rotnr6
      real*8   fac2, fac3, a1, a2, a3, rot(3)

c     Compute norm of rot

      rotnr2 = rot(1)*rot(1) + rot(2)*rot(2) + rot(3)*rot(3)

      if (rotnr2 .lt. 1.d-08) then

         fac2  = - 1.d0/12.d0 + rotnr2*(1.d0/180.d0
     &                        - rotnr2*(1.d0/6720.d0
     &                        - rotnr2/453600.d0))

         fac3  = - 1.d0/60.d0 + rotnr2*(1.d0/1260.d0
     &                        - rotnr2*(1.d0/60480
     &                        - rotnr2/4989600.d0))

         a1    = fac3 - fac2

         a2    = 1.d0/90.d0   - rotnr2*(1.d0/1680.d0
     &                        - rotnr2*(1.d0/75600.d0
     &                        - rotnr2/5987520.d0))

         a3    = 1.d0/630.d0  - rotnr2*(1.d0/15120.d0
     &                        - rotnr2*(1.d0/831600.d0
     &                        - rotnr2/77837760.d0))

      else

        rotnrm = sqrt ( rotnr2 )
        rotnr3 = rotnrm * rotnr2
        rotnr4 = rotnrm * rotnr3
        rotnr5 = rotnrm * rotnr4
        rotnr6 = rotnrm * rotnr5

        a1     = (3.0d0*sin(rotnrm)
     &         -  2.0d0*rotnrm-rotnrm*cos(rotnrm)               )/rotnr5
     &         - (rotnrm*sin(rotnrm) + 2.0d0*cos(rotnrm) - 2.0d0)/rotnr4

        a2     = (rotnrm*cos(rotnrm) - sin(rotnrm)              )/rotnr5
     &         -  4.0d0*(rotnrm*sin(rotnrm) + 2.0d0*cos(rotnrm)
     &                                                   - 2.0d0)/rotnr6

        a3     = ((rotnrm*sin(rotnrm) + 2.0d0*cos(rotnrm)
     &                                                   - 2.0d0)/rotnr4
     &         -  5.0d0*((3.0d0*sin(rotnrm) - 2.0d0*rotnrm
     &         -  rotnrm*cos(rotnrm)                   )/rotnr5))/rotnr2

      endif

      end

      subroutine bmshp(s,xl,shp,ndm,nel,ds)

      implicit none

      integer ndm,nel,i,j
      real*8  ds,s,ss, xl(ndm,nel),shp(2,nel),dsc(3)

c     Derivatives

      if (nel.eq.2) then
        shp(1,1) = -0.5d0
        shp(1,2) =  0.5d0
      elseif (nel.eq.3) then
        shp(1,1) = -0.5d0 + s
        shp(1,2) = -2.0d0 * s
        shp(1,3) =  0.5d0 + s
      end if

      do  i = 1,ndm
        dsc(i) = 0.0d0
      end do ! i
      do i = 1,nel,nel-1
        do j = 1,ndm
          dsc(j) = dsc(j) + shp(1,i)*xl(j,i)
        end do ! i
      end do ! j
      ds = dsc(1)*dsc(1) + dsc(2)*dsc(2) + dsc(3)*dsc(3)
      ds = sqrt(ds)
      if (ds.eq.0.0d0) stop ' Jacobian = 0.'
      do i = 1,nel
        shp(1,i) = shp(1,i)/ds
      end do ! i

c     Shape functions

      if (nel.eq.2) then
        shp(2,1) = 0.5d0*(1.d0 - s)
        shp(2,2) = 0.5d0*(1.d0 + s)
      elseif (nel.eq.3) then
        ss       = 1.d0 - s*s
        shp(2,1) = 0.5d0*((1.d0 - s) - ss)
        shp(2,2) = ss
        shp(2,3) = 0.5d0*((1.d0 + s) - ss)
      end if

      end

      subroutine extrfb(shp,tm,th,v1,a1)

c     Extract rotation parameters at Gauss points

c     tm(*) - Total rotation vector at t_n    theta(lamda_n)
c     th(*) - Relative rotation vector        theta_n
c     v1(*) - Spatial rotaional velocity      at t_n+1
c     a1(*) - Spatial rotational acceleration at t_n+1

      implicit none

      include 'erotas.h'

      integer  k
      real*8   shp(2,3),tm(3),th(3),a1(3),v1(3), tmi(3),tmj(3)

      call quarot(xln(1,1,2),tmi)
      call quarot(xln(1,2,2),tmj)

      do k = 1,3
        tm(k) = shp(2,1)*tmi(k)      + shp(2,2)*tmj(k)
        th(k) = shp(2,1)*rots(k,1,2) + shp(2,2)*rots(k,2,2)
        v1(k) = shp(2,1)*rvel(k,1,2) + shp(2,2)*rvel(k,2,2)
        a1(k) = shp(2,1)*racc(k,1,2) + shp(2,2)*racc(k,2,2)
      end do ! k

      end

      subroutine accefb(ul,d,shp,xl,utaro,tm,th,v1,a1,dya,tmat,
     &                  ula,urota, ndm,ndf,nen)

      implicit none

      integer  ndm,ndf,nen,i,j,k
      real*8   ul(ndf,nen,5),d(*),shp(2,3),ula(3),urota(3)
      real*8   ds0(3),xlam0(3,3),xl(ndm,*),tmps,dxsm,tmp(3),dya(3,3)
      real*8   th(3),xlam1(3,3),xlam(3,3), theta
      real*8   tm(3),tmat(3,3),utaro(3),v1(3),a1(3)

      do k = 1,3
        ds0(k) = shp(1,1)*xl(k,1)   + shp(1,2)*xl(k,2)
        ula(k) = shp(2,1)*ul(k,1,5) + shp(2,2)*ul(k,2,5)
      end do ! k

c     Set reference configuration triad

      if    (nint(d(96)).eq.1) then ! REFE,NODE option
        do j = 1,3
          tmp(j) = d(96+j) - xl(j,1)
        end do ! j
      elseif(nint(d(96)).eq.2) then ! REFE,VECT option
        do j = 1,3
          tmp(j) = d(96+j)
        end do ! j
      elseif(nint(d(96)).eq.3) then ! REFE,POLAr option
        utaro(1) = 0.5d0*(xl(1,1) + xl(1,2))
        utaro(2) = 0.5d0*(xl(2,1) + xl(2,2))
        theta    = atan2(utaro(2),utaro(1))
        tmp(1)   = cos(theta)
        tmp(2)   = sin(theta)
        tmp(3)   = 0.0d0
      else
        tmp(1)   = 0.0d0
        tmp(2)   = 0.0d0
        tmp(3)   = 1.0d0
      endif
      tmps = tmp(1)*tmp(1) + tmp(2)*tmp(2) + tmp(3)*tmp(3)

      if (tmps.eq.0.0d0) then

        dxsm = abs(ds0(1))
        j    = 1
        do i = 2,3
          if(abs(ds0(i)).lt.dxsm) then
            dxsm = abs(ds0(i))
            j    = i
          endif
        end do ! i
        tmp(j) = 1.0d0

      endif

      tmps = 1.d0/sqrt(ds0(1)*ds0(1) + ds0(2)*ds0(2) + ds0(3)*ds0(3))
      do i = 1,3
        xlam0(i,1) = ds0(i)*tmps
      end do ! i

c     Compute orthonormal initial transformation triad

      call vecp(tmp,xlam0(1,1), xlam0(1,2))
      tmps = 1.d0/sqrt(xlam0(1,2)**2 + xlam0(2,2)**2 + xlam0(3,2)**2)
      do i = 1,3
        xlam0(i,2) = xlam0(i,2)*tmps
      end do ! i
      call vecp(xlam0(1,1),xlam0(1,2),xlam0(1,3))

c     Spatial inertia tensor

      call lamrot(tm,xlam1)

      do j = 1,3
        do i = 1,3
          xlam(i,j) = xlam1(i,1)*xlam0(1,j)
     &              + xlam1(i,2)*xlam0(2,j)
     &              + xlam1(i,3)*xlam0(3,j)
        end do ! k
      end do ! j

      do j = 1,3
        do i = 1,3
          dya(i,j) =  d(4) *(d(36) * xlam(i,1)*xlam(j,1)
     &                     + d(33) * xlam(i,2)*xlam(j,2)
     &                     + d(34) * xlam(i,3)*xlam(j,3)
     &                     + d(35) *(xlam(i,2)*xlam(j,3)
     &                             + xlam(i,3)*xlam(j,2)))
        end do ! i
      end do ! j

      do i = 1,3
        tmp(i) = dya(i,1)*a1(1) + dya(i,2)*a1(2) + dya(i,3)*a1(3)
        ds0(i) = dya(i,1)*v1(1) + dya(i,2)*v1(2) + dya(i,3)*v1(3)
      end do ! i

c     v1 cross ( dya * v1 )

      call vecp(v1,ds0,utaro)

      do i = 1,3
        utaro(i) = utaro(i) + tmp(i)
      end do ! i

      call rottrq(th,tmat)

      do i = 1,3
        urota(i) = tmat(1,i)*utaro(1)
     &           + tmat(2,i)*utaro(2)
     &           + tmat(3,i)*utaro(3)
      end do ! i

      end

      subroutine massfb(v1,th,dya,tmat,utaro, matk,matc,matm)

c     Compute dynamic inertia terms for rotation part

      implicit none

      integer  i,j

      real*8   matk(3,3),matc(3,3),matm(3,3)
      real*8   v1(3)
      real*8   th(3)
      real*8   dya(3,3)
      real*8   utaro(3)
      real*8   tmat(3,3)
      real*8   tmpm(3,3),tnpn(3,3),tmp(3)
      real*8   xlam(3,3)

c     PSEUDO STIFFNESS TERM

      call pzero(matk,3*3)

      call settfb(th,utaro,matk)

      call pzero ( tmpm , 3*3 )

      tmpm(1,2) = - utaro(3)
      tmpm(2,1) =   utaro(3)
      tmpm(1,3) =   utaro(2)
      tmpm(3,1) = - utaro(2)
      tmpm(2,3) = - utaro(1)
      tmpm(3,2) =   utaro(1)

      do j = 1,3
        do i = 1,3
          tnpn(i,j) = - tmpm(i,1) * tmat(1,j)
     &                - tmpm(i,2) * tmat(2,j)
     &                - tmpm(i,3) * tmat(3,j)
        end do ! j
      end do ! i

      do j = 1,3
        do i = 1,3
          matk(i,j) = matk(i,j) + tmat(1,i)*tnpn(1,j)
     &                          + tmat(2,i)*tnpn(2,j)
     &                          + tmat(3,i)*tnpn(3,j)
        end do ! j
      end do ! i

c     MASS TERM

c     exp [ \bmtheta_{act,n+1} \times ]

      call lamrot(th,xlam)

      do j = 1,3
        do i = 1,3
          tmpm(i,j) = dya(i,1)*xlam(1,j)
     &              + dya(i,2)*xlam(2,j)
     &              + dya(i,3)*xlam(3,j)
        end do ! i
      end do ! j

      do j = 1,3
        do i = 1,3
          matm(i,j) = tmat(1,i)*tmpm(1,j)
     &              + tmat(2,i)*tmpm(2,j)
     &              + tmat(3,i)*tmpm(3,j)
        end do ! j
      end do ! i

c     PSEUDO DAMPING TERM

c---- [v1 x] dya - [dya v1 x]

      tmpm(1,1) =   0.0d0
      tmpm(1,2) = - v1(3)
      tmpm(1,3) =   v1(2)

      tmpm(2,1) =   v1(3)
      tmpm(2,2) =   0.0d0
      tmpm(2,3) = - v1(1)

      tmpm(3,1) = - v1(2)
      tmpm(3,2) =   v1(1)
      tmpm(3,3) =   0.0d0

      do j = 1,3
        do i = 1,3
          tnpn(i,j) = tmpm(i,1)*dya(1,j)
     &              + tmpm(i,2)*dya(2,j)
     &              + tmpm(i,3)*dya(3,j)
        end do ! i
      end do ! j

      do i = 1,3
        tmp(i) = dya(i,1)*v1(1) + dya(i,2)*v1(2) + dya(i,3)*v1(3)
      end do ! i

      tmpm(1,1) =   0.0d0
      tmpm(1,2) = - tmp(3)
      tmpm(1,3) =   tmp(2)

      tmpm(2,1) =   tmp(3)
      tmpm(2,2) =   0.0d0
      tmpm(2,3) = - tmp(1)

      tmpm(3,1) = - tmp(2)
      tmpm(3,2) =   tmp(1)
      tmpm(3,3) =   0.0d0

      do j = 1,3
        do i = 1,3
          tnpn(i,j) = tnpn(i,j) - tmpm(i,j)
        end do ! j
      end do ! i

      call pzero ( tmpm , 3*3 )

      do j = 1,3
        do i = 1,3
          tmpm(i,j) = tnpn(i,1)*xlam(1,j)
     &              + tnpn(i,2)*xlam(2,j)
     &              + tnpn(i,3)*xlam(3,j)
        end do ! i
      end do ! j

      do j = 1,3
        do i = 1,3
          matc(i,j) = tmat(1,i)*tmpm(1,j)
     &              + tmat(2,i)*tmpm(2,j)
     &              + tmat(3,i)*tmpm(3,j)
        end do ! j
      end do ! i

      end

      subroutine settfb(th,x,tmp)

c     Form big theta matrix

      implicit none

      integer  i,j
      real*8   c1,c2,c3,c4,c5,c6, th(3),x(3),tmp(3,3),a(3)

c     Initial step

      call vecp(th,x,a)

      call coeffc(th, c1, c2, c3, c4, c5)

      c6 = th(1)*x(1) + th(2)*x(2) + th(3)*x(3)

      do i = 1,3
        tmp(i,i)  = c5*c6
      end do ! i

      tmp(1,1) =  0.0d0
      tmp(2,2) =  0.0d0
      tmp(3,3) =  0.0d0

      tmp(1,2) = -x(3)*c4
      tmp(2,1) =  tmp(1,2)

      tmp(1,3) =  x(2)*c4
      tmp(3,1) =  tmp(1,3)

      tmp(2,3) = -x(1)*c4
      tmp(3,2) =  tmp(2,3)

      do j = 1,3
        do i = 1,3
          tmp(i,j) = tmp(i,j)
     &             + c1 * x(i)     * th(j)
     &             - c2 * a(i)     * th(j)
     &             + c3*c6 * th(i) * th(j)
     &             + c5 * th(i) * x(j)
        end do ! i
      end do ! j

      end
