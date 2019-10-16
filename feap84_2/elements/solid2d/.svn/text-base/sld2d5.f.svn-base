c$Id:$
      subroutine sld2d5(d,ul,xl,ix,tl,s,p,ndf,ndm,nst,isw)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Add computation of xref and xcur for modlsd call 25/07/2009
c       2. Compute v_avg and sig_33 for multiscale use      03/12/2010
c       3. Add 'l' to modlsd call                           05/01/2012
c       4. Add average of density for multiscale            10/05/2012
c       5. Add eps on call to slcn2d                        01/01/2013
c       6. Pass strains to stcn2z for z-zhu projections     01/01/2014
c-----[--.----+----.----+----.-----------------------------------------]
c     Purpose: Mixed-Enhanced formulation for small deformation
c              2-d plane problems.

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
      implicit   none

      include   'bdata.h'
      include   'cdata.h'
      include   'eldata.h'
      include   'elplot.h'
      include   'eltran.h'
      include   'fdata.h'
      include   'hdata.h'
      include   'iofile.h'
      include   'oelmt.h'
      include   'part0.h'
      include   'qudshp.h'
      include   'rdata.h'

      include   'comblk.h'

      logical    transp, noconv
      integer    ndf , ndm , nst , isw, tsw, ix(*)
      integer    i,i1,j,j1,k,l, nhv, nn,nenit, istrt
      real*8     dd(6,6,5,4),btd(2,4)
      real*8     bmat(4,4,2,4),bmate(4,2,4),eps(9,3),sig(10,4),dsig(4)
      real*8     ui(2),dui(2), bb(2),hh(2,2),gg(2,2,4),epsv(9),epsp(6,4)
      real*8     d(*),ul(ndf,nen,*),xl(ndm,*),tl(*),s(nst,*),p(ndf,*)
      real*8     alam,ha, ta,xx,yy, tol1,tolu, cfac,lfac,jac0, mas, dv

      save

      data       alam /   0.0d0 /, ha   / 0.0d0 /
      data       dsig / 4*0.0d0 /, jac0 / 0.0d0 /

c     Material data

      if(isw.eq.1) then

        if(nint(d(16)).eq.8 .and. nint(d(19)).ne.3) ix(3) = 0

c     Check data

      elseif(isw.eq.2) then

      else

c       Form constant parts to arrays

        tsw    =  max(1,min(4,nint(d(80))))
        transp = .not.(d(81).eq.0.0d0)
        call geom01(xl,ndm)
        call ttran1(tsw,transp)
        call atran1()
        call gstrn1()

        nhv   = nint(d(15))
        istrt = nint(d(84))

c       Get quadrature formula

        if(nint(d(182)).gt.0) then
          call int2dn(nel,lint,sg2)
        else
          l = 2
          call int2d(l,lint,sg2)
        endif

        ui(1)  = hr(nh2  )
        ui(2)  = hr(nh2+1)

c       Stiffness and residual computations

        if    (isw.eq.3 .or. isw.eq.6) then

c         Integration order set to static

          if((d(7).ge.0.0 .or. d(183).ne.0.0d0) .and.
     &               (ndfo(1).gt.0 .or. shflg)) then
            cfac = d(7)
            lfac = 1.d0 - cfac
          else
            cfac = 0.0d0
            lfac = 0.0d0
          endif

c         LOCAL ITERATION ON ENHANCED TERMS

          noconv = .true.
          nhv    = nint(d(15))
          tolu   = max(1.d-03*tol*rnmax/numel,tol)
          nenit  = 0

          do while (noconv)

            do i = 1,2
              bb(i) = 0.0d0
              do j = 1,2
                hh(j,i) = 0.0d0
              end do ! j
            end do ! i

c           LOOP OVER GAUSS POINTS

            nn = 2
            do l = 1,lint

c             Compute strain-displacement terms

              call mestrn(d,sg2(1,l),xl,ul,ui,tl,ndm,ndf, jac(l),
     &                    shp2(1,1,l), bmat(1,1,1,l),bmate(1,1,l),
     &                    ta,eps)
              jac(l) = jac(l)*sg2(3,l)

c             Compute Cauchy stresses & spatial tangent tensor at t-n+1

              call modlsd(l,d,ta,eps,hr(nh1+nn),hr(nh2+nn),nhv,istrt,
     &                    dd(1,1,1,l),sig(1,l),alam,ha,isw)

c             Multiply tangent moduli and stresses by volume element.

              k = 6*l-6
              do i = 1,4
                tt(i+k)  = sig(i,l)
                sig(i,l) = sig(i,l)*jac(l)
                do j = 1,4
                  dd(j,i,1,l) = dd(j,i,1,l)*jac(l)*ctan(1)
                end do ! j
              end do ! i

c             Compute parameters for multiscale plane strain

              v_avg  = v_avg  + jac(l)
              v_rho  = v_rho  + jac(l)*d(4)
              sig_33 = sig_33 + sig(3,l)

c             Enhanced mode only

              do j = 1,4
                do i = 1,2
                  btd(i,j) = bmate(1,i,l)*dd(1,j,1,l)
     &                     + bmate(2,i,l)*dd(2,j,1,l)
     &                     + bmate(4,i,l)*dd(4,j,1,l)
                end do ! i
              end do ! j

c             Compute residual and tangent stiffness

              do j = 1,2
                bb(j)     = bb(j) - bmate(1,j,l)*sig(1,l)
     &                            - bmate(2,j,l)*sig(2,l)
     &                            - bmate(4,j,l)*sig(4,l)
                do i = 1,2
                  hh(i,j) = hh(i,j) + btd(i,1)*bmate(1,j,l)
     &                              + btd(i,2)*bmate(2,j,l)
     &                              + btd(i,4)*bmate(4,j,l)
                end do ! i
              end do ! j

              nn = nn + nhv

            end do ! l

c           Invert diagonal array

            call invert(hh,2,2)

c           Compute K22-inverse * b

            dui(1) = hh(1,1)*bb(1) + hh(1,2)*bb(2)
            dui(2) = hh(2,1)*bb(1) + hh(2,2)*bb(2)

c           Check convergence

            tol1 = bb(1)*dui(1) + bb(2)*dui(2)

            if(abs(tol1).le.tolu .and. nenit.ge.1) then
              noconv = .false.
            elseif(tolu.eq.0.0d0  .or. nenit.ge.2) then
              noconv = .false.
            endif

c           Update enhanced parameters

c           if(noconv) then
              ui(1) = ui(1) + dui(1)
              ui(2) = ui(2) + dui(2)
c           endif

c           Increment loop counter

            nenit = nenit + 1

          end do ! while
          hr(nh2  ) = ui(1)
          hr(nh2+1) = ui(2)

          do l = 1,4
            do i = 1,2
              do j = 1,2
                gg(j,i,l) = 0.0d0
              end do ! j
            end do ! i
          end do ! l

c         Form mixed terms

          do l = 1,lint

c           Compute coordinates

            xx = shp2(3,1,l)*xl(1,1) + shp2(3,2,l)*xl(1,2)
     &         + shp2(3,3,l)*xl(1,3) + shp2(3,4,l)*xl(1,4)
            yy = shp2(3,1,l)*xl(2,1) + shp2(3,2,l)*xl(2,2)
     &         + shp2(3,3,l)*xl(2,3) + shp2(3,4,l)*xl(2,4)

c           Form residual

            call resid2d(cfac,lfac,jac(l),jac0,xx,yy,shp2(1,1,l),dsig,d,
     &                   ul(1,1,4),ul(1,1,5),p,ndf)

            do i = 1,4
              p(1,i) = p(1,i) - bmat(i,1,1,l)*sig(1,l)
     &                        - bmat(i,2,1,l)*sig(2,l)
     &                        - bmat(i,4,1,l)*sig(4,l)
              p(2,i) = p(2,i) - bmat(i,1,2,l)*sig(1,l)
     &                        - bmat(i,2,2,l)*sig(2,l)
     &                        - bmat(i,4,2,l)*sig(4,l)
            end do ! i

c           Form tangent matrix

            if(isw.eq.3) then

              dv = d(4)*ctan(3)*jac(l)

              i1 = 0
              do i = 1,4
                do j = 1,4
                  btd(1,j) = bmat(i,1,1,l)*dd(1,j,1,l)
     &                     + bmat(i,2,1,l)*dd(2,j,1,l)
     &                     + bmat(i,4,1,l)*dd(4,j,1,l)
                  btd(2,j) = bmat(i,1,2,l)*dd(1,j,1,l)
     &                     + bmat(i,2,2,l)*dd(2,j,1,l)
     &                     + bmat(i,4,2,l)*dd(4,j,1,l)
                end do ! j

c               Lumped mass tangent

                mas = shp2(3,i,l)*dv
                s(i1+1,i1+1) = s(i1+1,i1+1) + mas*lfac
                s(i1+2,i1+2) = s(i1+2,i1+2) + mas*lfac

c               Mixed stiffness

                j1 = 0
                do j = 1,4
                  s(i1+1,j1+1) = s(i1+1,j1+1) + btd(1,1)*bmat(j,1,1,l)
     &                                        + btd(1,2)*bmat(j,2,1,l)
     &                                        + btd(1,4)*bmat(j,4,1,l)
     &                                        + cfac*mas*shp2(3,j,l)
                  s(i1+1,j1+2) = s(i1+1,j1+2) + btd(1,1)*bmat(j,1,2,l)
     &                                        + btd(1,2)*bmat(j,2,2,l)
     &                                        + btd(1,4)*bmat(j,4,2,l)
                  s(i1+2,j1+1) = s(i1+2,j1+1) + btd(2,1)*bmat(j,1,1,l)
     &                                        + btd(2,2)*bmat(j,2,1,l)
     &                                        + btd(2,4)*bmat(j,4,1,l)
                  s(i1+2,j1+2) = s(i1+2,j1+2) + btd(2,1)*bmat(j,1,2,l)
     &                                        + btd(2,2)*bmat(j,2,2,l)
     &                                        + btd(2,4)*bmat(j,4,2,l)
     &                                        + cfac*mas*shp2(3,j,l)
                  j1 = j1 + ndf
                end do ! j

c               Coupling stiffness

                gg(1,1,i) = gg(1,1,i) + btd(1,1)*bmate(1,1,l)
     &                                + btd(1,2)*bmate(2,1,l)
     &                                + btd(1,4)*bmate(4,1,l)
                gg(2,1,i) = gg(2,1,i) + btd(1,1)*bmate(1,2,l)
     &                                + btd(1,2)*bmate(2,2,l)
     &                                + btd(1,4)*bmate(4,2,l)
                gg(1,2,i) = gg(1,2,i) + btd(2,1)*bmate(1,1,l)
     &                                + btd(2,2)*bmate(2,1,l)
     &                                + btd(2,4)*bmate(4,1,l)
                gg(2,2,i) = gg(2,2,i) + btd(2,1)*bmate(1,2,l)
     &                                + btd(2,2)*bmate(2,2,l)
     &                                + btd(2,4)*bmate(4,2,l)
                i1 = i1 + ndf
              end do ! i
            endif ! isw = 3

          end do ! l

c         Condense enhanced stiffness and residual parts

          call stcon2d(hh,gg,dui,ndf,nst,s,p)

c         Multiply by thickness if not unity

          if((isw.eq.3 .or. isw.eq.6) .and. d(14).ne.1.d0) then

            do j = 1,nst
              do i = 1,nst
                s(i,j) = s(i,j)*d(14)
              end do ! i
            end do ! j
            do j = 1,nel
              do i = 1,ndf
                p(i,j) = p(i,j)*d(14)
              end do ! i
            end do ! j

          endif

c       Output of element quantities

        elseif(isw.eq.4.or.isw.eq.8.or.isw.eq.16.or.isw.eq.25) then

          nn = 2
          do l = 1,lint

c           Form strain terms

            call mestrn(d,sg2(1,l),xl,ul,ui,tl,ndm,ndf, jac(l),
     &                  shp2(1,1,l),bmat,bmate, ta,eps)

            epsv(l) = eps(1,1) + eps(2,1) + eps(3,1)

            jac(l) = jac(l)*sg2(3,l)

            call modlsd(l,d,ta,eps,hr(nh1+nn),hr(nh2+nn),nhv,istrt,
     &                  dd,sig(1,l),alam,ha,isw)

            if(isw.eq.4) then

              call pstr2d(sig(1,l),sig(7,l))

c             Compute coordinates

              xx = shp2(3,1,l)*xl(1,1) + shp2(3,2,l)*xl(1,2)
     &           + shp2(3,3,l)*xl(1,3) + shp2(3,4,l)*xl(1,4)
              yy = shp2(3,1,l)*xl(2,1) + shp2(3,2,l)*xl(2,2)
     &           + shp2(3,3,l)*xl(2,3) + shp2(3,4,l)*xl(2,4)

c             Output stresses and strains

              mct = mct - 2
              if(mct.le.0) then
              write(iow,2001) o,head
                if(ior.lt.0 .and. pfr) then
                  write(*,2001) o,head
                endif
                mct = 50
              endif
              write(iow,2002)  n,ma,sig(9,l),(sig(i,l),i=1,4),sig(7,l),
     &                                 xx,yy,(eps(i,1),i=1,4),sig(8,l)
              if(ior.lt.0 .and. pfr) then
                write(*,2002)  n,ma,sig(9,l),(sig(i,l),i=1,4),sig(7,l),
     &                                 xx,yy,(eps(i,1),i=1,4),sig(8,l)
              endif

c           Store strains for plots

            elseif(isw.eq.8) then

              do j = 1,6
                epsp(j,l) = eps(j,1)
              end do ! j
            endif
            nn = nn + nhv
          end do ! l

c         Compute nodal stress values

          if(isw.eq.8) then

            call slcn2d(ix,sig,epsp,p,s,p(nen+1,1),nel,10)

c         Compute J-integrals and material forces

          elseif(isw.eq.16) then

            call pjint2d(d,ul,tl,shp2,jac,epsv,sig,p,ndf,ndm,lint,10)

c         Compute Z-Z projections

          elseif(isw.eq.25) then

            call stcn2z(xl,sig,epsp,shp2,jac,lint,ndm,nel,10)

          endif
        endif
      endif

c     Formats for input-output

2001  format(a1,20a4//5x,'Element Stresses'//'    Elmt Mat Angle',
     &   '   11-stress   22-stress   33-stress   12-stress',
     &   '    1-stress'/'  1-coord  2-coord   11-strain',
     &   '   22-strain   33-strain   12-strain    2-stress')
2002  format(i8,i4,0p,f6.1,1p,5e12.3/0p,2f9.3,1p,5e12.3/1x)

      end

      subroutine geom01(xl,ndm)

      implicit   none

      include   'egeom1.h'
      include   'egeom2.h'

      integer    ndm
      real*8     j0r, xl(ndm,*)

c     Shape function natural derivative constants

      x0  = 0.25d0*( xl(1,1) + xl(1,2) + xl(1,3) + xl(1,4))
      y0  = 0.25d0*( xl(2,1) + xl(2,2) + xl(2,3) + xl(2,4))

      x1  = 0.25d0*(-xl(1,1) + xl(1,2) + xl(1,3) - xl(1,4))
      y1  = 0.25d0*(-xl(2,1) + xl(2,2) + xl(2,3) - xl(2,4))

      x2  = 0.25d0*(-xl(1,1) - xl(1,2) + xl(1,3) + xl(1,4))
      y2  = 0.25d0*(-xl(2,1) - xl(2,2) + xl(2,3) + xl(2,4))

      x12 = 0.25d0*( xl(1,1) - xl(1,2) + xl(1,3) - xl(1,4))
      y12 = 0.25d0*( xl(2,1) - xl(2,2) + xl(2,3) - xl(2,4))

c     Jacobian constants
c       j = j0 + j1*xi + j2*eta

      j0  = x1*y2  - x2*y1
      j1  = x1*y12 - x12*y1
      j2  = x12*y2 - x2*y12

c     Shape fucntion constants
c       N(i),x = nx0(i) + (nx1(i)*xi + nx2(i)*eta)/j
c       N(i),y = ny0(i) + (ny1(i)*xi + ny2(i)*eta)/j

      j0r    = 0.25d0/j0
      nx0(1) = (-y2 + y1)*j0r
      nx0(2) = ( y2 + y1)*j0r
      nx0(3) = ( y2 - y1)*j0r
      nx0(4) = (-y2 - y1)*j0r

      ny0(1) = ( x2 - x1)*j0r
      ny0(2) = (-x2 - x1)*j0r
      ny0(3) = (-x2 + x1)*j0r
      ny0(4) = ( x2 + x1)*j0r

      nx1(1) = 0.25d0*(-y12 - y1)
      nx1(2) = 0.25d0*( y12 + y1)
      nx1(3) = 0.25d0*( y12 - y1)
      nx1(4) = 0.25d0*(-y12 + y1)

      ny1(1) = 0.25d0*( x12 + x1)
      ny1(2) = 0.25d0*(-x12 - x1)
      ny1(3) = 0.25d0*(-x12 + x1)
      ny1(4) = 0.25d0*( x12 - x1)

      nx2(1) = 0.25d0*( y12 + y2)
      nx2(2) = 0.25d0*( y12 - y2)
      nx2(3) = 0.25d0*(-y12 + y2)
      nx2(4) = 0.25d0*(-y12 - y2)

      ny2(1) = 0.25d0*(-x12 - x2)
      ny2(2) = 0.25d0*(-x12 + x2)
      ny2(3) = 0.25d0*( x12 - x2)
      ny2(4) = 0.25d0*( x12 + x2)

      end

      subroutine atran1()

      implicit   none

      include   'egeom3.h'
      include   'egeom4.h'

      integer    i

      do i = 1,2

c       Strain transformation terms

        ai(1,i) = ti(i,1)*ti(i,1)
        ai(2,i) = ti(i,2)*ti(i,2)
        ai(3,i) = 0.0d0
        ai(4,i) = ti(i,1)*ti(i,2)*2.d0

c       Stress transformation terms

        aa(1,i) = tb(1,i)*tb(1,i)
        aa(2,i) = tb(2,i)*tb(2,i)
        aa(3,i) = 0.0d0
        aa(4,i) = tb(1,i)*tb(2,i)

      end do ! i

      end

      subroutine ttran1(tsw,transp)

      implicit   none

      include   'egeom1.h'
      include   'egeom4.h'

      logical    transp
      integer    i,j, tsw
      real*8     j0r,jr1,jr2

c     Transpose check

      if(transp) then
        i = 2
        j = 1
      else
        i = 1
        j = 2
      endif ! transp

c     Average jacobian

      if(tsw.eq.1) then

        jr1     = 0.3333333333333333d0*j1/j0
        jr2     = 0.3333333333333333d0*j2/j0

        tb(1,1) = x1 + x12*jr2
        tb(i,j) = y1 + y12*jr2
        tb(j,i) = x2 + x12*jr1
        tb(2,2) = y2 + y12*jr1

c     Center jacobian

      elseif(tsw.eq.2) then

        tb(1,1) = x1
        tb(i,j) = y1
        tb(j,i) = x2
        tb(2,2) = y2

c     Center jacobian inverse

      elseif(tsw.eq.3) then

        j0r     =  1.d0/j0
        tb(1,1) =  y2*j0r
        tb(i,j) = -y1*j0r
        tb(j,i) = -x2*j0r
        tb(2,2) =  x1*j0r

c     Average jacobian inverse

      elseif(tsw.eq.4) then

        jr1     =  0.3333333333333333d0*j1/j0
        jr2     =  0.3333333333333333d0*j2/j0

        tb(2,2) =  x1 + x12*jr2
        tb(i,j) =  y1 + y12*jr2
        tb(j,i) =  x2 + x12*jr1
        tb(1,1) =  y2 + y12*jr1

        j0r     =  1.d0/(tb(1,1)*tb(2,2) - tb(1,2)*tb(2,1))
        tb(1,1) =  tb(1,1)*j0r
        tb(i,j) = -tb(i,j)*j0r
        tb(j,i) = -tb(j,i)*j0r
        tb(2,2) =  tb(2,2)*j0r

      endif ! tsw

c     Compute the inverse to tb

      j0r     =  tb(1,1)*tb(2,2) - tb(1,2)*tb(2,1)
      j0r     =  1.d0/j0r

      ti(1,1) =  tb(2,2)*j0r
      ti(1,2) = -tb(1,2)*j0r
      ti(2,1) = -tb(2,1)*j0r
      ti(2,2) =  tb(1,1)*j0r

      end

      subroutine gstrn1()

      implicit   none

      include   'egeom1.h'
      include   'egeom2.h'
      include   'egeom3.h'

      integer    i,j,k
      real*8     gg(2,2)

      do i = 1,4
        gg(1,1) = aa(1,1)*(nx2(i) - j2*nx0(i))
     &          + aa(4,1)*(ny2(i) - j2*ny0(i))
        gg(1,2) = aa(2,1)*(ny2(i) - j2*ny0(i))
     &          + aa(4,1)*(nx2(i) - j2*nx0(i))
        gg(2,1) = aa(1,2)*(nx1(i) - j1*nx0(i))
     &          + aa(4,2)*(ny1(i) - j1*ny0(i))
        gg(2,2) = aa(2,2)*(ny1(i) - j1*ny0(i))
     &          + aa(4,2)*(nx1(i) - j1*nx0(i))

        do k = 1,2
          do j = 1,4
            a1ig(i,j,k) = ai(j,1)*gg(1,k)
            a2ig(i,j,k) = ai(j,2)*gg(2,k)
          end do ! j
        end do ! k

      end do ! i

      end

      subroutine mestrn(d,sg,xl,ul,ui,tl,ndm,ndf, xsj,shp,
     &                  bmat,bmate,ta,eps)

      implicit   none

      include   'egeom1.h'
      include   'egeom2.h'
      include   'egeom3.h'
      include   'elcoor.h'

      integer    ndm,ndf,i
      real*8     d(*),sg(3),xl(ndm,*),ul(ndf,*),ui(2),tl(*)
      real*8     shp(3,4),xsj(*)
      real*8     bmat(4,4,2),bmate(4,2),eps(4), s1r,s2r, ta

c     Shape functions

      shp(3,1) = 0.25d0*(1.d0 - sg(1))*(1.d0 - sg(2))
      shp(3,2) = 0.25d0*(1.d0 + sg(1))*(1.d0 - sg(2))
      shp(3,3) = 0.25d0*(1.d0 + sg(1))*(1.d0 + sg(2))
      shp(3,4) = 0.25d0*(1.d0 - sg(1))*(1.d0 + sg(2))

c     Reference and current coordinates

      xref(1) = 0.0d0
      xref(2) = 0.0d0
      xcur(1) = 0.0d0
      xcur(2) = 0.0d0
      do i = 1,4
        xref(1) = xref(1) + shp(3,i)*xl(1,i)
        xref(2) = xref(2) + shp(3,i)*xl(2,i)
        xcur(1) = xcur(1) + shp(3,i)*ul(1,i)
        xcur(2) = xcur(2) + shp(3,i)*ul(2,i)
      end do
      xcur(1) = xcur(1) + xref(1)
      xcur(2) = xcur(2) + xref(1)

c     Jacobian factors

      xsj(1) = j0 + sg(1)*j1 + sg(2)*j2
      s1r = sg(1)/xsj(1)
      s2r = sg(2)/xsj(1)

c     Strain displacement matrix for mixed element

      do i = 1,4

        bmat(i,1,1) = nx0(i) + a1ig(i,1,1)*s2r + a2ig(i,1,1)*s1r
        bmat(i,2,1) =          a1ig(i,2,1)*s2r + a2ig(i,2,1)*s1r
        bmat(i,4,1) = ny0(i) + a1ig(i,4,1)*s2r + a2ig(i,4,1)*s1r

        bmat(i,1,2) =          a1ig(i,1,2)*s2r + a2ig(i,1,2)*s1r
        bmat(i,2,2) = ny0(i) + a1ig(i,2,2)*s2r + a2ig(i,2,2)*s1r
        bmat(i,4,2) = nx0(i) + a1ig(i,4,2)*s2r + a2ig(i,4,2)*s1r

      end do ! i

c     Enhanced strain matrix

      bmate(1,1) = ai(1,1)*s1r
      bmate(2,1) = ai(2,1)*s1r
      bmate(4,1) = ai(4,1)*s1r

      bmate(1,2) = ai(1,2)*s2r
      bmate(2,2) = ai(2,2)*s2r
      bmate(4,2) = ai(4,2)*s2r

c     Temperature

      ta = shp(3,1)*tl(1) + shp(3,2)*tl(2)
     &   + shp(3,3)*tl(3) + shp(3,4)*tl(4) - d(9)

c     Strain

      eps(1) = 0.0d0
      eps(2) = 0.0d0
      eps(3) = 0.0d0
      eps(4) = 0.0d0
      do i = 1,4
        eps(1) = eps(1) + bmat(i,1,1)*ul(1,i)
     &                  + bmat(i,1,2)*ul(2,i)
        eps(2) = eps(2) + bmat(i,2,1)*ul(1,i)
     &                  + bmat(i,2,2)*ul(2,i)
        eps(4) = eps(4) + bmat(i,4,1)*ul(1,i)
     &                  + bmat(i,4,2)*ul(2,i)
      end do ! i

c     Enhanced parts

      eps(1) = eps(1) + bmate(1,1)*ui(1) + bmate(1,2)*ui(2)
      eps(2) = eps(2) + bmate(2,1)*ui(1) + bmate(2,2)*ui(2)
      eps(4) = eps(4) + bmate(4,1)*ui(1) + bmate(4,2)*ui(2)

      end

      subroutine stcon2d(hh,gg,dui,ndf,nst,s,p)

c     Purpose: Perform static condensation of internal modes

      implicit none

      integer ndf,nst, i,j, i1,j1
      real*8  hh(2,2), gg(2,2,4), dui(2), s(nst,*), p(ndf,*), tt(2,2)

      save

      i1 = 1
      do i = 1,4

c       Reduce load vector

        p(1,i) = p(1,i) - gg(1,1,i)*dui(1) - gg(2,1,i)*dui(2)
        p(2,i) = p(2,i) - gg(1,2,i)*dui(1) - gg(2,2,i)*dui(2)

c       Reduce stiffness array

        tt(1,1) = gg(1,1,i)*hh(1,1) + gg(2,1,i)*hh(2,1)
        tt(1,2) = gg(1,1,i)*hh(1,2) + gg(2,1,i)*hh(2,2)
        tt(2,1) = gg(1,2,i)*hh(1,1) + gg(2,2,i)*hh(2,1)
        tt(2,2) = gg(1,2,i)*hh(1,2) + gg(2,2,i)*hh(2,2)
        j1 = 1
        do j = 1,4
          s(i1  ,j1  ) = s(i1  ,j1  ) - tt(1,1)*gg(1,1,j)
     &                                - tt(1,2)*gg(2,1,j)

          s(i1  ,j1+1) = s(i1  ,j1+1) - tt(1,1)*gg(1,2,j)
     &                                - tt(1,2)*gg(2,2,j)

          s(i1+1,j1  ) = s(i1+1,j1  ) - tt(2,1)*gg(1,1,j)
     &                                - tt(2,2)*gg(2,1,j)

          s(i1+1,j1+1) = s(i1+1,j1+1) - tt(2,1)*gg(1,2,j)
     &                                - tt(2,2)*gg(2,2,j)
          j1 = j1 + ndf
        end do ! j
        i1 = i1 + ndf
      end do ! i

      end
