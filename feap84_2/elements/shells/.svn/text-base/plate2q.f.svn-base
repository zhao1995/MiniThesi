c$Id:$
      subroutine plate2q(d,ul,xl,ix,s,p,ndf,ndm,nst,isw)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Add rotational lumped masses                     05/03/2010
c       2. Correct coding for mass matrix                   13/01/2011
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose:
c             Quadrilateral plate: 3 dofs per node (w, theta_x, theta_y)
c                                  4 internal bubble on rotation field

c             Mixed approach for shear stiffness.
c             Step 1: Condensation of bubble terms
c             Step 2: Condensation of shear terms

c**********************************************************************

      implicit  none

      include  'bdata.h'
      include  'cdata.h'
      include  'eldata.h'
      include  'elplot.h'
      include  'eltran.h'
      include  'iofile.h'
      include  'mdata.h'
      include  'part0.h'
      include  'prld1.h'
      include  'prstrs.h'
      include  'ptdat2.h'
      include  'rdata.h'
      include  'comblk.h'

      logical   shflag

      integer   ndm, ndf, nst, isw
      integer   i, ii, i1,i2, j, jj, j1,j2, k, l, lint
      integer   ix(*)

      real*8    thk, thk3, q0,qt, xx, yy , xsj, ctan1,ctan3
      real*8    jac(2,2)    , jac0(2,2)   , jinv(2,2) , td(12)

      real*8    k_tt(12,12) , k_wt(12,12)
      real*8    k_st(4,12)  , k_bt(4,12)  , k_bb(4,4)
      real*8    k_bs(4)     , k_ss(4,4)   , k_sw(4,12)

      real*8    bt_b(3,12)  , bt_b1(12,3) , bb_b(3,4) , bb_b1(4,3)
      real*8    bw_s(2,12)  , bwt_s(2,12) , ns(2,4)   , ns1(4,2)

      real*8    app1(4,12)  , app2(4,12)  , bpp1(4)   , bpp2(4)
      real*8    s_hat(4)    , b_hat(4)    , shear(2)

      real*8    d(*),ul(ndf,nen,*),xl(ndm,*), s(nst,*),p(ndf,*)
      real*8    shp(3,4), shpn(3), shpm(3,4), sg(3,25)
      real*8    eps(3), sig(7), dd(3,3), ds(2,2), alphg(3), co(4),si(4)
      real*8    ss(12,12),pp(12)

      save

c     Go to correct array processor

      go to(1,2,3,3,5,3,1,3), isw

c     Return

1     return

c     Check element for errors

2     call ckisop(ix,xl,bt_b,ndm)

      return

c     Compute stiffness

c     Zero all matrices

3     do i = 1,12
        do j = 1,4
          k_bt(j,i) = 0.0d0
          k_st(j,i) = 0.0d0
          k_sw(j,i) = 0.0d0
        end do ! j
        do j = 1,12
          k_tt(j,i) = 0.0d0
        end do ! j
      end do ! i
      do i = 1,4
        do j = 1,4
          k_bb(j,i) = 0.0d0
          k_ss(j,i) = 0.0d0
        end do ! j
        k_bs(i) = 0.0d0
      end do ! i

c     Set shear flag: true = include shear deformation

      shflag = d(79).eq.0.0d0

c     Set time varying pressure load

      if(int(d(74)).gt.0) then
        qt = d(10) + prldv(nint(d(74)))*d(71)
      else
        qt = d(10)*dm
      endif

c     Compute location Gauss points, weights and geometry

      l   = 3
      call int2d (l,lint,sg)
      call geompq(xl,ndm,co,si,jac0)

c     Compute integrals ( LOOP ON GAUSS POINTS )

      do l = 1, lint

        call shpspq(sg(1,l),xl,shp,shpn,shpm,xsj,jac,jinv,ndm)
        call bmatpq(sg(1,l),xl,shp,shpn,shpm,si,co,jac0,jinv,xsj,
     &              bw_s,bt_b,bb_b,bwt_s,ns,ndm)
        xsj = xsj*sg(3,l)

c       Define material constants x jacobian x weight

        call dmatpl(d,d(31),dd,ds,alphg)

c       Bending stiffness

        thk  = d(14)
        thk3 = thk**3/12.d0*xsj
        do j = 1,3
          do i = 1,3
            dd(i,j) = dd(i,j)*thk3
          end do ! i
        end do ! j

c       Shear compliance

        if(shflag) then
          thk3    =  xsj/((ds(1,1)*ds(2,2) - ds(1,2)*ds(2,1))*thk)
          bpp1(1) =  ds(2,2)*thk3
          ds(2,2) =  ds(1,1)*thk3
          ds(1,2) = -ds(1,2)*thk3
          ds(2,1) = -ds(2,1)*thk3
          ds(1,1) =  bpp1(1)
        else
          ds(1,1) =  0.0d0
          ds(1,2) =  0.0d0
          ds(2,1) =  0.0d0
          ds(2,2) =  0.0d0
        endif

        do i = 1,12

c         Construct  "bt_b1 = bt_b * Db"

          bt_b1(i,1) = bt_b(1,i)*dd(1,1)
     &               + bt_b(2,i)*dd(2,1)
     &               + bt_b(3,i)*dd(3,1)
          bt_b1(i,2) = bt_b(1,i)*dd(1,2)
     &               + bt_b(2,i)*dd(2,2)
     &               + bt_b(3,i)*dd(3,2)
          bt_b1(i,3) = bt_b(1,i)*dd(1,3)
     &               + bt_b(2,i)*dd(2,3)
     &               + bt_b(3,i)*dd(3,3)

        end do ! i

        do i = 1,4

c         Construct  "bb_b1 = bb_b * Db"

          bb_b1(i,1) = bb_b(1,i)*dd(1,1)
     &               + bb_b(2,i)*dd(2,1)
     &               + bb_b(3,i)*dd(3,1)
          bb_b1(i,2) = bb_b(1,i)*dd(1,2)
     &               + bb_b(2,i)*dd(2,2)
     &               + bb_b(3,i)*dd(3,2)
          bb_b1(i,3) = bb_b(1,i)*dd(1,3)
     &               + bb_b(2,i)*dd(2,3)
     &               + bb_b(3,i)*dd(3,3)

c         Construct   (N_s)^T * Ds^(-1)

          ns1(i,1)   = ns(1,i) * ds(1,1)
     &               + ns(2,i) * ds(2,1)
          ns1(i,2)   = ns(1,i) * ds(1,2)
     &               + ns(2,i) * ds(2,2)

        end do ! i

c       Construct stiffness matrices with dimension equal to NST

        do i = 1,12
          do j = 1,i
            k_tt(j,i) = k_tt(j,i) + bt_b1(i,1)*bt_b(1,j)
     &                            + bt_b1(i,2)*bt_b(2,j)
     &                            + bt_b1(i,3)*bt_b(3,j)
          end do ! j

          do j = 1,4
            k_sw(j,i) = k_sw(j,i) + bw_s(1,i) *ns(1,j)*xsj
     &                            + bw_s(2,i) *ns(2,j)*xsj
            k_st(j,i) = k_st(j,i) + bwt_s(1,i)*ns(1,j)*xsj
     &                            + bwt_s(2,i)*ns(2,j)*xsj
            k_bt(j,i) = k_bt(j,i) + bt_b1(i,1)*bb_b(1,j)
     &                            + bt_b1(i,2)*bb_b(2,j)
     &                            + bt_b1(i,3)*bb_b(3,j)
          end do ! j
        end do ! i

c       Construct stiffness matrices

        do i = 1, 4
          do j = 1,i

c           Bubble modes

            k_bb(i,j) = k_bb(i,j) + bb_b1(i,1)*bb_b(1,j)
     &                            + bb_b1(i,2)*bb_b(2,j)
     &                            + bb_b1(i,3)*bb_b(3,j)

c           Shear modes

            k_ss(i,j) = k_ss(i,j) - ns1(i,1)*ns(1,j)
     &                            - ns1(i,2)*ns(2,j)
          end do ! j
        end do ! i

c       Consistent load

        q0 = qt * xsj
        do i = 1,4
          k       = mod(i+2 ,4) + 1
          p(1,i) = p(1,i) + q0*shp(3,i)
          p(2,i) = p(2,i) + q0*(shpm(3,i)*co(i) - shpm(3,k)*co(k))
          p(3,i) = p(3,i) + q0*(shpm(3,i)*si(i) - shpm(3,k)*si(k))
        end do ! i

      end do ! l

c     END LOOP ON GAUSS POINTS

c     Generate specific (diagonal) form for K_sb

      k_bs(1) = 16.0d0/9.0d0*(jac0(1,1)*jac0(2,2)-jac0(1,2)*jac0(2,1))
      k_bs(2) = k_bs(1)
      k_bs(3) = k_bs(1)*0.2d0
      k_bs(4) = k_bs(3)

c     Make symmetric parts

      do i = 2,4
        do j = 1,i-1
          k_bb(j,i) = k_bb(i,j)
        end do ! j
      end do ! i

c     Condense stiffness matrices: First (bubble modes)

      call invert(k_bb,4,4)

      do j = 1,4
        do i = 1,j
          k_ss(i,j) = k_ss(i,j) - k_bs(i)*k_bb(i,j)*k_bs(j)
          k_ss(j,i) = k_ss(i,j)
        end do ! i
      end do ! j

      do j = 1,12
        do i = 1,4
          app1(i,j) = k_bb(i,1)*k_bt(1,j) + k_bb(i,2)*k_bt(2,j)
     &              + k_bb(i,3)*k_bt(3,j) + k_bb(i,4)*k_bt(4,j)
        end do ! i
      end do ! j

      do j = 1,12
        do i = 1,4
          k_st(i,j) = k_st(i,j) - k_bs(i)*app1(i,j)
        end do ! i
      end do ! j

      do j = 1,12
        do i = 1,j
          k_tt(i,j) = k_tt(i,j) - k_bt(1,i)*app1(1,j)
     &                          - k_bt(2,i)*app1(2,j)
     &                          - k_bt(3,i)*app1(3,j)
     &                          - k_bt(4,i)*app1(4,j)
        end do ! i
      end do ! j

c     Condense stiffness matrices: Second (shear modes)

      call invert(k_ss,4,4)

      do j = 1,12
        do i = 1,4
          app1(i,j) = k_ss(i,1)*k_sw(1,j) + k_ss(i,2)*k_sw(2,j)
     &              + k_ss(i,3)*k_sw(3,j) + k_ss(i,4)*k_sw(4,j)
          app2(i,j) = k_ss(i,1)*k_st(1,j) + k_ss(i,2)*k_st(2,j)
     &              + k_ss(i,3)*k_st(3,j) + k_ss(i,4)*k_st(4,j)
        end do ! i
      end do ! j

      do j = 1,12
        do i = 1,j
          k_tt(i,j) = k_tt(i,j)
     &              - k_sw(1,i)*app1(1,j) - k_sw(2,i)*app1(2,j)
     &              - k_sw(3,i)*app1(3,j) - k_sw(4,i)*app1(4,j)
     &              - k_st(1,i)*app2(1,j) - k_st(2,i)*app2(2,j)
     &              - k_st(3,i)*app2(3,j) - k_st(4,i)*app2(4,j)
        end do ! i

        do i = 1,12
          k_wt(i,j) = k_sw(1,i)*app2(1,j) + k_sw(2,i)*app2(2,j)
     &              + k_sw(3,i)*app2(3,j) + k_sw(4,i)*app2(4,j)
        end do ! i
      end do ! j

c     Form residual and assemble stiffness

      if (isw.eq.3 .or.isw.eq.6) then

        ctan1 = ctan(1) + d(78)*ctan(2)
        if(ndfo(1).gt.0 .or. shflg) then
          ctan3 = ctan(3) + d(77)*ctan(2)
        else
          ctan3 = 0.0d0
        endif

c       Accumulate stiffness parts

        do j = 1,12
          do i = 1,j
            k_tt(i,j) = k_tt(i,j) - k_wt(i,j) - k_wt(j,i)
            k_tt(j,i) = k_tt(i,j)
          end do ! i
        end do ! j

c       Assemble element array and residual

        ii = 0
        i1 = 0
        do i = 1,4
          jj = 0
          j1 = 0
          do j = 1,4
            do k = 1,3
              do l = 1,3
c               Compute residual
                p(k,i) = p(k,i) - k_tt(ii+k,jj+l)*(ul(l,j,1)
     &                                     + d(78)*ul(l,j,4))
c               Assemble element tangent matrix
                s(i1+k,j1+l) = k_tt(ii+k,jj+l)*ctan1
              end do ! l
            end do ! k
            jj = jj + 3
            j1 = j1 + ndf
          end do ! j
          ii = ii + 3
          i1 = i1 + ndf
        end do ! i

c       Add inertial parts if necessary

        if(ctan3.ne.0.0d0) then

          call massqp(d,xl,ndm,3,12, pp,ss)

          i1 = 0
          i2 = 0
          do i = 1,4
            j1 = 0
            j2 = 0
            do j = 1,4
              do k = 1,3
                do l = 1,3
                  p(k,i)       = p(k,i) - ss(i2+k,j2+l)*(ul(l,j,5)
     &                                           + d(77)*ul(l,j,4))
                  s(i1+k,j1+l) = s(i1+k,j1+l) + ctan3*ss(i2+k,j2+l)
                end do ! l
              end do ! k
              j1 = j1 + ndf
              j2 = j2 + 3
            end do ! j
            i1 = i1 + ndf
            i2 = i2 + 3
          end do ! i

        endif

c     Recover stress modes

      elseif ((isw.eq.4).or.(isw.eq.8)) then

c       Set local displacement order

        ii = 0
        do i=1,4
          do j = 1,3
            td(ii+j) = ul(j,i,1)
          end do ! j
          ii = ii + ndf
        end do ! i

c       Multiply stiffness order

        do i=1,4
          s_hat(i) = 0.0d0
          b_hat(i) = 0.0d0
          bpp1(i)  = 0.0d0
          bpp2(i)  = 0.0d0
        end do ! i

        do j=1,12
          do i=1,4
            bpp1(i) = bpp1(i) + (k_sw(i,j) + k_st(i,j))*td(j)
            bpp2(i) = bpp2(i) +  k_bt(i,j)             *td(j)
          end do ! i
        end do ! j

        do j=1,4
          do i=1,4
            s_hat(i) = s_hat(i) - k_ss(i,j)*bpp1(j)
          end do ! i
        end do ! j

        do j=1,4
          bpp2(j) = bpp2(j) + k_bs(j)*s_hat(j)
          do i=1,4
            b_hat(i) = b_hat(i) - k_bb(i,j)*bpp2(j)
          end do ! i
        end do ! j

c       Compute curvatures and moments

        if(isw.eq.4 .or. nsplts.gt.0) then

c         Compute Gauss points, weights and geometry

          l   = 1
          call int2d (l,lint,sg)
          call geompq(xl,ndm,co,si,jac0)

          do l=1,lint

            call shpspq(sg(1,l),xl,shp,shpn,shpm,xsj,
     &                  jac,jinv,ndm)
            call bmatpq(sg(1,l),xl,shp,shpn,shpm,si,co,jac0,jinv,
     &                  xsj,bw_s,bt_b,bb_b,bwt_s,ns,ndm)

            call dmatpl(d,d(31),dd,ds,alphg)

            thk3 = d(14)**3/12.d0
            do j = 1,3
              do i = 1,3
                dd(i,j) = dd(i,j)*thk3
              end do ! i
            end do ! j

            call strepq(dd,bt_b,bb_b,b_hat,ns,s_hat,ul,ndf,
     &                  eps,sig,shear)

            if(isw.eq.4) then

              xx = shp(3,1)*xl(1,1) + shp(3,2)*xl(1,2)
     &           + shp(3,3)*xl(1,3) + shp(3,4)*xl(1,4)
              yy = shp(3,1)*xl(2,1) + shp(3,2)*xl(2,2)
     &           + shp(3,3)*xl(2,3) + shp(3,4)*xl(2,4)

c             Output moments and curvatures

              mct = mct -2
              if(mct.le.0) then
                write(iow,2001) o,head
                if(ior.lt.0) write (*,2001) o,head
                mct = 50
              endif
              write(iow,2002) n,ma,(sig(j),j=1,5),xx,yy,eps,sig(6),
     &                        shear
              if(ior.lt.0) then
                write(*,2002) n,ma,(sig(j),j=1,5),xx,yy,eps,sig(6),
     &                        shear
              endif

c           Save for tplot outputs

            else
              ii       = 5*(l-1)
              tt(ii+1) = sig(1)
              tt(ii+2) = sig(2)
              tt(ii+3) = sig(3)
              tt(ii+4) = shear(1)
              tt(ii+5) = shear(2)
            endif

          end do ! l

c       Compute nodal stress values

        else
          call stcnpq(d,xl,ul,shp,p,s,s_hat,b_hat,ndf,ndm,nel)
        end if

      endif
      return

c     Compute consistent and lumped mass matrix

5     call massqp(d,xl,ndm,ndf,nst, p,s)
      return

c     Formats for input-output

2001  format(a1,20a4//5x,'Element Moments'//'  Element Material',
     &   3x,'11-Moment',3x,'22-Moment',3x,'12-Moment',4x,
     &   '1-Moment',4x,'2-Moment'/2x,'1-Coord',2x,'2-Coord',3x,
     &   '11-Strain',3x,'22-Strain',3x,'12-Strain',12x,'Angle'/
     &   21x,'Shear_x  ',3x,'Shear_y')

2002  format(2i9,1p,5e12.3/0p,2f9.3,1p,3e12.3,0p,f18.2/18x,1p,2e12.3/1x)

      end

      subroutine stcnpq(d,xl,ul,shp,dt,st,s_hat,b_hat,ndf,ndm,nel)

c     Lumped and consistent projection routine

      implicit  none

      include  'cdata.h'
      include  'strnum.h'
      include  'iofile.h'

      integer   ndf, ndm, nel
      integer   i, j, l, lint
      real*8    thk, xsj, xg

      real*8    dt(*),st(nen,*),xl(ndm,*),shp(3,4),d(*)
      real*8    ul(ndf,*),sg(3,25)
      real*8    shpn(3), shpm(3,4), jac0(2,2), jac(2,2), jinv(2,2)
      real*8    eps(3), sig(7), dd(3,3), ds(2,2), alphg(3), co(4),si(4)

      real*8    bw_s(2,12), bt_b(3,12), bb_b(3,4), bwt_s(2,12)
      real*8    ns(2,4), s_hat(4), b_hat(4), shear(2)

      save

c     Zero all matrices

      do i = 1,12
        do j = 1,3
          bt_b(j,i) = 0.0d0
        end do ! j
      end do ! i

c     Compute Gauss points, weights and geometry

      l   = 3
      call int2d (l,lint,sg)
      call geompq(xl,ndm,co,si,jac0)

      do j = 1,nel
        dt(j) = 0.0d0
      end do ! j

      do l=1,lint

        call shpspq(sg(1,l),xl,shp,shpn,shpm,xsj,jac,jinv,ndm)
        call bmatpq(sg(1,l),xl,shp,shpn,shpm,si,co,jac0,jinv,xsj,
     &              bw_s,bt_b,bb_b,bwt_s,ns,ndm)

        xsj = xsj*sg(3,l)

        call dmatpl(d,d(31),dd,ds,alphg)

        thk = d(14)**3/12.d0
        do j = 1,3
          do i = 1,3
            dd(i,j) = dd(i,j)*thk
          end do ! i
        end do ! j

        call strepq(dd,bt_b,bb_b,b_hat,ns,s_hat,ul,ndf,eps,sig,shear)

c       Compute lumped projection and assemble stress integrals

        do j = 1,nel
          xg      = xsj*shp(3,j)
          dt(j)   = dt(j)     +        xg
          st(j,1) = st(j,1)   + sig(1)*xg
          st(j,2) = st(j,2)   + sig(2)*xg
          st(j,4) = st(j,4)   + sig(3)*xg
          st(j,5) = st(j,5)   + shear(1)*xg
          st(j,6) = st(j,6)   + shear(2)*xg
        end do ! j

      end do ! l

      iste = 6

      end

      subroutine shpspq(xi,xl,shp,shpn,shpm,xsj,jac,jinv,ndm)

c     Shape functions with linked edges

      implicit  none

      include  'eldata.h'

      integer   i, j, k, ndm

      real*8    xl(ndm,4),shp(3,4),shpn(3),shpm(3,4),jac(2,2),jinv(2,2)
      real*8    xi(2),xi1p,xi1m,eta1m,eta1p,xi2,eta2,xsj,xsjinv

      save

c     Set up parameters

      xi1p  = 1.0d0 + xi(1)
      xi1m  = 1.0d0 - xi(1)
      eta1p = 1.0d0 + xi(2)
      eta1m = 1.0d0 - xi(2)
      xi2   = xi1p  * xi1m
      eta2  = eta1p * eta1m

c     Natural coordinate shape functions

      shp(3,1) =  0.25d0*xi1m*eta1m
      shp(3,2) =  0.25d0*xi1p*eta1m
      shp(3,3) =  0.25d0*xi1p*eta1p
      shp(3,4) =  0.25d0*xi1m*eta1p

c     Natural coordinate derivatives for mid-surface

      shp(1,1) = -0.25d0*eta1m
      shp(1,2) = -shp(1,1)
      shp(1,3) =  0.25d0*eta1p
      shp(1,4) = -shp(1,3)

      shp(2,1) = -0.25d0*xi1m
      shp(2,2) = -0.25d0*xi1p
      shp(2,3) = -shp(2,2)
      shp(2,4) = -shp(2,1)

c     Construct jacobian-transpose and its inverse

      do i = 1,2
        do j = 1,2
          jac(i,j) = 0.0d0
          do k = 1,nel
            jac(i,j) = jac(i,j) + xl(j,k)*shp(i,k)
          end do ! k
        end do ! j
      end do ! i

      xsj = jac(1,1)*jac(2,2)-jac(1,2)*jac(2,1)

      if(xsj.lt.0.0d0) then
        write(*,*) 'NEGATIVE JACOBIAN: Element = ',n
        xsj = -xsj
      elseif(xsj.eq.0.0d0) then
        write(*,*) 'ZERO JACOBIAN: Element = ',n
        call plstop()
      endif

      xsjinv    = 1.d0/xsj
      jinv(1,1) = jac(2,2)*xsjinv
      jinv(2,2) = jac(1,1)*xsjinv
      jinv(1,2) =-jac(1,2)*xsjinv
      jinv(2,1) =-jac(2,1)*xsjinv

c     Shape function 'N' and derivatives

      shpn(1)   = - 2.0d0*xi(1)*eta2
      shpn(2)   = - 2.0d0*xi(2)*xi2
      shpn(3)   =   xi2 * eta2

c     Shape function 'M' and derivatives (0.125=1/8; 0.0625=1/16)

      shpm(3,1) =   xi2*eta1m * 0.0625d0
      shpm(3,2) =   xi1p*eta2 * 0.0625d0
      shpm(3,3) =   xi2*eta1p * 0.0625d0
      shpm(3,4) =   xi1m*eta2 * 0.0625d0

      shpm(1,1) = - xi(1) * eta1m * 0.125d0
      shpm(1,2) =   eta2* 0.0625d0
      shpm(1,3) = - xi(1) * eta1p * 0.125d0
      shpm(1,4) = - eta2* 0.0625d0

      shpm(2,1) = - xi2* 0.0625d0
      shpm(2,2) = - xi(2) * xi1p * 0.125d0
      shpm(2,3) =   xi2* 0.0625d0
      shpm(2,4) = - xi(2) * xi1m * 0.125d0

      end

      subroutine geompq(xl,ndm,co,si,jac0)

c     Compute tangential directions
c     ATT.: co() and si() and cosine and sine times side length

      implicit  none

      integer   ndm, i, j
      real*8    xl(ndm,*), co(4), si(4), jac0(2,2)

      save

c     Compute side length * angles

      do i = 1,4
        j     =   mod(i,4) + 1
        co(i) = - xl(2,i)  + xl(2,j)
        si(i) =   xl(1,i)  - xl(1,j)
      end do ! i

c     Compute element center jacobian

      jac0(1,1) = 0.25d0*(-xl(1,1) + xl(1,2) + xl(1,3) - xl(1,4))
      jac0(2,1) = 0.25d0*(-xl(1,1) - xl(1,2) + xl(1,3) + xl(1,4))
      jac0(1,2) = 0.25d0*(-xl(2,1) + xl(2,2) + xl(2,3) - xl(2,4))
      jac0(2,2) = 0.25d0*(-xl(2,1) - xl(2,2) + xl(2,3) + xl(2,4))

      end

      subroutine bmatpq(xi,xl,shp,shpn,shpm,si,co,jac0,jinv,xsj,
     &                  bw_s,bt_b,bb_b,bwt_s,ns,ndm)

c     Compute all B-strain-type matrices for mixed formulation

      implicit  none

      integer   ndm

      real*8    xi(2)
      real*8    xl(ndm,*), shp(3,4), shpn(3), shpm(3,4), si(4), co(4)
      real*8    jac0(2,2), jinv(2,2), xsj, xsjinv, xsjinv2
      real*8    bw_s(2,3,4), bt_b(3,3,4), bb_b(3,*), bwt_s(2,3,4)
      real*8    ns(2,*)

      integer   i, c

      real*8    b1(4), b2(4), f1(4),f1c(4),f1s(4), f2(4),f2c(4),f2s(4)
      real*8    aa, bb, dj1, dj2
      real*8    Nj, Nj_xi, Nj_eta, Nj_x, Nj_y
      real*8    xiNj_xi, xiNj_eta, etaNj_xi, etaNj_eta
      real*8    N_x, N_y, xiNj_x, xiNj_y, etaNj_x, etaNj_y
      real*8    ax, ay, bx, by, cx, cy, dx, dy

      save

c     Construct parameters for B-strain-type matrices

c     Parameters for linear shape functions and for rotational

c     Contribution to transverse displacement

      do i = 1,4
        b1(i)   = jinv(1,1)*shp(1,i)  + jinv(1,2)*shp(2,i)
        b2(i)   = jinv(2,1)*shp(1,i)  + jinv(2,2)*shp(2,i)
        f1(i)   = jinv(1,1)*shpm(1,i) + jinv(1,2)*shpm(2,i)
        f1c(i)  = f1(i)*co(i)
        f1s(i)  = f1(i)*si(i)
        f2(i)   = jinv(2,1)*shpm(1,i) + jinv(2,2)*shpm(2,i)
        f2c(i)  = f2(i)*co(i)
        f2s(i)  = f2(i)*si(i)
      end do ! i

c     Parameters for bubble functions

      dj1       = (xl(1,1) - xl(1,2) + xl(1,3) - xl(1,4))*0.25d0
      dj2       = (xl(2,1) - xl(2,2) + xl(2,3) - xl(2,4))*0.25d0
      aa        = jac0(1,1)*dj2 - jac0(1,2)*dj1
      bb        = jac0(2,2)*dj1 - jac0(2,1)*dj2

      xsjinv    = 1.0d0 / xsj
      xsjinv2   = xsjinv * xsjinv

      N_x       = jinv(1,1)*shpn(1) + jinv(1,2)*shpn(2)
      N_y       = jinv(2,1)*shpn(1) + jinv(2,2)*shpn(2)

      Nj        = shpn(3) * xsjinv
      Nj_xi     = xsjinv2 * ( shpn(1)*xsj - shpn(3)*aa)
      Nj_eta    = xsjinv2 * ( shpn(2)*xsj - shpn(3)*bb)
      Nj_x      = jinv(1,1)*Nj_xi + jinv(1,2)*Nj_eta
      Nj_y      = jinv(2,1)*Nj_xi + jinv(2,2)*Nj_eta

      xiNj_xi   = Nj + xi(1)*Nj_xi
      xiNj_eta  =      xi(1)*Nj_eta
      etaNj_xi  =      xi(2)*Nj_xi
      etaNj_eta = Nj + xi(2)*Nj_eta

      xiNj_x    = jinv(1,1)*xiNj_xi + jinv(1,2)*xiNj_eta
      xiNj_y    = jinv(2,1)*xiNj_xi + jinv(2,2)*xiNj_eta
      etaNj_x   = jinv(1,1)*etaNj_xi + jinv(1,2)*etaNj_eta
      etaNj_y   = jinv(2,1)*etaNj_xi + jinv(2,2)*etaNj_eta

      ax        =  etaNj_x*jac0(2,1)
      ay        =  etaNj_y*jac0(2,1)
      bx        = -xiNj_x*jac0(1,1)
      by        = -xiNj_y*jac0(1,1)

      cx        =  etaNj_x*jac0(2,2)
      cy        =  etaNj_y*jac0(2,2)
      dx        = -xiNj_x*jac0(1,2)
      dy        = -xiNj_y*jac0(1,2)

c     Construct B-strain-type matrices

      do i = 1,4
        bw_s(1,1,i)  =   b1(i)                            !  B_w_s
        bw_s(2,1,i)  =   b2(i)

        bt_b(1,3,i)  =   b1(i)                            !  B_t_b
        bt_b(2,2,i)  = - b2(i)
        bt_b(3,2,i)  = - b1(i)
        bt_b(3,3,i)  =   b2(i)

        c            =   mod(i+2,4) + 1                   !  B_wt_s
        bwt_s(1,2,i) =           f1c(i) - f1c(c)
        bwt_s(1,3,i) =  shp(3,i)+f1s(i) - f1s(c)
        bwt_s(2,2,i) = -shp(3,i)+f2c(i) - f2c(c)
        bwt_s(2,3,i) =           f2s(i) - f2s(c)

      end do ! i

      bb_b(1,1) =   jac0(2,2) * Nj_x                     !    B_b_b
      bb_b(1,2) = - jac0(1,2) * Nj_x
      bb_b(2,1) = - jac0(2,1) * Nj_y
      bb_b(2,2) =   jac0(1,1) * Nj_y
      bb_b(3,1) = - jac0(2,1) * Nj_x + jac0(2,2) * Nj_y
      bb_b(3,2) =   jac0(1,1) * Nj_x - jac0(1,2) * Nj_y

      bb_b(1,3) =   cx
      bb_b(1,4) =   dx
      bb_b(2,3) = - ay
      bb_b(2,4) = - by
      bb_b(3,3) = - ax + cy
      bb_b(3,4) = - bx + dy

c     Construct N_s

      ns(1,1)   = jac0(1,1)
      ns(1,2)   = jac0(2,1)
      ns(2,1)   = jac0(1,2)
      ns(2,2)   = jac0(2,2)
      ns(1,3)   = jac0(1,1)*xi(2)
      ns(1,4)   = jac0(2,1)*xi(1)
      ns(2,3)   = jac0(1,2)*xi(2)
      ns(2,4)   = jac0(2,2)*xi(1)

      end

      subroutine strepq(dd,bt_b,bb_b,b_hat,ns,s_hat,ul,ndf,
     &                  eps,sig,shear)

c     Compute curvatures [eps(1)=eps_x,eps(2)=eps_y,eps(3)=gamma_xy]
c             stresses   [sig(1)=mom_x,sig(2)=mom_y,sig(3)=mom_xy]
c                         shear(1)=q_x,shear(2)=q_y]
      implicit  none

      integer   i,j,ndf
      real*8    dd(3,3)    , eps(3)   , sig(7)   , shear(2)
      real*8    bt_b(3,3,4), bb_b(3,4), ns(2,4)
      real*8    b_hat(4)   , s_hat(4) , ul(ndf,4)

      save

c     Initialize

      eps(1)   = 0.0d0
      eps(2)   = 0.0d0
      eps(3)   = 0.0d0
      shear(1) = 0.0d0
      shear(2) = 0.0d0

c     Compute bending strains

      do i=1,4
        do j = 1,3
          eps(1) = eps(1) + bt_b(1,j,i)*ul(j,i)
          eps(2) = eps(2) + bt_b(2,j,i)*ul(j,i)
          eps(3) = eps(3) + bt_b(3,j,i)*ul(j,i)
        end do ! j
      end do ! i
      do i=1,4
        eps(1) = eps(1) + bb_b(1,i)*b_hat(i)
        eps(2) = eps(2) + bb_b(2,i)*b_hat(i)
        eps(3) = eps(3) + bb_b(3,i)*b_hat(i)
      end do ! i

c     Compute moments

      sig(1) = dd(1,1)*eps(1) + dd(1,2)*eps(2) + dd(1,3)*eps(3)
      sig(2) = dd(2,1)*eps(1) + dd(2,2)*eps(2) + dd(2,3)*eps(3)
      sig(3) = 0.0d0
      sig(4) = dd(3,1)*eps(1) + dd(3,2)*eps(2) + dd(3,3)*eps(3)

      call pstr2d(sig,sig(5))
      do i = 3,6
        sig(i) = sig(i+1)
      end do ! i

c     Compute shears

      do i=1,4
        shear(1) = shear(1) + ns(1,i)*s_hat(i)
        shear(2) = shear(2) + ns(2,i)*s_hat(i)
      end do ! i

      end

      subroutine massqp(d,xl,ndm,ndf,nst, p,s)

      implicit  none

      integer   ndm,ndf,nst, i,i1,ii, j,j1,jj, l,lint
      real*8    d(*),xl(ndm,*), p(ndf,*),s(nst,nst)
      real*8    co(4),si(4),jac0(2,2),an(3,4),sg(3,25)
      real*8    jac(2,2),jinv(2,2),shp(3,4),shpm(3,4),shpn(3)
      real*8    xsj,ccm,clm, ar24, h12

      save

c     Initialize arrays

      do j = 1,nst
        do i = 1,nst
          s(i,j) = 0.0d0
        end do ! i
        p(j,1) = 0.0d0
      end do ! j

c     Compute location Gauss points, weights and geometry

      l   = 4
      call int2d (l,lint,sg)
      call geompq(xl,ndm,co,si,jac0)

c     Compute integrals ( LOOP ON GAUSS POINTS )

      do l = 1, lint

        call shpspq(sg(1,l),xl,shp,shpn,shpm,xsj,jac,jinv,ndm)
        xsj = xsj*sg(3,l)*d(4)*d(14)

        ccm = d(7)*xsj
        clm = xsj - ccm

c       Lumped mass matrix part

        i1 = 1
        do i=1,4
          p(1,i)   = p(1,i) + shp(3,i)*xsj
          s(i1,i1) = s(i1,i1) + clm*shp(3,i)
          i1       = i1 + ndf
        end do ! i

c       Consistent mass interpolation functions

        do i=1,4
          j       = mod(i,4) + 1
          i1      = ndf*(j-1)
          an(1,j) = shp(3,j)
          an(2,j) = shpm(3,j)*co(j) - shpm(3,i)*co(i)
          an(3,j) = shpm(3,j)*si(j) - shpm(3,i)*si(i)
        end do ! i

c       Consistent mass matrix part

        do i = 1,4
          i1 = ndf*(i-1)
          do ii = 1,3
            ar24 = an(ii,i)*ccm
            do j = 1,4
              j1 = ndf*(j-1)
              do jj = 1,3
                s(i1+ii,j1+jj) = s(i1+ii,j1+jj) +ar24*an(jj,j)
              end do ! jj
            end do ! ii
          end do ! j
        end do ! i

      end do ! l

c     Add rotational mass

      h12 = d(8)*d(14)*d(14)/12.0d0
      do i = 1,4
        p(2,i) = p(1,i)*h12
        p(3,i) = p(2,i)
      end do ! i


      end
