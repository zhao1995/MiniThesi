c$Id:$
      subroutine shl3di(d,ul,xl,ix,s,r,ndf,ndm,nst,isw)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Add computation of xref and xcur in elcoor.h     26/07/2009
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose:  Quadrilateral shell element for feap
c                incorporating membrane with normal drilling dof
c                modified to remove effects of constant strains and
c                including deep shell curvature corrections to the
c                discrete kirchhoff quadrilateral plate bending element

c         June 14, 1986  -  Basic developments with DKQ and JET mods
c         Oct. 22, 1989  -  Add rotation rank corrections
c         Dec  12, 1991  -  Add triangular element capability
c         Nov  10, 1995  -  Input from 'inmate', orthotropic
c         April 7, 1998  -  Through thickness integrations

c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c.... Input parameters set as follows:

c         ndm = 3 (x,y,z cartesian coordinates at nodes)
c         ndf = 6 (u-x,u-y,u-z,r-x,r-y,r-z at nodes)
c         nen = 3 or 4 nodes (counterclockwise around element)

c       Note: 1-direction bisects diagonals between 2-3 element and
c             2-direction bisects diagonals between 3-4 element nodes.
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      implicit  none

      include  'bdata.h'
      include  'cdata.h'
      include  'eldata.h'
      include  'elplot.h'
      include  'eltran.h'
      include  'evdata.h'
      include  'iofile.h'
      include  'part0.h'
      include  'prstrs.h'
      include  'rdata.h'
      include  'shlc16.h'
      include  'shld16.h'
      include  'shpf16.h'
      include  'sstr16.h'
      include  'comblk.h'

      integer   ndm,ndf,nst,isw
      integer   i,j,k,l,lint,ll,lt,ialph
      integer   i1,i2, j1,j2, ii,jj, nn

      real*8    a11,a13, a23,a33, xx,yy,zz, hh
      real*8    dv, dv1,dv2,dvm, pen,ggi,ggv,tc, xsj, thk
      real*8    xyn, x1n,x2n,y1n,y2n, ctan1,ctan3
      real*8    shp1i,shp2i,shp3i
      real*8    shp13, shp23

      integer   ix(*)

      real*8    d(*),xl(ndm,*),ul(ndf,nen,*),s(nst,*),r(ndf,*)
      real*8    dvl(9),eps(6),sig(10),ss(6,4),rr(6,4)
      real*8    yl(3,4),vl(6,4,3),wl(6,4,3),tr(3,3),bl(3),bg(3)
      real*8    gshp1(3,4),gshp2(3,4),gg(6,4),dd(6,6,5)
      real*8    bmat(4,6,4), dbmat(6,4), sc(24,24),rc(6,24)
      real*8    bf(3),bt(3,3),x3(3)

      save

      data      ialph /2/

c     Transfer to correct processor

      go to (1,2,3,3,3,3,1,3), isw

1     return

c     Check element for errors

2     call tran3d(xl,yl,tr,ndm,nel)
      call ckisop(ix,yl,shp,3)
      return

c     Compute element tangent array

c     Compute transformation and midsurface coords

3     call tran3d(xl,yl,tr,ndm,nel)

c     Compute local coordinates

      do i = 1,nel
        do j = 1,3
          vl(j  ,i,1) = 0.0d0
          vl(j+3,i,1) = 0.0d0
          do k = 1,3
            vl(j  ,i,1) = vl(j  ,i,1) + tr(j,k)*ul(k  ,i,1)
            vl(j+3,i,1) = vl(j+3,i,1) + tr(j,k)*ul(k+3,i,1)
          end do ! k
        end do ! j
        do j = 1,6
          wl(j,i,1) = ul(j,i,1)
          wl(j,i,2) = ul(j,i,4)
          wl(j,i,3) = ul(j,i,5)
        end do ! j
      end do ! i
      do i = 1,nel
        do j = 1,6
          rc(j,i) = 0.0d0
          rr(j,i) = 0.0d0
          ss(j,i) = 0.0d0
        end do ! j
      end do ! i
      do i = 1,24
        do j = 1,24
          sc(j,i) = 0.0d0
        end do ! j
      end do ! i

c     Get quadratrure data

      l   = nint(d(5))
      lt  = nint(d(102))
      call int2d (l ,lint,sg)
      call int1dl(lt,sgt)

c     Test for triangular element

      if(nel.eq.3) then
        qdflg = .false.
      elseif( ix(1) .eq. ix(2)  .or. ix(2) .eq. ix(3) .or.
     &        ix(3) .eq. ix(4)  .or. ix(4) .eq. ix(1) ) then
        qdflg = .false.
      else
        qdflg = .true.
      endif

      if(qdflg) then
        call pzero(gshp1,12)
        call pzero(gshp2,12)
        call pzero(gg   ,24)
        dv    = 0.0d0
        pen   = d(60)
        call jacq3d(yl)
      else
        call pzero(  bm, 162)
        call pzero(shp1, 12*lint)
        call pzero(shp2, 12*lint)
        shp13 = 0.0d0
        shp23 = 0.0d0
        a13   = 0.0d0
        a23   = 0.0d0
        a33   = 0.0d0
        x1n   = 0.0d0
        y1n   = 0.0d0
        x2n   = 0.0d0
        y2n   = 0.0d0
        xyn   = 0.0d0
        call jtri3d(yl,xsjt)
        yl(3,1) = 0.0d0
      endif

      if(yl(3,1).eq.0.0d0) then
        ii1 = 3
        ii2 = 5
      else
        ii1 = 1
        ii2 = 6
      endif

c     Construct integrals of drilling shape functions

      do l = 1,lint

c       Form shape functions and their integrals

        if(qdflg) then
          call rshp3d(sg(1,l),yl,shp(1,1,l),shp1(1,1,l),
     &                shp2(1,1,l),xsj,3)
          dvl(l) = xsj*sg(3,l)
          dv     = dv + dvl(l)
          do j = 1,4
            do i = 1,3
              gshp1(i,j) = gshp1(i,j) + shp1(i,j,l)*dvl(l)
              gshp2(i,j) = gshp2(i,j) + shp2(i,j,l)*dvl(l)
            end do ! i
          end do ! j
        else
          call shp2d (sg(1,l),yl,shp(1,1,l),xsj,3,nel,ix,.false.)
          dvl(l) = xsj*sg(3,l)
        endif
      end do ! l

      if(isw.eq.5) go to 5

c     Compute thickness for element

      thk  = d(14)

      if(isw.eq.4) go to 4

c     Construct modified drilling shape functions

      if(qdflg) then
        do j = 1,4
          do i = 1,3
            dv1 = gshp1(i,j)/dv
            dv2 = gshp2(i,j)/dv
            do l = 1,lint
              shp1(i,j,l) = shp1(i,j,l) - dv1
              shp2(i,j,l) = shp2(i,j,l) - dv2
            end do ! l
          end do ! j
        end do ! i
      endif

      if(isw.eq.8) go to 8

c     Compute element load vectors

      bl(1) = 0.0d0
      bl(2) = 0.0d0

c     Set body loading factors

      call sbodyf(d, bg)

      if(int(d(76)).gt.0) then
        bl(3) = d(10)
      else
        bl(3) = d(10)*dm
      endif

c     Multiply body loads by thickness

      do i = 1,3
        bg(i) = bg(i)*thk
      end do ! i

c     Transform body loads to local system

      do i = 1,3
        do j = 1,3
          bl(i) = bl(i) + tr(i,j)*bg(j)
        end do ! j
      end do ! i

      tc = 1.0d0

c     Transform displacements

      do i = 1,nel
        vl(1,i,1) = vl(1,i,1) - yl(3,i)*vl(5,i,1)
        vl(2,i,1) = vl(2,i,1) + yl(3,i)*vl(4,i,1)
      end do ! i

c     Transient factors

      ctan1 = ctan(1) + d(78)*ctan(2)
      if(ndfo(1).gt.0 .or. shflg) then
        ctan3 = ctan(3) + d(77)*ctan(2)
      else
        ctan3 = 0.0d0
      endif

c     Compute local velocity and acceleration

      do i = 1,nel
        do j = 1,3
          vl(j  ,i,2) = 0.0d0
          vl(j+3,i,2) = 0.0d0
          vl(j  ,i,3) = 0.0d0
          vl(j+3,i,3) = 0.0d0
          do k = 1,3
            vl(j  ,i,2) = vl(j  ,i,2) + tr(j,k)*wl(k  ,i,2)
            vl(j+3,i,2) = vl(j+3,i,2) + tr(j,k)*wl(k+3,i,2)
            vl(j  ,i,3) = vl(j  ,i,3) + tr(j,k)*wl(k  ,i,3)
            vl(j+3,i,3) = vl(j+3,i,3) + tr(j,k)*wl(k+3,i,3)
          end do ! k
        end do ! j
      end do ! i

c     Compute parts

      nn = 0
      do l = 1,lint

c       Compute loading term

        do i = 1,nel
          shp3i   = shp(3,i,l)
          rc(1,i) = rc(1,i) + shp3i*bl(1)*dvl(l)
          rc(2,i) = rc(2,i) + shp3i*bl(2)*dvl(l)
          rc(3,i) = rc(3,i) + shp3i*bl(3)*dvl(l)
          if(qdflg) then
            shp13   = shp1(3,i,l)
            shp23   = shp2(3,i,l)
            rc(4,i) = rc(4,i) + shp13*bl(3)*dvl(l)*tc
            rc(5,i) = rc(5,i) + shp23*bl(3)*dvl(l)*tc
            rc(6,i) = rc(6,i) -(shp13*bl(1) + shp23*bl(2))*dvl(l)*tc
          end if
        end do ! i

c       Membrane and bending stiffness part

        dv1 = 0.5d0*thk*dvl(l)*ctan1

        do ll = 1,lt

          hh = 0.5d0*thk*sgt(1,ll)
          call str3di(d,yl,vl,3,nel,l, xx,yy,zz,hh,eps,sig,dd,nn,isw)
          nn = nn + nint(d(15))

c         Store time history plot data for element

          i = 8*(l-1)
          do j = 1,4
            tt(j+i  ) = sig(j)
            tt(j+i+4) = eps(j)
          end do ! j

c         Multiply stress and tangent by volume and quadratrue weight

          dv2  = dv1*sgt(2,ll)
          do i = 1,4
            sig(i) = sig(i)*dv2
            do j = 1,4
              dd(j,i,1) = dd(j,i,1)*dv2
            end do ! j
          end do ! i

c         Recover previously computed shape functions

          do i = 1,nel
            shp1i       =  shp(1,i,l)
            shp2i       =  shp(2,i,l)

            bmat(1,1,i) =  bm(1,1,i)*hh + shp1i
            bmat(2,1,i) =  bm(2,1,i)*hh
            bmat(4,1,i) =  bm(3,1,i)*hh + shp2i
            bmat(1,2,i) =  bm(1,2,i)*hh
            bmat(2,2,i) =  bm(2,2,i)*hh + shp2i
            bmat(4,2,i) =  bm(3,2,i)*hh + shp1i
            bmat(1,3,i) =  bm(1,3,i)*hh
            bmat(2,3,i) =  bm(2,3,i)*hh
            bmat(4,3,i) =  bm(3,3,i)*hh
            bmat(1,4,i) =  bm(1,4,i)*hh
            bmat(2,4,i) =  bm(2,4,i)*hh
            bmat(4,4,i) =  bm(3,4,i)*hh
            bmat(1,5,i) =  bm(1,5,i)*hh
            bmat(2,5,i) =  bm(2,5,i)*hh
            bmat(4,5,i) =  bm(3,5,i)*hh

            if(qdflg) then
              bmat(1,6,i) = bm(1,6,i)*hh - shp1(1,i,l)
              bmat(2,6,i) = bm(2,6,i)*hh - shp2(2,i,l)
              bmat(4,6,i) = bm(3,6,i)*hh - shp1(2,i,l) - shp2(1,i,l)
            else
              bmat(1,6,i) = 0.0d0
              bmat(2,6,i) = 0.0d0
              bmat(4,6,i) = 0.0d0
            endif

c           Compute stress divergence terms

            do ii = 1,6
              rc(ii,i) = rc(ii,i) - bmat(1,ii,i)*sig(1)
     &                            - bmat(2,ii,i)*sig(2)
     &                            - bmat(4,ii,i)*sig(4)
            end do ! ii
          end do ! i

c         Form stiffness

          i1 = 0
          do i = 1,nel

c           Form stress-displacement matrix (Bi-trans * D)

            do ii = 1,6
              dbmat(ii,1) = bmat(1,ii,i)*dd(1,1,1)
     &                    + bmat(2,ii,i)*dd(2,1,1)
     &                    + bmat(4,ii,i)*dd(4,1,1)

              dbmat(ii,2) = bmat(1,ii,i)*dd(1,2,1)
     &                    + bmat(2,ii,i)*dd(2,2,1)
     &                    + bmat(4,ii,i)*dd(4,2,1)

              dbmat(ii,4) = bmat(1,ii,i)*dd(1,4,1)
     &                    + bmat(2,ii,i)*dd(2,4,1)
     &                    + bmat(4,ii,i)*dd(4,4,1)

            end do ! ii

c           Loop on columns

            j1 = 0
            do j = 1,nel

c             Compute stiffness contribution

              do jj = 1,6
                do ii = 1,6
                  sc(ii+i1,jj+j1) = sc(ii+i1,jj+j1)
     +                            + dbmat(ii,1)*bmat(1,jj,j)
     +                            + dbmat(ii,2)*bmat(2,jj,j)
     +                            + dbmat(ii,4)*bmat(4,jj,j)
                end do ! ii
              end do ! jj
              j1 = j1 + 6
            end do ! j
            i1 = i1 + 6
          end do ! i

        end do ! ll

c       Inertial parts

        do i = 1,nel
          dvm = dvl(l)*d(4)*d(14)*shp(3,i,l)
          do ii = 1,3
            rr(ii,i) = rr(ii,i) - dvm*( wl(ii,i,3) + d(77)*wl(ii,i,2))
            ss(ii,i) = ss(ii,i) + dvm*ctan3
          end do ! ii
          dvm = max(0.0d0,d(8))*dvm*d(14)**2/12.d0
          do ii = 4,6
            rr(ii,i) = rr(ii,i) - dvm*( wl(ii,i,3) + d(77)*wl(ii,i,2))
            ss(ii,i) = ss(ii,i) + dvm*ctan3
          end do ! ii
        end do ! i

c       Compute Hughes/Brezzi rotation matrix

        if(ialph.ne.0 .and. qdflg) then
          do j = 1,4
            gg(1,j) = gg(1,j) - shp(2,j,l)*dvl(l)
            gg(2,j) = gg(2,j) + shp(1,j,l)*dvl(l)
            gg(6,j) = gg(6,j) - 2.d0*shp(3,j,l)*dvl(l)
            if(ialph.gt.1) then
              gg(6,j) = gg(6,j) + (shp1(2,j,l) - shp2(1,j,l))*dvl(l)
            endif
          end do ! j
        endif
      end do ! l

c     Perform rank one update for Huges/Brezzi term

      if(ialph.gt.0 .and. qdflg) then

c       Compute H/B strain

        xx = 0.0d0
        do i = 1,4
          do j = 1,6
            xx = xx + gg(j,i)*vl(j,i,1)
          end do ! j
        end do ! i

c       Compute H/B residual and tangent

        ggv = pen*d(27)*thk*ctan1/dv

        i1 = 0
        do i = 1,4
          do ii = 1,6
            ggi      = ggv*gg(ii,i)
            rc(ii,i) = rc(ii,i) - ggi*xx
            j1 = 0
            do j = 1,4
              do jj = 1,6
                sc(ii+i1,jj+j1) = sc(ii+i1,jj+j1) + ggi*gg(jj,j)
              end do ! jj
              j1 = j1 + 6
            end do ! j
          end do ! ii
          i1 = i1 + 6
        end do ! i
      endif

c     Correct residual and tangent matrix for element warpage

      if(yl(3,1).ne.0.0d0) then
        i1 = 1
        do i = 1,4
          rc(4,i) = rc(4,i) + yl(3,i)*rc(2,i)
          rc(5,i) = rc(5,i) - yl(3,i)*rc(1,i)
          j1 = 1
          do j = 1,4
            call proj3d(sc(i1,j1),yl(3,i),yl(3,j),24)
            j1 = j1 + 6
          end do ! j
          i1 = i1 + 6
        end do ! i
      endif

c     Move to final locations

      i1 = 0
      i2 = 0
      do i = 1,nel
        do j = 1,6
          r(j,i) = rc(j,i)
        end do ! j
        j1 = 0
        j2 = 0
        do j = 1,nel
          do ii = 1,6
            do jj = 1,6
              s(j1+jj,i1+ii) = sc(j2+jj,i2+ii)
            end do ! jj
          end do ! ii
          j1 = j1 + ndf
          j2 = j2 + 6
        end do ! j
        i1 = i1 + ndf
        i2 = i2 + 6
      end do ! i

c     Rotate to global frame

      call rots3d(s,r,tr,nst,ndf,nel)

c     Rotational body force terms

      if(d(4).gt.0.0d0 .and. d(65).gt.0.0d0) then

        do l = 1,lint ! {

c         Angular velocity: d(4) = rho; d(65) = omega

          do i = 1,3
            bf(i) = 0.0d0
            x3(i) = 0.0d0
            do j = 1,4
              x3(i) = x3(i) + shp(3,j,l)*xl(i,j)
            end do ! j
          end do ! i
          call sbodyw(d(4)*thk,d(65),x3, bf,bt, .false.)

          do j = 1,4
            shp3i = shp(3,j,l)*dvl(l)
            do i = 1,3
              r(i,j) = r(i,j) + shp3i*bf(i)
            end do ! i
          end do ! j
        end do ! l
      endif

c     Add the inertia parts

      i1 = 0
      do i = 1,nel
        do j = 1,6
          r(j   ,i   ) = r(j   ,i   ) + rr(j,i)
          s(j+i1,j+i1) = s(j+i1,j+i1) + ss(j,i)
        end do ! j
        i1 = i1 + ndf
      end do ! i

      return

c     Compute and output element variables

4     nn = 0
      do l = 1,lint

c       Form shape functions

        if(qdflg) then
          call rshp3d(sg(1,l),yl,shp,shp1,shp2,xsj,3)

c         Modify rotational shape functions

          do j = 1,4
            do i = 1,3
              shp1(i,j,1) = shp1(i,j,1) - gshp1(i,j)/dv
              shp2(i,j,1) = shp2(i,j,1) - gshp2(i,j)/dv
            end do ! i
          end do ! j
        else
          call shp2d (sg(1,l),yl,shp,xsj,3,nel,ix,.false.)
        endif

        do ll = 1,lt
          hh  = 0.5d0*thk*sgt(1,ll)
          call str3di(d,xl,vl,ndm,nel,1, xx,yy,zz,hh,eps,sig,dd,nn,isw)
          nn  = nn + nint(d(15))
          mct = mct - 3
          if(mct.le.0) then
            write(iow,2002) o,head
            if(ior.lt.0) then
            write(*,2002) o,head
            endif
            mct = 50
          end if
          xx = xx + tr(1,3)*hh
          yy = yy + tr(2,3)*hh
          zz = zz + tr(3,3)*hh
          write(iow,2003) n,xx,(sig(i),i=1,6),ma,yy,eps,zz
          if(ior.lt.0) then
            write(*,2003) n,xx,(sig(i),i=1,6),ma,yy,eps,zz
          endif
        end do ! ll
      end do ! l
      return

c     Compute element mass or geometric stifness arrays

5     if(imtyp.eq.1) then

c       Compute mass

        do l = 1,lint
          dv1 = dvl(l)*d(4)*d(14)
          dv2 = dv1*d(14)**2/12.d0*max(0.0d0,d(8))
          i1  = 0
          do j = 1,nel
            do i = 1,3
              rc(i  ,j)         = rc(i  ,j) + shp(3,j,l)*dv1
              rc(i+3,j)         = rc(i+3,j) + shp(3,j,l)*dv2
              sc(i1+i  ,i1+i  ) = rc(i  ,j)
              sc(i1+i+3,i1+i+3) = rc(i+3,j)
            end do ! i
            i1 = i1 + 6
          end do ! j
        end do ! l

c       Move to tangent locations

        i1 = 0
        j1 = 0
        do i = 1,nel
          do ii = 1,ndf
            r(ii,i) = rc(ii,i)
            s(i1+ii,i1+ii) = sc(j1+ii,j1+ii)
          end do ! ii
          i1 = i1 + ndf
          j1 = j1 + 6
        end do !i

      elseif(imtyp.eq.2) then

c       Compute geometric stiffness

        do i = 1,nel
          do j = 1,3
            vl(j  ,i,1) = 0.0d0
            vl(j+3,i,1) = 0.0d0
            do k = 1,3
              vl(j  ,i,1) = vl(j  ,i,1) + tr(j,k)*wl(k  ,i,1)
              vl(j+3,i,1) = vl(j+3,i,1) + tr(j,k)*wl(k+3,i,1)
            end do ! k
          end do ! j
        end do ! i

        nn = 0
        do l = 1,lint

c         Modify rotational shape functions

          if(qdflg) then
            do j = 1,4
              do i = 1,3
                shp1(i,j,l) = shp1(i,j,l) - gshp1(i,j)/dv
                shp2(i,j,l) = shp2(i,j,l) - gshp2(i,j)/dv
              end do ! i
            end do ! j
          endif

          hh = 0.0d0
          call str3di(d,xl,vl,ndm,nel,l, xx,yy,zz,hh,eps,sig,dd,nn,isw)
          nn = nn + nint(d(15))

          do i = 1,4
            sig(i) = sig(i)*thk
          end do ! i

          i1 = 0
          do i = 1,nel !{
            j1 = 0
            dv1 = (shp(1,i,l)*sig(1) + shp(2,i,l)*sig(4))*dvl(l)
            dv2 = (shp(1,i,l)*sig(4) + shp(2,i,l)*sig(2))*dvl(l)
            do j = 1,nel !{
              a11 = dv1*shp(1,j,l) + dv2*shp(2,j,l)
              do k = 1,3 !{
                sc(i1+k,j1+k) = sc(i1+k,j1+k) - a11
              end do ! k
              j1 = j1 + ndf
            end do ! j
            i1 = i1 + ndf
          end do ! i

        end do ! l

        if(yl(3,1).ne.0.0d0) then
          i1 = 1
          do i = 1,4 !{
            j1 = 1
            do j = 1,4 !{
              call proj3d(sc(i1,j1),yl(3,i),yl(3,j),24)
              j1 = j1 + 6
            end do ! j
            i1 = i1 + 6
          end do ! i
        endif

c       Move to tangent and residual locations

        i1 = 0
        i2 = 0
        do i = 1,nel

c         Residual

          do j = 1,6
            r(j,i) = rc(j,i)
          end do ! j

c         Tangent

          j1 = 0
          j2 = 0
          do j = 1,nel
            do ii = 1,6
              do jj = 1,6
                s(j1+jj,i1+ii) = sc(j2+jj,i2+ii)
              end do ! jj
            end do ! ii
            j1 = j1 + ndf
            j2 = j2 + 6
          end do ! j
          i1 = i1 + ndf
          i2 = i2 + 6
        end do ! i

        call rots3d(s,r,tr,nst,ndf,nel)

      endif

      return

c     Compute nodal output quantities

8     call stcn3si(d,yl,ul,tr,r,s,dvl,thk,ndf,nel,isw)
      return

c     Formats

2002  format(a1,20a4//5x,'S h e l l   S t r e s s e s'//
     & ' elmt x-coord  xx-stress  yy-stress  zz-stress  xy-stress',
     & ' matl y-coord  xx-strain  yy-strain  zz-strain  xy-strain',
     & '      z-coord'/38(' -'))

2003  format(/i5,0p,1f8.3,1p,5e11.3,0p,1f8.2/i5,0p,1f8.3,1p,5e11.3,
     &       0p,1f8.2/5x,0p,1f8.3,1p,6e11.3/(13x,1p,6e11.3))

      end

      subroutine stcn3si(d,yl,ul,tr,dt,st,dvl,thk,ndf,nel,isw)

      implicit  none

      include  'cdata.h'
      include  'strnum.h'
      include  'shpf16.h'
      include  'shlc16.h'
      include  'sstr16.h'

      integer   ndf,nel, i,j,l, ll,lt, nn,isw
      real*8    xsji, xx,yy,zz, hh, dh, thk

      real*8    dt(*),st(nen,*),yl(3,*),tr(3,3),d(*),ul(ndf,*)
      real*8    vl(6,4),dvl(4)
      real*8    eps(6),sig(10),momt(6),norm(6),dd(6,6,5),eps0(6)

      save

c     Compute membrane and bending stresses for projection.

      do i = 1,nel
        do j = 1,3
          vl(j  ,i) = 0.0d0
          vl(j+3,i) = 0.0d0
          do l = 1,3
            vl(j  ,i) = vl(j  ,i) + tr(j,l)*ul(l  ,i)
            vl(j+3,i) = vl(j+3,i) + tr(j,l)*ul(l+3,i)
          end do ! l
        end do ! j
      end do ! i

      lt = nint(d(102))
      nn = 0
      do l = 1,4
        do j = 1,6
          norm(j) = 0.0d0
          momt(j) = 0.0d0
        end do ! j
        do ll = 1,lt
          hh = 0.5d0*thk*sgt(1,ll)
          dh = 0.5d0*thk*sgt(2,ll)
          call str3di(d,yl,vl,3,nel,l, xx,yy,zz,hh,eps,sig,dd,nn,isw)
          norm(1) = norm(1) + sig(1)*dh
          norm(2) = norm(2) + sig(2)*dh
          norm(3) = norm(3) + sig(4)*dh
          momt(1) = momt(1) + sig(1)*hh*dh
          momt(2) = momt(2) + sig(2)*hh*dh
          momt(3) = momt(3) + sig(4)*hh*dh

          if(ll.eq.1) then
            do j = 1,nel
              xsji     = dvl(l)*shp(3,j,l)
              st(j,20) = st(j,20) + sig(1)*xsji
              st(j,21) = st(j,21) + sig(2)*xsji
              st(j,22) = st(j,22) + sig(4)*xsji
            end do ! j
            do j = 1,3
              eps0(j) = eps(j)
            end do ! j

          elseif(ll.eq.lt) then

            do j = 1,nel
              xsji     = dvl(l)*shp(3,j,l)
              st(j,17) = st(j,17) + sig(1)*xsji
              st(j,18) = st(j,18) + sig(2)*xsji
              st(j,19) = st(j,19) + sig(4)*xsji
            end do ! j
            do j = 1,3
              eps0(j+3) = (eps(j) - eps0(j))/thk
              eps0(j  ) = (eps(j) + eps0(j))*0.5d0
            end do ! j

          endif

          nn = nn + nint(d(15))
        end do ! ll
        do j = 1,nel
          xsji     = dvl(l)*shp(3,j,l)
          dt(j)    = dt(j)    + xsji

          st(j,1)  = st(j,1)  + norm(1)*xsji
          st(j,2)  = st(j,2)  + norm(2)*xsji
          st(j,3)  = st(j,3)  + norm(3)*xsji
          st(j,4)  = st(j,4)  + norm(4)*xsji
          st(j,5)  = st(j,5)  + norm(5)*xsji

          st(j,6)  = st(j,6)  + momt(1)*xsji
          st(j,7)  = st(j,7)  + momt(2)*xsji
          st(j,8)  = st(j,8)  + momt(3)*xsji
          st(j,9)  = st(j,9)  + momt(4)*xsji
          st(j,10) = st(j,10) + momt(5)*xsji

          st(j,11) = st(j,11) + eps0(1)*xsji
          st(j,12) = st(j,12) + eps0(2)*xsji
          st(j,13) = st(j,13) + eps0(3)*xsji
          st(j,14) = st(j,14) + eps0(4)*xsji
          st(j,15) = st(j,15) + eps0(5)*xsji
          st(j,16) = st(j,16) + eps0(6)*xsji

        end do ! j
      end do ! l

      iste = 22

      end

      subroutine str3di(d,xl,vl,ndm,nel,l, xx,yy,zz,hh, eps,sig,dd,
     &                  nn,isw)

      implicit  none

      include  'elcoor.h'
      include  'hdata.h'
      include  'shlc16.h'
      include  'shld16.h'
      include  'shpf16.h'
      include  'sstr16.h'

      include  'comblk.h'

      integer   ndm,nel,l, i,j, nn,isw, nhv, istrt

      real*8    alam,ha, xx,yy,zz, hh, ta

      real*8    d(*), xl(ndm,*),vl(6,*), eps0(6),eps(6),sig(10)
      real*8    dd(6,6,5), uu(3)

      save

      data      alam,ha / 2*0.0d0 /

c     Compute membrane and bending strains

      if(qdflg) then
        call dktq3d(shp(1,5,l),shp(1,1,l))
      else
        call dktb3d(sg(1,l),xsjt)
      endif

      do i = 1,6
        eps0(i) = 0.0
        eps (i) = 0.0
      end do ! i
      xx = 0.0d0
      yy = 0.0d0
      zz = 0.0d0
      do i = 1,3
        uu(i) = 0.0d0
      end do ! i
      do j = 1,nel
        xx = xx + shp(3,j,l)*xl(1,j)
        yy = yy + shp(3,j,l)*xl(2,j)
        zz = zz + shp(3,j,l)*xl(3,j)
        eps0(1) = eps0(1) + shp(1,j,l)*vl(1,j)
     &                    - shp1(1,j,l)*vl(6,j)
        eps0(2) = eps0(2) + shp(2,j,l)*vl(2,j)
     &                    - shp2(2,j,l)*vl(6,j)
        eps0(3) = eps0(3) + shp(1,j,l)*vl(2,j) + shp(2,j,l)*vl(1,j)
     &                    - (shp1(2,j,l) + shp2(1,j,l))*vl(6,j)
        do i = ii1,ii2
          eps0(4) = eps0(4) + bm(1,i,j)*vl(i,j)
          eps0(5) = eps0(5) + bm(2,i,j)*vl(i,j)
          eps0(6) = eps0(6) + bm(3,i,j)*vl(i,j)
        end do ! i
        do i = 1,3
          uu(i) = uu(i) + shp(3,j,l)*vl(i,j)
        end do ! i
      end do ! j

c     Compute strain at current elevation

      eps(1) = eps0(1) + hh*eps0(4)
      eps(2) = eps0(2) + hh*eps0(5)
      eps(4) = eps0(3) + hh*eps0(6)

c     Save coordinates

      xref(1) = xx
      xref(2) = yy
      xref(3) = zz
      do i = 1,3
        xcur(i) = xref(i) + uu(i)
      end do ! i

      ta     = 0.0d0

c     Compute stress at point

      nhv   = nint(d(15))
      istrt = nint(d(84))
      call modlsd(l,d,ta,eps,hr(nh1+nn),hr(nh2+nn),nhv,istrt,dd,sig,
     &            alam,ha,isw)

      end
