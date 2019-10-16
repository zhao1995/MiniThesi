c$Id:$
      subroutine mem3ds(d,ul,xl,ix,s,p,ndf,ndm,nst,isw)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Revise set of quadrature on d(5)                 27/03/2009
c       2. Add principal values to nodal projections        15/04/2009
c       3. Compute xref and xcur in 'elcoor.h'              26/03/2009
c       4. Dimension dd(6,6,5) in stcn3mb                   18/06/2011
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose:  Quadrilateral membrane element for feap

c         12 January  1998 -  Basic developments
c         05 February 2003 -  Add inertial properties/mass

c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c.... Input parameters set as follows:

c         ndm = 3 (x,y,z cartesian coordinates at nodes)
c         ndf = 3 (u-x,u-y,u-z at nodes)
c         nen = 4 nodes (counterclockwise around element)

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
      integer   i,j,k,l,lint, i1,j1, i2,j2

      real*8    a11,a12, a21,a22, a31,a32, xx,yy,zz
      real*8    dv, dv1,dv2, xsj, thk
      real*8    xn,yn, shp1i,shp2i,shp3i

      integer   ix(*)

      real*8    d(*),xl(ndm,*),ul(ndf,nen,*),s(nst,*),p(nst)
      real*8    norm(6),dvl(9),eps(9,3)
      real*8    yl(3,9),vl(6,9),tr(3,3),bl(3),bf(3),bg(3),bt(3,3)
      real*8    pm(27),sm(27,27),el(3),x3(3), dd(6,6,5)

      save

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

      do i = 1,nel ! {
        do j = 1,3 ! {
          vl(j  ,i) = 0.0d0
          do k = 1,3
            vl(j  ,i) = vl(j  ,i) + tr(j,k)*ul(k,i,1)
          end do ! k  }
        end do ! j  }
      end do ! i  }

c     Test for triangular element

      if(nel.eq.3) then
        qdflg = .false.
      elseif( ix(1) .eq. ix(2)  .or. ix(2) .eq. ix(3) .or.
     &        ix(3) .eq. ix(4)  .or. ix(4) .eq. ix(1) ) then
        qdflg = .false.
      else
        qdflg = .true.
      endif

      if(qdflg)then
        call jacq3d(yl)
        l = min(3,nint(d(5)))
        if(l.eq.0) then
          l = 2
        endif
        call int2d(l ,lint,sg)
      else
        call jtri3d(yl,xsjt)
        el(1)   = 1.d0/3.d0
        el(2)   = el(1)
        el(3)   = el(1)
        lint    = 1
        sg(3,1) = 1.d0
      endif

      dv = 0.0d0
      do l = 1,lint ! {

c       Form shape functions and their integrals

        if(qdflg) then
          call shp2d(sg(1,l),yl,shp(1,1,l),xsj,3,nel,ix,.false.)
        else
          call trishp(el,yl,3,1, xsj,shp(1,1,l) )
        endif
        dvl(l) = xsj*sg(3,l)
        dv     = dv + dvl(l)
      end do ! l  }

c     Compute thickness for element

      thk  = d(14)

      if(isw.eq.4) go to 4
      if(isw.eq.5) go to 5
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

c     Transfer body loads to local frame

      do i = 1,3 ! {
        do j = 1,3 ! {
          bl(i) = bl(i) + tr(i,j)*bg(j)
        end do ! j  }
      end do ! i  }

c     Compute membrane/load parts

      do l = 1,lint ! {

        dv1 = thk *dvl(l)*ctan(1)

        call stre3m(d,yl,vl,3,nel,l, xx,yy,zz,eps,norm,dd,isw)

c       Store time history plot data for element

        i = 3*(l-1)
        do j = 1,3
          tt(j+i  ) = norm(j)
        end do ! j

c       Recover previously computed shape functions

        i1 = 1
        do i = 1,nel ! {
          shp1i = shp(1,i,l)
          shp2i = shp(2,i,l)
          shp3i = shp(3,i,l)

c         Compute loading term

          p(i1  ) = p(i1  ) + shp3i*bl(1)*dvl(l)
          p(i1+1) = p(i1+1) + shp3i*bl(2)*dvl(l)
          p(i1+2) = p(i1+2) + shp3i*bl(3)*dvl(l)

c         Compute stress divergence terms

          p(i1  ) = p(i1  ) - (shp1i*norm(1) + shp2i*norm(3))*dvl(l)

          p(i1+1) = p(i1+1) - (shp2i*norm(2) + shp1i*norm(3))*dvl(l)

c         Form stress-displacement matrix (Bi-trans * D)

          a11 = (dd(1,1,1)*shp1i + dd(1,4,1)*shp2i)*dv1
          a12 = (dd(1,2,1)*shp2i + dd(1,4,1)*shp1i)*dv1
          a21 = (dd(2,1,1)*shp1i + dd(2,4,1)*shp2i)*dv1
          a22 = (dd(2,2,1)*shp2i + dd(2,4,1)*shp1i)*dv1
          a31 = (dd(4,1,1)*shp1i + dd(4,4,1)*shp2i)*dv1
          a32 = (dd(4,2,1)*shp2i + dd(4,4,1)*shp1i)*dv1

c         Loop on columns

          j1 = i1
          do j = i,nel ! {
            xn = shp(1,j,l)
            yn = shp(2,j,l)

c           Compute membrane part

            s(i1  ,j1  ) = s(i1  ,j1  ) + (a11*xn + a31*yn)
            s(i1+1,j1  ) = s(i1+1,j1  ) + (a12*xn + a32*yn)
            s(i1  ,j1+1) = s(i1  ,j1+1) + (a21*yn + a31*xn)
            s(i1+1,j1+1) = s(i1+1,j1+1) + (a22*yn + a32*xn)

            j1 = j1 + ndf
          end do ! j  }
          i1 = i1 + ndf
        end do ! i }
      end do ! l }

c     Rotate to global frame

      call rotm3d(s,p,tr,nel,nst,ndf)

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

          j1 = 0
          do j = 1,4
            shp3i = shp(3,j,l)*dvl(l)
            do i = 1,3
              p(j1+i) = p(j1+i) + shp3i*bf(i)
            end do ! i
            j1 = j1 + ndf
          end do ! j
        end do ! l
      endif

c     Inertial terms

      if(ndfo(1).gt.0 .or. shflg) then
        do i = 1,3*nel ! {
          pm(i) = 0.0d0
          do j = 1,3*nel ! {
            sm(j,i) = 0.0d0
          end do ! j }
        end do ! i }
        call me3smas(d,dvl,pm,sm,3,12,lint)
        i1 = 0
        i2 = 0
        do i = 1,nel ! {
          j1 = 0
          j2 = 0
          do j = 1,nel ! {
            do k = 1,3 ! {
              p(i1+k)      = p(i1+k) - sm(i2+k,j2+k)*ul(k,j,4)
              s(i1+k,j1+k) = s(i1+k,j1+k) + sm(i2+k,j2+k)*ctan(3)
            end do ! k }
            j1 = j1 + ndf
            j2 = j2 + 3
          end do ! j }
          i1 = i1 + ndf
          i2 = i2 + 3
        end do ! i }
      endif

      return

c     Compute and output element variables

4     do l = 1,lint ! {

        call stre3m(d,xl,vl,ndm,nel,l, xx,yy,zz,eps,norm,dd,isw)
        mct = mct - 3
        if(mct.le.0) then
          write(iow,2002) o,head
          if(ior.lt.0) then
            write(*,2002) o,head
          endif
          mct = 50
        end if
        write(iow,2003) n,xx,norm,ma,yy,sigt,zz,(eps(i,1),i=1,3)
        if(ior.lt.0) then
          write(*,2003) n,xx,norm,ma,yy,sigt,zz,(eps(i,1),i=1,3)
        endif
      end do ! l  }
      return

c     Compute element mass or geometric stifness arrays

5     if(imtyp.eq.1) then

c       Compute mass

        call me3smas(d,dvl,p,s,ndf,nst,lint)

      elseif(imtyp.eq.2) then

c       Compute geometric stiffness

        do i = 1,nel ! {
          do j = 1,3 ! {
            vl(j  ,i) = 0.0d0
            do k = 1,3 ! {
              vl(j  ,i) = vl(j  ,i) + tr(j,k)*ul(k,i,1)
            end do ! k  }
          end do ! j  }
        end do ! i  }

        do l = 1,lint ! {

          call stre3m(d,xl,vl,ndm,nel,l, xx,yy,zz,eps,norm,dd,isw)

          i1 = 0
          do i = 1,nel !{
            j1 = 0
            dv1 = (shp(1,i,l)*norm(1) + shp(2,i,l)*norm(3))*dvl(l)
            dv2 = (shp(1,i,l)*norm(3) + shp(2,i,l)*norm(2))*dvl(l)
            do j = 1,nel !{
              a11 = dv1*shp(1,j,l) + dv2*shp(2,j,l)
              do k = 1,3 !{
                s(i1+k,j1+k) = s(i1+k,j1+k) - a11
              end do ! k  }
              j1 = j1 + ndf
            end do ! j  }
            i1 = i1 + ndf
          end do ! i  }
        end do ! l  }

        call rotm3d(s,p,tr,nel,nst,ndf)
      endif

      return

c     Compute nodal output quantities

8     call stcn3mb(d,yl,ul,tr,p,s,dvl,ndf,nel,lint,isw)
      return

c     Formats

2002  format(a1,20a4//5x,'M e m b r a n e   S t r e s s e s'//
     & 4x,'Elmt x-Coord  11-Force   22-Force   12-Force    1-Force ',
     & 3x,'2-Force    Angle'/
     & 4x,'Matl y-Coord  11-Stress  22-Stress  12-Stress   1-Stress',
     & 3x,'2-Stress   Angle'/
     & 9x,'z-Coord  11-Strain  22-Strain  12-Strain')

2003  format(/i8,0p,1f8.3,1p,5e11.3,0p,1f8.2
     &       /i8,0p,1f8.3,1p,5e11.3,0p,1f8.2
     &       /8x,0p,1f8.3,1p,3e11.3)

      end

      subroutine rotm3d(s,p,t,nel,nst,ndf)

c     Transform loads and stiffness to global coords.

      implicit  none

      integer   nel,nst,ndf, i,i0, ir,ii, j,j0, jc,jj

      real*8    s(nst,nst),p(nst),t(3,3),a(3,3),b(3)

      real*8    dot

      i0 = 0
      do ir = 1,nel ! {
        do ii = 1,3 ! {
          b(ii  ) = dot(t(1,ii),p(i0+1),3)
        end do ! ii }
        do ii = 1,3 ! {
          p(i0+ii) = b(ii)
        end do ! ii }
        j0 = i0
        do jc = ir,nel ! {
          do ii = 1,3 ! {
            do jj = 1,3 ! {
              a(jj,ii) = dot(t(1,ii),s(i0+1,jj+j0),3)
            end do ! jj }
          end do ! ii }
          do jj = 1,3 ! {
            do ii = 1,3 ! {
              s(ii+i0,jj+j0) = dot(a(1,ii),t(1,jj),3)
            end do ! ii }
          end do ! jj }

c         Compute symmetric block

          if(ir.ne.jc) then
            do i = 1,3 ! {
              do j = 1,3 ! {
                s(j0+j,i0+i) = s(i0+i,j0+j)
              end do ! j  }
            end do ! i  }
          endif
          j0 = j0 + ndf
        end do ! jc  }
        i0 = i0 + ndf
      end do ! ir  }

      end

      subroutine stcn3mb(d,yl,ul,tr,dt,st,dvl,ndf,nel,lint,isw)

      implicit  none

      include  'cdata.h'
      include  'strnum.h'
      include  'shpf16.h'
      include  'sstr16.h'

      integer   ndf,nel,lint,isw, i,j,k
      real*8    xsji, xx,yy,zz

      real*8    dt(*),st(nen,*),yl(3,*),tr(3,3),d(*),ul(ndf,*)
      real*8    vl(6,4),eps(6),norm(6),dd(6,6,5),dvl(4)
      real*8    ps(3),pe(3),pn(3)

      save

c     Compute membrane and bending stresses for projection.

      do i = 1,nel ! {
        do j = 1,3 ! {
          vl(j  ,i) = 0.0d0
          do k = 1,3 ! {
            vl(j  ,i) = vl(j  ,i) + tr(j,k)*ul(k,i)
          end do ! k  }
        end do ! j  }
      end do ! i  }

      do k = 1,lint ! {
        call stre3m(d,yl,vl,3,nel,k, xx,yy,zz,eps,norm,dd,isw)
        do j = 1,nel ! {
          xsji = dvl(k)*shp(3,j,k)
          dt(j)    = dt(j)    + xsji

          st(j, 1) = st(j, 1) + norm(1)*xsji
          st(j, 2) = st(j, 2) + norm(2)*xsji
          st(j, 4) = st(j, 4) + norm(3)*xsji

          call pstr2d(norm,pn)
          st(j, 5) = st(j, 5) + pn(1)*xsji
          st(j, 6) = st(j, 6) + pn(2)*xsji

          st(j, 7) = st(j, 7) + eps(1)*xsji
          st(j, 8) = st(j, 8) + eps(2)*xsji
          st(j,10) = st(j,10) + eps(3)*xsji

          call pstr2d(eps,pe)
          st(j,11) = st(j,11) + pe(1)*xsji
          st(j,12) = st(j,12) + pe(2)*xsji

          st(j,13) = st(j,13) + sigt(1)*xsji
          st(j,14) = st(j,14) + sigt(2)*xsji
          st(j,15) = st(j,15) + sigt(3)*xsji

          call pstr2d(sigt,ps)
          st(j,16) = st(j,16) + ps(1)*xsji
          st(j,17) = st(j,17) + ps(2)*xsji

        end do ! j  }
      end do ! k  }

      iste = 17
      end

      subroutine stre3m(d,xl,vl,ndm,nel,l, xx,yy,zz,eps,norm,dd,isw)

      implicit  none

      include  'elcoor.h'
      include  'hdata.h'
      include  'shlc16.h'
      include  'shld16.h'
      include  'comblk.h'
      include  'shpf16.h'
      include  'sstr16.h'

      integer   ndm,nel,l,isw, i,j, nhv,istrt

      real*8    alam,ha, xx,yy,zz,thk,ta, uu(3)

      real*8    d(*), xl(ndm,*),vl(6,*), eps(6),norm(6),temp(4)
      real*8    dd(6,6,5),alp(6)

      save

      data      alam,ha / 2*0.0d0 /

c     Compute membrane strains

      do i = 1,6 ! {
        eps(i) = 0.0
      end do ! i  }
      xx = 0.0d0
      yy = 0.0d0
      zz = 0.0d0
      do j = 1,3
        uu(j) = 0.0d0
      end do ! j
      do i = 1,nel ! {
        xx = xx + shp(3,i,l)*xl(1,i)
        yy = yy + shp(3,i,l)*xl(2,i)
        zz = zz + shp(3,i,l)*xl(3,i)
        eps(1) = eps(1) + shp(1,i,l)*vl(1,i)
        eps(2) = eps(2) + shp(2,i,l)*vl(2,i)
        eps(3) = eps(3) + shp(1,i,l)*vl(2,i) + shp(2,i,l)*vl(1,i)
        do j = 1,3
          uu(j) = uu(j) + shp(3,i,l)*vl(j,i)
        end do ! j
      end do ! i  }

c     Save coordinates

      xref(1) = xx
      xref(2) = yy
      xref(3) = zz
      do j = 1,3
        xcur(j) = xref(j) + uu(j)
      end do ! j

c     Elastic behavior

      if(nint(d(40)).eq.0) then
        call dmat2d(d,d(31),dd,alp)
      else
        ta    = 0.0d0
        nhv   = nint(d(15))
        istrt = nint(d(84))
        call modlsd(l,d,ta,eps,hr(nh1),hr(nh2),nhv,istrt,dd,sigt,
     &              alam,ha,isw)
      endif

c     Compute surface stresses

      sigt(1) = dd(1,1,1)*eps(1) + dd(1,2,1)*eps(2) + dd(1,4,1)*eps(3)
      sigt(2) = dd(2,1,1)*eps(1) + dd(2,2,1)*eps(2) + dd(2,4,1)*eps(3)
      sigt(3) = dd(4,1,1)*eps(1) + dd(4,2,1)*eps(2) + dd(4,4,1)*eps(3)
      temp(1) = sigt(1)
      temp(2) = sigt(2)
      temp(3) = 0.0d0
      temp(4) = sigt(3)

      call pstr2d(temp,sigt(4))

c     Compute in-plane loading

      thk     = d(14)
      do i = 1,5 ! {
        norm(i) = sigt(i)*thk
      end do ! i   }
      norm(6) = sigt(6)

      end

      subroutine me3smas(d,dvl,p,s,ndf,nst,lint)

      implicit   none

      include   'eldata.h'
      include   'shpf16.h'

      integer    ndf,nst,lint
      real*8     d(*),dvl(*),p(*),s(nst,nst)

      integer    i,j,k,l, i1,j1
      real*8     dv1,dv2,a11,lmas,cmas

c     Compute mass

      do l = 1,lint ! {
        dv1 = dvl(l)*d(4)*d(14)
        j1  = 0
        do j = 1,nel ! {
          dv2 = shp(3,j,l)*dv1
          do k = 1,3 ! {
            p(j1+k)      = p(j1+k) + dv2
          end do ! i }
          i1 = 0
          do i = 1,nel ! {
            a11 = shp(3,i,l)*dv2
            do k = 1,3 ! {
              s(i1+k,j1+k) = s(i1+k,j1+k) + a11
            end do ! i }
            i1 = i1 + ndf
          end do ! i }
          j1 = j1 + ndf
        end do ! j }
      end do ! l }

c     Interpolate the consistent mass

      if(d(7).eq.0.0d0) then
        lmas = d(7)
        cmas = 1.d0 - d(7)
        do j = 1,nst ! {
          do i = 1,nst ! {
            s(i,j) = cmas*s(i,j)
          end do ! i }
          s(j,j) = s(j,j) + lmas*s(j,j)
        end do ! j }
      endif

      end
