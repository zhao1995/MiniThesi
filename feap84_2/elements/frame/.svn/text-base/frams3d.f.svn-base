c$Id:$
      subroutine frams3d(d,ul,xl,s,r,ndf,ndm,nst,isw)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Use constants from 'constant.h'                  14/11/2006
c       2. Remove mass3s to separate function               06/09/2007
c       3. Remove mprint from debug print                   04/10/2008
c       4. Increase quadrature for 'rect' to 10             15/11/2008
c       5. Subroutine framtr: correct division by dl in t_2 10/01/2010
c       6. Replace d(21) by d(1) to get E                   11/01/2010
c-----[--.----+----.----+----.-----------------------------------------]
c     Small Deformation Three dimensional frame element

c     Control data:
c         ndm - 3 (x,y,z)
c         ndf - 6 (u,v,w, theta_x,theta_y,theta_z)
c         nen - 3 or more (see below)

c      Beam end nodes 1 and 2
c      Plane defined by nodes 1, 2, 3 contains z-axis
c                            (perpendicular to x-axis)

c      Vector products: e_1 =  (x_2 - x_1)/|x_2 - x_1|
c                       v_2 = - e_1  x ( x_3 - x_1)
c                       e_2 =   v_2/|v_2|
c                       e_3 =   e_1 x e_2

c                       z (e_3)  x (e_1)
c                     3 o- - - /
c                       |     o 2
c                       |    /
c                       |   / <--- Frame axis
c                       |  /
c                       | /
c                       |/
c     (e_2) y ----------o 1

c     Displacement:     u_x = u_0 + z * theta_y - y * theta_z

c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'bdata.h'
      include  'cdata.h'
      include  'eldata.h'
      include  'eltran.h'
      include  'evdata.h'
      include  'hdata.h'
      include  'iofile.h'
      include  'part0.h'
      include  'prstrs.h'
      include  'rdata.h'
      include  'comblk.h'

      integer   ndf,ndm,nst,isw, i,ii,i1,i2,j,jj,j1,j2,k,l,lint,nh,nn
      real*8    le,dl,dl2,dl3,ctan1,ctan3
      real*8    d(*),xl(ndm,*),ul(ndf,nen,*),s(nst,*),r(ndf,*)
      real*8    t(3,3),fi(6,2),sm(12,12),pm(12),eps(6,2),sig(6),dd(6,6)
      real*8    sg(2,5),shpw(4,2),shpt(4,2),shpu(2,2),nxi(3)

      save

c     Set values for history storage, etc.

      if(isw.eq.1) then

c     Compute direction cosine terms and member length

      else

        call framtr(d,xl,ndm, le,t)
        dl  = 1.0d0/le
        dl2 = dl*dl
        dl3 = dl*dl2

c       Check

        if(isw.eq.2) then

c       Mass and Geometric stiffness computation for eigen problem

        elseif(isw.eq.5) then

          if(imtyp.eq.1) then !  Mass
            call mass3s(s,r,d(7),d,le,nst,ndm,ndf)

c         Geometric stiffness computation

          elseif(imtyp.eq.2) then

            le     = sqrt((xl(1,2) - xl(1,1))**2
     &                  + (xl(2,2) - xl(2,1))**2
     &                  + (xl(3,2) - xl(3,1))**2)
            sig(1) = d(1)*d(32)*(ul(1,2,1) - ul(1,1,1))*0.5d0
            lint = 3
            call int1d(lint, sg)
            do l = 1,lint

c             Shape functions

              call shp1dh(sg(1,l),le,shpw,shpt)
              shpu(1,1) = -dl
              shpu(1,2) = -shpu(1,1)
              shpu(2,1) = 0.5d0 - 0.5d0*sg(1,l)
              shpu(2,2) = 0.5d0 + 0.5d0*sg(1,l)

              sig(2) = sig(1)*sg(2,l)

              i1 = 0
              do ii = 1,2

                nxi(1) = shpu(1,ii)*sig(2)
                nxi(2) = shpw(1,ii)*sig(2)
                nxi(3) = shpt(1,ii)*sig(2)
                j1 = 0
                do jj = 1,2

                  s(i1+1,j1+1) = s(i1+1,j1+1) - nxi(1)*shpu(1,jj)
                  s(i1+2,j1+2) = s(i1+2,j1+2) - nxi(2)*shpw(1,jj)
                  s(i1+2,j1+3) = s(i1+2,j1+3) - nxi(2)*shpt(1,jj)
                  s(i1+3,j1+2) = s(i1+3,j1+2) - nxi(3)*shpw(1,jj)
                  s(i1+3,j1+3) = s(i1+3,j1+3) - nxi(3)*shpt(1,jj)

                  j1 = j1 + ndf
                end do ! jj
                i1 = i1 + ndf
              end do ! ii
            end do ! ll
          endif
          call bm3trn(s,t,nst,ndf,1)

c       Compute tangent and residual

        elseif(isw.eq.3 .or. isw.eq.6) then

c         Resultant elastic integrated form

          if(nint(d(100)).eq.0) then

c           Compute axial stiffness terms

            i2       =  ndf + 1
            s(1,1)   =  d(1)*d(32)*dl
            s(i2,1)  = -s(1,1)
            s(1,i2)  = -s(1,1)
            s(i2,i2) =  s(1,1)

c           Compute torsional stiffness terms

            i2       =  ndf + 4
            s(4,4)   =  d(1)*d(36)*dl/(1.d0 + d(2))*0.5d0
            s(i2,4)  = -s(4,4)
            s(4,i2)  = -s(4,4)
            s(i2,i2) =  s(4,4)

c         Compute bending stiffness terms for z-displacements

            i1       = ndf + 3
            i2       = ndf + 5

            s(3,3)   =  12.0d0*d(1)*d(33)*dl3
            s(3,i1)  = -s(3,3)
            s(i1,3)  = -s(3,3)
            s(i1,i1) =  s(3,3)

            s(5,i2)  = 2.d0*d(1)*d(33)*dl
            s(5,5)   =  s(5,i2) + s(5,i2)
            s(i2,5)  =  s(5,i2)
            s(i2,i2) =  s(5,5)

            s(3,5)   = -6.d0*d(1)*d(33)*dl2
            s(5,3)   =  s(3,5)
            s(3,i2)  =  s(3,5)
            s(i2,3)  =  s(3,i2)
            s(5,i1)  = -s(3,5)
            s(i1,5)  =  s(5,i1)
            s(i1,i2) = -s(3,5)
            s(i2,i1) = s(i1,i2)

c           Compute bending stiffness terms for y-displacement

            i1       = ndf + 2
            i2       = ndf + 6

            s(2,2)   =  12.d0*d(1)*d(34)*dl3
            s(2,i1)  = -s(2,2)
            s(i1,2)  = -s(2,2)
            s(i1,i1) =  s(2,2)

            s(6,i2)  = 2.d0*d(1)*d(34)*dl
            s(6,6)   =  s(6,i2) + s(6,i2)
            s(i2,6)  =  s(6,i2)
            s(i2,i2) =  s(6,6)

            s(2,6)   =  6.d0*d(1)*d(34)*dl2
            s(6,2)   =  s(2,6)
            s(2,i2)  =  s(2,6)
            s(i2,2)  =  s(2,i2)
            s(6,i1)  = -s(2,6)
            s(i1,6)  =  s(6,i1)
            s(i1,i2) = -s(2,6)
            s(i2,i1) =  s(i1,i2)

c           Compute bending stiffness terms for yz-displacement

            if(d(35).ne.0.0d0) then

              i1       = ndf + 3
              i2       = ndf + 5
              j1       = ndf + 2
              j2       = ndf + 6

              s(2,3)   =  12.d0*d(1)*d(35)*dl3
              s(3,2)   =  s(2,3)

              s(2,i1)  = -s(2,3)
              s(i1,2)  = -s(2,3)

              s(j1,3)  = -s(2,3)
              s(3,j1)  = -s(2,3)

              s(j1,i1) =  s(2,3)
              s(i1,j1) =  s(2,3)

              s(6,i2)  = 2.d0*d(1)*d(35)*dl
              s(i2,6)  =  s(6,i2)

              s(6,5)   =  s(6,i2) + s(6,i2)
              s(5,6)   =  s(6,5)

              s(j2,5)  =  s(6,i2)
              s(5,j2)  =  s(6,i2)

              s(j2,i2) =  s(6,5)
              s(i2,j2) =  s(6,5)

              s(2,5)   =  6.d0*d(1)*d(35)*dl2
              s(5,2)   =  s(2,5)

              s(6,3)   =  s(2,5)
              s(3,6)   =  s(2,5)

              s(2,i2)  =  s(2,5)
              s(i2,2)  =  s(2,5)

              s(j2,3)  =  s(2,i2)
              s(3,j2)  =  s(2,i2)

              s(6,i1)  = -s(2,5)
              s(i1,6)  = -s(2,5)

              s(j1,5)  =  s(6,i1)
              s(5,j1)  =  s(6,i1)

              s(j1,i2) = -s(2,5)
              s(i2,j1) =  s(j1,i2)

              s(i1,j2) = -s(2,5)
              s(j2,i1) =  s(j1,i2)

            endif

c           Transform to global coordinate displacements

            call bm3trn(s,t,nst,ndf,1)

            do i = 1,6

c             Do deformation state

              do j = 1,6
                r(i,1) = r(i,1)
     &                 - s(i    ,j    )*(ul(j,1,1) + d(78)*ul(j,1,4))
     &                 - s(i    ,j+ndf)*(ul(j,2,1) + d(78)*ul(j,2,4))
                r(i,2) = r(i,2)
     &                 - s(i+ndf,j    )*(ul(j,1,1) + d(78)*ul(j,1,4))
     &                 - s(i+ndf,j+ndf)*(ul(j,2,1) + d(78)*ul(j,2,4))
              end do ! j
            end do ! i

            ctan1 = ctan(1) + d(78)*ctan(2)
            do j = 1,nst
              do i = 1,nst
                s(i,j) = s(i,j)*ctan1
              end do ! i
            end do ! j

c         Integrated cross section formulation

          else

c           Transform displacements

            do i = 1,6
              ul(i,1,2) = ul(i,1,1)
              ul(i,2,2) = ul(i,2,1)
            end do ! i
            do i = 1,3
              ul(i  ,1,1) = t(i,1)*ul(1,1,2)
     &                    + t(i,2)*ul(2,1,2)
     &                    + t(i,3)*ul(3,1,2)
              ul(i+3,1,1) = t(i,1)*ul(4,1,2)
     &                    + t(i,2)*ul(5,1,2)
     &                    + t(i,3)*ul(6,1,2)
              ul(i  ,2,1) = t(i,1)*ul(1,2,2)
     &                    + t(i,2)*ul(2,2,2)
     &                    + t(i,3)*ul(3,2,2)
              ul(i+3,2,1) = t(i,1)*ul(4,2,2)
     &                    + t(i,2)*ul(5,2,2)
     &                    + t(i,3)*ul(6,2,2)
            end do ! i

c           Set quadrature order

            lint = nel
            if(nint(d(182)).gt.0) then
              call int1dn(lint, sg)
            else
              call int1d(lint, sg)
            endif

c           Quadrature loop

            nh = nint(d(15))
            nn = 0
            do l = 1,lint

              call shp1dh(sg,le,shpw,shpt)

c             Compute residual for stresses

              call bm3res(d,hr(nh1+nn),hr(nh2+nn),nh, eps, sig, dd, isw)

c             Compute tangent array

              nn = nn + nh
            end do ! l

          endif

c         Set body loading factors

          call fbody3d(d,xl, r, ndm,ndf, isw)

c         Inertia computation

          if(ndfo(1).gt.0 .or. shflg) then
            ctan3 = ctan(3) + d(77)*ctan(2)
            do i = 1,12
              do j = 1,12
                sm(j,i) = 0.0d0
              end do ! i
              pm(i) = 0.0d0
            end do ! i
            call mass3s(sm,pm,d(7),d,le,12,ndm,6)
            call bm3trn(sm,t,12,6,1)
            i1 = 0
            i2 = 0
            do i = 1,2
              do k = 1,6
                j1 = 0
                j2 = 0
                do j = 1,2
                  do l = 1,6
                    r(k,i)  = r(k,i)
     &                      - sm(k+i2,l+j2)*(ul(l,j,5)+d(77)*ul(l,j,4))
                    s(k+i1,l+j1) = s(k+i1,l+j1) + sm(k+i2,l+j2)*ctan3
                  end do ! l
                  j1 = j1 + ndf
                  j2 = j2 + 6
                end do ! j
              end do ! k
              i1 = i1 + ndf
              i2 = i2 + 6
            end do ! i
          endif

        endif

      endif

c     Output member forces

      if(isw.eq.4 .or. isw.eq.8) then

c       Transform displacements

        do i = 1,ndf
          ul(i,1,2) = ul(i,1,1)
          ul(i,2,2) = ul(i,2,1)
        end do ! i
        do i = 1,3
          ul(i  ,1,1) = t(i,1)*ul(1,1,2)
     &                + t(i,2)*ul(2,1,2)
     &                + t(i,3)*ul(3,1,2)
          ul(i+3,1,1) = t(i,1)*ul(4,1,2)
     &                + t(i,2)*ul(5,1,2)
     &                + t(i,3)*ul(6,1,2)
          ul(i  ,2,1) = t(i,1)*ul(1,2,2)
     &                + t(i,2)*ul(2,2,2)
     &                + t(i,3)*ul(3,2,2)
          ul(i+3,2,1) = t(i,1)*ul(4,2,2)
     &                + t(i,2)*ul(5,2,2)
     &                + t(i,3)*ul(6,2,2)
        end do ! i

c       Compute member strains

        eps(1,1) = (ul(1,2,1) - ul(1,1,1))*dl
        eps(1,2) =  eps(1,1)

        eps(2,1) = 12.d0*(ul(2,2,1) - ul(2,1,1))*dl3
     &           -  6.d0*(ul(6,2,1) + ul(6,1,1))*dl2
        eps(2,2) =  eps(2,1)

        eps(3,1) = 12.d0*(ul(3,2,1) - ul(3,1,1))*dl3
     &           +  6.d0*(ul(5,2,1) + ul(5,1,1))*dl2
        eps(3,2) =  eps(3,1)

        eps(4,1) = (ul(4,2,1) - ul(4,1,1))*dl
        eps(4,2) =  eps(4,1)

        eps(5,1) =-(2.d0* ul(5,2,1) + 4.d0*ul(5,1,1))*dl
     &           -  6.d0*(ul(3,2,1) - ul(3,1,1))*dl2
        eps(5,2) = (4.d0* ul(5,2,1) + 2.d0*ul(5,1,1))*dl
     &           +  6.d0*(ul(3,2,1) - ul(3,1,1))*dl2

        eps(6,1) =-(2.d0* ul(6,2,1) + 4.d0*ul(6,1,1))*dl
     &           +  6.d0*(ul(2,2,1) - ul(2,1,1))*dl2
        eps(6,2) = (4.d0* ul(6,2,1) + 2.d0*ul(6,1,1))*dl
     &           -  6.d0*(ul(2,2,1) - ul(2,1,1))*dl2

c       Compute elastic (resultant) member forces

        do i = 1,2
          fi(1,i) = d(1)*d(32)*eps(1,i)
          fi(2,i) = d(1)*d(34)*eps(2,i)
          fi(3,i) = d(1)*d(33)*eps(3,i)
          fi(4,i) = d(1)*d(36)*eps(4,i)/(1.d0 + d(2))*0.5d0
          fi(5,i) = d(1)*d(33)*eps(5,i)
          fi(6,i) = d(1)*d(34)*eps(6,i)
        end do ! i

c       Output member forces

        if(isw.eq.4) then
          mct = mct - 1
          if(mct.le.0) then
            write(iow,2001) o,head
            if(ior.lt.0) then
              write(*,2001) o,head
            endif
            mct = 50
          endif

          write(iow,2002) n,ma,(xl(i,1),i=1,3),(fi(i,1),i=1,6),
     &                         (xl(i,2),i=1,3),(fi(i,2),i=1,6)
          if(ior.lt.0) then
            write(*,2002) n,ma,(xl(i,1),i=1,3),(fi(i,1),i=1,6),
     &                         (xl(i,2),i=1,3),(fi(i,2),i=1,6)
          endif

c       Project member forces resultants to nodes

        else

          call frcn3d(fi,r,s)

        endif

      endif

c     Format statements

2001  format(a1,20a4//5x,'3-D Frame Element Forces'//
     &   '    Elmt  Mat     x-Coor     y-Coor     z-Coor'/
     & 7x,'I-end:      Force    1-Shear    2-Shear   1-Torque',
     &     '   1-Moment   2-Moment'/
     &   '                  x-Coor     y-Coor     z-Coor'/
     & 7x,'J-end:      Force    1-Shear    2-Shear   1-Torque',
     &     '   1-Moment   2-Moment'/1x,78('-'))

2002  format(i8,i5,1p,3e11.3/13x,1p,6e11.3/
     &       13x,  1p,3e11.3/13x,1p,6e11.3/1x)

      end

      subroutine framtr(d,xl,ndm, le,t)

c-----[--.----+----.----+----.-----------------------------------------]
c     Purpose: Three dimensional frame element

c     Inputs:
c         d(*)      - Material parameters
c         xl(ndm,*) - Element coordinates
c         ndm       - Dimension for 'xl'

c     Outputs:
c         le        - Element length
c         t(3,3)    - Transformation array
c-----[--.----+----.----+----.-----------------------------------------]
      implicit   none

      include   'debugs.h'
      include   'eldata.h'
      include   'iofile.h'
      include   'refnd.h'

      character  wd(4)*6
      integer    ndm, nmax, nmin, nnit, i
      real*8     d(*),xl(ndm,*), le, t(3,3), dl,theta, tol

      data       tol / 1.d-08 /

      data       wd  /'Node','Vector','Polar','Axial'/

      lref    = nint(d(96))
      refx(1) = d(97)
      refx(2) = d(98)
      refx(3) = d(99)

      t(1,1) = xl(1,2) - xl(1,1)
      t(1,2) = xl(2,2) - xl(2,1)
      t(1,3) = xl(3,2) - xl(3,1)
      le  = sqrt(t(1,1)*t(1,1)+t(1,2)*t(1,2)+t(1,3)*t(1,3))
      dl  = 1.0d0/le
      t(1,1) = t(1,1)*dl
      t(1,2) = t(1,2)*dl
      t(1,3) = t(1,3)*dl

c     Reference Node

      if    (lref.eq.1) then    ! Nodal reference vector
        t(3,1) = refx(1) - xl(1,1)
        t(3,2) = refx(2) - xl(2,1)
        t(3,3) = refx(3) - xl(3,1)

c     Reference vector

      elseif(lref.eq.2) then    ! Specified reference vector
        t(3,1) = refx(1)
        t(3,2) = refx(2)
        t(3,3) = refx(3)

c     Reference polar

      elseif(lref.eq.3) then    ! Polar reference vector
        t(3,1) = 0.5d0*(xl(1,1) + xl(1,2))
        t(3,2) = 0.5d0*(xl(2,1) + xl(2,2))
        theta  = atan2(t(3,2),t(3,1))
        dl     = sqrt(t(3,1)**2 + t(3,2)**2)
        t(3,1) = dl*cos(theta)
        t(3,2) = dl*sin(theta)
        t(3,3) = 0.0d0
      elseif(lref.eq.4) then    ! Axial reference vector
        nmax = 1
        nmin = 1
        do i = 2,3
          if(t(1,i).gt.t(1,nmax)) then
            nmax = i
          endif
          if(t(1,i).lt.t(1,nmin)) then
            nmin = i
          endif
        end do ! i
        nnit      =  6 - nmax - nmin
        t(3,nmax) = -t(1,nmin)
        t(3,nnit) =  t(1,nnit)
        t(3,nmin) =  t(1,nmax)
      else
        write(iow,3000)
        if(ior.lt.0) then
          write(*,3000)
        endif
        call plstop()
      endif

      t(2,1) = (t(3,2)*t(1,3) - t(3,3)*t(1,2))
      t(2,2) = (t(3,3)*t(1,1) - t(3,1)*t(1,3))
      t(2,3) = (t(3,1)*t(1,2) - t(3,2)*t(1,1))
      dl  = sqrt(t(2,1)*t(2,1)+t(2,2)*t(2,2)+t(2,3)*t(2,3))
      if(dl.lt.tol*le) then
        write(  *,3001) n,wd(lref),(i,refx(i),i=1,3)
        write(iow,3001) n,wd(lref),(i,refx(i),i=1,3)
        call mprint(t,3,3,3,'T_frame')
        call plstop()
      else
        dl     = 1.0d0/dl
        t(2,1) = t(2,1)*dl
        t(2,2) = t(2,2)*dl
        t(2,3) = t(2,3)*dl
        t(3,1) = t(1,2)*t(2,3) - t(1,3)*t(2,2)
        t(3,2) = t(1,3)*t(2,1) - t(1,1)*t(2,3)
        t(3,3) = t(1,1)*t(2,2) - t(1,2)*t(2,1)
      endif

c     Formats

3000  format('    *ERROR*  No Reference type specified for frame',
     &       ' element'/ 13x,'Add to MATErial data')

3001  format('    *ERROR*  Bad reference data in frame element:',i9/
     &       '             Reference type specified: ',a:/
     &      (13x,'V(',i1,') = ',1p,1e12.5))

      end

      subroutine frcn3d(fi,dt,st)

      implicit  none

      include  'cdata.h'
      include  'strnum.h'

      integer   i,j
      real*8    dt(*),st(nen,*),fi(6,*)

      save

c     Stress projections

      do i = 1,2
        dt(i) = 1.d0
        do j = 1,6
          st(i,j) = fi(j,i)
        end do ! j
      end do ! i
      iste = 6

      end

      subroutine b3tubm(d, in)

c-----[--.----+----.----+----.-----------------------------------------]
c     Purpose: Compute beam mass for tube section:

c     Inputs:
c        d(*)       - Material parameters
c        stype      - Section type

c     Outputs:
c        in(6,6)    - Mass-inertia array
c-----[--.----+----.----+----.-----------------------------------------]
      implicit none

      include 'pconstant.h'

      real*8   rad,thick,area, d(*),in(6,6)

      save

c     Compute cross section radius, thickness, sector area

      rad   = d(103)
      thick = d(104)
      area  = pi*rad*thick*d(4)

c     Set cross sectional  properties

      in(1,1) = 2.d0*area
      in(2,2) = in(1,1)
      in(3,3) = in(1,1)
      in(5,5) = area*rad*rad
      in(4,4) = in(5,5)*2.0d0
      in(5,5) = in(5,5)*d(8)
      in(6,6) = in(5,5)

      end

      subroutine b3rctm(d, in)

c-----[--.----+----.----+----.-----------------------------------------]
c     Purpose: Compute beam mass for rectangular sections:

c     Inputs:
c        d(*)       - Material parameters

c     Outputs:
c        in(6,6)    - Cross section mass-inertia array
c-----[--.----+----.----+----.-----------------------------------------]
      implicit none

      integer  nrct, nr,nqy,nqz,ny,nz
      real*8   yy,zz,yl,zl,yr,zr,da,ww
      real*8   d(*),in(6,6), sy(2,11),sz(2,11)

      save

c     Compute constitution using Gauss-Lobbato quadrature in layers

      nrct  = int(d(101))
      do nr = 1, nrct
        nqy = nint(d(101+5*nr))/100
        nqz = mod(nint(d(101+5*nr)),100)
        call int1dl(nqy,sy)
        call int1dl(nqz,sz)
        yl = d(97 +5*nr)
        zl = d(98 +5*nr)
        yr = d(99 +5*nr)
        zr = d(100+5*nr)
        da = 0.25d0*(yr-yl)*(zr-zl)*d(4)
        do nz = 1, nqz
          zz = 0.5d0*((1.d0 - sz(1,nz))*zl + (1.d0 + sz(1,nz))*zr)
          do ny = 1, nqy
            yy      = 0.5d0*((1.d0-sy(1,ny))*yl + (1.d0+sy(1,ny))*yr)
            ww      = sy(2,ny)*sz(2,nz)*da
            in(1,1) = in(1,1) + ww
            in(5,5) = in(5,5) + ww*zz*zz
            in(6,6) = in(6,6) + ww*yy*yy
            in(5,6) = in(5,6) + ww*yy*zz
          end do !nz
        end do !ny
      end do !nr
      in(2,2) = in(1,1)
      in(3,3) = in(1,1)
      in(4,4) = in(5,5) + in(6,6)
      in(6,5) = in(5,6)*d(8)
      in(5,6) = in(6,5)
      in(5,5) = in(5,5)*d(8)
      in(6,6) = in(6,6)*d(8)

      end

      subroutine b3secm(d,stype, in)

c-----[--.----+----.----+----.-----------------------------------------]
c     Purpose: Compute beam mass for shaped sections:
c              stype = 1: Wide flange
c              stype = 2: Channel
c              stype = 3: Angle
c              stype = 4: Solid circular

c     Inputs:
c        d(*)       - Material parameters
c        stype      - Section type

c     Outputs:
c        in(6,6)    - Mass-inertia array
c-----[--.----+----.----+----.-----------------------------------------]
      implicit none

      include 'counts.h'
      include 'eldata.h'
      include 'elplot.h'
      include 'pconstant.h'
      include 'tdata.h'

      integer  ii,stype,nqudr
      real*8   xy(3,17),sw(3,17)
      real*8   d(*),in(6,6)
      real*8   yy,zz,hh,wt,wb,tu,tb,tw,a1,a2,a3,ay,az, ww

      save

c     Wide flange section

      if(stype.eq.1) then

        hh = d(101)
        wt = d(102)
        wb = d(103)
        tu = d(104)
        tb = d(105)
        tw = d(106)
        a1 = (hh-tu-tb)*tw
        a2 = wt*tu
        a3 = wb*tb
        ay = a1*(hh-tu+tb)*0.5d0 + a2*(hh-0.5d0*tu) + a3*tb*0.5d0
        ay = ay/(a1+a2+a3)

        xy(1, 1) = -0.50d0*wb
        xy(2, 1) = -ay
        xy(3, 1) =  0.25d0*a3

        xy(1, 2) =  0.50d0*wb
        xy(2, 2) = -ay
        xy(3, 2) =  0.25d0*a3

        xy(1, 3) =  0.50d0*wb
        xy(2, 3) =  tb-ay
        xy(3, 3) =  0.25d0*a3

        xy(1, 4) =  0.50d0*tw
        xy(2, 4) =  tb-ay
        xy(3, 4) =  0.25d0*a1

        xy(1, 5) =  0.50d0*tw
        xy(2, 5) =  hh-tu-ay
        xy(3, 5) =  0.25d0*a1

        xy(1, 6) =  0.50d0*wt
        xy(2, 6) =  hh-tu-ay
        xy(3, 6) =  0.25d0*a2

        xy(1, 7) =  0.50d0*wt
        xy(2, 7) =  hh-ay
        xy(3, 7) =  0.25d0*a2

        xy(1, 8) = -0.50d0*wt
        xy(2, 8) =  hh-ay
        xy(3, 8) =  0.25d0*a2

        xy(1, 9) = -0.50d0*wt
        xy(2, 9) =  hh-tu-ay
        xy(3, 9) =  0.25d0*a2

        xy(1,10) = -0.50d0*tw
        xy(2,10) =  hh-tu-ay
        xy(3,10) =  0.25d0*a1

        xy(1,11) = -0.50d0*tw
        xy(2,11) =  tb-ay
        xy(3,11) =  0.25d0*a1

        xy(1,12) = -0.50d0*wb
        xy(2,12) =  tb-ay
        xy(3,12) =  0.25d0*a3

        nqudr    =  12

c     Channel

      elseif(stype.eq.2) then

        hh = d(101)
        wt = d(102)
        wb = d(103)
        tu = d(104)
        tb = d(105)
        tw = d(106)
        a1 = (hh-tu-tb)*tw
        a2 = wt*tu
        a3 = wb*tb
        ay = a1*(hh-tu+tb)*0.5d0 + a2*(hh-0.5d0*tu) + a3*tb*0.5d0
        ay = ay/(a1+a2+a3)
        az = (a1*tw + a2*wt + a3*wb)*0.5d0
        az = az/(a1+a2+a3)

        xy(1, 1) = -az
        xy(2, 1) = -ay
        xy(3, 1) =  0.25d0*a3

        xy(1, 2) =  wb - az
        xy(2, 2) = -ay
        xy(3, 2) =  0.25d0*a3

        xy(1, 3) =  wb - az
        xy(2, 3) =  tb - ay
        xy(3, 3) =  0.25d0*a3

        xy(1, 4) =  tw - az
        xy(2, 4) =  tb - ay
        xy(3, 4) =  0.25d0*a1

        xy(1, 5) =  tw - az
        xy(2, 5) =  hh - tu - ay
        xy(3, 5) =  0.25d0*a1

        xy(1, 6) =  wt - az
        xy(2, 6) =  hh - tu - ay
        xy(3, 6) =  0.25d0*a2

        xy(1, 7) =  wt - az
        xy(2, 7) =  hh - ay
        xy(3, 7) =  0.25d0*a2

        xy(1, 8) = -az
        xy(2, 8) =  hh - ay
        xy(3, 8) =  0.25d0*a2

        xy(1, 9) = -az
        xy(2, 9) =  hh - tu - ay
        xy(3, 9) =  0.25d0*(a1 + a2)

        xy(1,10) = -az
        xy(2,10) =  tb - ay
        xy(3,10) =  0.25d0*(a1 + a3)

        nqudr    =  10

c     Angle

      elseif(stype.eq.3) then

        hh = d(101)
        wb = d(102)
        tw = d(103)
        tb = d(104)
        a1 = (hh-tu-tb)*tw
        a3 = wb*tb
        ay = (a1*(hh-tb) + a3*tb)*0.5d0
        ay = ay/(a1+a3)
        az = (a1*tw + a3*wb)*0.5d0
        az = az/(a1+a3)

        xy(1, 1) = -az
        xy(2, 1) = -ay
        xy(3, 1) =  0.25d0*a3

        xy(1, 2) =  wb - az
        xy(2, 2) = -ay
        xy(3, 2) =  0.25d0*a3

        xy(1, 3) =  wb - az
        xy(2, 3) =  tb - ay
        xy(3, 3) =  0.25d0*a3

        xy(1, 4) =  tw - az
        xy(2, 4) =  tb - ay
        xy(3, 4) =  0.25d0*a1

        xy(1, 5) =  tw - az
        xy(2, 5) =  hh - ay
        xy(3, 5) =  0.25d0*a1

        xy(1, 6) = -az
        xy(2, 6) =  hh - ay
        xy(3, 6) =  0.25d0*a1

        xy(1, 7) = -az
        xy(2, 7) =  tb - ay
        xy(3, 7) =  0.25d0*(a1 + a3)

        nqudr    =  7

c     Solid circular

      elseif(stype.eq.4) then

        hh    = d(101)             ! radius
        ii    = nint(d(102))       ! quadrature order
        a3    = pi*hh*hh  ! area weight

        call int2dc(ii,nqudr,sw)

        do ii = 1,nqudr
          xy(1,ii) = hh*sw(1,ii)
          xy(2,ii) = hh*sw(2,ii)
          xy(3,ii) = a3*sw(3,ii)
        end do ! ii

      else
        write(*,*) ' FRAMS3D: Section ',stype,' NOT CODED'
      endif

c     Compute mass effects using Gauss-Lobbato quadrature in layers

      do ii = 1, nqudr
        yy = xy(1,ii)
        zz = xy(2,ii)
        ww = xy(3,ii)*d(4)
        in(1,1) = in(1,1) + ww
        in(1,5) = in(1,5) + ww*zz
        in(1,6) = in(1,6) - ww*yy
        in(5,5) = in(5,5) + ww*zz*zz
        in(6,6) = in(6,6) + ww*yy*yy
        in(5,6) = in(5,6) + ww*yy*zz
      end do !ii
      in(2,2) =  in(1,1)
      in(2,4) = -in(1,5)
      in(3,3) =  in(1,1)
      in(3,4) = -in(1,6)
      in(4,2) = -in(1,5)
      in(4,3) = -in(1,6)
      in(4,4) =  in(5,5) + in(6,6)
      in(5,1) =  in(1,5)
      in(6,1) =  in(1,6)
      in(6,5) =  in(5,6)

      end
