c$Id:$
      subroutine frans3d(d,ul,xl,s,r,ndf,ndm,nst,isw)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c     Linear Three dimensional frame element

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
      include  'bm2str.h'
      include  'cdata.h'
      include  'debugs.h'
      include  'eldata.h'
      include  'eltran.h'
      include  'evdata.h'
      include  'hdata.h'
      include  'iofile.h'
      include  'part0.h'
      include  'pmod2d.h'
      include  'prstrs.h'
      include  'ptdat6.h'
      include  'rdata.h'
      include  'comblk.h'

      logical   noconv
      integer   ndf,ndm,nst,isw, i,ii,i1, j,jj,j1, k,l,ll,lint
      integer   itmax, nh,nn, nqud
      real*8    le,hle,dl,dl2,dl3,ctan1,ctan3
      real*8    itol, ben,hen,uen,duen, energy,dva,dvi
      real*8    d(*),xl(ndm,*),ul(ndf,nen,*),s(nst,*),r(ndf,*)
      real*8    t(3,3),fi(6,4),sm(12,12),pm(12)
      real*8    xx(3), aa(6,6,4), eps(6,4), sig(6,4)
      real*8    bmat(6,6,2),baii(6,6),nxi(3), pp(6,4),gen(6,4)
      real*8    dx(4),dudx(4),sg(2,4)
      real*8    shpu(2,2,4),shpt(4,2,4),shpw(4,2,4)

      save

      data      itol / 1.d-8 /, itmax / 5 /

c     Compute direction cosine terms and member length

      if(isw.eq.1) then

c       Increment history storage if necessary

        nh1 = nh1 + 1   ! Enhanced strain parameter

        do i = 1,6
          do j = 1,6
            bmat(j,i,1) = 0.0d0
            bmat(j,i,2) = 0.0d0
          end do ! j
        end do ! i

c       Set effective modulus to E

        d(21) = d(1)

      elseif(isw.eq.7) then

      else

        call framtr(d,xl,ndm, le,t)
        dl  = 1.0d0/le
        dl2 = dl*dl
        dl3 = dl*dl2

c       Compute stiffness and residual

        if(isw.eq.3 .or. isw.eq.6) then

c         Quadrature terms

          nqud = nint(d(82))
c         lint = nel + 1
          lint = nel
          if(nint(d(182)).gt.0) then
            call int1dn(lint, sg)
          else
            call int1d(lint, sg)
          endif

c         Transform to local coordinates

          call bm3trn (ul(1,1,1),t,ndf,ndf,2)
          call bm3trn (ul(1,1,2),t,ndf,ndf,2)
          call bm3trn (ul(1,1,4),t,ndf,ndf,2)
          call bm3trn (ul(1,1,5),t,ndf,ndf,2)
          hle = le*0.5d0

          do ll = 1,lint

c           Shape functions

            call shp1dh(sg(1,ll),le,shpw(1,1,ll),shpt(1,1,ll))
            shpu(1,1,ll) = -dl
            shpu(1,2,ll) = -shpu(1,1,ll)
            shpu(2,1,ll) = 0.5d0 - 0.5d0*sg(1,ll)
            shpu(2,2,ll) = 0.5d0 + 0.5d0*sg(1,ll)

            dx(ll)      = sg(2,ll)*hle

c           Form displacement derivatives from nodal displacements

            dudx(ll)  = (ul(1,2,1) - ul(1,1,1))*dl

c           Strains: V_1, V_2, N, M_1, M_2, T

            eps(1,ll) = shpw(3,1,ll)*ul(2,1,1) + shpw(3,2,ll)*ul(2,2,1)
     &                + shpt(3,1,ll)*ul(6,1,1) + shpt(3,2,ll)*ul(6,2,1)

            eps(2,ll) = shpw(3,1,ll)*ul(3,1,1) + shpw(3,2,ll)*ul(3,2,1)
     &                - shpt(3,1,ll)*ul(5,1,1) - shpt(3,2,ll)*ul(5,2,1)

            eps(3,ll) = dudx(ll)

            eps(4,ll) = shpw(2,1,ll)*ul(3,1,1) + shpw(2,2,ll)*ul(3,2,1)
     &                - shpt(2,1,ll)*ul(5,1,1) - shpt(2,2,ll)*ul(5,2,1)
            eps(4,ll) = -eps(4,ll)

            eps(5,ll) = shpw(2,1,ll)*ul(2,1,1) + shpw(2,2,ll)*ul(2,2,1)
     &                + shpt(2,1,ll)*ul(6,1,1) + shpt(2,2,ll)*ul(6,2,1)

            eps(6,ll) = (ul(4,2,1) - ul(4,1,1))*dl
          end do ! ll

c         Enhanced strain computation

          if(etype.eq.3) then
            uen = hr(nh1)
            ii  = 0
            noconv = .true.
            do while(noconv)

              ii  = ii + 1

c             Zero enhanced terms

              ben = 0.0d0
              hen = 0.0d0

              nn    = 1
              nh    = nint(d(15))
              do ll = 1,lint

c               Update axial strain for enhanced term

                eps(3,ll) = eps(3,ll) + sg(1,ll)*uen

c               Call beam 3d model

                call bm3mat(d,hr(nh1+nn),hr(nh2+nn),nh,
     &                      eps(1,ll),sig(1,ll),aa(1,1,ll),isw)

                hen = hen + sg(1,ll)*aa(3,3,ll)*sg(1,ll)*dx(ll)*ctan(1)
                ben = ben - sg(1,ll)*sig(3,ll)*dx(ll)
                nn  = nn + nh*nqud

              end do ! ll

              hen  = 1.d0/ hen
              duen = ben * hen
              uen  = uen + duen

              if(abs(duen).le.itol*abs(uen) .or. ben.eq.0.0d0) then
                noconv = .false.
              elseif(ii.gt.itmax) then
                noconv = .false.
              endif

            end do ! while

c           Save enhance mode parameter

            hr(nh2) = uen

          else

            hen = 0.0d0

            nn    = 0
            nh    = nint(d(15))
            do ll = 1,lint

c             Call beam 3d model

              call bm3mat(d,hr(nh1+nn),hr(nh2+nn),nh,
     &                    eps(1,ll),sig(1,ll),aa(1,1,ll),isw)
              nn = nn + nh*nqud
            end do ! ll

          end if ! etype

c         Stiffness and residual computation

          if(isw.eq.3 .or. isw.eq.6) then

c           Zero enhanced coupling array

            do i = 1,6
              gen(i,1) = 0.0d0
              gen(i,2) = 0.0d0
            end do ! i

c           Final tangent form

            do ll = 1,lint

              ctan1 = ctan(1)*dx(ll)
              do j = 1,6
                do i = 1,6
                  aa(i,j,ll) = aa(i,j,ll)*ctan1
                end do ! i
                sig(j,ll) = sig(j,ll)*dx(ll)
              end do ! j

c             Compute strain-displacement matrices for two nodes

              do i = 1,2
                bmat(3,1,i) =  shpu(1,i,ll)

                bmat(4,3,i) = -shpw(2,i,ll)
                bmat(4,5,i) =  shpt(2,i,ll)

                bmat(5,2,i) =  shpw(2,i,ll)
                bmat(5,6,i) =  shpt(2,i,ll)

                bmat(6,4,i) =  shpu(1,i,ll)
              end do ! i

c             Mechanical tangent terms

              i1 = 0
              do ii = 1,nel


c               B^T * AA

                do i = 3,6
                  baii(1,i) =  bmat(3,1,ii)*aa(3,i,ll)
                  baii(2,i) =  bmat(5,2,ii)*aa(5,i,ll)
                  baii(3,i) =  bmat(4,3,ii)*aa(4,i,ll)
                  baii(4,i) =  bmat(6,4,ii)*aa(6,i,ll)
                  baii(5,i) =  bmat(4,5,ii)*aa(4,i,ll)
                  baii(6,i) =  bmat(5,6,ii)*aa(5,i,ll)

c                 Enhanced stiffness

                  gen(i,ii) = gen(i,ii) + baii(4,i)*sg(1,ll)

                end do ! i

c               Residual

                r(1,ii) = r(1,ii) - bmat(3,1,ii)*sig(3,ll)
                r(2,ii) = r(2,ii) - bmat(5,2,ii)*sig(5,ll)
                r(3,ii) = r(3,ii) - bmat(4,3,ii)*sig(4,ll)
                r(4,ii) = r(4,ii) - bmat(6,4,ii)*sig(6,ll)
                r(5,ii) = r(5,ii) - bmat(4,5,ii)*sig(4,ll)
                r(6,ii) = r(6,ii) - bmat(5,6,ii)*sig(5,ll)

c               Tangent

                if(isw.eq.3) then

                  j1 = 0
                  do jj = 1,nel

c                   Material part

                    do j = 1,6
                      do i = 1,6
                        s(i1+i,j1+j) = s(i1+i,j1+j)
     &                               + baii(i,3)*bmat(3,j,jj)
     &                               + baii(i,4)*bmat(4,j,jj)
     &                               + baii(i,5)*bmat(5,j,jj)
     &                               + baii(i,6)*bmat(6,j,jj)
                      end do ! i
                    end do ! j

                    j1 = j1 + ndf
                  end do ! jj
                endif
                i1 = i1 + ndf
              end do ! ii
            end do ! ll

c           Transform stiffness and residual to global coordinates

            if(isw.eq.3) then

c             Static condensation

              j1 = 0
              do jj = 1,nel
                do j = 1,6
                  duen = gen(j,jj)*hen
                  i1 = 0
                  do ii = 1,nel
                    do i = 1,6
                      s(i+i1,j+j1) = s(i+i1,j+j1) - gen(i,ii)*duen
                    end do ! i
                    i1 = i1 + ndf
                  end do ! ii
                end do ! j
                j1 = j1 + ndf
              end do ! jj

              call bm3trn (s,t,nst,ndf,1)
            endif
            call bm3trn (r,t,ndf,ndf,3)

c           Set body loading factors

            call fbody3d(d,xl, r, ndm,ndf, isw)

c           Inertia contributions

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
              ii = 0
              do i = 1,2
                do k = 1,6
                  j1 = 0
                  jj = 0
                  do j = 1,2
                    do l = 1,6
                      r(k,i) = r(k,i)
     &                       - sm(k+ii,l+jj)*(ul(l,j,5)+d(77)*ul(l,j,4))
                      s(k+i1,l+j1) = s(k+i1,l+j1) + sm(k+ii,l+jj)*ctan3
                    end do ! l
                    j1 = j1 + ndf
                    jj = jj + 6
                  end do ! j
                end do ! k
                i1 = i1 + ndf
                ii = ii + 6
              end do ! i
            endif

c         Geometric stiffness computation

          elseif(isw.eq.5 .and. imtyp.eq.2) then

            do ll = 1,lint

              ctan1      = dx(ll)*ctan(1)
              sig(3,ll) = d(1)*d(32)*(eps(3,ll)+sg(1,ll)*uen)*ctan1

              i1 = 0
              do ii = 1,nel

                nxi(1) = shpu(1,ii,ll)*sig(3,ll)
                nxi(2) = shpw(1,ii,ll)*sig(3,ll)
                nxi(3) = shpt(1,ii,ll)*sig(3,ll)
                j1 = 0
                do jj = 1,nel

c                 s(i1+1,j1+1) = s(i1+1,j1+1) - nxi(1)*shpu(1,jj,ll)
                  s(i1+2,j1+2) = s(i1+2,j1+2) - nxi(2)*shpw(1,jj,ll)
                  s(i1+2,j1+3) = s(i1+2,j1+3) - nxi(2)*shpt(1,jj,ll)
                  s(i1+3,j1+2) = s(i1+3,j1+2) - nxi(3)*shpw(1,jj,ll)
                  s(i1+3,j1+3) = s(i1+3,j1+3) - nxi(3)*shpt(1,jj,ll)

                  j1 = j1 + ndf
                end do ! jj
                i1 = i1 + ndf
              end do ! ii
            end do ! ll

c         Compute energy

          elseif(isw.eq.13) then

            dva = hle*d(4)*d(32)
            dvi  =hle*d(4)*d(33)*d(8)

c           Compute internal energy

            nn = 1
            nh = nint(d(15))
            energy = 0.0d0
            do ll = 1,lint

c             Compute energy density from stress and deformation

              call shp1d(sg(1,ll),xl,shpu,ndm,nel,dx(ll))
              dx(ll) = sg(2,ll)*dx(ll)
              call b2mods (d,ul,sig,aa,energy,shpu,
     &                     hr(nh1+nn),hr(nh2+nn),ndf,nen,isw,ll)

c             Accumulate energy

              epl(8) = epl(8) + 0.5d0*energy*dx(ll)

              nn = nn + nh
            end do ! ll

c           Compute kinetic energy for lumped mass

            epl(7) = epl(7) + 0.5d0*dva*(ul(1,1,4)**2 + ul(1,2,4)**2
     &                                 + ul(2,1,4)**2 + ul(2,2,4)**2)
     &                      + 0.5d0*dvi*(ul(3,1,4)**2 + ul(3,2,4)**2)

c         Initialize history variables

          elseif(isw.eq.14) then

            call modl1d(d,le,xx,hr(nh1),hr(nh2),nint(d(15)),1,0,
     &                  aa,sig,isw)

          endif

c       Output member forces

        elseif(isw.eq.4 .or. isw.eq.8) then

c         Output forces

          elseif(isw.eq.4 .or. isw.eq.8) then

c           Member forces

            do i = 1,6
              pp(i,1) = 0.0d0
              pp(i,2) = 0.0d0
            end do ! i
            do ll = 1,lint

c           Stress output

            if(isw.eq.4) then

              do i = 1,ndm
                xx(i) = 0.0d0
                do ii = 1,nel
                  xx(i) = xx(i) + xl(i,ii)*shpu(2,ii,ll)
                end do ! ii
              end do ! i
              mct = mct - 3
              if (mct.le.0) then
                write(iow,2001) o,head
                if(ior.lt.0) write(*,2001) o,head
                mct = 50
              endif
              write(iow,2002) n,ma,(xx(i),i=1,3),
     &                        (sig(i,ll),i=1,6),(eps(i,1),i=1,6)
              if(nout.gt.0) then
                write(iow,2003) (siglr(i),i=1,nout)
                write(iow,2004) (epslr(i),i=1,nout)
              endif
              if(ior.lt.0) then
                write(*,2002) n,ma,(xx(i),i=1,3),
     &                        (sig(i,ll),i=1,6),(eps(i,1),i=1,6)
                if(nout.gt.0) then
                  write(*,2003) (siglr(i),i=1,nout)
                  write(*,2004) (epslr(i),i=1,nout)
                endif
              endif

c           Stress projections save

            else

c             Compute strain-displacement matrices for two nodes

              do i = 1,2
                bmat(1,1,i) =  shpu(1,i,ll)
                bmat(2,1,i) =  0.0d0
                bmat(2,2,i) =  shpw(2,i,ll)
                bmat(2,3,i) =  shpt(2,i,ll)
              end do ! i

c             End forces

              do ii = 1,nel

                do i = 1,3
                  pp(i,ii) = pp(i,ii) - (bmat(1,i,ii)*sig(1,ll)
     &                                +  bmat(2,i,ii)*sig(3,ll))*dx(ll)
                end do ! i

              end do ! ii
            end if
          end do ! ll

c         Projection on end nodes (uses reactions)

          if(isw.eq.8) then
            do i = 1,3
              pp(i,2) = -pp(i,2)
            end do ! i
            call frcn2d(pp,r,s)
          endif

          do i = 1,ndf
            ul(i,1,2) = ul(i,1,1)
            ul(i,2,2) = ul(i,2,1)
          end do ! i
          do i = 1,3
            ul(i  ,1,1) = t(i,1)*ul(1,1,2)
     &                  + t(i,2)*ul(2,1,2)
     &                  + t(i,3)*ul(3,1,2)
            ul(i+3,1,1) = t(i,1)*ul(4,1,2)
     &                  + t(i,2)*ul(5,1,2)
     &                  + t(i,3)*ul(6,1,2)
            ul(i  ,2,1) = t(i,1)*ul(1,2,2)
     &                  + t(i,2)*ul(2,2,2)
     &                  + t(i,3)*ul(3,2,2)
            ul(i+3,2,1) = t(i,1)*ul(4,2,2)
     &                  + t(i,2)*ul(5,2,2)
     &                  + t(i,3)*ul(6,2,2)
          end do ! i

c         Compute member strains

          eps(1,1) = (ul(1,2,1) - ul(1,1,1))*dl
          eps(1,2) =  eps(1,1)

          eps(2,1) = 12.d0*(ul(2,2,1) - ul(2,1,1))*dl3
     &             -  6.d0*(ul(6,2,1) + ul(6,1,1))*dl2
          eps(2,2) =  eps(2,1)

          eps(3,1) = 12.d0*(ul(3,2,1) - ul(3,1,1))*dl3
     &             +  6.d0*(ul(5,2,1) + ul(5,1,1))*dl2
          eps(3,2) =  eps(3,1)

          eps(4,1) = (ul(4,2,1) - ul(4,1,1))*dl
          eps(4,2) =  eps(4,1)

          eps(5,1) =-(2.d0* ul(5,2,1) + 4.d0*ul(5,1,1))*dl
     &             -  6.d0*(ul(3,2,1) - ul(3,1,1))*dl2
          eps(5,2) = (4.d0* ul(5,2,1) + 2.d0*ul(5,1,1))*dl
     &             +  6.d0*(ul(3,2,1) - ul(3,1,1))*dl2

          eps(6,1) =-(2.d0* ul(6,2,1) + 4.d0*ul(6,1,1))*dl
     &             +  6.d0*(ul(2,2,1) - ul(2,1,1))*dl2
          eps(6,2) = (4.d0* ul(6,2,1) + 2.d0*ul(6,1,1))*dl
     &             -  6.d0*(ul(2,2,1) - ul(2,1,1))*dl2

c         Compute elastic (resultant) member forces

          do i = 1,2
            fi(1,i) = d(21)*d(32)*eps(1,i)
            fi(2,i) = d(21)*d(34)*eps(2,i)
            fi(3,i) = d(21)*d(33)*eps(3,i)
            fi(4,i) = d(21)*d(36)*eps(4,i)/(1.d0 + d(2))*0.5d0
            fi(5,i) = d(21)*d(33)*eps(5,i)
            fi(6,i) = d(21)*d(34)*eps(6,i)
          end do ! i

c         Output member forces

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
     &                           (xl(i,2),i=1,3),(fi(i,2),i=1,6)
            if(ior.lt.0) then
              write(*,2002) n,ma,(xl(i,1),i=1,3),(fi(i,1),i=1,6),
     &                           (xl(i,2),i=1,3),(fi(i,2),i=1,6)
            endif

c         Project member forces resultants to nodes

          else

            call frcn3d(fi,r,s)

          endif

c       Compute mass matrix (lumped only)

        elseif(isw.eq.5) then

          call mass3s(s,r,d(7),d,le,nst,ndm,ndf)
          call bm3trn(s,t,nst,ndf,1)

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
     &       13x,  1p,6e11.3/13x,1p,6e11.3/1x)

2003  format('  Stress_Layer',1p,5e13.4)
2004  format('  Strain_Layer',1p,5e13.4)

      end
