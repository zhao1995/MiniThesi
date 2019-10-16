c$Id:$
      subroutine franf3d(d,ul,xl,s,r,ndf,ndm,nst,isw)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Remove transformation of v and a to local frame  12/12/2008
c-----[--.----+----.----+----.-----------------------------------------]
c     Second order strain Three dimensional frame element

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

      logical   noconv, small
      integer   ndf,ndm,nst,isw, i,ii,i1, j,jj,j1, k,l,ll,lint
      integer   itmax, nh,nn, nhv
      real*8    le,hle,dl,dl2,dl3,ctan1,ctan3,ltan3, cfac,lfac
      real*8    itol, ben,hen,uen,duen, energy
      real*8    d(*),xl(ndm,*),ul(ndf,nen,*),s(nst,*),r(ndf,*)
      real*8    t(3,3),sm(12,12),pm(12), cc(6,6,2)
      real*8    xx(3), aa(6,6,4), eps(6,4),sig(6,4), epsl(9,2),sigl(9,2)
      real*8    bmat(6,6,2),baii(6,6),nxi(3), pp(6,2),gen(6,2)
      real*8    dx(4),dudx(4),dvdx(4),dwdx(4),sg(2,4)
      real*8    shpu(2,2,4),shpt(4,2,4),shpw(4,2,4)

      save

      data      itol / 1.d-8 /, itmax / 5 /

c     Set small/large flag

      small = d(18).gt.0.0d0

c     Number quadrature/section

      if(isw.ge.2) then
        if(nint(d(100)).gt.0) then
          nh   = nint(d(15))
          nhv  = nint(d(82))
          nhv  = nh*nhv
        else
          nh   = 0
          nhv  = 0
        endif
      endif

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

c     Compute mass matrix

      elseif(isw.eq.5 .and. imtyp.eq.1) then

        call framtr(d,xl,ndm, le,t)
        call mass3s(s,r,d(7),d,le,nst,ndm,ndf)
        call bm3trn(s,t,nst,ndf,1)
        call bm3trn(r,t,ndf,ndf,3)

      elseif(isw.eq.7) then

c     Initialize history variables

      elseif(isw.eq.14) then

        nn = 1
        ii = nint(d(82))
        do ll = 1,nel+1
          do i = 1,ii
            call modl1d(d,le,xx,hr(nh1+nn),hr(nh2+nn),nh,1,0,
     &                  aa,sig,isw)
            nn = nn + nh
          end do ! i
        end do ! ll

      elseif(isw.ne.12) then

        if(nint(d(182)).gt.0) then
          lint = nel
          call int1dn(lint, sg)
        else
          lint = nel + 1
          call int1d(lint, sg)
        endif

        call framtr(d,xl,ndm, le,t)
        dl  = 1.0d0/le
        dl2 = dl*dl
        dl3 = dl*dl2

c       Quadrature pointer increment

        nn = 1

c       Transform displacement and increment to local coordinates

        call bm3trn (ul(1,1,1),t,ndf,ndf,2)
        call bm3trn (ul(1,1,2),t,ndf,ndf,2)
        hle = le*0.5d0

        do ll = 1,lint

c         Shape functions

          call shp1dh(sg(1,ll),le,shpw(1,1,ll),shpt(1,1,ll))
          shpu(1,1,ll) = -dl
          shpu(1,2,ll) = -shpu(1,1,ll)
          shpu(2,1,ll) = 0.5d0 - 0.5d0*sg(1,ll)
          shpu(2,2,ll) = 0.5d0 + 0.5d0*sg(1,ll)

          dx(ll)      = sg(2,ll)*hle

c         Form displacement derivatives from nodal displacements

          dudx(ll)  = (ul(1,2,1) - ul(1,1,1))*dl

          dvdx(ll)  = shpw(1,1,ll)*ul(2,1,1) + shpw(1,2,ll)*ul(2,2,1)
     &              + shpt(1,1,ll)*ul(6,1,1) + shpt(1,2,ll)*ul(6,2,1)

          dwdx(ll)  = shpw(1,1,ll)*ul(3,1,1) + shpw(1,2,ll)*ul(3,2,1)
     &              - shpt(1,1,ll)*ul(5,1,1) - shpt(1,2,ll)*ul(5,2,1)

c         Strains: V_1, V_2, N, M_1, M_2, T

          eps(1,ll) = shpw(3,1,ll)*ul(2,1,1) + shpw(3,2,ll)*ul(2,2,1)
     &              + shpt(3,1,ll)*ul(6,1,1) + shpt(3,2,ll)*ul(6,2,1)

          eps(2,ll) = shpw(3,1,ll)*ul(3,1,1) + shpw(3,2,ll)*ul(3,2,1)
     &              - shpt(3,1,ll)*ul(5,1,1) - shpt(3,2,ll)*ul(5,2,1)

          if(small) then
            eps(3,ll) = dudx(ll)
            dudx(ll)  = 0.0d0
            dvdx(ll)  = 0.0d0
            dwdx(ll)  = 0.0d0
          else
            eps(3,ll) = dudx(ll) + 0.5d0*(dudx(ll)*dudx(ll)
     &                                  + dvdx(ll)*dvdx(ll)
     &                                  + dwdx(ll)*dwdx(ll))
          endif

          eps(5,ll) = shpw(2,1,ll)*ul(2,1,1) + shpw(2,2,ll)*ul(2,2,1)
     &              + shpt(2,1,ll)*ul(6,1,1) + shpt(2,2,ll)*ul(6,2,1)

          eps(4,ll) = shpw(2,1,ll)*ul(3,1,1) + shpw(2,2,ll)*ul(3,2,1)
     &              - shpt(2,1,ll)*ul(5,1,1) - shpt(2,2,ll)*ul(5,2,1)
          eps(4,ll) = -eps(4,ll)

          eps(6,ll) = (ul(4,2,1) - ul(4,1,1))*dl
        end do ! ll

c       Enhanced strain computation

        if(etype.eq.3) then
          uen = hr(nh1)
          ii  = 0
          noconv = .true.
          do while(noconv)

            ii  = ii + 1

c           Zero enhanced terms

            ben = 0.0d0
            hen = 0.0d0

            nn  = 1
            do ll = 1,lint

c             Update axial strain for enhanced term

              do i = 1,6
                epsl(i,1) = eps(i,ll)
              end do ! i
              epsl(3,1) = epsl(3,1) + sg(1,ll)*uen

c             Call beam 3d model

              call bm3mat(d,hr(nh1+nn),hr(nh2+nn),nh,epsl,sigl,cc,isw)

              hen = hen + sg(1,ll)*cc(3,3,1)*sg(1,ll)*dx(ll)*ctan(1)
              ben = ben - sg(1,ll)*sigl(3,1)*dx(ll)

              do i = 1,6
                do j = 1,6
                  aa(j,i,ll) = cc(j,i,1)
                end do ! j
                sig(i,ll) = sigl(i,1)
              end do ! i

              nn  = nn + nhv
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

c         Save enhance mode parameter

          hr(nh2) = uen

        else

          hen = 0.0d0
          uen = 0.0d0

          nn = 1
          do ll = 1,lint

c           Call beam 3d model

            call bm3mat(d,hr(nh1+nn),hr(nh2+nn),nh,
     &                  eps(1,ll),sig(1,ll),aa(1,1,ll),isw)
            nn = nn + nhv
          end do ! ll

        end if ! etype

c       Stiffness and residual computation

        if(isw.eq.3 .or. isw.eq.6) then

c         Zero enhanced coupling array

          do i = 1,6
            gen(i,1) = 0.0d0
            gen(i,2) = 0.0d0
          end do ! i

c         Final tangent form

          do ll = 1,lint

            ctan1 = ctan(1)*dx(ll)
            do j = 1,6
              do i = 1,6
                aa(i,j,ll) = aa(i,j,ll)*ctan1
              end do ! i
              sig(j,ll) = sig(j,ll)*dx(ll)
            end do ! j

c           Compute strain-displacement matrices for two nodes

            do i = 1,2
              bmat(3,1,i) =  shpu(1,i,ll)*(1.d0 + dudx(ll))
              bmat(3,2,i) =  shpw(1,i,ll)*dvdx(ll)
              bmat(3,3,i) =  shpw(1,i,ll)*dwdx(ll)
              bmat(3,5,i) = -shpt(1,i,ll)*dwdx(ll)
              bmat(3,6,i) =  shpt(1,i,ll)*dvdx(ll)

              bmat(5,2,i) =  shpw(2,i,ll)
              bmat(5,6,i) =  shpt(2,i,ll)

              bmat(4,3,i) = -shpw(2,i,ll)
              bmat(4,5,i) =  shpt(2,i,ll)

              bmat(6,4,i) =  shpu(1,i,ll)
            end do ! i

c           Mechanical tangent terms

            i1 = 0
            do ii = 1,nel

c             B^T * AA

              do i = 3,6
                baii(1,i) =  bmat(3,1,ii)*aa(3,i,ll)
                baii(2,i) =  bmat(3,2,ii)*aa(3,i,ll)
     &                    +  bmat(5,2,ii)*aa(5,i,ll)
                baii(3,i) =  bmat(3,3,ii)*aa(3,i,ll)
     &                    +  bmat(4,3,ii)*aa(4,i,ll)
                baii(4,i) =  bmat(6,4,ii)*aa(6,i,ll)
                baii(5,i) =  bmat(3,5,ii)*aa(3,i,ll)
     &                    +  bmat(4,5,ii)*aa(4,i,ll)
                baii(6,i) =  bmat(3,6,ii)*aa(3,i,ll)
     &                    +  bmat(5,6,ii)*aa(5,i,ll)
              end do ! i

c             Enhanced stiffness

              do i = 1,6
                gen(i,ii) = gen(i,ii) + baii(i,3)*sg(1,ll)
              end do ! i

c             Residual

              r(1,ii) = r(1,ii) - bmat(3,1,ii)*sig(3,ll)
              r(2,ii) = r(2,ii) - bmat(3,2,ii)*sig(3,ll)
     &                          - bmat(5,2,ii)*sig(5,ll)
              r(3,ii) = r(3,ii) - bmat(3,3,ii)*sig(3,ll)
     &                          - bmat(4,3,ii)*sig(4,ll)
              r(4,ii) = r(4,ii) - bmat(6,4,ii)*sig(6,ll)
              r(5,ii) = r(5,ii) - bmat(3,5,ii)*sig(3,ll)
     &                          - bmat(4,5,ii)*sig(4,ll)
              r(6,ii) = r(6,ii) - bmat(3,6,ii)*sig(3,ll)
     &                          - bmat(5,6,ii)*sig(5,ll)

c             Tangent

              if(isw.eq.3) then

                nxi(1) = shpu(1,ii,ll)*sig(3,ll)*ctan(1)
                nxi(2) = shpw(1,ii,ll)*sig(3,ll)*ctan(1)
                nxi(3) = shpt(1,ii,ll)*sig(3,ll)*ctan(1)

                j1 = 0
                do jj = 1,nel

c                 Material part

                  do j = 1,6
                    do i = 1,6
                      s(i1+i,j1+j) = s(i1+i,j1+j)
     &                             + baii(i,3)*bmat(3,j,jj)
     &                             + baii(i,4)*bmat(4,j,jj)
     &                             + baii(i,5)*bmat(5,j,jj)
     &                             + baii(i,6)*bmat(6,j,jj)
                    end do ! i
                  end do ! j

c                 Geometric part

                  if(gflag .and. .not.small) then
                    s(i1+1,j1+1) = s(i1+1,j1+1) + nxi(1)*shpu(1,jj,ll)

                    s(i1+2,j1+2) = s(i1+2,j1+2) + nxi(2)*shpw(1,jj,ll)
                    s(i1+2,j1+6) = s(i1+2,j1+6) + nxi(2)*shpt(1,jj,ll)

                    s(i1+6,j1+2) = s(i1+6,j1+2) + nxi(3)*shpw(1,jj,ll)
                    s(i1+6,j1+6) = s(i1+6,j1+6) + nxi(3)*shpt(1,jj,ll)

                    s(i1+3,j1+3) = s(i1+3,j1+3) + nxi(2)*shpw(1,jj,ll)
                    s(i1+3,j1+5) = s(i1+3,j1+5) - nxi(2)*shpt(1,jj,ll)

                    s(i1+5,j1+3) = s(i1+5,j1+3) - nxi(3)*shpw(1,jj,ll)
                    s(i1+5,j1+5) = s(i1+5,j1+5) + nxi(3)*shpt(1,jj,ll)
                  endif

                  j1 = j1 + ndf
                end do ! jj
              endif
              i1 = i1 + ndf
            end do ! ii
          end do ! ll

c         Transform stiffness and residual to global coordinates

          if(isw.eq.3) then

c           Static condensation

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

c         Set body loading factors

          call fbody3d(d,xl, r, ndm,ndf, isw)

c         Inertia contributions

          if(ndfo(1).gt.0 .or. shflg) then
            ctan3 = ctan(3) + d(77)*ctan(2)
            if(d(7).ge.0.0d0) then
              cfac = d(7)
              lfac = 1.0d0 - cfac
            else
              cfac = 0.0d0
              lfac = 0.0d0
            endif
            ctan3 = cfac*ctan3
            ltan3 = lfac*ctan3
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
                s(k+i1,k+i1) = s(k+i1,k+i1) + pm(k+ii)*ltan3
                r(k,i)       = r(k,i) - pm(k+ii)
     &                       * (ul(k,i,5)+d(77)*ul(k,i,4))*lfac
                do j = 1,2
                  do l = 1,6
                    r(k,i) = r(k,i) - sm(k+ii,l+jj)
     &                     * (ul(l,j,5)+d(77)*ul(l,j,4))*cfac
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

c       Output member forces

        elseif(isw.eq.4 .or. isw.eq.8) then

c         Member forces

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

              eps(3,ll) = eps(3,ll) + sg(1,ll)*uen

              write(iow,2002) n,ma,(xx(i),i=1,3),
     &                        (sig(i,ll),i=3,6),(eps(i,ll),i=3,6)
              if(nout.gt.0) then
                write(iow,2003) (siglr(i),i=1,nout)
                write(iow,2004) (epslr(i),i=1,nout)
              endif
              if(ior.lt.0) then
                write(*,2002) n,ma,(xx(i),i=1,3),
     &                        (sig(i,ll),i=3,6),(eps(i,ll),i=3,6)
                if(nout.gt.0) then
                  write(*,2003) (siglr(i),i=1,nout)
                  write(*,2004) (epslr(i),i=1,nout)
                endif
              endif

c           Stress projections save

            else

c             Compute strain-displacement matrices for two nodes

              do i = 1,2
                bmat(3,1,i) =  shpu(1,i,ll)*(1.d0 + dudx(ll))
                bmat(3,2,i) =  shpw(1,i,ll)*dvdx(ll)
                bmat(3,3,i) =  shpw(1,i,ll)*dwdx(ll)
                bmat(3,5,i) = -shpt(1,i,ll)*dwdx(ll)
                bmat(3,6,i) =  shpt(1,i,ll)*dvdx(ll)
                bmat(5,2,i) =  shpw(2,i,ll)
                bmat(5,6,i) =  shpt(2,i,ll)

                bmat(4,3,i) = -shpw(2,i,ll)
                bmat(4,5,i) =  shpt(2,i,ll)

                bmat(6,4,i) =  shpu(1,i,ll)

c               End forces

                pp(1,i) = pp(1,i) - bmat(3,1,i)*sig(3,ll)*dx(ll)
                pp(2,i) = pp(2,i) - bmat(3,2,i)*sig(3,ll)*dx(ll)
     &                            - bmat(5,2,i)*sig(5,ll)*dx(ll)
                pp(3,i) = pp(3,i) - bmat(3,3,i)*sig(3,ll)*dx(ll)
     &                            - bmat(4,3,i)*sig(4,ll)*dx(ll)
                pp(4,i) = pp(4,i) - bmat(6,4,i)*sig(6,ll)*dx(ll)
                pp(5,i) = pp(5,i) - bmat(3,5,i)*sig(3,ll)*dx(ll)
     &                            - bmat(4,5,i)*sig(4,ll)*dx(ll)
                pp(6,i) = pp(6,i) - bmat(3,6,i)*sig(3,ll)*dx(ll)
     &                            - bmat(5,6,i)*sig(5,ll)*dx(ll)
              end do ! i
            end if
          end do ! ll

c         Projection on end nodes (uses reactions)

          if(isw.eq.8) then
            do i = 1,6
              pp(i,2) = -pp(i,2)
            end do ! i
            call frcn3d(pp,r,s)
          endif

c       Geometric stiffness computation

        elseif(isw.eq.5 .and. imtyp.eq.2) then

          do ll = 1,lint

            sig(3,ll) = sig(3,ll)*dx(ll)

            i1 = 0
            do ii = 1,nel

              nxi(1) = shpu(1,ii,ll)*sig(3,ll)
              nxi(2) = shpw(1,ii,ll)*sig(3,ll)
              nxi(3) = shpt(1,ii,ll)*sig(3,ll)
              j1 = 0
              do jj = 1,nel

                s(i1+1,j1+1) = s(i1+1,j1+1) - nxi(1)*shpu(1,jj,ll)
                s(i1+2,j1+2) = s(i1+2,j1+2) - nxi(2)*shpw(1,jj,ll)
                s(i1+2,j1+3) = s(i1+2,j1+3) - nxi(2)*shpt(1,jj,ll)
                s(i1+3,j1+2) = s(i1+3,j1+2) - nxi(3)*shpw(1,jj,ll)
                s(i1+3,j1+3) = s(i1+3,j1+3) - nxi(3)*shpt(1,jj,ll)

                j1 = j1 + ndf
              end do ! jj
              i1 = i1 + ndf
            end do ! ii
          end do ! ll

c       Compute elastic stored and kinetic energy

        elseif(isw.eq.13) then
          if(ndfo(1).gt.0 .or. shflg .and. d(7).ge.0.0d0) then
            cfac = d(7)
            lfac = 1.0d0 - cfac
          else
            cfac = 0.0d0
            lfac = 0.0d0
          endif

          do ll = 1,lint

c           Compute elastic energy density from stress and deformation

            call shp1d(sg(1,ll),xl,shpu,ndm,nel,dx(ll))
            dx(ll) = sg(2,ll)*dx(ll)
            energy = sig(1,ll)*eps(1,ll) + sig(2,ll)*eps(2,ll)
     &             + sig(3,ll)*eps(3,ll) + sig(4,ll)*eps(4,ll)
     &             + sig(5,ll)*eps(3,ll) + sig(5,ll)*eps(4,ll)
     &             + sig(3,ll)*sg(1,ll)*uen

            epl(8) = epl(8) + 0.5d0*energy*dx(ll)

          end do ! ll

c         Compute kinetic energy

          call framtr(d,xl,ndm, le,t)
          call mass3s(s,r,d(7),d,le,nst,ndm,ndf)
          i1 = 0
          do ii = 1,2
            do i = 1,6
              pp(i,ii) = 0.0d0
              j1       = 0
              do jj = 1,2
                do j = 1,6
                  pp(i,ii) = pp(i,ii) + s(i+i1,j+j1)*ul(j,jj,4)
                end do ! j

                epl(7) = epl(7) + 0.5d0*r(i,ii)*ul(i,ii,4)**2*lfac

                j1 = j1 + ndf
              end do ! jj
            end do ! i
            i1 = i1 + ndf
          end do ! ii
          epl(7) = epl(7)
     &           + 0.5d0*(pp(1,1)*ul(1,1,4) + pp(1,2)*ul(1,2,4)
     &                  + pp(2,1)*ul(2,1,4) + pp(2,2)*ul(2,2,4)
     &                  + pp(3,1)*ul(3,1,4) + pp(3,2)*ul(3,2,4)
     &                  + pp(4,1)*ul(4,1,4) + pp(4,2)*ul(4,2,4)
     &                  + pp(5,1)*ul(5,1,4) + pp(5,2)*ul(5,2,4)
     &                  + pp(6,1)*ul(6,1,4) + pp(6,2)*ul(6,2,4))*cfac

        endif

      endif

c     Format statements

2001  format(a1,20a4//5x,'3-D Frame Element Forces'//
     &   '    Elmt  Mat     x-Coor     y-Coor     z-Coor'/
     & 7x,'            Force   1-Torque   1-Moment   2-Moment'/
     & 1x,78('-'))

2002  format(i8,i5,1p,3e11.3/13x,1p,4e11.3/
     &       13x,  1p,6e11.3/13x,1p,4e11.3/1x)

2003  format('  Stress_Layer',1p,5e13.4)
2004  format('  Strain_Layer',1p,5e13.4)

      end

      subroutine bm3trn(s,t,nst,ndf,isw)

      implicit  none

      integer   nst,ndf,isw, i,j,j1,j2,k, nsiz
      real*8    s(nst,nst),ss(6),t(3,3)

      save

c     Transform if not identity

      if(t(1,1)+t(2,2)+t(3,3).lt.2.999999d+00) then

c       Stiffness transformation

        if(isw.eq.1) then

c         Postmultiply local stiffness by transformation array

          j1 = 0
          nsiz = ndf + ndf
          do k = 1,2
            do i = 1,nsiz
              do j = 1,6
                ss(j) = s(i,j+j1)
              end do ! i
              j2 = j1 + 3
              do j = 1,3
                s(i,j+j1) = ss(1)*t(1,j) + ss(2)*t(2,j) + ss(3)*t(3,j)
                s(i,j+j2) = ss(4)*t(1,j) + ss(5)*t(2,j) + ss(6)*t(3,j)
              end do ! j
            end do ! i
            j1 = j1 + ndf
          end do ! k

c         Premultiply result by transpose of transformation array

          j1 = 0
          do k = 1,2
            do i = 1,nsiz
              do j = 1,6
                ss(j) = s(j+j1,i)
              end do ! j
              j2 = j1 + 3
              do j = 1,3
                s(j+j1,i) = t(1,j)*ss(1) + t(2,j)*ss(2) + t(3,j)*ss(3)
                s(j+j2,i) = t(1,j)*ss(4) + t(2,j)*ss(5) + t(3,j)*ss(6)
              end do ! j
            end do ! i
            j1 = j1 + ndf
          end do ! k

        elseif(isw.eq.2) then

c         Premultiply result by transformation array

          do i = 1,2
            do j = 1,ndf
              ss(j) = s(j,i)
            end do ! j
            do j = 1,3
              s(j  ,i) = t(j,1)*ss(1) + t(j,2)*ss(2) + t(j,3)*ss(3)
            end do ! j
            if(ndf.ge.6) then
              do j = 1,3
                s(j+3,i) = t(j,1)*ss(4) + t(j,2)*ss(5) + t(j,3)*ss(6)
              end do ! j
            endif
          end do ! i

        elseif(isw.eq.3) then

c         Premultiply result by transpose of transformation array

          do i = 1,2
            do j = 1,ndf
              ss(j) = s(j,i)
            end do ! j
            do j = 1,3
              s(j  ,i) = t(1,j)*ss(1) + t(2,j)*ss(2) + t(3,j)*ss(3)
            end do ! j
            if(ndf.ge.6) then
              do j = 1,3
                s(j+3,i) = t(1,j)*ss(4) + t(2,j)*ss(5) + t(3,j)*ss(6)
              end do ! j
            endif
          end do ! i

        endif

      endif

      end
