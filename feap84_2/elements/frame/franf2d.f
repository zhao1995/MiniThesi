c$Id:$
      subroutine franf2d(d,ul,xl,s,r,ndf,ndm,nst,isw)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Add 'b1' to compute consistent mass              06/09/2007
c       2. Remove duplicate lumped inertia effects          11/09/2007
c       3. Correct inertia terms for transient analysis     23/07/2008
c       4. Use 'nv' to increment history variables          17/04/2011
c       5. Add u geometric stiffness with tolerance         10/10/2013
c-----[--.----+----.----+----.-----------------------------------------]
c     Purpose: Two dimensional Euler-Bernoulli frame element: Second
c              order theory. Includes one enhanced mode for axial and
c              bending strain match.
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'bdata.h'
      include  'bm2com.h'
      include  'bm2str.h'
      include  'cdata.h'
      include  'cdat1.h'
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
      include  'tdata.h'
      include  'comblk.h'

      logical   noconv, small
      integer   ndf,ndm,nst,isw, nlay,nlob
      integer   i,ii,itmax, j,jj, ll,lint, mm,nn,nh,nhv, nv
      real*8    cs,sn,b1,b2,b3,b4,dva,dvi,xjac, energy,len,hle,utol
      real*8    rhoa,rhoi, ctan1,ctan3, duen,uen,ben,hen
      real*8    cfac,lfac, cfac3,lfac3
      real*8    aa(3,3,4),shpw(4,2,4),shpt(4,2,4),shpu(2,2,4),dx(4)
      real*8    cc(3,3,2),ac(3)
      real*8    sg(2,4), dudx(4), dwdx(4), eps(4), kap(4), gam(4)
      real*8    bmat(2,3,2),baii(4,2),forc(3,5),xx(2),nxi(3)
      real*8    xl(ndm,*),ul(ndf,nen,*)
      real*8    d(*),r(ndf,*),s(nst,nst), gen(3,2), pp(3,2)

      save

      data      itmax / 5     /
      data      utol  / 1.d-10 /

c     Second order (nonlinear) deformation BEAM element.

c     d(1)*d(32)       = EA
c     d(1)*d(33)       = EI
c     d(4)*d(32)       = rho*A
c     d(4)*d(33)       = rho*I

c     Check for small deformation

      small = d(18).gt.0.0d0

c     Compute element length and direction cosines

      if(isw.ge.2) then
        cs  = xl(1,nel) - xl(1,1)
        sn  = xl(2,nel) - xl(2,1)
        len = sqrt(cs*cs + sn*sn)
        cs  = cs/len
        sn  = sn/len
        hle = 0.5d0*len

        nlay = nint(d(101))
        if(nlay.eq.0) then
          nout = 0
          nh   = 0
          nhv  = 0
          rhoa = d(4)*d(32)
          rhoi = d(4)*d(33)
        else
          nlob = nint(d(102))
          nout = (nlay - 1)*(nlob - 1) + 1
          nh   = nint(d(15))
          nhv  = nh*nout
          call int1dl(nlob,sl)
          call bm2rho(nlay,d,d(103), rhoa,rhoi)
        endif
      endif

c     Read data

      if(isw.eq.1) then

c       Increment history storage if necessary

        nh1 = nh1 + 1   ! Enhanced strain parameter

c     Compute mass array

      elseif(isw.eq.5 .and. imtyp.eq.1) then

c       Inertia factors

        cfac = d(7)
        lfac = 1.d0 - cfac

        if(nint(d(182)).gt.0) then
          lint = nel
          call int1dn(lint, sg)
        else
          lint = nel + 2
          call int1d(lint, sg)
        endif

c       Compute lumped and consistent mass matrix

        call bm2trn (xl,cs,sn,ndm*nel,ndm,2)

        do ll = 1,lint

c         Compute shape functions

          call shp1d (sg(1,ll),xl,shpu(1,1,ll),ndm,nel,xjac)
          call shp1dh(sg(1,ll),len,shpw(1,1,ll),shpt(1,1,ll))

          dva  = sg(2,ll)*xjac*rhoa
          dvi  = sg(2,ll)*xjac*rhoi*d(8)

c         For each node j compute db = rho*shape*dv

          jj = 0
          do j = 1,nel
            b1 = shpu(2,j,ll)*dva*lfac
            b2 = shpu(2,j,ll)*dvi*lfac

c           Compute a lumped mass

            r(1,j) = r(1,j) + b1
            r(2,j) = r(2,j) + b1
            r(3,j) = r(3,j) + b2

            b1 = shpu(2,j,ll)*dva*cfac
            b2 = shpw(4,j,ll)*dva*cfac
            b3 = shpt(4,j,ll)*dva*cfac
            ii = 0
            do i = 1,nel
              s(ii+1,jj+1) = s(ii+1,jj+1) + shpu(2,i,ll)*b1
              s(ii+2,jj+2) = s(ii+2,jj+2) + shpw(4,i,ll)*b2
              s(ii+2,jj+3) = s(ii+2,jj+3) + shpw(4,i,ll)*b3
              s(ii+3,jj+2) = s(ii+3,jj+2) + shpt(4,i,ll)*b2
              s(ii+3,jj+3) = s(ii+3,jj+3) + shpt(4,i,ll)*b3
              ii = ii + ndf
            end do ! i
            jj = jj + ndf
          end do ! j
        end do ! ll

c       Place in consistent mass

        jj = 0
        do j = 1,nel
          s(jj+1,jj+1) = s(jj+1,jj+1) + r(1,j)
          s(jj+2,jj+2) = s(jj+2,jj+2) + r(2,j)
          s(jj+3,jj+3) = s(jj+3,jj+3) + r(3,j)
          jj = jj + ndf
        end do ! j

c       Transform to global frame

        call bm2trn (s,cs,sn,nst,ndf,1)

      elseif(isw.ne.12) then

c       Quadrature terms

        if(nint(d(182)).gt.0) then
          lint = nel
          call int1dn(lint, sg)
        else
          lint = nel + 1
          call int1d(lint, sg)
        endif

c       Transform to local coordinates

        call bm2trn (ul(1,1,1),cs,sn,ndf*nel,ndf,2)
        call bm2trn (ul(1,1,2),cs,sn,ndf*nel,ndf,2)
        call bm2trn (ul(1,1,4),cs,sn,ndf*nel,ndf,2)
        call bm2trn (ul(1,1,5),cs,sn,ndf*nel,ndf,2)
        call bm2trn (xl       ,cs,sn,ndm*nel,ndm,2)

        do ll = 1,lint

c         Shape functions

          call shp1dh(sg(1,ll),len,shpw(1,1,ll),shpt(1,1,ll))
          shpu(1,1,ll) = -1.d0/len
          shpu(1,2,ll) = -shpu(1,1,ll)
          shpu(2,1,ll) = 0.5d0 - 0.5d0*sg(1,ll)
          shpu(2,2,ll) = 0.5d0 + 0.5d0*sg(1,ll)

          dx(ll)      = sg(2,ll)*hle

c         Form displacement derivatives from nodal displacements

          dwdx(ll) = shpw(1,1,ll)*ul(2,1,1) + shpw(1,2,ll)*ul(2,2,1)
     &             + shpt(1,1,ll)*ul(3,1,1) + shpt(1,2,ll)*ul(3,2,1)

          dudx(ll) = 0.5d0*(ul(1,2,1) - ul(1,1,1))/hle

          if(small) then
            eps(ll)  = dudx(ll)
          else
            eps(ll)  = dudx(ll) + 0.5d0*(dudx(ll)*dudx(ll)
     &                                 + dwdx(ll)*dwdx(ll))
          endif

          kap(ll)  = shpw(2,1,ll)*ul(2,1,1) + shpw(2,2,ll)*ul(2,2,1)
     &             + shpt(2,1,ll)*ul(3,1,1) + shpt(2,2,ll)*ul(3,2,1)

          gam(ll)  = shpw(3,1,ll)*ul(2,1,1) + shpw(3,2,ll)*ul(2,2,1)
     &             + shpt(3,1,ll)*ul(3,1,1) + shpt(3,2,ll)*ul(3,2,1)
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

            nv  = 1
            do ll = 1,lint

              defa(1,1) = eps(ll) + sg(1,ll)*uen
              defa(2,1) = 0.0d0
              defa(3,1) = kap(ll)
              call bm2con (d,hr(nh1+nv),hr(nh2+nv),nh,cc,strs,
     &                     defa,defa,isw)

              do j = 1,3
                forc(j,ll) = strs(j,1)
                do i = 1,3
                  aa(i,j,ll) = cc(i,j,1)*dx(ll)*ctan(1)
                end do ! i
              end do ! j

              hen = hen + sg(1,ll)*aa(1,1,ll)*sg(1,ll)
              ben = ben - sg(1,ll)*forc(1,ll)*dx(ll)

              nv = nv + nhv
            end do ! ll

            hen  = 1.d0/ hen
            duen = ben * hen
            uen  = uen + duen

            if(abs(duen).le.utol*abs(uen) .or. ben.eq.0.0d0) then
              noconv = .false.
            elseif(ii.gt.itmax) then
              noconv = .false.
c              write(*,*) 'WARNING - No convergence in FRANF2D'
c              write(*,*) '          ISW =',isw,duen,uen,hen
            endif

          end do ! while

c         Save enhance mode parameter

          hr(nh2) = uen

        else

          hen = 0.0d0
          nv  = 1
          do ll = 1,lint

            defa(1,1) = eps(ll)
            defa(2,1) = 0.0d0
            defa(3,1) = kap(ll)
            call bm2con (d,hr(nh1+nv),hr(nh2+nv),nh,cc,strs,
     &                   defa,defa,isw)
            do j = 1,3
              do i = 1,3
                forc(j,ll) = strs(j,1)
                aa(i,j,ll) = cc(i,j,1)*dx(ll)*ctan(1)
              end do ! i
            end do ! j

            nv = nv + nhv
          end do ! ll

        end if ! etype

c       Stiffness and residual computation

        if(isw.eq.3 .or. isw.eq.6) then

c         Zero enhanced coupling array

          do i = 1,3
            gen(i,1)   = 0.0d0
            gen(i,2)   = 0.0d0
          end do ! i

c         Final tangent form

          do ll = 1,lint

            do j = 1,3
              forc(j,ll) = forc(j,ll)*dx(ll)
            end do ! j

c           Compute strain-displacement matrices for two nodes

            if(small) then
              do i = 1,2
                bmat(1,1,i) =  shpu(1,i,ll)
                bmat(2,1,i) =  0.0d0
                bmat(1,2,i) =  0.0d0
                bmat(2,2,i) =  shpw(2,i,ll)
                bmat(1,3,i) =  0.0d0
                bmat(2,3,i) =  shpt(2,i,ll)
              end do ! i
            else
              do i = 1,2
                bmat(1,1,i) =  shpu(1,i,ll)*(1.d0 + dudx(ll))
                bmat(2,1,i) =  0.0d0
                bmat(1,2,i) =  shpw(1,i,ll)*dwdx(ll)
                bmat(2,2,i) =  shpw(2,i,ll)
                bmat(1,3,i) =  shpt(1,i,ll)*dwdx(ll)
                bmat(2,3,i) =  shpt(2,i,ll)
              end do ! i
            endif

c           Mechanical tangent terms

            mm = 0
            do ii = 1,nel

              do i = 1,3

c               B^T * AA

                baii(i,1) = (bmat(1,i,ii)*aa(1,1,ll) +
     &                       bmat(2,i,ii)*aa(3,1,ll))
                baii(i,2) = (bmat(1,i,ii)*aa(1,3,ll) +
     &                       bmat(2,i,ii)*aa(3,3,ll))

c               Residual

                r(i,ii) = r(i,ii) - bmat(1,i,ii)*forc(1,ll)
     &                            - bmat(2,i,ii)*forc(3,ll)

c               Enhanced stiffness

                gen(i,ii) = gen(i,ii) + baii(i,1)*sg(1,ll)

              end do ! i

c             Tangent

              if(isw.eq.3) then

                nxi(1) = shpu(1,ii,ll)*forc(1,ll)*ctan(1)
                nxi(2) = shpw(1,ii,ll)*forc(1,ll)*ctan(1)
                nxi(3) = shpt(1,ii,ll)*forc(1,ll)*ctan(1)
                nn = 0
                do jj = 1,nel

c                 Material part

                  do j = 1,3
                    do i = 1,3
                      s(mm+i,nn+j) = s(mm+i,nn+j)
     &                             + baii(i,1)*bmat(1,j,jj)
     &                             + baii(i,2)*bmat(2,j,jj)
                    end do ! i
                  end do ! j

c                 Geometric part

                  if(gflag. and. .not.small) then
                    s(mm+1,nn+1) = s(mm+1,nn+1) + nxi(1)*shpu(1,jj,ll)
                    s(mm+2,nn+2) = s(mm+2,nn+2) + nxi(2)*shpw(1,jj,ll)
                    s(mm+2,nn+3) = s(mm+2,nn+3) + nxi(2)*shpt(1,jj,ll)
                    s(mm+3,nn+2) = s(mm+3,nn+2) + nxi(3)*shpw(1,jj,ll)
                    s(mm+3,nn+3) = s(mm+3,nn+3) + nxi(3)*shpt(1,jj,ll)
                  endif

                  nn = nn + ndf
                end do ! jj
              endif
              mm = mm + ndf
            end do ! ii
          end do ! ll

c         Compute lumped and consistent inertia arrays

          if(d(7).ge.0.0 .and. (ndfo(1).gt.0 .or. shflg)) then
            ctan3 = ctan(3) + d(77)*ctan(2)
            cfac  = d(7)
            cfac3 = cfac*ctan3
            lfac  = 1.0d0 - cfac
            lfac3 = lfac*ctan3

            if(nint(d(182)).gt.0) then
              lint = nel
              call int1dn(lint, sg)
            else
              lint = nel + 2
              call int1d(lint, sg)
            endif

c           Compute lumped and consistent mass matrix

            do ll = 1,lint

c             Compute shape functions

              call shp1d (sg(1,ll),xl,shpu(1,1,ll),ndm,nel,xjac)
              call shp1dh(sg(1,ll),len,shpw(1,1,ll),shpt(1,1,ll))

              dva  = sg(2,ll)*xjac*rhoa
              dvi  = sg(2,ll)*xjac*rhoi*d(8)

c             Compute acceleration

              ac(1) = 0.0d0
              ac(2) = 0.0d0
              do j = 1,nel
                ac(1) = ac(1) + shpu(2,j,ll)*ul(1,j,5)  ! axial accel.
                ac(2) = ac(2) + shpw(4,j,ll)*ul(2,j,5)  ! transv. accel.
     &                        + shpt(4,j,ll)*ul(3,j,5)
              end do
              ac(1) = ac(1)*cfac
              ac(2) = ac(2)*cfac

c             For each node j compute db = rho*shape*dv

              jj = 0
              do j = 1,nel
                b1 = shpu(2,j,ll)*dva
                b2 = shpu(2,j,ll)*dvi
                b3 = shpw(4,j,ll)*dva
                b4 = shpt(4,j,ll)*dva
                r(1,j) = r(1,j)
     &                  - b1*(ac(1) + ul(1,j,5)*lfac + d(77)*ul(1,j,4))
                r(2,j) = r(2,j) - b3*ac(2)
     &                  - b1*(ul(2,j,5)*lfac + d(77)*ul(2,j,4))
                r(3,j) = r(3,j) - b4*ac(2)
     &                  - b2*(ul(3,j,5)*lfac + d(77)*ul(3,j,4))

c               Compute lumped part

                s(jj+1,jj+1) = s(jj+1,jj+1) + b1*lfac3
                s(jj+2,jj+2) = s(jj+2,jj+2) + b1*lfac3
                s(jj+3,jj+3) = s(jj+3,jj+3) + b2*lfac3

c               Compute consistent part

                b1 = shpu(2,j,ll)*dva*cfac3
                b2 = shpw(4,j,ll)*dva*cfac3
                b3 = shpt(4,j,ll)*dva*cfac3
                ii = 0
                do i = 1,nel
                  s(ii+1,jj+1) = s(ii+1,jj+1) + shpu(2,i,ll)*b1
                  s(ii+2,jj+2) = s(ii+2,jj+2) + shpw(4,i,ll)*b2
                  s(ii+2,jj+3) = s(ii+2,jj+3) + shpw(4,i,ll)*b3
                  s(ii+3,jj+2) = s(ii+3,jj+2) + shpt(4,i,ll)*b2
                  s(ii+3,jj+3) = s(ii+3,jj+3) + shpt(4,i,ll)*b3
                  ii = ii + ndf
                end do ! i
                jj = jj + ndf
              end do ! j
            end do ! ll
          endif

c         Transform stiffness and residual to global coordinates

          if(isw.eq.3) then

c           Static condensation

            nn = 0
            do jj = 1,nel
              do j = 1,3
                duen = gen(j,jj)*hen
                mm = 0
                do ii = 1,nel
                  do i = 1,3
                    s(i+mm,j+nn) = s(i+mm,j+nn) - gen(i,ii)*duen
                  end do ! i
                  mm = mm + ndf
                end do ! ii
              end do ! j
              nn = nn + ndf
            end do ! jj

            call bm2trn (s,cs,sn,nst,ndf,1)
          endif
          call bm2trn ( r,cs,-sn,nst,ndf,2)

c         Set body loading factors and follower forces

          call fbody2d(d,xl,ul, r,s, ndm,ndf,nst, isw)

c       Output forces

        elseif(isw.eq.4 .or. isw.eq.8) then

c         Member forces

          do i = 1,3
            pp(i,1) = 0.0d0
            pp(i,2) = 0.0d0
          end do ! i
          do ll = 1,lint

            defa(1,1)  = eps(ll) + sg(1,ll)*uen
            defa(2,1)  = gam(ll)
            defa(3,1)  = kap(ll)

c           Stress output

            if(isw.eq.4) then

              do i = 1,ndm
                xx(i) = 0.
                do ii = 1,nel
                  xx(i) = xx(i) + xl(i,ii)*shpu(2,ii,ll)
                end do ! ii
              end do ! i
              mct = mct - 3
              if (mct.le.0) then
                write(iow,2001) o,head,ttim
                if(ior.lt.0) write(*,2001) o,head,ttim
                mct = 50
              endif
              write(iow,2002) n,ma,(xx(i),i=1,2),
     &                        (forc(i,ll),i=1,3),(defa(i,1),i=1,3)
              if(nout.gt.0) then
                write(iow,2003) (siglr(i),i=1,nout)
                write(iow,2004) (epslr(i),i=1,nout)
              endif
              if(ior.lt.0) then
                write(*,2002) n,ma,(xx(i),i=1,2),
     &                        (forc(i,ll),i=1,3),(defa(i,1),i=1,3)
                if(nout.gt.0) then
                  write(*,2003) (siglr(i),i=1,nout)
                  write(*,2004) (epslr(i),i=1,nout)
                endif
              endif

c           Stress projections save

            else

c             Compute strain-displacement matrices for two nodes

              if(small) then
                do i = 1,2
                  bmat(1,1,i) =  shpu(1,i,ll)
                  bmat(2,1,i) =  0.0d0
                  bmat(1,2,i) =  0.0d0
                  bmat(2,2,i) =  shpw(2,i,ll)
                  bmat(1,3,i) =  0.0d0
                  bmat(2,3,i) =  shpt(2,i,ll)
                end do ! i
              else
                do i = 1,2
                  bmat(1,1,i) =  shpu(1,i,ll)*(1.d0 + dudx(ll))
                  bmat(2,1,i) =  0.0d0
                  bmat(1,2,i) =  shpw(1,i,ll)*dwdx(ll)
                  bmat(2,2,i) =  shpw(2,i,ll)
                  bmat(1,3,i) =  shpt(1,i,ll)*dwdx(ll)
                  bmat(2,3,i) =  shpt(2,i,ll)
                end do ! i
              endif

c             End forces

              do ii = 1,nel

                do i = 1,3
                  pp(i,ii) = pp(i,ii) - (bmat(1,i,ii)*forc(1,ll)
     &                                +  bmat(2,i,ii)*forc(3,ll))*dx(ll)
                end do ! i

              end do ! ii
            end if
          end do ! ll

c         Projection on end notes (uses reactions)

          if(isw.eq.8) then
            do i = 1,3
              pp(i,2) = -pp(i,2)
            end do ! i
            call frcn2d(pp,r,s)
          endif

c       Geometric stiffness computation

        elseif(isw.eq.5 .and. imtyp.eq.2) then

          do ll = 1,lint

            ctan1      = dx(ll)*ctan(1)
            forc(1,ll) = forc(1,ll)*ctan1

            mm = 0
            do ii = 1,nel

              nxi(1) = shpu(1,ii,ll)*forc(1,ll)*1.e-2 ! Add to avoid 0
              nxi(2) = shpw(1,ii,ll)*forc(1,ll)
              nxi(3) = shpt(1,ii,ll)*forc(1,ll)
              nn = 0
              do jj = 1,nel

                s(mm+1,nn+1) = s(mm+1,nn+1) - nxi(1)*shpu(1,jj,ll)
                s(mm+2,nn+2) = s(mm+2,nn+2) - nxi(2)*shpw(1,jj,ll)
                s(mm+2,nn+3) = s(mm+2,nn+3) - nxi(2)*shpt(1,jj,ll)
                s(mm+3,nn+2) = s(mm+3,nn+2) - nxi(3)*shpw(1,jj,ll)
                s(mm+3,nn+3) = s(mm+3,nn+3) - nxi(3)*shpt(1,jj,ll)

                nn = nn + ndf
              end do ! jj
              mm = mm + ndf
            end do ! ii
          end do ! ll

c         Transform to global coordinates

          call bm2trn (s,cs,sn,nst,ndf,1)

c       Compute energy

        elseif(isw.eq.13) then

          dva = hle*rhoa
          dvi  =hle*rhoi*d(8)

c         Compute internal energy

          do ll = 1,lint

c           Compute energy density from stress and deformation

            call shp1d(sg(1,ll),xl,shpu,ndm,nel,xjac)
            dx(ll) = sg(2,ll)*xjac

            energy = forc(1,ll)*(eps(ll) + sg(1,ll)*uen)
     &             + forc(3,ll)* kap(ll)

c           Accumulate energy

            epl(8) = epl(8) + 0.5d0*energy*dx(ll)

          end do ! ll

c         Compute kinetic energy for lumped mass

          epl(7) = epl(7) + 0.5d0*dva*(ul(1,1,4)**2 + ul(1,2,4)**2
     &                               + ul(2,1,4)**2 + ul(2,2,4)**2)
     &                    + 0.5d0*dvi*(ul(3,1,4)**2 + ul(3,2,4)**2)

c       Initialize history variables

        elseif(isw.eq.14) then

          call bm2init(d,hr(nh1+1),hr(nh2+1),nh)

        endif

      endif

c     Formats

2001  format(a1,20a4/5x,'time',e13.5,5x,' element forces '//
     &  43x,'*********  FORCE / STRAIN  *********'/
     &   3x,'element  material',
     &  3x,'1-coord',3x,'2-coord',6x,'n-dir',8x,'s-dir',8x,'m-dir'/)

2002  format(2i10,1p,2f10.3,1p,3e13.4/40x,1p,3e13.4)

2003  format('  Stress_Layer',1p,5e13.4)
2004  format('  Strain_Layer',1p,5e13.4)

      end
