c$Id:$
      subroutine frams2e(d,ul,xl,s,r,ndf,ndm,nst,isw)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Use constants from 'pconstant.h'                 14/11/2006
c       2. Add geometric stiffness for eigen solutions      09/06/2009
c       3. Set d(21) to d(1) in isw = 1 to get E            14/01/2010
c-----[--.----+----.----+----.-----------------------------------------]
c     Purpose: Two dimensional small displacement frame element
c              Enhanced formulation
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
      include  'pconstant.h'
      include  'part0.h'
      include  'prstrs.h'
      include  'ptdat6.h'
      include  'rdata.h'
      include  'tdata.h'
      include  'comblk.h'

      logical   enflg
      integer   ndf,ndm,nst,isw, nlay,nlob
      integer   i,ii,istrt, j,jj, ll,lint, nn,nhv,niter
      integer   is(12)
      real*8    cs,sn,a1,a2,a3,le,lr,dx,dva,dvi,cfac,dfac,xjac
      real*8    energy, rhoa, rhoi, ctan1, ctan3, shpen,hdet, itol
      real*8    xn,ux,wx
      real*8    sg(2,4),xx(2),dma(3), shp(2,3)
      real*8    xl(ndm,*),ul(ndf,nen,*)
      real*8    d(*),r(ndf,*),s(nst,nst),pp(3,3),alp(2),dalp(2)
      real*8    eps0(3),eps(3,2),bc(2,3),bb(2),hh(2,2),gg(6,2),ghi(2)
      real*8    cc(3,3,2),aa(3,3),ae(2,3),sig(3,2),forc(3),mass(6,6)

      save

      data      itol / 1.d-8 /

c     Compute element length and direction cosines

      if(isw.ge.2) then
        cs = xl(1,nel) - xl(1,1)
        sn = xl(2,nel) - xl(2,1)
        le = sqrt(cs*cs + sn*sn)
        lr = 1.d0/le
        cs = cs*lr
        sn = sn*lr

        nlay = nint(d(101))
        if(nlay.eq.0) then
          nout = 0
          nhv  = 0
          rhoa = d(4)*d(32)
          rhoi = d(4)*d(33)
        else
          nlob = nint(d(102))
          nout = (nlay - 1)*(nlob - 1) + 1
          nhv  = nint(d(15))*nout
          call int1dl(nlob,sl)
          call bm2rho(nlay,d,d(103), rhoa,rhoi)
        endif
      endif

c     Read data

      if(isw.eq.1) then

c       Set effective modulus to E

        d(21) = d(1)

c       Increment history storage if necessary

        nh1 = nh1 + 2  ! Add for enhanced parameters

        do i = 1,3
          is(i  ) = i
          is(i+3) = i + ndf
        end do ! i

c     Compute mass or geometric stiffness array

      elseif(isw.eq.5) then

c       Mass matrix

        if(imtyp.eq.1) then

c         Lumped mass factors

          dma(1) = 0.5d0*rhoa*le
          dma(2) = dma(1)
          dma(3) = 0.5d0*rhoi*d(8)*le

          do i = 1,3
            r(i,1) = dma(i)
            r(i,2) = dma(i)
          end do !i

c         Set mass factors

          cfac      = d(7)
          dfac      = 1.d0 - d(7)
          mass(1,1) =  rhoa*le*(dfac*0.5d0 + cfac/3.d0)
          mass(1,4) =  mass(1,1)*0.5d0*cfac
          mass(2,2) =  mass(1,1)
          mass(2,3) =  rhoa*le*le/24.d0*cfac
          mass(2,5) =  mass(1,4)
          mass(2,6) = -mass(2,3)
          mass(3,2) =  mass(2,3)
          mass(3,3) = (rhoa*le*le*le/120.d0 + rhoi*le/3.0d0)*cfac
          mass(3,5) =  mass(2,3)
          mass(3,6) = -mass(3,3)            - rhoi*le/6.0d0*cfac
          mass(3,3) =  mass(3,3)            + rhoi*le*0.5d0*dfac
          mass(4,1) =  mass(1,4)
          mass(4,4) =  mass(1,1)
          mass(5,2) =  mass(2,5)
          mass(5,3) =  mass(3,5)
          mass(5,5) =  mass(2,2)
          mass(5,6) = -mass(2,3)
          mass(6,2) =  mass(2,6)
          mass(6,3) =  mass(3,6)
          mass(6,5) =  mass(5,6)
          mass(6,6) =  mass(3,3)

          do j = 1,6
            do i = 1,6
              s(is(i),is(j)) = s(is(i),is(j)) + mass(i,j)
            end do ! i
          end do ! j

c       Geometric stiffness

        elseif(imtyp.eq.2) then

          ux  = (cs*(ul(1,2,1)-ul(1,1,1))
     &        +  sn*(ul(2,2,1)-ul(2,1,1)))/le
          wx  = (cs*(ul(2,2,1)-ul(2,1,1))
     &        -  sn*(ul(1,2,1)-ul(1,1,1)))/le
          xn  = -d(21)*d(32)*(ux + 0.5d0*(ux*ux + wx*wx))*ctan(1)

          s(1    ,1    ) =  xn/le
          s(2    ,2    ) =  s(1,1)

          s(1    ,1+ndf) = -s(1,1)
          s(2    ,2+ndf) = -s(2,2)

          s(1+ndf,1    ) = -s(1,1)
          s(2+ndf,2    ) = -s(2,2)

          s(1+ndf,1+ndf) =  s(1,1)
          s(2+ndf,2+ndf) =  s(2,2)

        endif

      elseif(isw.eq.12) then

        lint = nel
        if(nint(d(182)).gt.0) then
          call int1dn(lint, sg)
        else
          call int1d(lint, sg)
        endif
        call bm2trn (ul(1,1,1),cs,sn,ndf*nel,ndf,2)
        call bm2trn (ul(1,1,2),cs,sn,ndf*nel,ndf,2)
        call bm2trn (xl       ,cs,sn,ndm*nel,ndm,2)
        nn = 2
        do ll = 1,lint
          call b2mode(d,sig,eps,cc,energy,hr(nh1+nn),hr(nh2+nn),isw,ll)
          nn = nn + nhv
        end do ! ll

c     Compute energy

      elseif(isw.eq.13) then

        if(nint(d(182)).gt.0) then
          lint = nel
          call int1dn(lint, sg)
        else
          lint = nel - 1
          call int1d(lint, sg)
        endif
        call bm2trn (ul(1,1,1),cs,sn,ndf*nel,ndf,2)
        call bm2trn (ul(1,1,2),cs,sn,ndf*nel,ndf,2)
        call bm2trn (ul(1,1,4),cs,sn,ndf*nel,ndf,2)
        call bm2trn (xl       ,cs,sn,ndm*nel,ndm,2)
        dva = 0.5d0*rhoa*le
        dvi  =0.5d0*rhoi*d(8)*le

c       Compute internal energy

        nn = 2
        do ll = 1,lint

c         Compute energy density from stress and deformation

          call shp1d(sg(1,ll),xl,shp,ndm,nel,xjac)
          dx = sg(2,ll)*xjac
          call b2mode(d,sig,eps,cc,energy,hr(nh1+nn),hr(nh2+nn),isw,ll)

c         Accumulate energy

          epl(8) = epl(8) + 0.5d0*energy*dx

          nn = nn + nhv
        end do ! ll

c       Compute kinetic energy for lumped mass

        epl(7) = epl(7) + 0.5d0*dva*(ul(1,1,4)**2 + ul(1,2,4)**2
     &                             + ul(2,1,4)**2 + ul(2,2,4)**2)
     &                  + 0.5d0*dvi*(ul(3,1,4)**2 + ul(3,2,4)**2)

c     Initialize history variables

      elseif(isw.eq.14) then

        istrt  = nint(d(84))
        call modl1d(d,le,xx,hr(nh1),hr(nh2),nint(d(15)),1,istrt,
     &              forc,aa,isw)

c     Residual, tangent, stress, projection options

      elseif(isw.eq.3 .or. isw.eq.4 .or. isw.eq.6 .or. isw.eq.8) then

c       Transform nodal parameters to local frame

        call bm2trn (ul(1,1,1),cs,sn,ndf*nel,ndf,2)
        call bm2trn (ul(1,1,2),cs,sn,ndf*nel,ndf,2)
        call bm2trn (ul(1,1,5),cs,sn,ndf*nel,ndf,2)
        call bm2trn (xl       ,cs,sn,ndm*nel,ndm,2)

c       Extract enhanced parameters

        alp(1) = hr(nh2)
        alp(2) = hr(nh2+1)

c       Set quadrature

        lint   = nel
        if(nint(d(182)).gt.0) then
          call int1dn(lint, sg)
        else
          call int1d(lint, sg)
        endif

c       Compute enhanced terms

        eps0(1) = (ul(1,2,1) - ul(1,1,1))*lr     ! Axial
        eps0(2) = (ul(2,2,1) - ul(2,1,1))*lr     ! Shear
     &          - (ul(3,2,1) + ul(3,1,1))*0.5d0
        eps0(3) = (ul(3,2,1) - ul(3,1,1))*lr     ! Bending

        niter = 0
        enflg = .true.
        do while (enflg)

          niter   = niter + 1

          hh(1,1) = 0.0d0
          hh(1,2) = 0.0d0
          hh(2,1) = 0.0d0
          hh(2,2) = 0.0d0
          bb(1)   = 0.0d0
          bb(2)   = 0.0d0

          do j = 1,3
            do i = 1,3
              aa(i,j) = 0.0d0
            end do ! i
            ae(1,j) = 0.0d0
            ae(2,j) = 0.0d0
            forc(j) = 0.0d0
          end do ! j


          nn = 2
          do ll = 1,lint

            dx     = sg(2,ll)*le*0.5d0
            shpen  = 4.d0*sg(1,ll)*lr
            eps(1,1) = eps0(1) + shpen*alp(1)
            eps(2,1) = eps0(2) + two3 *alp(2)
            eps(3,1) = eps0(3) + shpen*alp(2)

            call b2mode(d,sig,eps,cc,energy,hr(nh1+nn),hr(nh2+nn),
     &                  isw,ll)

c           Multiply moduli by solution parameter: ctan(1)

c           ctan1 = ctan(1) + d(78)*ctan(2)
            ctan1 = ctan(1)*dx
            do jj = 1,3
              do ii = 1,3
                cc(ii,jj,1) = cc(ii,jj,1)*ctan1
              end do ! ii
              sig(jj,1) = sig(jj,1)*dx
            end do ! jj

c           Compute enhanced residual and tangent

            bb(1) = bb(1) - shpen*sig(1,1)
            bb(2) = bb(2) -  two3*sig(2,1) - shpen*sig(3,1)
            do i = 1,3
              bc(1,i) = shpen*cc(1,i,1)
              bc(2,i) = two3*cc(2,i,1) + shpen*cc(3,i,1)
            end do ! i
            hh(1,1) = hh(1,1) + bc(1,1)*shpen
            hh(1,2) = hh(1,2) + bc(1,2)*two3 + bc(1,3)*shpen
            hh(2,1) = hh(2,1) + bc(2,1)*shpen
            hh(2,2) = hh(2,2) + bc(2,2)*two3 + bc(2,3)*shpen

c           Accumulate enhanced integrals for stiffness

            do j = 1,3
              do i = 1,3
                aa(i,j) = aa(i,j) + cc(i,j,1)
              end do ! i
              ae(1,j) = ae(1,j) + bc(1,j)
              ae(2,j) = ae(2,j) + bc(2,j)
              forc(j) = forc(j) + sig(j,1)
            end do ! j

c           Increment history counter

            nn = nn + nhv

          end do ! ll

c         Invert enhanced stiffness and compute incremental alpha

          hdet    =  1.d0/(hh(1,1)*hh(2,2) - hh(1,2)*hh(1,2))
          dx      =  hh(1,1)*hdet
          hh(1,1) =  hh(2,2)*hdet
          hh(1,2) = -hh(1,2)*hdet
          hh(2,1) = -hh(2,1)*hdet
          hh(2,2) =  dx
          dalp(1) =  hh(1,1)*bb(1) + hh(1,2)*bb(2)
          dalp(2) =  hh(2,1)*bb(1) + hh(2,2)*bb(2)
          alp(1)  =  alp(1) + dalp(1)
          alp(2)  =  alp(2) + dalp(2)

c         Check convergence

          if(max(abs(dalp(1)),abs(dalp(2))) .le.
     &       max(abs( alp(1)),abs( alp(2)))*itol) then
            enflg = .false.
          elseif(niter.gt.3) then
            enflg = .false.
          endif
        end do ! while

c       Save history variables

        hr(nh2  ) = alp(1)
        hr(nh2+1) = alp(2)

c       Residual

        do i = 1,3
          r(i,1) = r(i,1) + forc(i)*lr
          r(i,2) = r(i,2) - forc(i)*lr
        end do ! i
        r(3,1) = r(3,1) + 0.5d0*forc(2)
        r(3,2) = r(3,2) + 0.5d0*forc(2)

c       Stiffness for main terms

        if(isw.eq.3) then
          a1 = lr*lr
          a2 = lr*0.5d0
          a3 = 0.25d0*aa(2,2)
          do j = 1,3
            do i = 1,3
              s(i    ,j    ) = s(i    ,j    ) + a1*aa(i,j)
              s(i+ndf,j    ) = s(i+ndf,j    ) - a1*aa(i,j)
              s(i    ,j+ndf) = s(i    ,j+ndf) - a1*aa(i,j)
              s(i+ndf,j+ndf) = s(i+ndf,j+ndf) + a1*aa(i,j)
            end do ! i

            s(j    ,3    ) = s(j    ,3    ) +a2*aa(j,2)
            s(j    ,3+ndf) = s(j    ,3+ndf) +a2*aa(j,2)
            s(j+ndf,3    ) = s(j+ndf,3    ) -a2*aa(j,2)
            s(j+ndf,3+ndf) = s(j+ndf,3+ndf) -a2*aa(j,2)

            s(3    ,j    ) = s(3    ,j    ) +a2*aa(2,j)
            s(3+ndf,j    ) = s(3+ndf,j    ) +a2*aa(2,j)
            s(3    ,j+ndf) = s(3    ,j+ndf) -a2*aa(2,j)
            s(3+ndf,j+ndf) = s(3+ndf,j+ndf) -a2*aa(2,j)
          end do ! j
          s(3    ,3    ) = s(3    ,3    ) + a3
          s(3    ,3+ndf) = s(3    ,3+ndf) + a3
          s(3+ndf,3    ) = s(3+ndf,3    ) + a3
          s(3+ndf,3+ndf) = s(3+ndf,3+ndf) + a3

c         Coupling and static condensation

          do i = 1,3
            gg(i  ,1) = -lr*ae(i,1)
            gg(i  ,2) = -lr*ae(i,2)
            gg(i+3,1) =  lr*ae(i,1)
            gg(i+3,2) =  lr*ae(i,2)
          end do ! i
          do i = 1,2
            gg(3,i) = gg(3,i) - 0.5d0*ae(2,i)
            gg(6,i) = gg(6,i) - 0.5d0*ae(2,i)
          end do ! i

          do i = 1,6
            ghi(1) = gg(i,1)*hh(1,1) + gg(i,2)*hh(2,1)
            ghi(2) = gg(i,1)*hh(1,2) + gg(i,2)*hh(2,2)
            do j = 1,6
              s(is(i),is(j)) = s(is(i),is(j)) - ghi(1)*gg(j,1)
     &                                        - ghi(2)*gg(j,2)
            end do ! j
          end do ! i
        endif

c       Lumped and consistent inertia contributions

        if((ndfo(1).gt.0 .or. shflg)  .and.
     &          isw.eq.3 .or. isw.eq.6) then
          cfac = d(7)
          dfac = 1.d0 - d(7)

          dma(1) = 0.5d0*rhoa*le
          dma(2) = dma(1)
          dma(3) = 0.5d0*rhoi*d(8)*le

c         Set mass factors

          mass(1,1) =  rhoa*le*(dfac*0.5d0 + cfac/3.d0)
          mass(1,4) =  mass(1,1)*0.5d0*cfac
          mass(2,2) =  mass(1,1)
          mass(2,3) =  rhoa*le*le/24.d0*cfac
          mass(2,5) =  mass(1,4)
          mass(2,6) = -mass(2,3)
          mass(3,2) =  mass(2,3)
          mass(3,3) = (rhoa*le*le*le/120.d0 + rhoi*le/3.0d0)*cfac
          mass(3,5) =  mass(2,3)
          mass(3,6) = -mass(3,3)            - rhoi*le/6.0d0*cfac
          mass(3,3) =  mass(3,3)            + rhoi*le*0.5d0*dfac
          mass(4,1) =  mass(1,4)
          mass(4,4) =  mass(1,1)
          mass(5,2) =  mass(2,5)
          mass(5,3) =  mass(3,5)
          mass(5,5) =  mass(2,2)
          mass(5,6) = -mass(2,3)
          mass(6,2) =  mass(2,6)
          mass(6,3) =  mass(3,6)
          mass(6,5) =  mass(5,6)
          mass(6,6) =  mass(3,3)

          do i = 1,3
            do j = 1,3
              r(i,1) = r(i,1)- mass(i  ,j  )*(ul(i,1,5)+d(77)*ul(i,1,4))
     &                       - mass(i  ,j+3)*(ul(i,2,5)+d(77)*ul(i,2,4))
              r(i,2) = r(i,2)- mass(i+3,j  )*(ul(i,1,5)+d(77)*ul(i,1,4))
     &                       - mass(i+3,j+3)*(ul(i,2,5)+d(77)*ul(i,2,4))
            end do ! j
          end do !i

          if(isw.eq.3) then
            ctan3  = ctan(3) + d(77)*ctan(2)
            do j = 1,6
              do i = 1,6
                s(is(i),is(j)) = s(is(i),is(j)) + mass(i,j)*ctan3
              end do ! i
            end do ! j
          endif
        endif

c       Transform stiffness and residual to global coordinates

        if(isw.eq.3) then
          call bm2trn (s,cs,sn,nst,ndf,1)
        endif
        if(isw.eq.3 .or. isw.eq.6) then
          call bm2trn ( r,cs,-sn,nst,ndf,2)

c         Set body loads

          call fbody2d(d,xl,ul, r,s, ndm,ndf,nst, isw)

c       Output forces

        elseif(isw.eq.4 .or. isw.eq.8) then

          call bm2trn (xl       ,cs,-sn,ndm*nel,ndm,2)

          do i = 1,3
            pp(i,1) =  r(i,1)
            pp(i,2) = -r(i,2)
            r(i,1)  =  0.0d0
            r(i,2)  =  0.0d0
          end do ! i

          if(isw.eq.4) then
            mct = mct - 3
            if (mct.le.0) then
              write(iow,2001) o,head,ttim
              if(ior.lt.0) write(*,2001) o,head,ttim
              mct = 50
            endif
            write(iow,2002) n,ma,((xl(i,j),i=1,2),(pp(i,j),i=1,3),j=1,2)
            if(ior.lt.0) then
              write(*,2002) n,ma,((xl(i,j),i=1,2),(pp(i,j),i=1,3),j=1,2)
            endif

c         Stress projections save

          else

            call frcn2d(pp,r,s)

          endif

        endif

      endif

c     Formats

2001  format(a1,20a4/5x,'Time',1p,e12.4,6x,' Element Forces '//
     & 43x,'*************  FORCE  *************'/
     &  2x,'Element Material',
     &  5x,'1-Coord',5x,'2-Coord',5x,'n-dir',7x,'s-dir',7x,'m-dir'/)

2002  format(2i9,0p,2e12.3,1p,3e12.3/18x,0p,2e12.3,1p,3e12.3)

      end

      subroutine b2mode(d,sig,eps,cc,engy,hn,h1,isw,ll)

      implicit none

      include 'bm2com.h'
      include 'ddata.h'
      include 'eldata.h'
      include 'elplot.h'
      include 'tdata.h'

      integer  ii,isw, ll, nh
      real*8   d(*),cc(3,3,2), eps(3,2),sig(3,2), engy,hn(*),h1(*)

      save

c     Compute constitutive model type

      nh   = nint(d(15))

c     Compute forces

      call bm2con (d,hn(1),h1(1),nh,cc,sig,eps,eps,isw)

c     Save forces and deformation for outputs

      ii       = 6*(ll-1)
      tt(ii+1) = sig(1,1)
      tt(ii+2) = eps(1,1)
      tt(ii+3) = sig(2,1)
      tt(ii+4) = eps(2,1)
      tt(ii+5) = sig(3,1)
      tt(ii+6) = eps(3,1)

c     Compute stored energy density

      if(isw.eq.13) then

        engy = sig(1,1)*eps(1,1)+sig(2,1)*eps(2,1)+sig(3,1)*eps(3,1)

      endif

      end
