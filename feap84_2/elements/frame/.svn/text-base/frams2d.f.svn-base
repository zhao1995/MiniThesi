c$Id:$
      subroutine frams2d(d,ul,xl,s,r,ndf,ndm,nst,isw)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Add geometric stiffness for eigen solutions      09/06/2009
c       2. Set d(21) to d(1) in isw = 1 to get E            14/01/2010
c-----[--.----+----.----+----.-----------------------------------------]
c     Purpose: Two dimensional small displacement frame element
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
      include  'prstrs.h'
      include  'ptdat6.h'
      include  'rdata.h'
      include  'tdata.h'
      include  'comblk.h'

      integer   ndf,ndm,nst,isw, nlay,nlob
      integer   i,ii,i1,istrt, j,jj,j1, ll,lint, mm,nn,nhv, nv
      real*8    cs,sn,a1,a2,a3,a4,b1,b2,le,dx,dva,dvi,cfac,dfac,xjac
      real*8    energy, rhoa, rhoi, ctan1, ctan3, xn,ux,wx
      real*8    aa(4,4,2),shp(2,3),sg(2,3),forc(4,2),xx(2),dma(3)
      real*8    xl(ndm,*),ul(ndf,nen,*)
      real*8    d(*),r(ndf,*),s(nst,nst),pp(3,3)

      save

c     Small deformation element.

c     d(21)*d(32)       = EA
c     d(37)*d(27)*d(32) = kappa*GA
c     d(21)*d(33)       = EI
c  ?  d(x)              = gamma
c  ?  d(x)              = No
c  ?  d(x)              = Vo
c  ?  d(x)              = Mo
c     d( 4)*d(32)       = rho*A
c     d( 4)*d(33)       = rho*I

c     Compute element length and direction cosines

      if(isw.ge.2) then
        cs = xl(1,nel) - xl(1,1)
        sn = xl(2,nel) - xl(2,1)
        le = sqrt(cs*cs + sn*sn)
        cs = cs/le
        sn = sn/le

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

      elseif(isw.eq.3 .or. isw.eq.6) then

        nv   = 0
        if(nint(d(182)).gt.0) then
          lint = nel
          call int1dn(lint, sg)
        else
          lint = nel - 1
          call int1d(lint, sg)
        endif
        call bm2trn (ul(1,1,1),cs,sn,ndf*nel,ndf,2)
        call bm2trn (ul(1,1,2),cs,sn,ndf*nel,ndf,2)
        call bm2trn (ul(1,1,5),cs,sn,ndf*nel,ndf,2)
        call bm2trn (xl       ,cs,sn,ndm*nel,ndm,2)
        dma(1) = 0.5d0*rhoa*le
        dma(2) = dma(1)
        dma(3) = 0.5d0*rhoi*d(8)*le
        cfac   = d(7)/3.d0
        dfac   = 1.d0 - cfac
        do ll = 1,lint
          call shp1d(sg(1,ll),xl,shp,ndm,nel,xjac)
          dx = sg(2,ll)*xjac
          call b2mods(d,ul,forc,aa,energy,shp,hr(nh1+nv),hr(nh2+nv),
     &                ndf,nen,isw,ll)

c         Multiply moduli by solution parameter: ctan(1)

c         ctan1 = ctan(1) + d(78)*ctan(2)
          ctan1 = ctan(1)
          do jj = 1,4
            do ii = 1,4
              aa(ii,jj,1) = aa(ii,jj,1)*ctan1
            end do ! ii
          end do ! jj

c         Compute residual and tangent

c         Mechanical tangent terms

          mm = 0
          do ii = 1,nel
            b1 = shp(1,ii)*dx
            b2 = shp(2,ii)*dx
            if(isw.eq.3) then
              nn = 0
              do jj = 1,nel
                a1 = b1*shp(1,jj)
                a2 = b1*shp(2,jj)
                a3 = b2*shp(1,jj)
                a4 = b2*shp(2,jj)
                do i = 1,3
                  s(i+mm,3+nn)   = s(i+mm,3+nn) + a2*aa(i,4,1)
                  s(3+mm,i+nn)   = s(3+mm,i+nn) + a3*aa(4,i,1)
                  do j = 1,3
                    s(i+mm,j+nn) = s(i+mm,j+nn) + a1*aa(i,j,1)
                  end do ! j
                end do ! i
                s(3+mm,3+nn) = s(3+mm,3+nn) + a4*aa(4,4,1)
                nn = nn + ndf
              end do ! jj
            endif
            do i = 1,3
              r(i,ii)   = r(i,ii)   - forc(i,1)*b1
            end do ! i
            r(3,ii)   = r(3,ii) - forc(4,1)*b2
            mm = mm + ndf
          end do ! ii
          nv = nv + nhv
        end do ! ll

c       Lumped and consistent inertia contributions

        if(ndfo(1).gt.0 .or. shflg) then
          ctan3 = ctan(3) + d(77)*ctan(2)
          do i = 1,3
            r(i,1) = r(i,1) - dma(i)*(dfac*(ul(i,1,5)+d(77)*ul(i,1,4))
     &                              + cfac*(ul(i,2,5)+d(77)*ul(i,2,4)))
            r(i,2) = r(i,2) - dma(i)*(cfac*(ul(i,1,5)+d(77)*ul(i,1,4))
     &                              + dfac*(ul(i,2,5)+d(77)*ul(i,2,4)))

            s(i    ,i    ) = s(i    ,i    ) + dma(i)*ctan3*dfac
            s(i    ,i+ndf) = s(i    ,i+ndf) + dma(i)*ctan3*cfac
            s(i+ndf,i    ) = s(i+ndf,i    ) + dma(i)*ctan3*cfac
            s(i+ndf,i+ndf) = s(i+ndf,i+ndf) + dma(i)*ctan3*dfac
          end do ! i
        endif

c       Transform stiffness and residual to global coordinates

        if(isw.eq.3) then
          call bm2trn (s,cs,sn,nst,ndf,1)
        endif
        call bm2trn ( r,cs,-sn,nst,ndf,2)

c       Set body loads

        call fbody2d(d,xl,ul, r,s, ndm,ndf,nst, isw)

c     Output forces

      elseif(isw.eq.4 .or. isw.eq.8) then
        nv   = 0
        if(nint(d(182)).gt.0) then
          lint = nel
          call int1dn(lint, sg)
        else
          lint = nel - 1
          call int1d(lint, sg)
        endif

c       Loop over quadrature points

        call bm2trn (ul(1,1,1),cs,sn,ndf*nel,ndf,2)
        call bm2trn (ul(1,1,2),cs,sn,ndf*nel,ndf,2)
        call bm2trn (ul(1,1,5),cs,sn,ndf*nel,ndf,2)
        call bm2trn (xl,cs,sn,ndm*nel,ndm,2)
        call pzero  (pp,9)
        do ll = 1,lint
          call shp1d(sg(1,ll),xl,shp,ndm,nel,xjac)

c         Output forces

          call b2mods(d,ul,forc,aa,energy,shp,hr(nh1+nv),hr(nh2+nv),
     &                ndf,nen,isw,ll)
          if(isw.eq.4) then
            do i = 1,ndm
              xx(i) = 0.
              do ii = 1,nel
                xx(i) = xx(i) + xl(i,ii)*shp(2,ii)
              end do ! ii
            end do ! i
            mct = mct - 3
            if (mct.le.0) then
              write(iow,2001) o,head,ttim
              if(ior.lt.0) write(*,2001) o,head,ttim
              mct = 50
            endif
            write(iow,2002) n,ma,(xx(i),i=1,2),
     &                      (strs(i,1),i=1,3),(defa(i,1),i=1,3)
            if(nout.gt.0) then
              write(iow,2003) (siglr(i),i=1,nout)
              write(iow,2004) (epslr(i),i=1,nout)
            endif
            if(ior.lt.0) then
              write(*,2002) n,ma,(xx(i),i=1,2),
     &                      (strs(i,1),i=1,3),(defa(i,1),i=1,3)
              if(nout.gt.0) then
                write(*,2003) (siglr(i),i=1,nout)
                write(*,2004) (epslr(i),i=1,nout)
              endif
            endif

c         Stress projections save

          else

            dx = sg(2,ll)*xjac
            do ii = 1,nel
              b1 = shp(1,ii)*dx
              do i = 1,3
                pp(i,ii)   = pp(i,ii)   - forc(i,1)*b1
              end do ! i
              pp(3,ii)   = pp(3,ii) - forc(4,1)*shp(2,ii)*dx
            end do ! ii

          endif

          nv = nv + nhv
        end do ! ll

        if(isw.eq.8) then
          do i = 1,3
            pp(i,2) = -pp(i,2)
          end do ! i
          call frcn2d(pp,r,s)
        endif

c     Compute mass or geometric stiffness array

      elseif(isw.eq.5) then

c       Mass matrix

        if(imtyp.eq.1) then

          cfac = d(7)
          dfac = 1.d0 - d(7)

          lint = nel
          if(nint(d(182)).gt.0) then
            call int1dn(lint, sg)
          else
            call int1d(lint, sg)
          endif

c         Compute mass matrix

          call bm2trn (xl,cs,sn,ndm*nel,ndm,2)
          do ll = 1,lint

c           Compute shape functions

            call shp1d(sg(1,ll),xl,shp,ndm,nel,xjac)

            dma(1) = sg(2,ll)*xjac*rhoa
            dma(2) = dma(1)
            dma(3) = sg(2,ll)*xjac*rhoi*d(8)

c           Compute db = rho*shape*dv

            j1 = 0
            do jj = 1,nel

c             Compute lumped mass

              do i = 1,3
                r(i,jj)      = r(i,jj)      + dma(i)*shp(2,jj)
                s(j1+i,j1+i) = s(j1+i,j1+i) + dma(i)*shp(2,jj)*dfac
              end do ! i

c             Compute consistent mass

              i1 = 0
              do ii = 1,nel
                ctan3 = shp(2,ii)*shp(2,jj)*cfac
                do i = 1,3
                  s(i1+i,j1+i) = s(i1+i,j1+i) +dma(i)*ctan3
                end do ! i
                i1 = i1 + ndf
              end do ! ii
              j1 = j1 + ndf
            end do ! j
          end do ! ll

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

        nv   = 0
        if(nint(d(182)).gt.0) then
          lint = nel
          call int1dn(lint, sg)
        else
          lint = nel - 1
          call int1d(lint, sg)
        endif
        call bm2trn (ul(1,1,1),cs,sn,ndf*nel,ndf,2)
        call bm2trn (ul(1,1,2),cs,sn,ndf*nel,ndf,2)
        call bm2trn (xl       ,cs,sn,ndm*nel,ndm,2)
        do ll = 1,lint
          call shp1d(sg(1,ll),xl,shp,ndm,nel,xjac)
          call b2mods(d,ul,forc,aa,energy,shp,hr(nh1+nv),hr(nh2+nv),
     &                ndf,nen,isw,ll)
          nv = nv + nhv
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

        nv   = 0
        do ll = 1,lint

c         Compute energy density from stress and deformation

          call shp1d(sg(1,ll),xl,shp,ndm,nel,xjac)
          dx = sg(2,ll)*xjac
          call b2mods(d,ul,forc,aa,energy,shp,hr(nh1+nv),hr(nh2+nv),
     &                ndf,nen,isw,ll)

c         Accumulate energy

          epl(8) = epl(8) + 0.5d0*energy*dx

          nv = nv + nhv
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

      endif

c     Formats

2001  format(a1,20a4/5x,'Time',1p,e12.4,6x,' Element Forces '//
     & 43x,'*********  FORCE / STRAIN  *********'/
     &  2x,'Element Material',
     &  5x,'1-Coord',5x,'2-Coord',5x,'n-dir',7x,'s-dir',7x,'m-dir'/)

2002  format(2i9,0p,2e12.3,1p,3e12.3/42x,1p,3e12.3)

2003  format('  Stress_Layer',1p,5e13.4)
2004  format('  Strain_Layer',1p,5e13.4)

      end

      subroutine b2mods(d,ul,forca,aa,engy,shp,hn,h1,ndf,nen,isw,ll)

      implicit none

      include 'bm2com.h'
      include 'ddata.h'
      include 'eldata.h'
      include 'elplot.h'
      include 'tdata.h'

      integer  ndf,nen, i,ii,isw, ll, nh
      real*8   d(*),ul(ndf,nen,*),cc(3,3,2),hn(*),h1(*)
      real*8   forca(4,2),aa(4,4,2),def1(3),shp(2,3), engy

      save

c     Compute constitutive model type

      nh   = nint(d(15))

c     Compute beam strains

      do i = 1,3
        defa(i,1) = 0.0d0
        defa(i,2) = 0.0d0
        do ii = 1,nel
          defa(i,1) = defa(i,1) + ul(i,ii,1)*shp(1,ii)
          defa(i,2) = defa(i,2) + ul(i,ii,4)*shp(1,ii)
        end do ! ii
      end do ! i

      do ii = 1,nel
        defa(2,1) = defa(2,1) - ul(3,ii,1)*shp(2,ii)
        defa(2,2) = defa(2,2) - ul(3,ii,4)*shp(2,ii)
      end do ! ii

c     Compute forces

      call bm2con (d,hn(1),h1(1),nh,cc,strs,defa,def1,isw)

c     Save forces and deformation for outputs

      ii       = 6*(ll-1)
      tt(ii+1) = strs(1,1)
      tt(ii+2) = defa(1,1)
      tt(ii+3) = strs(2,1)
      tt(ii+4) = defa(2,1)
      tt(ii+5) = strs(3,1)
      tt(ii+6) = defa(3,1)

c     Compute stored energy density

      if(isw.eq.13) then

        engy = strs(1,1)*defa(1,1)
     &       + strs(2,1)*defa(2,1)
     &       + strs(3,1)*defa(3,1)

      elseif(isw.ne.12) then

        do ii = 1,2

c         Compute first Piola-material frame

          forca(1,ii) =  strs(1,ii)
          forca(2,ii) =  strs(2,ii)
          forca(3,ii) =  strs(3,ii)
          forca(4,ii) = -forca(2,ii)

c         Compute tangent tensor

          aa(1,1,ii) = cc(1,1,ii)
          aa(1,2,ii) = cc(1,2,ii)
          aa(1,3,ii) = cc(1,3,ii)

          aa(2,1,ii) = cc(1,2,ii)
          aa(2,2,ii) = cc(2,2,ii)
          aa(2,3,ii) = cc(2,3,ii)

          aa(3,1,ii) = cc(1,3,ii)
          aa(3,2,ii) = cc(2,3,ii)
          aa(3,3,ii) = cc(3,3,ii)

          aa(1,4,ii) = -cc(1,2,ii)
          aa(2,4,ii) = -cc(2,2,ii)
          aa(3,4,ii) = -cc(2,3,ii)

          aa(4,1,ii) = aa(1,4,ii)
          aa(4,2,ii) = aa(2,4,ii)
          aa(4,3,ii) = aa(3,4,ii)

          aa(4,4,ii) = cc(2,2,ii)
        end do ! ii

      endif

      end
