c$Id:$
      subroutine framf2d(d,ul,xl,s,r,ndf,ndm,nst,isw)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Add geometric stiffness for eigen solutions      09/06/2009
c       2. Set d(21) to d(1) in isw = 1 to get E            14/01/2010
c       3. Use 'nv' as incrment to nh1,nh2                  17/04/2011
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: 2-d finite deformation frame element

c      Inputs:

c      Outputs:
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'augdat.h'
      include  'bdata.h'
      include  'bm2com.h'
      include  'bm2str.h'
      include  'cdata.h'
      include  'cdat1.h'
      include  'debugs.h'
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
      integer   i,ii, j,jj, ll,lint, mm,nn,nhv,nhi, nv
      real*8    cs,sn,a1,a2,a3,a4,b1,b2,le,dx,dva,dvi,cfac,dfac,xjac
      real*8    energy, xn,ux,wx
      real*8    rhoa, rhoi,a(4,4),shp(2,3),sg(2,3),forc(4),xx(2)
      real*8    xl(ndm,*),ul(ndf,nen,*)
      real*8    d(*),r(ndf,*),s(nst,nst),pp(3,3)

      save

c     Finite deformation visco-plastic BEAM element.

c     d(1)*d(32)                  = EA
c     d(37)*d(2)*d(32)/2/(1+d(2)) = kappa*GA
c     d(1)*d(33)                  = EI
c     d(4)*d(32)                  = rho*A
c     d(4)*d(33)                  = rho*I

      if(isw.eq.1) then

c       Set effective modulus to E

        d(21) = d(1)

        if(nint(d(160)).eq.2) then
          nh1    = nh1 + 1 ! Augmented storage
          d(166) = 1
        else
          d(166) = 0
        endif
        if(nint(d(20)).eq.6) then
          write(iow,2000)
        endif

c     Compute element length and direction cosines

      elseif(isw.ge.2) then
        cs = xl(1,nel) - xl(1,1)
        sn = xl(2,nel) - xl(2,1)
        le = sqrt(cs*cs + sn*sn)
        cs = cs/le
        sn = sn/le

        nlay = int(d(101))
        if(nlay.eq.0) then
          nout = 0
          nhv  = 2
          rhoa = d(4)*d(32)
          rhoi = d(4)*d(33)
        else
          nlob = int(d(102))
          nout = (nlay - 1)*(nlob - 1) + 1
          nhv  = int(d(15))*nout + 2
          call int1dl(nlob,sl)
          call bm2rho(nlay,d,d(103), rhoa,rhoi)
        endif
      endif

c     Read data

      if(isw.eq.1) then

c       History storage for rotation parameters

        nh1 = nh1 + 2

      elseif(isw.eq.3 .or. isw.eq.6) then

        nhi  = nint(d(166))
        nv   = nhi

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
        dva  = 0.5d0*rhoa*le
        dvi  = 0.5d0*rhoi*d(8)*le
        cfac = d(7)/3.d0
        dfac = 1.d0 - cfac
        do ll = 1,lint
          call shp1d(sg(1,ll),xl,shp,ndm,nel,xjac)
          dx  = sg(2,ll)*xjac
          call bm2modl(d,ul,forc,a,energy,shp,hr(nh1+nv),hr(nh2+nv),
     &                ndf,nen,isw)

c         Multiply moduli by solution parameter: ctan(1)
          do jj = 1,4
            do ii = 1,4
              a(ii,jj) = a(ii,jj)*ctan(1)
            end do ! ii
          end do ! jj

c         Compute residual and tangent

c         Mechanical and Geometric tangent terms

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
                  s(i+mm,3+nn)   = s(i+mm,3+nn) + a2*a(i,4)
                  s(3+mm,i+nn)   = s(3+mm,i+nn) + a3*a(4,i)
                  do j = 1,3
                    s(i+mm,j+nn) = s(i+mm,j+nn) + a1*a(i,j)
                  end do ! j
                end do ! i
                s(3+mm,3+nn) = s(3+mm,3+nn) + a4*a(4,4)
                nn = nn + ndf
              end do ! jj
            endif
            do i = 1,3
              r(i,ii)   = r(i,ii)   - forc(i)*b1
            end do ! i
            r(3,ii)   = r(3,ii) - forc(4)*b2
            mm = mm + ndf
          end do ! ii
          nv = nv + nhv
        end do ! ll

c       Lumped and consistent inertia contributions

        if(ndfo(1).gt.0 .or. shflg) then
          r(1,1) = r(1,1) - dva*(cfac*(ul(1,2,5)-ul(1,1,5)) + ul(1,1,5))
          r(2,1) = r(2,1) - dva*(cfac*(ul(2,2,5)-ul(2,1,5)) + ul(2,1,5))
          r(3,1) = r(3,1) - dvi*(cfac*(ul(3,2,5)-ul(3,1,5)) + ul(3,1,5))
          r(1,2) = r(1,2) + dva*(cfac*(ul(1,2,5)-ul(1,1,5)) - ul(1,2,5))
          r(2,2) = r(2,2) + dva*(cfac*(ul(2,2,5)-ul(2,1,5)) - ul(2,2,5))
          r(3,2) = r(3,2) + dvi*(cfac*(ul(3,2,5)-ul(3,1,5)) - ul(3,2,5))

          s(1    ,1    ) = s(1    ,1    ) + dva*ctan(3)*dfac
          s(2    ,2    ) = s(2    ,2    ) + dva*ctan(3)*dfac
          s(3    ,3    ) = s(3    ,3    ) + dvi*ctan(3)*dfac
          s(1    ,1+ndf) = s(1    ,1+ndf) + dva*ctan(3)*cfac
          s(2    ,2+ndf) = s(2    ,2+ndf) + dva*ctan(3)*cfac
          s(3    ,3+ndf) = s(3    ,3+ndf) + dvi*ctan(3)*cfac
          s(1+ndf,1    ) = s(1+ndf,1    ) + dva*ctan(3)*cfac
          s(2+ndf,2    ) = s(2+ndf,2    ) + dva*ctan(3)*cfac
          s(3+ndf,3    ) = s(3+ndf,3    ) + dvi*ctan(3)*cfac
          s(1+ndf,1+ndf) = s(1+ndf,1+ndf) + dva*ctan(3)*dfac
          s(2+ndf,2+ndf) = s(2+ndf,2+ndf) + dva*ctan(3)*dfac
          s(3+ndf,3+ndf) = s(3+ndf,3+ndf) + dvi*ctan(3)*dfac
        endif

c       Transform stiffness and residual to global coordinates

        if(isw.eq.3) then
          call bm2trn (s,cs,sn,nst,ndf,1)
        endif
        call bm2trn ( r,cs,-sn,nst,ndf,2)

c       Set body loads and follower forces

        call fbody2d(d,xl,ul, r,s, ndm,ndf,nst, isw)

c     Output forces

      elseif(isw.eq.4 .or. isw.eq.8) then
        nhi  = nint(d(166))
        nv   = nhi

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

          call bm2modl(d,ul,forc,a,energy,shp,hr(nh1+nv),hr(nh2+nv),
     &                ndf,nen,isw)
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
            write(iow,2002) n,ma,(xx(i),i=1,ndm),
     &                     (strs(i,1),i=1,3),(defa(i,1),i=1,3)
            if(nout.gt.0) then
              write(iow,2003) (siglr(i),i=1,nout)
              write(iow,2004) (epslr(i),i=1,nout)
            endif
            if(ior.lt.0) then
              write(*,2002) n,ma,(xx(i),i=1,ndm),
     &                     (strs(i,1),i=1,3),(defa(i,1),i=1,3)
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
                pp(i,ii)   = pp(i,ii)   - forc(i)*b1
              end do ! i
              pp(3,ii)   = pp(3,ii) - forc(4)*shp(2,ii)*dx
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

c       Mass arrays

        if(imtyp.eq.1) then

          call bm2trn (xl,cs,sn,ndm*nel,ndm,2)
          dva  = 0.5d0*rhoa*le
          dvi  = 0.5d0*rhoi*d(8)*le
          cfac = d(7)/3.d0
          dfac = 1.d0 - cfac

c         Compute lumped mass matrix

          r(1,1)         = r(1,1) + dva
          r(2,1)         = r(2,1) + dva
          r(3,1)         = r(3,1) + dvi
          r(1,2)         = r(1,2) + dva
          r(2,2)         = r(2,2) + dva
          r(3,2)         = r(3,2) + dvi

c         Compute consistent mass matrix

          s(1    ,1    ) = s(1    ,1    ) + dfac*dva
          s(2    ,2    ) = s(2    ,2    ) + dfac*dva
          s(3    ,3    ) = s(3    ,3    ) + dfac*dvi
          s(1    ,1+ndf) = s(1    ,1+ndf) + cfac*dva
          s(2    ,2+ndf) = s(2    ,2+ndf) + cfac*dva
          s(3    ,3+ndf) = s(3    ,3+ndf) + cfac*dvi
          s(1+ndf,1    ) = s(1+ndf,1    ) + cfac*dva
          s(2+ndf,2    ) = s(2+ndf,2    ) + cfac*dva
          s(3+ndf,3    ) = s(3+ndf,3    ) + cfac*dvi
          s(1+ndf,1+ndf) = s(1+ndf,1+ndf) + dfac*dva
          s(2+ndf,2+ndf) = s(2+ndf,2+ndf) + dfac*dva
          s(3+ndf,3+ndf) = s(3+ndf,3+ndf) + dfac*dvi

c       Geometric stiffness term

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

c     Augmented update or time update of history variables

      elseif(isw.eq.10 .or. isw.eq.12) then

        if(isw.eq.10 .and. nint(d(160)).eq.2 .and. debug) then
          write(  *,*) ' AUGMENT FRAME ELEMENT:',n,nint(d(160)),hr(nh2)
          write(iow,*) ' AUGMENT FRAME ELEMENT:',n,nint(d(160)),hr(nh2)
        endif
        nhi  = nint(d(166))
        nv   = nhi
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
          call bm2modl(d,ul,forc,a,energy,shp,hr(nh1+nv),hr(nh2+nv),
     &                ndf,nen,isw)
          nv = nv + nhv
        end do ! ll

c     Compute energy

      elseif(isw.eq.13) then

        nhi  = nint(d(166))
        nv   = nhi

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
        dva = 0.5d0*rhoa*a1
        dvi  =0.5d0*rhoi*d(8)*a1

c       Compute internal energy

        do ll = 1,lint

c         Compute energy density from stress and deformation

          call shp1d(sg(1,ll),xl,shp,ndm,nel,xjac)
          dx = sg(2,ll)*xjac
          nhi = nint(d(166))
          call bm2modl(d,ul,forc,a,energy,shp,hr(nh1+nv),hr(nh2+nv),
     &                ndf,nen,isw)

c         Accumulate energy

          epl(8) = epl(8) + 0.5d0*energy*dx

          nv = nv + nhv
        end do ! ll

c       Compute kinetic energy for inertia

        cfac   = d(7)/3.d0
        epl(7) = epl(7) + 0.5d0*dva*(ul(1,1,4)**2 + ul(1,2,4)**2
     &                             + ul(2,1,4)**2 + ul(2,2,4)**2
     &                       - cfac*((ul(1,1,4) - ul(1,2,4))**2
     &                             + (ul(2,1,4) - ul(2,2,4))**2))
     &                  + 0.5d0*dvi*(ul(3,1,4)**2 + ul(3,2,4)**2
     &                       - cfac*(ul(3,1,4) - ul(3,2,4))**2)

c     Initialize history variables

      elseif(isw.eq.14) then

        nhv    = nint(d(15))
        nhi    = nint(d(166))
        hr(nhi+nh1) = 1.d0
        hr(nhi+nh2) = 1.d0
        call bm2init(d,hr(nh1+nhi+2),hr(nh2+nhi+2),nhv)
      endif

c     Formats

2000  format(10x,'Conserving form element')

2001  format(a1,20a4/5x,'time',e13.5,5x,' element forces '//
     &  43x,'*********  FORCE / STRAIN  *********'/
     &   3x,'element  material',
     &  3x,'1-coord',3x,'2-coord',6x,'n-dir',8x,'s-dir',8x,'m-dir'/)

2002  format(2i10,1p,2f10.3,1p,3e13.4/40x,1p,3e13.4)

2003  format('  Stress_Layer',1p,5e13.4)
2004  format('  Strain_Layer',1p,5e13.4)

      end
