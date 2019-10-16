c$Id:$
      subroutine trussnd(d,ul,xl,ix,tl,s,p,ndf,ndm,nst,isw)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Use constants from 'pconstant.h'                 14/11/2006
c       2. Add implicit/explicit integration option         04/04/2007
c       3. Add 'jj' to 'nh3' for isw.eq.17                  17/04/2011
c-----[--.----+----.----+----.-----------------------------------------]
c     1, 2 and 3 dimensional thermo-mechanical truss element routine

c     Outputs: (isw = 4)

c         xx  - Coordinates at mid-length of truss bar
c         sig - Force on truss bar
c         eps - Strain on truss bar
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'augdat.h'
      include  'bdata.h'
      include  'cdata.h'
      include  'eldata.h'
      include  'elplot.h'
      include  'eltran.h'
      include  'hdata.h'
      include  'iofile.h'
      include  'pconstant.h'
      include  'part0.h'
      include  'part1.h'
      include  'pmod2d.h'
      include  'prstrs.h'
      include  'ptdat6.h'
      include  'rdata.h'
      include  'strnum.h'
      include  'comblk.h'

      logical   consrv,nonlin,ther,mech
      integer   i,j,ii,jj,i1,j1,istrt,nhv,nhi,tdof,ndf,ndm,nst,isw,ix(*)
      real*8    xlen, xlen0, xlen2, sig(2), eps(2), ta, body, dd(2)
      real*8    lm,cmd,cmo,hmd,hmo,geo,forc, alphar, flux, ctan1
      real*8    db(3,2),bl(3,2),br(3,2),xx(3),th(2), bf(3)
      real*8    d(*),ul(ndf,nen,*),xl(ndm,nen),tl(*),s(nst,nst),p(ndf,*)

      save

c     Set Analysis Type

      if(nint(d(39)).eq.1) then
        dtype  =  1
        consrv =  nint(d(17)).eq.4
        nonlin = .true.
      else
        dtype  =  nint(d(18))
        nonlin =  dtype.lt.0
        consrv = .false.
      endif
      hflag  = nint(d(67)).eq.1
      nhv    = nint(d(15))

c     Set nodal temperatures: Can be specified or computed

      if(isw.gt.1) then
        tdof = nint(d(19))
        if(tdof.le.0) then
          do i = 1,nel ! {
            th(i) = tl(i)
          end do ! i     }
        else
          do i = 1,nel ! {
            th(i) = ul(tdof,i,1)
          end do ! i     }
        endif
      endif

c     INPUT MATERIAL PROPERTIES

      if(isw.eq.1) then

        write(iow,2000)
        if(ior.lt.0) write(*,2000)

c       Input material parameters

        call inmate(d,tdof,ndm*2,2)
        if(nint(d(160)).eq.2) then
          nh1    = nh1 + 1 ! Augmented storage
          d(166) = 1
        else
          d(166) = 0
        endif

c       Set tdof to zero if 1, 2, or larger than ndf

        if(tdof.gt.ndf) then
          write(iow,3003)
          if(ior.lt.0) write(*,3003)
          tdof = 0
        elseif(tdof.gt.0 .and. tdof.le.ndm) then
          write(iow,3004)
          if(ior.lt.0) write(*,3004)
          tdof = 0
        endif

c       Deactivate dof in element for dof > ndm

        do i = ndm+1,ndf
          ix(i) = 0
        end do ! i

c       If temperature dof is specified activate dof

        if(tdof.gt.0) then
          ix(tdof) = 1
        endif

c       Set plot sequence

        pstyp = 1

c       Maximum plot projections

        istv  = max(2,istv)

c     CHECK FOR ZERO LENGTH ELEMENTS

      elseif(isw.eq.2) then

        if(ix(1).eq.0 .or. ix(2).eq.0) then
          write(iow,4000) n,ix(1),ix(2)
          if(ior.lt.0) then
            write(*,4000) n,ix(1),ix(2)
          endif
        else
          xlen = 0.0d0
          eps(1)  = 0.0d0
          do i = 1,ndm
            eps(1)  = max(eps(1) , abs(xl(i,1)), abs(xl(i,2)))
            xlen = max(xlen, abs(xl(i,2) - xl(i,1)))
          end do ! i
          if(xlen.eq.1.0d-10*eps(1)) then
            write(iow,4001) n
            if(ior.lt.0) then
              write(*,4001) n
            endif
          endif
        endif

c     COMPUTE ELEMENT STIFFNESS AND RESIDUAL ARRAYS

      elseif(isw.eq.3 .or. isw.eq.6) then

c       Compute residuals and tangents for mechanical/thermal parts

        mech =  npart.eq.ndfp(1)
        if(tdof.gt.ndm .and. hflag) then
          ther = npart.eq.ndfp(tdof)
        else
          ther = .false.
        endif

c       Compute stress-divergence vector (p) and stiffness matrix (s)

        if(mech) then

c         Explicit/Implicit element solutions

          if(isw.eq.3) then
            ctan1 = ctan(1)
            if(d(187).gt.0.0d0 .and.
     &         min(ctan(2),ctan(3)).gt.0.0d0) then
              ctan(1) = 0.0d0
            endif
          endif

c         Compute element stress, strain, and tangent modulus

          call strn1d(d,xl,ul,th,ndm,ndf,nen,xlen,xlen0,xx,eps,ta,
     &                consrv,nonlin,isw)
          nhi    = nint(d(166))
          if(nint(d(160)).eq.2) then ! Set augmented value
            sig(1) = hr(nh2)
          else
            sig(1) = 0.0d0
          endif
          istrt = nint(d(84))
          call modl1d(d,ta,eps,hr(nhi+nh1),hr(nhi+nh2),nhv,
     &                1,istrt,dd,sig,isw)

c         Save stress and strain for tplot

          tt(1) = sig(1)
          tt(2) = eps(1)

c         Multiply tangent modulus by length and area

          if(dtype.lt.0) then
            dd(1) = dd(1) - sig(1)
          endif

          dd(1) = (dd(1)*ctan(1) + dd(2)*d(78)*ctan(2))*d(32)*xlen

c         Compute strain-displacement matrix

          xlen2  = 1.d0/(xlen*xlen)
          alphar = xlen2/ctan(1) - xlen2
          do i = 1,ndm

            bl(i,1) = (xl(i,1) - xl(i,2))*xlen2
            br(i,1) =  bl(i,1)

c           Set non-linear terms in strain-displacement matrix

            if(nonlin) then
              bl(i,1) = bl(i,1) + (ul(i,1,1) - ul(i,2,1))*xlen2
              if(consrv) then
                br(i,1) = bl(i,1) + (ul(i,1,2) - ul(i,2,2))*alphar
              else
                br(i,1) = bl(i,1)
              endif
            endif

c           Set remaining terms and multiply by elastic modulus

            bl(i,2) = -bl(i,1)
            br(i,2) = -br(i,1)

            db(i,1) =  bl(i,1)*dd(1)
            db(i,2) = -db(i,1)

          end do ! i

c         Compute mass terms

          cmd = d(4)*d(32)*xlen0*one3
          cmo = cmd*0.5d0
          lm  = cmd*1.5d0
          hmd = d(7)*cmd + (1.d0 - d(7))*lm
          hmo = d(7)*cmo

c         Form internal, body  and inertia force vector
c         Include stiffness proportional Rayleigh damping in force
c         Include mass proportional Rayleigh damping in residual

          forc = (sig(1) + d(78)*sig(2))*xlen*d(32)

          call sbodyf(d, bf)

          do i = 1,ndm

c           Set body loading factors

            body = 0.5d0*xlen0*bf(i)

            p(i,1) = body   - bl(i,1)*forc
     &             - hmd*(ul(i,1,5) + d(77)*ul(i,1,4))
     &             - hmo*(ul(i,2,5) + d(77)*ul(i,2,4))

            p(i,2) = body   - bl(i,2)*forc
     &             - hmo*(ul(i,1,5) + d(77)*ul(i,1,4))
     &             - hmd*(ul(i,2,5) + d(77)*ul(i,2,4))

          end do ! i

c         Compute stiffness terms

          if(isw.eq.3) then

            i1 = 0
            do ii = 1,2
              j1 = 0
              do jj = 1,2
                do i = 1,ndm
                  do j = 1,ndm
                    s(i+i1,j+j1) = db(i,ii)*br(j,jj)
                  end do ! j
                end do ! i
                j1 = j1 + ndf
              end do ! jj
              i1 = i1 + ndf
            end do ! ii

c           Correct for geometric and inertial tangent effects

            if(nonlin .and. gflag) then
              geo = sig(1)*d(32)/xlen*ctan(1)
            else
              geo = 0.0d0
            endif

c           Set diagonal and off-diagonal terms for geometric stiffness
c           Include mass proportional rayleigh damping

            if(ndfo(1).gt.0 .or. shflg) then
              hmd =  hmd*(ctan(3) + d(77)*ctan(2)) + geo
              hmo =  hmo*(ctan(3) + d(77)*ctan(2)) - geo
            else
              hmd =  geo
              hmo = -geo
            endif

            do i = 1,ndm
              j      = i + ndf
              s(i,i) = s(i,i) + hmd
              s(i,j) = s(i,j) + hmo
              s(j,i) = s(j,i) + hmo
              s(j,j) = s(j,j) + hmd
            end do ! i

          endif
          if(isw.eq.3) then
            ctan(1) = ctan1
          endif

        endif

c       Compute thermal vector (p) and matrix (s)

        if(ther) then

          call thertrs(d,ul,xl,s(tdof,tdof),p(tdof,1),
     &                 ndf,ndm,nen,nst,tdof,isw)

        endif

c     OUTPUT STRESS AND STRAIN IN ELEMENT

      elseif(isw.eq.4 .or. isw.eq.8) then

c       Check for thermal parts

        if(tdof.gt.ndm .and. hflag) then
          ther = .true.
        else
          ther = .false.
        endif

c       Form strain and stress

        call strn1d(d,xl,ul,th,ndm,ndf,nen,xlen,xlen0,xx,eps,ta,
     &              consrv,nonlin,isw)
        nhi  = nint(d(166))
        if(nint(d(160)).eq.2) then ! Set augmented value
          sig(1) = hr(nh2)
        else
          sig(1) = 0.0d0
        endif
        istrt = nint(d(84))
        call modl1d(d,ta,eps,hr(nhi+nh1),hr(nhi+nh2),nhv,
     &              1,istrt,dd,sig,isw)

c       Form truss force: multiply stress by area

        sig(1) = d(32)*sig(1)

        if(ther) then
          flux = d(61)*(ul(tdof,2,1) - ul(tdof,1,1))/xlen0
        else
          flux = 0.0d0
        endif

c       Output element force/strains

        if(isw.eq.4) then
          mct = mct - 1
          if(mct.le.0) then
            write(iow,2001) o,head
            if(ior.lt.0) write(*,2001) o,head
            mct = 50
          endif
          write(iow,2002) n,ma,xx,sig(1),eps(1),flux
          if(ior.lt.0) then
            write(*,2002) n,ma,xx,sig(1),eps(1),flux
          endif
        elseif(nint(d(171)).eq.0) then
          if(tdof.gt.ndm .and. hflag) then
            ther = npart.eq.ndfp(tdof)
          else
            ther = .false.
          endif
          call trcnnd(sig(1),flux,p,s,nen,ther)
        endif

c     COMPUTE ELEMENT MASS MATRICES

      elseif(isw.eq.5) then

        xlen0 = 0.0d0
        do i = 1,ndm
          xlen0 = xlen0 + (xl(i,2)-xl(i,1))**2
        end do ! i

c       Compute mass matrix terms

        cmd = d(4)*d(32)*sqrt(xlen0)*one3
        lm  = 1.5d0*cmd
        cmo = 0.5d0*cmd
        hmd = d(7)*cmd + (1.d0 - d(7))*lm
        hmo = d(7)*cmo

        do i = 1,ndm

c         Higher order mass

          j      = i + ndf
          s(i,i) = hmd
          s(j,j) = hmd
          s(j,i) = hmo
          s(i,j) = hmo

c         Lumped mass

          p(i,1) = lm
          p(i,2) = p(i,1)

        end do ! i

c     Augmented update

      elseif(isw.eq.10) then

c       Form strain and stress

        if(nint(d(160)).eq.2) then ! Set augmented value
          call strn1d(d,xl,ul,th,ndm,ndf,nen,xlen,xlen0,xx,eps,ta,
     &                consrv,nonlin,isw)
          nhi    = nint(d(166))
          sig(1) = hr(nh2)
          istrt  = nint(d(84))
          call modl1d(d,ta,eps,hr(nhi+nh1),hr(nhi+nh2),nhv,
     &                1,istrt,dd,sig,isw)

          hr(nh2) = hr(nh2) + augf*(sig(1) - hr(nh2))
        endif

c     COMPUTE ELEMENT ENERGY

      elseif(isw.eq.13) then

c       Compute element stress, strain, and tangent modulus

        call strn1d(d,xl,ul,th,ndm,ndf,nen,xlen,xlen0,xx,eps,ta,
     &              consrv,nonlin,isw)
        nhi   = nint(d(166))
        istrt = nint(d(84))
        call modl1d(d,ta,eps,hr(nhi+nh1),hr(nhi+nh2),nhv,
     &              1,istrt,dd,sig,isw)

c       Stored energy

        epl(8) = epl(8) + 0.5d0*eps(1)*sig(1)*d(32)*xlen

c       Compute mass terms

        cmd = d(4)*d(32)*xlen0*one3
        cmo = cmd*0.5d0
        lm  = cmd*1.5d0
        hmd =(d(7)*cmd + (1.d0 - d(7))*lm)*0.5d0
        hmo = d(7)*cmo

c       Kinetic energy

        do i = 1,ndm
          epl(7) = epl(7) + hmd*(ul(i,1,4)**2 + ul(i,2,4)**2)
     &                    + hmo* ul(i,1,4)*ul(i,2,4)
        end do ! i

c     Initialize history variables in constitution

      elseif(isw.eq.14) then

        nhi   = nint(d(166))
        istrt = nint(d(84))
        call modl1d(d,ta,eps,hr(nhi+nh1),hr(nhi+nh2),nhv,
     &              1,istrt,dd,sig,isw)

c     Initialize element strains for activation

      elseif(isw.eq.17) then

        jj = 0
        do i = 1,2
          do j = 1,ndm
            hr(nh3+jj) = ul(j,i,1)
            jj     = jj + 1
          end do ! j
        end do ! i

c     Initialize element strains for deactivation

      elseif(isw.eq.18) then

        do i = 0,2*ndm-1
          hr(nh3+i) = 0.0d0
        end do ! i

c     Critical time step computation

      elseif(isw.eq.21) then

        call tcritnd(d,xl,ul,ndm,ndf,nel)

c     External node determination

      elseif(isw.eq.26) then

      endif

c     FORMATS

2000  format(9x,'T r u s s    E l e m e n t'/1x)
2001  format(a1,20a4//9x,'Truss Element'//' Elmt Matl    ',
     & '1-coord    2-coord    3-coord     Force',9x,'Strain',10x,'Flux')
2002  format(2i5,1p,3e11.3,1p,3e14.5)

3003  format(' *WARNING* Thermal d.o.f. > active d.o.f.s : Set to 0')
3004  format(' *WARNING* Thermal d.o.f. can not be <= ndm: Set to 0')

4000  format(' *ERROR* Element',i7,' has nodes',2i8)
4001  format(' *ERROR* Element',i7,' has zero length')

      end

      subroutine strn1d(d,xl,ul,tl,ndm,ndf,nen,xlen,xlen0,xx,eps,ta,
     &                  consrv,nonlin,isw)

c     Compute constitutive equation

      implicit  none

      include  'iofile.h'
      include  'elcoor.h'
      include  'eltran.h'
      include  'pmod2d.h'

      logical   consrv,nonlin
      integer   i,ndm,ndf,nen, isw
      real*8    xlen,xlen0,xlenn,xlen1,xlena,dx,du,dd, ta, alpha,alphar
      real*8    xx(3),d(*),xl(ndm,*),ul(ndf,nen,*),tl(*),eps(*)

      save

c     Set integration parameter: For energy = 1.0

      if(isw.eq.13) then
        alpha = 1.d0
      else
        alpha = ctan(1)
      endif
      alphar = 1.d0/alpha - 1.d0

c     Compute length and strain terms

      xlen0  = 0.0d0
      xlenn  = 0.0d0
      xlen1  = 0.0d0
      xlena  = 0.0d0
      eps(1) = 0.0d0
      eps(2) = 0.0d0
      do i = 1,ndm
        dx    = xl(i,2)   - xl(i,1)
        du    = ul(i,2,1) - ul(i,1,1)
        dd    = ul(i,2,2) - ul(i,1,2)
        xlen0 = xlen0 +  dx**2
        xlen1 = xlen1 + (dx + du + dd*alphar)**2
        xlenn = xlenn + (dx + du - dd       )**2
        xlena = xlena + (dx + du            )**2
        eps(1)= eps(1) +  dx*du
        eps(2)= eps(2) +  dx*(ul(i,2,4) - ul(i,1,4))
        xx(i) = (xl(i,2) + xl(i,1))*0.5d0
      end do ! i

      do i = 1,ndm
        xref(i) = xx(i)
        xcur(i) = xx(i) + (ul(i,1,1) + ul(i,2,1))*0.5d0
      end do ! i
      do i = ndm+1,3
        xref(i) = 0.0d0
        xcur(i) = 0.0d0
      end do ! i

c     Compute temperature change

      ta   = 0.5d0*(tl(1) + tl(2)) - d(9)

c     Compute strain forms

      if(dtype.lt.0) then

c       Logarithmic stretch strain and deformed length

        eps(1) = 0.5d0*log(xlena/xlen0)
        xlen0 = sqrt(xlen0)
        xlen  = sqrt(xlena)

      else

c       Green or linear strain

        if(nonlin) then
          if(consrv) then
            eps(1)= 0.5d0*(alpha*(xlen1 - xlenn) + xlenn)/xlen0 - 0.5d0
          else
            eps(1)= 0.5d0*xlena/xlen0 - 0.5d0
          endif
        else
          eps(1) = eps(1)/xlen0
          eps(2) = eps(2)/xlen0
        endif
        xlen0 = sqrt(xlen0)
        xlen  = xlen0
      endif

      end

      subroutine thertrs(d,ul,xl,s,p,ndf,ndm,nen,nst,tdof,isw)

c     Thermal tangent/residual computation

      implicit  none

      include  'eltran.h'
      include  'pconstant.h'

      integer   ndf,ndm,nen,nst,tdof,isw,i
      real*8    k, rhoc,xlen
      real*8    d(*),ul(ndf,nen,*),xl(ndm,*),s(nst,*),p(ndf,*)

      save

c     Compute length

      xlen = 0.0d0
      do i = 1,ndm
        xlen  = xlen + (xl(i,2) - xl(i,1))**2
      end do ! i
      xlen = sqrt(xlen)
      k    = d(61)/xlen
      rhoc = d(4)*d(64)*xlen*0.5d0

c     Form transient (mass) part of tangent

      i      = ndf + 1
      s(1,i) = rhoc*d(7)*one3
      s(1,1) = rhoc - s(1,i)
      s(i,1) = s(1,i)
      s(i,i) = s(1,1)

      if(isw.eq.3 .or. isw.eq.6) then

c       Form steady state part of residual

        p(1,1) =  k*(ul(tdof,2,1) - ul(tdof,1,1))
        p(1,2) = -p(1,1)

        p(1,1) = p(1,1) - s(1,1)*ul(tdof,1,4) - s(1,i)*ul(tdof,2,4)
        p(1,2) = p(1,2) - s(i,1)*ul(tdof,1,4) - s(i,i)*ul(tdof,2,4)

c       Form final tangent

        s(1,1) =  ctan(1)*k + ctan(2)*s(1,1)
        s(i,1) = -ctan(1)*k + ctan(2)*s(i,1)
        s(1,i) = -ctan(1)*k + ctan(2)*s(1,i)
        s(i,i) =  ctan(1)*k + ctan(2)*s(i,i)

      endif

      end

      subroutine trcnnd(sig,flux,dt,st,nen,ther)

      implicit  none

      include  'strnum.h'

      logical   ther
      integer   nen, i
      real*8    dt(*),st(nen,*),sig(*),flux

c     Stress projections

      do i = 1,2
c       dt(i)   = dt(i)   + 1.d0
c       st(i,1) = st(i,1) + sig(1)
        dt(i)   = 1.d0
        st(i,1) = sig(1)
      end do ! i

c     Thermal part

      if(ther) then
        do i = 1,2
          st(i,2) = flux
c         st(i,2) = st(i,2) + flux
        end do ! i
        iste  = 2
      else
        iste  = 1
      endif

      end
