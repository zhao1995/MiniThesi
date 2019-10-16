c$Id:$
      subroutine stvtor(d,ul,xl,ix,s,p,ndf,ndm,nst,isw)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All Rights Reserved

c-----[--+---------+---------+---------+---------+---------+---------+-]
c     Purpose:  Two dimensional elastic/plastic St. Venant torsion
c               element

c-----[--+---------+---------+---------+---------+---------+---------+-]
c     1. Parameters input by this routine when isw = 1 (i.e., mate)

c        Record 1. (type,mu,tau_0,H_iso)

c           type  - elastic or mises (plastic)
c           mu    - elastic shear modulus
c           tau_0 - initial yield stress in shear
c           H_iso - plastic hardening modulus in shear
c-----[--+---------+---------+---------+---------+---------+---------+-]
c     2. Control parameters

c        This is a two dimensional element which can analyze plane
c        geometries.  Set control parameters as follows:

c           ndm - set to 2     (x,y or r,z-coords)
c           ndf - set > or = 1 (nodal warping)
c           nel - set > or = 4  and < or = 9

c                    A eta
c                    |
c             4      7      3
c              o-----o-----o
c              |     |     |
c              |     |     |
c             8o     o-----o6---> xi
c              |    9      |
c              |           |
c              o-----o-----o
c             1      5      2

c               Node numbering
c-----[--+---------+---------+---------+---------+---------+---------+-]
      implicit   none

      include   'bdata.h'
      include   'cdata.h'
      include   'debugs.h'
      include   'eldata.h'
      include   'eltran.h'
      include   'fdata.h'
      include   'hdata.h'
      include   'iofile.h'
      include   'prstrs.h'

      include   'comblk.h'

      character  type*15, region*7
      logical    errck, pcomp, tinput, start

      integer    ndf,ndm,nst,isw
      integer    i,j, i1,j1, l,lint, nhv,nn
      real*8     xx,yy, xsj, a1,a2,bd1,bd2, shj, torque, jiner

      integer    ix(*)
      real*8     d(*),ul(ndf,*),xl(ndm,*),s(nst,*),p(*)
      real*8     dd(2,2),gam(2),tau(5),shp(3,9),sg(3,9),td(3)

      save

c     Output Element Type

      if(isw.eq.0) then

        if(ior.lt.0) then
          write(*,*) '   Elmt  2: 2-d Linear Heat Conduction.'
        endif

c     Input material properties

      elseif(isw.eq.1) then

        d(2)  = 1.d+20
        nh1   = 0
        start = .true.
        do while (start)

          errck = tinput(type,1,td,3)

c         Elastic case

          if(pcomp(type,'elas',4)) then
            d(1) = td(1)

c         Plastic case (Mises)

          elseif(pcomp(type,'plas',4)) then
            d(2) = td(1)
            d(3) = td(2)
            nh1  = 3*nen
            d(6) = 3
            d(7) = 1.0d0

c         Error values

          elseif(pcomp(type,'erro',4)) then
            d(22) = td(1)
          elseif(pcomp(type,'    ',4)) then
            start = .false.
          endif

        end do ! while

c       Set number of quadrature points/direction

        if(nen.le.4) then
          d(5) = 2
        else
          d(5) = 3
        endif

c       Output current parameters

        if(nint(d(7)).eq.0) then
          write(iow,2000) d(1)
          if(ior.lt.0) then
            write(iow,2000) d(1)
          endif
        else
          write(iow,2001) d(1),d(2),d(3)
          if(ior.lt.0) then
            write(iow,2001) d(1),d(2),d(3)
          endif
        endif

c       Error value

        if(d(22).gt.0.0d0) then
          write(iow,2002) d(22)
          if(ior.lt.0) then
            write(iow,2002) d(22)
          endif
        endif

c       Initialize torque and inertia

        torque = 0.0d0
        jiner  = 0.0d0

c     Check of mesh if desired (chec)

      elseif(isw.eq.2) then

        call ckisop(ix,xl,shp,ndm)

c     Compute residual and tangent arrays

      elseif(isw.eq.3 .or. isw.eq.6) then

c       Set the total torque to zero

        if(n.le.1) then
          torque = 0.0d0
          jiner  = 0.0d0
        endif

c       Quadrature loop

        l   = nint(d(5))
        call int2d(l,lint,sg)
        nn  = 0
        nhv = nint(d(6))
        do l = 1,lint

          call shp2d(sg(1,l),xl,shp,xsj,ndm,nel,ix,.false.)

c         Compute strains

          call strator(xl,ul, xx,yy,shp,gam, ndm,ndf)

c         Compute stresses

          call stretor(d,gam,hr(nh1+nn),hr(nh2+nn),tau,dd, region)

          do i = 1,2
            do j = 1,2
              dd(i,j) = dd(i,j)*ctan(1)
            end do
          end do

          xsj = xsj*sg(3,l)

c         Accumulate the inertia

          jiner  = jiner + (gam(1)**2 + gam(2)**2)*xsj

c         Accumulate the torque

          torque = torque + (xx*tau(2) - yy*tau(1))*xsj

c         Compute the other arrays

          j1 = 1
          do j = 1,nel
            a1 = shp(1,j)*xsj
            a2 = shp(2,j)*xsj

c           Compute residual

            p(j1) = p(j1) - a1*tau(1) - a2*tau(2)

c           Compute tangent

            bd1 = a1*dd(1,1) + a2*dd(2,1)
            bd2 = a1*dd(1,2) + a2*dd(2,2)

            i1 = 1
            do i = 1,nel
              s(i1,j1) = s(i1,j1) + bd1*shp(1,i) + bd2*shp(2,i)
              i1 = i1 + ndf
            end do
            j1 = j1 + ndf
          end do

          nn = nn + nhv

        end do

c     Compute tau stress and print at center of element

      elseif(isw.eq.4) then

        l   = nint(d(5))
        call int2d(l,lint,sg)
        nn  = 0
        nhv = nint(d(6))
        do l = 1,lint

          call shp2d(sg(1,l),xl,shp,xsj,ndm,nel,ix,.false.)

c         Compute and output tau and gradients

c         Compute strains

          call strator(xl,ul, xx,yy,shp,gam, ndm,ndf)

c         Compute stresses

          call stretor(d,gam,hr(nh1+nn),hr(nh2+nn),tau,dd, region)

          mct = mct - 1
          if(mct.le.0) then
            write(iow,2003) o,head
            if(ior.lt.0 .and. pfr) then
              write(*,2003) o,head
            endif
            mct = 50
          endif
          write(iow,2004) n,ma,xx,yy,(tau(i),i=1,2),gam,
     &                    region,(tau(i),i=3,5),torque,jiner
          if(ior.lt.0 .and. pfr) then
            write(*,2004) n,ma,xx,yy,(tau(i),i=1,2),gam,
     &                    region,(tau(i),i=3,5),torque,jiner
          endif

          nn = nn + nhv

        end do

c     Compute tau capacity (mass) matrix for eigencomputations

      elseif(isw.eq.5) then
        l = nint(d(5))
        call int2d(l,lint,sg)
        do l = 1,lint
          call shp2d(sg(1,l),xl,shp,xsj,ndm,nel,ix,.false.)
          xsj = xsj*sg(3,l)
          j1 = 1
          do j = 1,nel
            shj = d(2)*shp(3,j)*xsj

c           Lumped capacity (lmas)

            p(j1) = p(j1) + shj
            i1 = 1

c           Consistent capacity (cmas)

            do i = 1,nel
              s(i1,j1) = s(i1,j1) + shj*shp(3,i)
              i1 = i1 + ndf
            end do
            j1 = j1 + ndf
          end do
        end do

c     Compute surface tau loading (not implemented)

c     elseif(isw.eq.7) then

c     Compute nodal tau for print/plots

      elseif(isw.eq.8) then

        call stcntor(ix,d,xl,ul,s,shp,hr(nph),hr(nph+numnp),hr(ner),
     &               erav,ndf,ndm,nel,numnp)

c     Compute error data for tau

      elseif(isw.eq.9 .or. isw.eq.11) then

        call stertor(ix,d,xl,ul,shp,hr(nph+numnp),ndf,ndm,nel,numnp)

      endif

c     Formats

2000  format(5x,'Elastic St. Venant Torsion Element'//
     &  5x,'Shear Modulus ',1p,e12.5/)

2001  format(5x,'Elastic-Plastic (Mises) St. Venant Torsion Element'//
     &  5x,'Shear Modulus ',1p,e12.5/5x,'Shear Yield   ',1p,e12.5/
     &  5x,'Isotropic H   ',1p,e12.5/)

2002  format( 5x,'Error Factor',1p,e13.5/)

2003  format(a1,20a4//5x,'Element Stress'//' Elmt Matl 1-coord  2-coord'

     &      ,'    1-stress    2-stress    1-strain    2-strain'/
     &   5x,'State    Ratio    Yield  e_p-strain      Torque   Inertia')

2004  format(2i5,0p,2f9.3,1p,4e12.3/3x,a7,1p,2e9.2,1p,3e12.3)

      end

      subroutine stcntor(ix,d,xl,ul,s,shp,dt,st,ser,erav,
     &                  ndf,ndm,nel,numnp)

      implicit   none

      include   'hdata.h'
      include   'iodata.h'

      include   'comblk.h'

      character  region*7
      integer    ndf,ndm,nel,numnp
      integer    i,j,l,ll,lint, nhv,nn
      real*8     rr,xx,yy,xsj,xg,erav

      integer    ix(*)
      real*8     dt(numnp),st(numnp,*),xl(ndm,*),ser(*),d(*),shp(3,*)
      real*8     gam(2),tau(5),ul(ndf,*),s(nel,*),sg(3,9)
      real*8     dd(2,2)

c     Lumped and consistent projection routine

      nn  = 0
      nhv = nint(d(6))
      l   = nint(d(5))
      call int2d(l,lint,sg)
      call pzero(s,nel*nel)
      do l = 1,lint
        call shp2d(sg(1,l),xl,shp,xsj,ndm,nel,ix,.false.)
        xsj = xsj*sg(3,l)

c       Compute strains

        call strator(xl,ul, xx,yy,shp,gam, ndm,ndf)

c       Compute reciprocal radius

        rr = sqrt(xx*xx + yy*yy)
        if(rr.gt.1.d-8) then
          rr = 1.d0/rr
        else
          rr = 1.d0
        endif

c       Compute stresses

        call stretor(d,gam,hr(nh1+nn),hr(nh2+nn),tau,dd, region)

c       Compute consistent projection matrix

        do i = 1,nel
          xg     = shp(3,i)*xsj
          do j = 1,nel
            s(i,j) = s(i,j) + xg*shp(3,j)
          end do
        end do

c       Compute lumped projection and assemble stress integrals

        do j = 1,nel
          ll = abs(ix(j))
          if(ll.gt.0) then
            xg       = xsj*shp(3,j)
            dt(ll)   = dt(ll)   + xg
            st(ll,1) = st(ll,1) + xg*tau(1)
            st(ll,2) = st(ll,2) + xg*tau(2)
            st(ll,3) = st(ll,3) + xg*sqrt(tau(1)**2 + tau(2)**2)
            st(ll,4) = st(ll,4) + xg*(yy*tau(2) + xx*tau(1))*rr
            st(ll,5) = st(ll,5) + xg*(xx*tau(2) - yy*tau(1))*rr
            st(ll,6) = st(ll,6) + xg*tau(3)
            st(ll,11)= st(ll,11)+ xg*tau(5)
            ser(ll)  = ser(ll)  + xg*erav
          endif
        end do
        nn = nn + nhv
      end do

      end

      subroutine strator(xl,ul, xx,yy,shp,gam, ndm,ndf)

c     Compute gamma

      implicit   none

      include   'eldata.h'

      integer    ndm,ndf, i
      real*8     xx,yy
      real*8     xl(ndm,*),ul(ndf,*), shp(3,*),gam(2)

      gam(1) = 0.0d0
      gam(2) = 0.0d0
      xx     = 0.0d0
      yy     = 0.0d0
      do i = 1,nel
        gam(1) = gam(1) + shp(1,i)*ul(1,i)
        gam(2) = gam(2) + shp(2,i)*ul(1,i)
        xx     = xx     + shp(3,i)*xl(1,i)
        yy     = yy     + shp(3,i)*xl(2,i)
      end do

c     Compute gamma

      gam(1) = gam(1) - yy*dm
      gam(2) = gam(2) + xx*dm

      end

      subroutine stertor(ix,d,xl,ul,shp,st,ndf,ndm,nel,numnp)

      implicit   none

      include   'adapt1.h'
      include   'adapt2.h'
      include   'errind.h'
      include   'hdata.h'

      include   'comblk.h'

      character  region*7
      integer    ndf,ndm,nel,numnp, i,j,l,ll,lint, nhv,nn
      real*8     xx,yy,xsj

      integer    ix(*)
      real*8     st(numnp,*),xl(ndm,*),shp(3,*)
      real*8     d(*),gam(2),tau(5),ul(ndf,*),sg(3,9)
      real*8     gradp(2),taup(2),dd(2,2)

c     Bad routine!

      vfem   = 0.0d0
      vproj  = 0.0d0
      verror = 0.0d0
      vener  = 0.0d0
      venere = 0.0d0
      heta   = 0.0d0

c     Quadrature loop

      nn  = 0
      nhv = nint(d(6))
      l   = nint(d(5))
      call int2d(l,lint,sg)
      do l = 1,lint

        call shp2d(sg(1,l),xl,shp,xsj,ndm,nel,ix,.false.)
        xsj = xsj*sg(3,l)

c       Compute strains

        call strator(xl,ul, xx,yy,shp,gam, ndm,ndf)

c       Compute stresses

        call stretor(d,gam,hr(nh1+nn),hr(nh2+nn),tau,dd, region)

        nn = nn + nhv

        do i = 1,2
          taup(i) = 0.0d0
        end do
        do i = 1,nel
          ll = iabs(ix(i))
          if(ll.ne.0) then
            do j = 1,2
              taup(j) = taup(j) + shp(3,i)*st(ll,j)
            end do
          endif
        end do

c       Compute integral of stress squares for error indicator use

        gradp(1) = taup(1)/d(1)
        gradp(2) = taup(2)/d(1)

        heta = heta + xsj
        do i = 1,2
         vfem    = vfem   + tau(i)*tau(i)*xsj
         vproj   = vproj  + taup(i)*taup(i)*xsj
         verror  = verror + ((taup(i)-tau(i))**2)*xsj
         vener   = vener  + tau(i)*gam(i)*xsj
         venere  = venere + (taup(i)-tau(i))*(gradp(i)-gam(i))*xsj
        end do

      end do

      arsq  =  arsq  + heta
      efem  =  efem  + vfem
      eproj =  eproj + vproj
      eerror=  eerror+ verror
      eener =  eener + vener
      eenere=  eenere+ venere
      areai = heta

c     Check for triangle

      if(nel.lt.4 .or. ix(1).eq.ix(2) .or. ix(2).eq.ix(3)
     &            .or. ix(3).eq.ix(4) .or. ix(4).eq.ix(1) ) then
        heta = heta*2.d0
      endif
      heta  =  d(22)*sqrt(heta)

      end

      subroutine stretor(d,gam,epn,ep1,tau,dd, region)

c     Compute tau from strains

c     Model: Elasticity if d(7) .eq. 0.0d0

c       Linear isotropic elastic

c          d(1) = mu

c     Model: Plasticity if d(7) .ne. 0.0d0

c          d(1) = mu
c          d(2) = tau_0
c          d(3) = H_iso

c          epn(*) - plastic internal variables for time t_n
c          ep1(*) - plastic internal variables for time t_n+1
c-----[-+---------+---------+---------+---------+---------+---------+-]

      implicit   none

      include   'counts.h'
      include   'debugs.h'

      character  region*7
      logical    elast
      integer    i,j

      real*8     d(*), gam(*), epn(*), ep1(*), tau(5), dd(2,2)

c     Set elasticity flag

      elast  =  d(7).eq.0.d0 .or. niter.eq.0
      region = 'Elastic'

c     Set isotropic modulus matrix

      do i = 1,2
        do j = 1,2
          dd(i,j) = 0.0d0
        end do ! j
        dd(i,i) = d(1)
      end do ! i

      if(debug) then
        call mprint(dd,2,2,2,'DD_ELAST')
      endif

c     Set yield function and accumulated plastic strain zero

      tau(3) = 0.0d0
      tau(4) = 0.0d0
      tau(5) = 0.0d0

c     Compute elastic stress

      if(elast) then

        do i = 1,2
          tau(i) = d(1)*gam(i)
        end do ! i

c     Compute elastic trial stress and plastic corrector

      else

        do i = 1,2
          tau(i) = d(1)*(gam(i) - epn(i))
        end do ! i

c       Check for plasticity and apply corrector if necessary

        call plastor(d(1),epn,ep1, tau, dd, region)

c       Print elastic plastic tangent

        if(debug) then
          call mprint(dd,2,2,2,'DD_PLAST')
        endif

      endif

      end

      subroutine plastor(d,epn,ep1, sig, dd, region)

      implicit   none

      character  region*7
      integer    i,j
      real*8     mu,ahat,f1,fac1,fac2,fac3,lam
      real*8     tnorm,told,yield
      real*8     d(*),epn(*),ep1(*),sig(5),dd(2,2),tau(2)

c     Yield:  f = || tau || - Y(q)

c             d(1) = mu
c             d(2) = tau_0
c             d(3) = H_iso

c             epn(i), ep1(i) - plastic strains for i=1,2
c             epn(3), ep1(3) - isotropic hardening variable

c     Compute deviator stress

      do i = 1,2
        tau(i) = sig(i)
        ep1(i) = epn(i)
      end do
      ep1(3) = epn(3)

c     Compute ||t|| trial

      tnorm = sqrt(tau(1)**2 + tau(2)**2)

c     Compute the yield function for trial stresses

      yield = d(2)  + d(3)*epn(3)
      f1    = tnorm - yield

c     Set output quantities for elastic state

      sig(3) = tnorm/d(2)
      sig(4) = yield
      sig(5) = ep1(3)

c     Check for yield

      if(f1.gt.0.0d0) then

c       Set elastic modulus

        mu    = d(1)

c       Compute plastic consistency parameter

        ahat  = 1.d0/(mu + d(3))
        lam   = f1*ahat

c       Compute norm of deviator stresses

        told  = tnorm
        tnorm = tnorm - mu*lam

c       Single surface case

        if(tnorm.gt.0.0d0) then

          region = 'Plastic'

c         Compute plastic stresses, update internal variables
c         N.B.,  n(i) is stored in tau(i), q in ep1(3)

          told  = 1.0d0/told
          fac1  = tnorm*told

          do i = 1,2
            sig(i) = tau(i)*fac1
            tau(i) = tau(i)*told
            ep1(i) = ep1(i) + lam*tau(i)
          end do
          ep1(3) = ep1(3) + lam

c         Compute current yield function for outputs

          sig(3) = tnorm/d(2)
          sig(4) = d(2) + d(3)*ep1(3)
          sig(5) = ep1(3)

c         Compute Elastic-Plastic tangent matrix

c         Stiffness factors

          fac1 = mu*lam*told
          fac2 = mu*(mu*ahat - fac1)
          fac1 = mu*(1.d0 - fac1)

c         Plastic tangent

          do i = 1,2
            do j = 1,2
              dd(i,j) = 0.0d0
            end do
            dd(i,i) = fac1
          end do

c         Add rank one update to tangent

          do i = 1,2
            fac3 = fac2*tau(i)
            do j = 1,2
              dd(i,j) = dd(i,j) - fac3*tau(j)
            end do
          end do

        endif
      endif

      end
