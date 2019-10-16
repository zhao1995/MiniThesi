c$Id:$
      subroutine plast_apl(d,eps, hn,nh,istrt, sig,dm,isw)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2011
c       1. Add use of 'istrt' to control iterations         20/05/2013
c-----[--.----+----.----+----.-----------------------------------------]
c     Purpose: Hill48 platicity model (small strain)

c     Input:
c          d(*)    -  Program material parameters   (ndd)
c          eps(*)  -  Current strains at point      (small deformation)
c          hn(nh)  -  History terms at point: t_n
c          nh      -  Number of history terms
c          istrt   -  Start state: 0 = elastic; 1 = last solution
c          isw     -  Solution option from element

c     Output:
c          hn(nh)  -  History terms at point: t_n+1
c          sig(*)  -  Stresses at point.
c                     N.B. 1-d models use only sig(1)
c          dm(6,*) -  Current material tangent moduli
c                     N.B. 1-d models use only dm(1,1) and dm(2,1)
c-----[--.----+----.----+----.-----------------------------------------]
      implicit   none

      include   'comblk.h'
      include   'counts.h'
      include   'debugs.h'
      include   'hdata.h'
      include   'iofile.h'
      include   'pconstant.h'
      include   'setups.h'
      include   'tdata.h'

      logical    cnv,plas, state

      integer    nh,istrt,isw
      integer    i,j,icp

      real*8     aln,al,cl,cln
      real*8     eps(6),d(*),hn(nh)
      real*8     sig(6),dm(6,6)

      real*8     cn(6)
      real*8     sigd(6)
      real*8     ee(6),ss(3,3),spr(3)
      real*8     xi(6,6)
      real*8     epp(6),epn(6),dsig(6),r(6),gf(6)
      real*8     sigy,nsig,fnc,dg,ncn,ddg,nr,nr0,rcn,tol,toln
      real*8     sigyo,sigyi,beta, Ksw,epo,n, fca, hiso, dyld

      real*8     ro(3,3),tt(6,6),dd(6,6),aa(6,6),am(6,6)

      data       tol/1.d-10/, toln /1.d-08/

c     Initialize history during initialization pass

      if(isw.eq.14) then

        if(debug) then
          write(*,*) 'Initializing History'
        endif

c     Anisotropic elasto-plastic computations

      else

c       Check state for iterations

        if(niter.eq.0) then         ! First iteration in step
          if(istrt.eq.0) then       ! Elastic state requested
            state = .false.
            dyld  =  0.0d0
          else                      ! Last state requested
            state = .true.
            dyld  =  1.0d-08*d(41)
          endif
        else                        ! Not first iteration in step
          state = .true.
          dyld  =  0.0d0
          if(rank.gt.0) dyld = -1.0d-08*d(41)
        endif

c       Extract rotation array for orthotropic vectors

        if(nint(d(242)).eq.1) then
          do i = 1,3
            ro(i,1) = d(242+i)
            ro(i,2) = d(245+i)
          end do ! i
          call triad(ro)

c         Set elastic moduli

          dd      = 0.0d0
          dd(1,1) = d(21)
          dd(2,2) = d(22)
          dd(3,3) = d(23)
          dd(1,2) = d(24)
          dd(2,1) = d(24)
          dd(2,3) = d(25)
          dd(3,2) = d(25)
          dd(3,1) = d(26)
          dd(1,3) = d(26)
          dd(4,4) = d(27)
          dd(5,5) = d(28)
          dd(6,6) = d(29)
          call tranr4(ro,ro,tt,.false.)
          call pushr4(tt,tt,dd,dm,1.d0)

c         Set plastic yield function: d(51:56) = F,G,H,L,M,N
c                |  G+H  -H    -G  |
c            A = |  -H   H+F   -F  |
c                |  -G   -F    F+G |

          aa      =  0.0d0
          aa(1,1) =  d(52) + d(53)      ! G + H
          aa(2,2) =  d(53) + d(51)      ! H + F
          aa(3,3) =  d(51) + d(52)      ! F + G

          aa(1,2) = -d(53)              ! H
          aa(2,1) =  aa(1,2)
          aa(2,3) = -d(51)              ! F
          aa(3,2) =  aa(2,3)
          aa(3,1) = -d(52)              ! G
          aa(1,3) =  aa(3,1)

          aa(4,4) =  2.d0*d(56)         ! 2N
          aa(5,5) =  2.d0*d(54)         ! 2L
          aa(6,6) =  2.d0*d(55)         ! 2M

          call trany4(ro,tt)
          call pushr4(tt,tt,aa,am,1.0d0)

        else
          write(*,*) ' ERROR no vectors specified'
        endif

c       Extract plastic strain at t_n
        epn(:) = hn(1:6) ! Grab plastic strains at time t_n
        aln    = hn(7)
        cln    = hn(8)

c       Compute trial elastic strain
        ee(:) = eps(:) - epn(:)

c       Compute Trial Stress

        sig(:) = dm(:,1)*ee(1)+dm(:,2)*ee(2)+dm(:,3)*ee(3)
     &         + dm(:,4)*ee(4)+dm(:,5)*ee(5)+dm(:,6)*ee(6)

c       Trial plastic strain
        epp(:) = epn(:)
        al     = aln

        call invert(dm,6,6)   ! Flip back to get compliances

c       Compute A:sig
        sigd(:) = am(:,1)*sig(1) + am(:,2)*sig(2) + am(:,3)*sig(3)
     &          + am(:,4)*sig(4) + am(:,5)*sig(5) + am(:,6)*sig(6)
        nsig    = sqrt(sig(1)*sigd(1) + sig(2)*sigd(2) + sig(3)*sigd(3)
     &          +      sig(4)*sigd(4) + sig(5)*sigd(5) + sig(6)*sigd(6))

c       Compute yield function and plastic strain residual

        sigyo = d(41)
        if    (nint(d(46)).eq.5) then
          sigyi = d(42)
          beta  = d(43)
          hiso  = d(44)
          sigy  = sigyi - (sigyi-sigyo)*exp(-beta*al) + hiso*al
          Ksw   = 0.0d0
          epo   = 0.0d0
          n     = 0.0d0

c       Swift hardening

        else
          Ksw   = d(42)
          epo   = d(43)
          n     = d(44)
          hiso  = d(45)
          sigy  = Ksw*(epo+al)**n + hiso*al
          sigyi = 0.0d0
          beta  = 0.0d0
        end if

        fnc = nsig - sigy

c       Check for yield
c       cnv = fnc.le.0.d0
        cnv = fnc.le.dyld
        plas = .not.cnv

        if(plas .and. state) then

c         Set up for iterations

          gf(:) = sigd(:)/nsig               ! grad(f)

          dg  = 0.d0
          nr0 = 0.d0
          do i = 1,6
            r(i) = 0.d0                      ! initial strain residual
            nr0  = nr0 + epn(i)*epn(i)
          end do ! i
          if(nr0.ne.0.d0) then
            nr0 = sqrt(nr0)                  ! value to iterate against
          else
            nr0 = 1.d0
          end if

          nr  = 0.d0
          icp = 0

c         Closest Point Projection
          do while(.not.cnv)

            icp = icp + 1

            if(icp.gt.100) then
              write(*,*) 'Local CP projection iterations exceeded 100'
              write(*,*) 'TOL',tol,'NR',nr,nr0,'F',abs(fnc),sigy
              call plstop()
            end if

c           Compute algorithmic moduli
            do j = 1,6
              xi(:,j) = dm(:,j)
     &                + (dg/nsig)*(am(:,j) - gf(:)*gf(j))
            end do ! j

            call invert(xi,6,6)     ! Algorithmic moduli Xi

            cn(:) = xi(:,1)*gf(1) + xi(:,2)*gf(2) + xi(:,3)*gf(3) ! Xi:n
     &            + xi(:,4)*gf(4) + xi(:,5)*gf(5) + xi(:,6)*gf(6)

            ncn = gf(1)*cn(1) + gf(2)*cn(2) + gf(3)*cn(3)       ! n:Xi:n
     &          + gf(4)*cn(4) + gf(5)*cn(5) + gf(6)*cn(6)

            rcn = r(1)*cn(1) + r(2)*cn(2) + r(3)*cn(3)          ! r:Xi:n
     &          + r(4)*cn(4) + r(5)*cn(5) + r(6)*cn(6)

            if    (nint(d(46)).eq.5) then
              fca = -beta*(sigyi-sigyo)*exp(-beta*al) - hiso
            else
              fca = -n*Ksw*(epo+al)**(n-1.d0) - hiso
            end if
            ddg = (fnc-rcn)/(ncn-fca*sqt23)
            dg  =  dg + ddg
            al  =  al + sqt23*ddg

            dsig(:) = -xi(:,1)*(r(1) + ddg*gf(1))
     &                -xi(:,2)*(r(2) + ddg*gf(2))
     &                -xi(:,3)*(r(3) + ddg*gf(3))
     &                -xi(:,4)*(r(4) + ddg*gf(4))
     &                -xi(:,5)*(r(5) + ddg*gf(5))
     &                -xi(:,6)*(r(6) + ddg*gf(6))

            sig(:) = sig(:) + dsig(:)

            epp(:) = epp(:) -dm(:,1)*dsig(1) -dm(:,2)*dsig(2)
     &                      -dm(:,3)*dsig(3) -dm(:,4)*dsig(4)
     &                      -dm(:,5)*dsig(5) -dm(:,6)*dsig(6)

c           Get new residuals

c           Compute A:sig
            sigd(:) = am(:,1)*sig(1) + am(:,2)*sig(2) + am(:,3)*sig(3)
     &              + am(:,4)*sig(4) + am(:,5)*sig(5) + am(:,6)*sig(6)

            nsig = sqrt(sig(1)*sigd(1) + sig(2)*sigd(2) +sig(3)*sigd(3)
     &                + sig(4)*sigd(4) + sig(5)*sigd(5) +sig(6)*sigd(6))

c           Compute yield function and plastic strain residual
            if    (nint(d(46)).eq.5) then
              sigy = sigyi - (sigyi-sigyo)*exp(-beta*al) + hiso*al
            else
              sigy = Ksw*(epo+al)**n + hiso*al
            end if
            fnc = nsig - sigy

c           Compute surface normal and strain residual
            gf(:) =  sigd(:)/nsig
            r(:)  = -epp(:) + epn(:) + dg*gf(:)
            nr = 0.d0
            do i=1,6
              nr    =  nr + r(i)*r(i)
            end do ! i
            nr = sqrt(nr)
c           write(*,*) icp,nr,tol*nr0,abs(fnc),tol*sigy
            cnv = abs(fnc).lt.(tol*sigy) .and. nr.lt.(toln*nr0)

          end do ! while
          if(icp.gt.1) then
          write(*,*) 'RANK',rank,'ICP',icp
          write(*,*) 'TOL',tol,'NR',nr,nr0,'F',abs(fnc),sigy
          endif

c         Compute algorithmic moduli
          do j = 1,6
            xi(:,j) = dm(:,j) ! Xi^{-1} = C^{-1} + dg* grad[grad[f]]
     &              + (dg/nsig)*(am(:,j) - gf(:)*gf(j))
          end do ! j

          call invert(xi,6,6)     ! Algorithmic moduli Xi

          cn(:) = xi(:,1)*gf(1)+xi(:,2)*gf(2)+xi(:,3)*gf(3) ! Xi:n
     &          + xi(:,4)*gf(4)+xi(:,5)*gf(5)+xi(:,6)*gf(6)

          ncn = gf(1)*cn(1)+gf(2)*cn(2)+gf(3)*cn(3)         ! n:Xi:n
     &        + gf(4)*cn(4)+gf(5)*cn(5)+gf(6)*cn(6)

          do j = 1,6 ! Cep = Xi - (Xi:n x Xi:n)/[n:Xi:c-sqrt(2/3)*fca]
            dm(:,j) = xi(:,j)-cn(:)*cn(j)/(ncn-fca*sqt23)
          end do ! j

        else   ! End Plasticity

          call invert(dm,6,6)   ! Elastic case
          dg = 0.0d0

        end if

c       Compute Cockcraft-Lantham parameter

        ss(1,1) = sig(1)
        ss(2,2) = sig(2)
        ss(3,3) = sig(3)
        ss(1,2) = sig(4)
        ss(2,1) = sig(4)
        ss(2,3) = sig(5)
        ss(3,2) = sig(5)
        ss(3,1) = sig(6)
        ss(1,3) = sig(6)
        call eig3(ss,spr,i)

        cl = cln + max(spr(1),spr(2),spr(3),0.d0)*dg*sqt23

c       Update history

        hn(1:6) = epp(:)
        hn(7)   = al
        hn(8)   = cl

        if (debug) call mprint(sig,6,1,6,'sig_n+1')

      endif

      end
